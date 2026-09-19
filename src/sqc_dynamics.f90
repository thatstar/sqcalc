! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Dynamic structure factor S(q,w) along a line in reciprocal space.
!!
!! The method follows the recipe used by dynasor and MDANSE: evaluate the
!! density amplitude rho_a(q,t) of every species at a few q points, correlate
!! it in time and Fourier transform the correlation.  Nothing here uses the
!! reciprocal grid:
!!
!!   * the q line is arbitrary (off lattice q included), so the amplitudes are
!!     summed directly, O(N * nq) per frame - for the ~100 q points of a line
!!     scan that is far cheaper than transforming the whole grid;
!!   * the correlation is accumulated in a ring buffer, so the memory does not
!!     depend on the length of the trajectory.
!!
!! Definitions (W(q) is the same normalization the static table uses):
!!
!!   C_ab(q,l) = <rho_a(q,t+l) rho_b*(q,t)>            l = 0..maxframes
!!   F(q,l)    = sum_ab w_a(q) w_b(q) C_ab(q,l) / W(q)
!!   S(q,w_n)  = (dt/2pi) sum_l F(q,l) exp(i w_n l dt)      (even extension)
!!
!! The sum over time origins carries its own count per lag (`C_cnt`), which is
!! the "nanmean" bookkeeping: a trajectory whose length is not a multiple of
!! the window needs no padding and its ragged tail cannot bias a lag.  With
!! the one-sided folding the zeroth moment is exact,
!!
!!   sum_n S(q,w_n) dw = F(q,0) = S(q) = the OUTPUT table of the same run,
!!
!! which is the primary consistency test of the method.
module sqc_dynamics
   use sqc_kinds
   use sqc_dump, only: frame_t
   use sqc_weights, only: weight_scheme_t
   use sqc_structure_factor, only: structure_factor_t, sf_shared_setup, sf_alloc_partials, &
                                   sf_prepare_species, mode_denominator, norm_self, &
                                   norm_natom, no_unit
   use, intrinsic :: iso_c_binding, only: c_double_complex
   use, intrinsic :: iso_fortran_env, only: output_unit
   use omp_lib, only: omp_get_max_threads
#ifdef SQC_HAVE_HDF5
   use sqc_hdf5, only: hdf5_write_dynamics
#endif
   implicit none
   private

   public :: dynamics_structure_factor_t, dyn_format_text, dyn_format_hdf5

   !> Output flavours of the two optional files.
   integer, parameter :: dyn_format_text = 0
   integer, parameter :: dyn_format_hdf5 = 1

   !> Which quantity a writer is asked for.
   integer, parameter :: dyn_sqw = 1
   integer, parameter :: dyn_fsq = 2

   !> Direct summation over a q line plus a multi-origin time correlation.
   type, extends(structure_factor_t) :: dynamics_structure_factor_t
      !> Time between two dumped frames [same unit as `--dt`].
      real(rk) :: frame_dt = 0.0_rk
      !> Correlation window in frames (the largest lag kept).
      integer :: maxframes = 0
      !> Frames between two consecutive time origins.
      integer :: lag_stride = 1
      !> The q line: `nintervals` intervals from `s0` to `s1` [1/A] along
      !! `direction` (unnormalized), through the Gamma point.
      integer :: nintervals = 0
      real(rk) :: s0 = 0.0_rk, s1 = 0.0_rk
      real(rk) :: direction(3) = 0.0_rk
      !> Optional output files (empty = not written) and their format.
      character(len=:), allocatable :: sqw_path, fsq_path
      integer :: sqw_format = dyn_format_text
      !> Species present and their grouping by contiguous atom ranges.
      integer(ik), allocatable :: species_of(:), species_type(:)
      integer(ik), allocatable :: atom_of(:), species_first(:)
      integer :: nspecies = 0
      !> Amplitude w_a(q) of each species at each q (weights applied at output).
      real(rk), allocatable :: amp(:, :)
      logical :: q_dependent = .false.
      real(rk), allocatable :: amp_const(:), amp_table(:, :)
      !> Ring buffer rho(nmodes, nspecies, 0:maxframes).
      complex(c_double_complex), allocatable :: rho(:, :, :)
      !> Correlation sums and counts per (mode, lag, species, species).
      complex(c_double_complex), allocatable :: csum(:, :, :, :)
      integer(lk), allocatable :: ccnt(:, :)
      !> Zero lag sums over the whole trajectory (the static OUTPUT table).
      complex(c_double_complex), allocatable :: c0sum(:, :, :)
      !> Results: the spectra and the intermediate scattering function.
      real(rk), allocatable :: omega(:), sqw(:, :), sqw_partial(:, :, :)
      real(rk), allocatable :: tau(:), ftau(:, :), ftau_partial(:, :, :)
      !> Frequency grid spacing, pi/(maxframes*frame_dt).
      real(rk) :: dw = 0.0_rk
   contains
      procedure :: configure => dyn_configure
      procedure :: method_setup => dyn_setup
      procedure :: accumulate_frame => dyn_accumulate
      procedure :: prepare_output => dyn_prepare_output
      procedure :: write_sqw => dyn_write_sqw
      procedure :: write_fsq => dyn_write_fsq
   end type dynamics_structure_factor_t

contains

   !> Build the q line and the shared accumulators from the first frame.
   subroutine dyn_configure(self, frame, scheme, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: u(3), norm_u, s_i
      integer :: i

      ierr = 0
      message = ''
      call sf_shared_setup(self, frame, scheme, ierr, message)
      if (ierr /= 0) return

      if (self%nintervals < 1) then
         ierr = 1
         message = '--dyn needs at least one interval on the q line'
         return
      end if
      norm_u = sqrt(sum(self%direction**2))
      if (norm_u <= 0.0_rk) then
         ierr = 1
         message = '--dyn needs a non-zero direction, e.g. 1,1,0'
         return
      end if
      u = self%direction/norm_u

      ! q_i = s_i u, anchored at the Gamma point.
      self%nmodes = int(self%nintervals + 1, lk)
      self%nq = int(self%nmodes)
      allocate (self%qvec(3, self%nmodes), self%qlen(self%nmodes), &
                self%shell(self%nmodes), self%gidx(self%nmodes))
      do i = 1, int(self%nmodes)
         s_i = self%s0 + real(i - 1, rk)*(self%s1 - self%s0)/real(self%nintervals, rk)
         self%qvec(:, i) = s_i*u
         self%qlen(i) = s_i
         self%shell(i) = i
         self%gidx(i) = int(i, lk)
      end do
      ! The shared writer prints qmin + (shell-1/2)*shell_dq, so these two make
      ! it print the q line itself.
      self%shell_dq = (self%s1 - self%s0)/real(self%nintervals, rk)
      self%qmin = self%s0 - 0.5_rk*self%shell_dq
      self%qmax = self%s1 + 0.5_rk*self%shell_dq

      allocate (self%num(self%nq), self%den(self%nq), self%shell_count(self%nq))
      allocate (self%mode_values(self%nmodes))
      self%num = 0.0_rk
      self%den = 0.0_rk
      self%mode_values = 0.0_rk
      self%shell_count = 1
      do i = 1, int(self%nmodes)
         self%den(i) = mode_denominator(self%norm, scheme, self%natoms, self%type_counts, &
                                        self%qlen(i))
      end do
      call sf_alloc_partials(self, scheme)

      call self%method_setup(frame, scheme, ierr, message)
   end subroutine dyn_configure

   !> Species bookkeeping, weights per q and the correlation buffers.
   subroutine dyn_setup(self, frame, scheme, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(ik), allocatable :: cursor(:)
      integer :: i, t, isp

      ierr = 0
      message = ''
      call sf_prepare_species(self, frame, scheme, self%species_of, self%species_type, &
                              self%amp_const, self%amp_table, self%q_dependent, ierr, message)
      if (ierr /= 0) return
      self%nspecies = size(self%species_type)

      ! Contiguous atom ranges per species, as in the direct method.
      allocate (self%species_first(self%nspecies + 1))
      self%species_first = 1
      do i = 1, frame%natoms
         t = int(self%species_of(i))
         self%species_first(t + 1) = self%species_first(t + 1) + 1
      end do
      do t = 1, self%nspecies
         self%species_first(t + 1) = self%species_first(t + 1) + self%species_first(t) - 1
      end do
      allocate (self%atom_of(frame%natoms))
      allocate (cursor(self%nspecies))
      cursor = self%species_first(1:self%nspecies)
      do i = 1, frame%natoms
         t = int(self%species_of(i))
         self%atom_of(cursor(t)) = int(i, ik)
         cursor(t) = cursor(t) + 1
      end do
      deallocate (cursor)

      ! Weight of every species at every q, applied when the results are
      ! assembled (never inside the per frame sum).
      allocate (self%amp(self%nmodes, self%nspecies))
      if (self%q_dependent) then
         self%amp = self%amp_table
      else
         do isp = 1, self%nspecies
            self%amp(:, isp) = self%amp_const(isp)
         end do
      end if

      allocate (self%rho(self%nmodes, self%nspecies, 0:self%maxframes))
      self%rho = (0.0_rk, 0.0_rk)
      allocate (self%csum(self%nmodes, 0:self%maxframes, self%nspecies, self%nspecies))
      self%csum = (0.0_rk, 0.0_rk)
      allocate (self%c0sum(self%nmodes, self%nspecies, self%nspecies))
      self%c0sum = (0.0_rk, 0.0_rk)
      allocate (self%ccnt(self%nmodes, 0:self%maxframes))
      self%ccnt = 0
   end subroutine dyn_setup

   !> Add one frame: the density amplitudes and their time correlations.
   subroutine dyn_accumulate(self, frame, scheme, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      complex(c_double_complex) :: acc
      real(rk) :: qx, qy, qz, phase
      integer :: t_frame, slot, l, oslot, im, isp, jsp, i

      ierr = 0
      message = ''
      if (.not. frame%unwrapped) then
         ierr = 1
         message = 'the dynamic method needs unwrapped coordinates'
         return
      end if
      if (frame%natoms /= self%natoms) then
         ierr = 1
         message = 'the number of atoms changed between frames'
         return
      end if

      t_frame = int(self%nframes)
      slot = mod(t_frame, self%maxframes + 1)

      ! --- rho_a(q,t) by direct summation over the q line -------------------
      !$omp parallel do schedule(static) private(im, isp, i, acc, phase, qx, qy, qz)
      do im = 1, int(self%nmodes)
         qx = self%qvec(1, im)
         qy = self%qvec(2, im)
         qz = self%qvec(3, im)
         do isp = 1, self%nspecies
            acc = (0.0_rk, 0.0_rk)
            do i = self%species_first(isp), self%species_first(isp + 1) - 1
               phase = qx*frame%pos(1, self%atom_of(i)) + qy*frame%pos(2, self%atom_of(i)) &
                       + qz*frame%pos(3, self%atom_of(i))
               acc = acc + cmplx(cos(phase), sin(phase), c_double_complex)
            end do
            self%rho(im, isp, slot) = acc
         end do
      end do
      !$omp end parallel do

      ! --- zero lag sums over the whole trajectory --------------------------
      do isp = 1, self%nspecies
         do jsp = 1, self%nspecies
            do im = 1, int(self%nmodes)
               self%c0sum(im, isp, jsp) = self%c0sum(im, isp, jsp) &
                  + self%rho(im, isp, slot)*conjg(self%rho(im, jsp, slot))
            end do
         end do
      end do

      ! --- multi-origin correlation of the lags -----------------------------
      !$omp parallel do schedule(static) private(l, oslot, im, isp, jsp)
      do l = 0, min(t_frame, self%maxframes)
         if (mod(t_frame - l, self%lag_stride) /= 0) cycle
         oslot = mod(t_frame - l, self%maxframes + 1)
         do im = 1, int(self%nmodes)
            do isp = 1, self%nspecies
               do jsp = 1, self%nspecies
                  self%csum(im, l, isp, jsp) = self%csum(im, l, isp, jsp) &
                     + self%rho(im, isp, slot)*conjg(self%rho(im, jsp, oslot))
               end do
            end do
         end do
      end do
      !$omp end parallel do
      do l = 0, min(t_frame, self%maxframes)
         if (mod(t_frame - l, self%lag_stride) /= 0) cycle
         self%ccnt(:, l) = self%ccnt(:, l) + 1
      end do

      self%nframes = self%nframes + 1
   end subroutine dyn_accumulate

   !> Assemble the static table, F(q,tau) and the spectra.
   subroutine dyn_prepare_output(self, scheme, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: f(:)
      real(rk) :: frames, value
      integer :: im, l, isp, jsp, t, u, p, n_pairs

      ierr = 0
      message = ''
      if (self%nframes <= self%maxframes) then
         write (message, '(a,i0,a,i0,a)') 'the trajectory has ', self%nframes, &
            ' frames but --maxframes is ', self%maxframes, &
            '; at least maxframes+1 frames are needed'
         ierr = 1
         return
      end if
      frames = real(self%nframes, rk)
      self%dw = acos(-1.0_rk)/(real(self%maxframes, rk)*self%frame_dt)

      allocate (self%omega(0:self%maxframes), self%tau(0:self%maxframes))
      do l = 0, self%maxframes
         self%tau(l) = real(l, rk)*self%frame_dt
         self%omega(l) = real(l, rk)*self%dw
      end do
      allocate (self%sqw(self%nmodes, 0:self%maxframes))
      allocate (self%ftau(self%nmodes, 0:self%maxframes))
      self%sqw = 0.0_rk
      self%ftau = 0.0_rk

      ! Static table: |rho(q)|^2 averaged over every frame, with the same
      ! normalization as the other methods.
      do im = 1, int(self%nmodes)
         value = 0.0_rk
         do isp = 1, self%nspecies
            do jsp = 1, self%nspecies
               value = value + self%amp(im, isp)*self%amp(im, jsp) &
                       *real(self%c0sum(im, isp, jsp), rk)
            end do
         end do
         self%num(im) = value
      end do
      if (self%partials) then
         self%partial_num = 0.0_rk
         do isp = 1, self%nspecies
            do jsp = 1, self%nspecies
               t = int(self%species_type(isp))
               u = int(self%species_type(jsp))
               do im = 1, int(self%nmodes)
                  self%partial_num(t, u, im) = real(self%c0sum(im, isp, jsp), rk)
               end do
            end do
         end do
      end if

      ! F(q,tau), and the partial intermediate functions.
      n_pairs = 0
      do t = 1, self%ntypes
         do u = t, self%ntypes
            if (.not. pair_present(self, t, u)) cycle
            n_pairs = n_pairs + 1
         end do
      end do
      if (self%partials) then
         allocate (self%ftau_partial(self%nmodes, 0:self%maxframes, max(n_pairs, 1)))
         allocate (self%sqw_partial(self%nmodes, 0:self%maxframes, max(n_pairs, 1)))
         self%ftau_partial = 0.0_rk
         self%sqw_partial = 0.0_rk
      else
         allocate (self%ftau_partial(0, 0, 0), self%sqw_partial(0, 0, 0))
      end if

      allocate (f(0:self%maxframes))
      do im = 1, int(self%nmodes)
         ! total
         do l = 0, self%maxframes
            if (l == 0) then
               value = 0.0_rk
               do isp = 1, self%nspecies
                  do jsp = 1, self%nspecies
                     value = value + self%amp(im, isp)*self%amp(im, jsp) &
                             *real(self%c0sum(im, isp, jsp), rk)/frames
                  end do
               end do
            else
               value = 0.0_rk
               if (self%ccnt(im, l) > 0) then
                  do isp = 1, self%nspecies
                     do jsp = 1, self%nspecies
                        value = value + self%amp(im, isp)*self%amp(im, jsp) &
                                *real(self%csum(im, l, isp, jsp), rk)/real(self%ccnt(im, l), rk)
                     end do
                  end do
               end if
            end if
            if (self%den(im) > 0.0_rk) then
               self%ftau(im, l) = value/self%den(im)
            else
               self%ftau(im, l) = 0.0_rk
            end if
         end do
         f = self%ftau(im, :)
         call dyn_spectrum(f, self%frame_dt, self%sqw(im, :))

         ! partials, in the OVITO (1/N) convention of the static columns
         if (self%partials) then
            p = 0
            do t = 1, self%ntypes
               do u = t, self%ntypes
                  if (.not. pair_present(self, t, u)) cycle
                  p = p + 1
                  isp = species_index(self, t)
                  jsp = species_index(self, u)
                  do l = 0, self%maxframes
                     if (l == 0) then
                        value = real(self%c0sum(im, isp, jsp), rk)/frames
                     else
                        value = 0.0_rk
                        if (self%ccnt(im, l) > 0) then
                           value = real(self%csum(im, l, isp, jsp), rk)/real(self%ccnt(im, l), rk)
                           ! The cross partial has to be symmetrized, otherwise
                           ! the OVITO sum rule S = S_aa + 2 S_ab + S_bb would
                           ! hold only at zero lag (C_ab and C_ba differ by the
                           ! finite-trajectory noise).
                           if (isp /= jsp) value = 0.5_rk*(value &
                              + real(self%csum(im, l, jsp, isp), rk)/real(self%ccnt(im, l), rk))
                        end if
                     end if
                     self%ftau_partial(im, l, p) = value/real(self%natoms, rk)
                  end do
                  f = self%ftau_partial(im, :, p)
                  call dyn_spectrum(f, self%frame_dt, self%sqw_partial(im, :, p))
               end do
            end do
         end if
      end do
      deallocate (f)
   end subroutine dyn_prepare_output

   !> One-sided spectrum of a real, even correlation function.
   !!
   !! `f(l)` is F(q, l*dt) for l = 0..L; the even extension of length 2L is
   !! transformed with the trapezoid rule and the interior points are doubled
   !! (S(q,-w) = S(q,w)), so that sum_n spec(n)*dw = f(0) holds exactly.
   pure subroutine dyn_spectrum(f, dt, spec)
      real(rk), intent(in) :: f(0:)
      real(rk), intent(in) :: dt
      real(rk), intent(out) :: spec(0:)
      real(rk) :: acc, pi
      integer :: l_max, n, l

      pi = acos(-1.0_rk)
      l_max = ubound(f, 1)
      do n = 0, l_max
         acc = f(0)
         do l = 1, l_max - 1
            acc = acc + 2.0_rk*f(l)*cos(pi*real(n*l, rk)/real(l_max, rk))
         end do
         if (l_max > 0) then
            if (mod(n, 2) == 0) then
               acc = acc + f(l_max)
            else
               acc = acc - f(l_max)
            end if
         end if
         spec(n) = acc*dt/(2.0_rk*pi)
      end do
      do n = 1, l_max - 1
         spec(n) = 2.0_rk*spec(n)
      end do
   end subroutine dyn_spectrum

   !> Number of atoms of a LAMMPS type id (0 when the type is absent).
   pure integer function count_pair(self, type_id) result(n)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer, intent(in) :: type_id
      n = 0
      if (allocated(self%type_counts)) then
         if (type_id >= 1 .and. type_id <= size(self%type_counts)) n = int(self%type_counts(type_id))
      end if
   end function count_pair

   !> True when a pair of type ids has the partial column the base writer uses.
   pure logical function pair_present(self, ia, ib) result(present)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer, intent(in) :: ia, ib
      present = count_pair(self, ia) > 0 .and. count_pair(self, ib) > 0
   end function pair_present

   !> Number of type pairs that get a column (`npair` of the base writer).
   pure integer function npairs_out(self) result(n)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer :: t, u
      n = 0
      do t = 1, self%ntypes
         do u = t, self%ntypes
            if (.not. pair_present(self, t, u)) cycle
            n = n + 1
         end do
      end do
   end function npairs_out

   !> Species index of a LAMMPS type id (0 when the type is absent).
   pure integer function species_index(self, type_id) result(idx)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer, intent(in) :: type_id
      idx = 0
      if (allocated(self%species_type)) then
         idx = findloc(self%species_type, int(type_id, ik), dim=1)
      end if
   end function species_index

   !> The S(q,w) table (text or HDF5).
   subroutine dyn_write_sqw(self, path, format, scheme, input, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input
      integer, intent(in) :: format
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message

      call dyn_write_table(self, path, format, scheme, input, dyn_sqw, 'omega', &
                           self%omega, self%sqw, self%sqw_partial, ierr, message)
   end subroutine dyn_write_sqw

   !> The F(q,tau) table (text or HDF5).
   subroutine dyn_write_fsq(self, path, format, scheme, input, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input
      integer, intent(in) :: format
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message

      call dyn_write_table(self, path, format, scheme, input, dyn_fsq, 'tau', &
                           self%tau, self%ftau, self%ftau_partial, ierr, message)
   end subroutine dyn_write_fsq

   !> Text or HDF5 writer shared by the spectra and F(q,tau).
   subroutine dyn_write_table(self, path, format, scheme, input, kind, axis_name, &
                              axis, spec, part, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input, axis_name
      integer, intent(in) :: format, kind
      type(weight_scheme_t), intent(in) :: scheme
      real(rk), intent(in) :: axis(0:), spec(:, 0:), part(:, 0:, :)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer :: unit, im, n, p, t, u

      ierr = 0
      message = ''
      if (format == dyn_format_hdf5) then
#ifdef SQC_HAVE_HDF5
         call dyn_write_hdf5(self, scheme, path, kind, axis_name, axis, spec, part, ierr, message)
#else
         ierr = 1
         message = 'this build has no HDF5 support; use --sqw-format text'
#endif
         return
      end if

      if (trim(path) == '-') then
         unit = output_unit
      else
         open (newunit=unit, file=trim(path), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            message = 'cannot write the dynamic output "'//trim(path)//'"'
            return
         end if
      end if

      call dyn_header(self, scheme, input, kind, axis_name, unit)
      do im = 1, int(self%nmodes)
         do n = 0, self%maxframes
            write (unit, '(3(f14.8,2x),f16.8,2x,es20.12)', advance='no') &
               self%qvec(1, im), self%qvec(2, im), self%qvec(3, im), axis(n), spec(im, n)
            if (self%partials) then
               p = 0
               do t = 1, self%ntypes
                  do u = t, self%ntypes
                     if (.not. pair_present(self, t, u)) cycle
                     p = p + 1
                     write (unit, '(2x,es20.12)', advance='no') part(im, n, p)
                  end do
               end do
            end if
            write (unit, '(a)') ''
         end do
      end do
      if (unit /= output_unit) close (unit)
   end subroutine dyn_write_table

   !> Metadata preamble of the text tables.
   subroutine dyn_header(self, scheme, input, kind, axis_name, unit)
      class(dynamics_structure_factor_t), intent(in) :: self
      type(weight_scheme_t), intent(in) :: scheme
      character(len=*), intent(in) :: input, axis_name
      integer, intent(in) :: kind, unit
      real(rk) :: unorm
      integer :: t, u
      character(len=24) :: label

      unorm = sqrt(sum(self%direction**2))
      if (kind == dyn_sqw) then
         write (unit, '(a)') '# sqcalc 0.1.0 dynamic structure factor S(q,w)'
      else
         write (unit, '(a)') '# sqcalc 0.1.0 intermediate scattering function F(q,t)'
      end if
      write (unit, '(a)') '# input '//trim(input)//'  weight '//trim(scheme%label())// &
         '  norm '//trim(norm_name(self%norm))
      if (self%q_dependent) then
         write (unit, '(a)') '# note the amplitudes depend on q (X-ray form factors)'
      end if
      write (unit, '(a,f12.6,a,i0,a,i0,a,i0)') '# frame_dt ', self%frame_dt, &
         '  maxframes ', self%maxframes, '  lag ', self%lag_stride, '  nframes ', self%nframes
      write (unit, '(a,i0,a,f10.6,a,f10.6,a)') '# q line ', self%nintervals, &
         ' intervals: |q| ', self%s0, ' .. ', self%s1, ' 1/A (through Gamma)'
      write (unit, '(a,3(f12.8,1x))') '# direction ', self%direction/unorm
      if (kind == dyn_sqw) then
         write (unit, '(a,i0,a,f0.6,a,f0.6,a)') '# nomega ', self%maxframes + 1, &
            '  omega 0 .. ', real(self%maxframes, rk)*self%dw, ' dw ', self%dw, &
            ' (1/time unit of dt, one-sided)'
         write (unit, '(a)') '# S(q,-w) = S(q,w);  sum_n S(q,w_n)*dw = F(q,0) '// &
            '= the OUTPUT table of this run'
      else
         write (unit, '(a,i0,a,f0.6,a)') '# nlag ', self%maxframes + 1, &
            '  tau 0 .. ', real(self%maxframes, rk)*self%frame_dt, ' (time unit of dt)'
         write (unit, '(a)') '# F(q,0) = the OUTPUT table of this run; the spectra are '// &
            'its Fourier transform'
      end if
      write (unit, '(a)', advance='no') '# qx qy qz '//trim(axis_name)//' '
      if (kind == dyn_sqw) then
         write (unit, '(a)', advance='no') 'S(q,w)'
      else
         write (unit, '(a)', advance='no') 'F(q,t)'
      end if
      if (self%partials) then
         do t = 1, self%ntypes
            do u = t, self%ntypes
               if (.not. pair_present(self, t, u)) cycle
               label = self%partial_column(t, u)
               write (unit, '(a)', advance='no') ' '//trim(label)
            end do
         end do
      end if
      write (unit, '(a)') ''
   end subroutine dyn_header

#ifdef SQC_HAVE_HDF5
   !> HDF5 flavour: one group per quantity, mirroring the shell/rdf writers.
   subroutine dyn_write_hdf5(self, scheme, path, kind, axis_name, axis, spec, part, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      type(weight_scheme_t), intent(in) :: scheme
      character(len=*), intent(in) :: path, axis_name
      integer, intent(in) :: kind
      real(rk), intent(in) :: axis(0:), spec(:, 0:), part(:, 0:, :)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: q4(:, :)
      character(len=18), allocatable :: labels(:)
      character(len=:), allocatable :: group
      integer :: im, p, t, u, npair

      ierr = 0
      message = ''
      group = 'sqw'
      if (kind == dyn_fsq) group = 'fsq'
      allocate (q4(self%nmodes, 4))
      do im = 1, int(self%nmodes)
         q4(im, 1:3) = self%qvec(:, im)
         q4(im, 4) = self%qlen(im)
      end do
      npair = npairs_out(self)
      allocate (labels(max(npair, 1)))
      labels = ' '
      p = 0
      do t = 1, self%ntypes
         do u = t, self%ntypes
            if (.not. pair_present(self, t, u)) cycle
            p = p + 1
            labels(p) = scheme%pair_label(t, u)
         end do
      end do
      call hdf5_write_dynamics(path, group, axis_name, q4, axis, spec, part, labels, self%ccnt, &
                               self%nframes, self%frame_dt, self%maxframes, self%lag_stride, &
                               trim(scheme%label()), trim(norm_name(self%norm)), ierr, message)
      deallocate (q4, labels)
   end subroutine dyn_write_hdf5
#endif

   !> Text label of a normalization code.
   pure function norm_name(norm) result(label)
      integer, intent(in) :: norm
      character(len=8) :: label
      select case (norm)
      case (norm_self)
         label = 'self'
      case (norm_natom)
         label = 'n'
      case default
         label = 'mean'
      end select
   end function norm_name

end module sqc_dynamics
