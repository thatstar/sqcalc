! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Structure factor accumulation.
!!
!! All methods evaluate the scattering amplitude on the reciprocal lattice of
!! the simulation cell,
!!
!!    rho(q) = sum_j w_j exp(i q . r_j),
!!
!! where w_j is the per-atom weight (unit, neutron or X-ray) and q = h b1 +
!! k b2 + l b3.  |rho(q)|^2 is averaged over the trajectory and normalized to
!! give S(q):
!!
!!    norm_mean : S = <|rho|^2> / (N <w>^2)      (Faber-Ziman total S(q))
!!    norm_self : S = <|rho|^2> / sum_j w_j^2    (S -> 1 at large q)
!!    norm_natom: S = <|rho|^2> / N              (debyer style)
!!
!! The 1D output is the average of S over the q shells [qmin, qmax]; the
!! optional grid output writes S at every reciprocal lattice point of that
!! range.
module sqc_structure_factor
   use sqc_kinds
   use sqc_cell, only: cell_t, two_pi
   use sqc_dump, only: frame_t
   use sqc_weights, only: weight_scheme_t, weight_xray
   use sqc_finufft, only: finufft_opts_t, finufft_default_opts, finufft_makeplan, &
                          finufft_setpts, finufft_execute, finufft_destroy, &
                          mode_cmcl, type1
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_double, c_double_complex, &
                                          c_int, c_int64_t, c_associated
   use omp_lib, only: omp_get_max_threads, omp_get_thread_num
   implicit none
   private

   public :: structure_factor_t, nufft_structure_factor_t, direct_structure_factor_t, &
             method_nufft, method_direct, norm_mean, norm_self, norm_natom, no_unit, &
             sf_prepare_species, mode_denominator, method_debye

   !> Sentinel for "do not write this table".  A plain negative test would be
   !! wrong because OPEN(NEWUNIT=) may hand out negative unit numbers.
   integer, parameter :: no_unit = -huge(1)

   !> Available evaluation methods.
   integer, parameter :: method_nufft = 1
   integer, parameter :: method_direct = 2
   integer, parameter :: method_debye = 3

   !> Normalization conventions.
   integer, parameter :: norm_mean = 1
   integer, parameter :: norm_self = 2
   integer, parameter :: norm_natom = 3

   !> Safety limit on the number of grid modes (memory: 16 bytes per mode).
   integer(lk), parameter :: max_grid_modes = 400000000_lk

   !> Common configuration, reciprocal grid bookkeeping and accumulators.
   type, abstract :: structure_factor_t
      !> Normalization convention (norm_mean, norm_self, norm_natom).
      integer :: norm = norm_mean
      !> Number of q shells in the 1D output.
      integer :: nq = 500
      !> Shell range.
      real(rk) :: qmin = 0.0_rk
      real(rk) :: qmax = 20.0_rk
      !> Requested FINUFFT tolerance.
      real(rk) :: eps = 1.0e-9_rk
      !> Number of OpenMP threads (0 = inherit the environment).
      integer :: nthreads = 0
      !> Also accumulate S(q) on every reciprocal lattice point.
      logical :: want_grid = .false.
      !> Reciprocal grid dimensions (odd, index h = p - (m+1)/2).
      integer :: modes(3) = [1, 1, 1]
      !> Total number of grid points (product of modes).
      integer(lk) :: gridpoints = 1
      !> Number of modes kept in the output range.
      integer(lk) :: nmodes = 0
      !> Shell width, (qmax - qmin)/nq.
      real(rk) :: shell_dq = 1.0_rk
      !> Miller indices of the kept modes.
      integer, allocatable :: hkl(:, :)
      !> Cartesian q vector and its length for each kept mode.
      real(rk), allocatable :: qvec(:, :), qlen(:)
      !> Shell index of each kept mode (1..nq).
      integer(ik), allocatable :: shell(:)
      !> Position of each kept mode in the flat grid arrays.
      integer(lk), allocatable :: gidx(:)
      !> Accumulated numerator and denominator of S(q) per shell.
      real(rk), allocatable :: num(:), den(:)
      !> Scratch: |rho(q)|^2 of each kept mode for the current frame.
      real(rk), allocatable :: mode_values(:)
      !> Write partial structure factors S_ab(q) next to the total.
      logical :: partials = .false.
      !> Report the partials in the Faber-Ziman normalization.
      logical :: faber_ziman = .false.
      !> Accumulated, unweighted partial sums Re(rho_a rho_b*) per shell.
      real(rk), allocatable :: partial_num(:, :, :)
      !> Human readable label of each type pair (for the partial columns).
      character(len=18), allocatable :: pair_label(:, :)
      !> Number of (frame, mode) contributions per shell.
      integer(lk), allocatable :: shell_count(:)
      !> Accumulated numerator and denominator per grid point.
      real(rk), allocatable :: gnum(:), gden(:)
      !> Frame bookkeeping.
      integer(ik) :: natoms = 0
      integer(ik) :: ntypes = 0
      !> Number of atoms per type id.
      integer(ik), allocatable :: type_counts(:)
      integer(lk) :: nframes = 0
   contains
      procedure :: configure => sf_configure
      procedure :: write_results => sf_write_results
      procedure :: accumulate_modes => sf_accumulate_modes
      procedure :: accumulate_values => sf_accumulate_values
      procedure :: accumulate_partials => sf_accumulate_partials
      procedure :: prepare_output => sf_prepare_output
      procedure :: shell_value => sf_shell_value
      procedure :: partial_value => sf_partial_value
      procedure :: partial_column => sf_partial_column
      procedure :: grid_value => sf_grid_value
      procedure(sf_setup_iface), deferred :: method_setup
      procedure(sf_accumulate_iface), deferred :: accumulate_frame
   end type structure_factor_t

   abstract interface
      !> Per-method set up driven by the first frame of the trajectory.
      subroutine sf_setup_iface(self, frame, scheme, ierr, message)
         import :: structure_factor_t, frame_t, weight_scheme_t
         class(structure_factor_t), intent(inout) :: self
         type(frame_t), intent(in) :: frame
         type(weight_scheme_t), intent(in) :: scheme
         integer, intent(out) :: ierr
         character(len=*), intent(out) :: message
      end subroutine sf_setup_iface

      !> Add one frame to the accumulators.
      subroutine sf_accumulate_iface(self, frame, scheme, ierr, message)
         import :: structure_factor_t, frame_t, weight_scheme_t
         class(structure_factor_t), intent(inout) :: self
         type(frame_t), intent(in) :: frame
         type(weight_scheme_t), intent(in) :: scheme
         integer, intent(out) :: ierr
         character(len=*), intent(out) :: message
      end subroutine sf_accumulate_iface
   end interface

   !> NUFFT based evaluation: one type-1 transform per species and frame.
   type, extends(structure_factor_t) :: nufft_structure_factor_t
      !> FINUFFT plan handle.
      type(c_ptr) :: plan = c_null_ptr
      !> FINUFFT options used for the plan.
      type(finufft_opts_t) :: opts
      !> Species index of each atom (0 when its type has no atoms).
      integer(ik), allocatable :: species_of(:)
      !> LAMMPS type id of each species.
      integer(ik), allocatable :: species_type(:)
      !> Scaled coordinates x_j in [-pi, pi) for all atoms, one array per axis:
      !! the FINUFFT C interface needs contiguous data, so a (3, natoms) array
      !! sliced by rows must not be passed directly.
      real(rk), allocatable :: xa(:), ya(:), za(:)
      !> True when the amplitude depends on q (X-ray).
      logical :: q_dependent = .false.
      !> Precomputed amplitude of each species at each kept mode.
      real(rk), allocatable :: amp_table(:, :)
      !> Constant amplitude per species for q-independent weights.
      real(rk), allocatable :: amp_const(:)
      !> Per-atom strengths and transform scratch space.
      complex(c_double_complex), allocatable :: strengths(:), fk(:), total(:)
      !> Per-species amplitudes on the flat grid (only with --partials).
      complex(c_double_complex), allocatable :: rho_species(:, :)
   contains
      procedure :: method_setup => nufft_setup
      procedure :: accumulate_frame => nufft_accumulate
      final :: nufft_finalize
   end type nufft_structure_factor_t

   !> Direct O(N * nmodes) evaluation, used as reference and offline fallback.
   type, extends(structure_factor_t) :: direct_structure_factor_t
      !> LAMMPS type id of each atom.
      integer(ik), allocatable :: type_of(:)
      !> Atoms grouped by type (species_first(t) .. species_first(t+1)-1).
      integer(ik), allocatable :: atom_of(:)
      integer(ik), allocatable :: species_first(:)
      !> Cartesian coordinates of the current frame.
      real(rk), allocatable :: cart(:, :)
      !> Scratch space for |rho|^2 per kept mode.
      real(rk), allocatable :: intensity(:)
   contains
      procedure :: method_setup => direct_setup
      procedure :: accumulate_frame => direct_accumulate
   end type direct_structure_factor_t

contains

   ! ---------------------------------------------------------------------
   ! Base class
   ! ---------------------------------------------------------------------

   !> Build the reciprocal grid and the accumulators from the first frame.
   subroutine sf_configure(self, frame, scheme, ierr, message)
      class(structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(ik), allocatable :: type_counts(:)
      real(rk) :: q(3), ql, dmode, length_a
      integer :: i, p1, p2, p3, h(3), n_keep, t, maxtype, pass
      integer(lk) :: g

      ierr = 0
      message = ''
      self%natoms = frame%natoms
      maxtype = 1
      do i = 1, frame%natoms
         maxtype = max(maxtype, int(frame%type_id(i)))
      end do
      self%ntypes = max(maxtype, scheme%mapped_types())

      ! Count atoms of each type id in the reference frame.
      allocate (type_counts(max(self%ntypes, 1)))
      type_counts = 0
      do i = 1, frame%natoms
         t = int(frame%type_id(i))
         if (t < 1 .or. t > size(type_counts)) then
            ierr = 1
            write (message, '(a,i0,a)') 'atom type id ', t, ' is outside the range given by -m'
            return
         end if
         type_counts(t) = type_counts(t) + 1
      end do
      if (sum(type_counts) /= frame%natoms) then
         ierr = 1
         message = 'could not classify all atoms by type'
         return
      end if
      self%type_counts = type_counts

      ! Grid extent: enough reciprocal lattice points to reach qmax.
      do i = 1, 3
         length_a = sqrt(sum(frame%cell%a(:, i)**2))
         if (length_a <= 0.0_rk) then
            ierr = 1
            message = 'degenerate simulation cell'
            return
         end if
         self%modes(i) = 2*ceiling(self%qmax*length_a/two_pi) + 1
         self%modes(i) = max(self%modes(i), 3)
      end do
      self%gridpoints = int(self%modes(1), lk)*int(self%modes(2), lk)*int(self%modes(3), lk)
      if (self%gridpoints > max_grid_modes) then
         ierr = 1
         write (message, '(a,i0,a,i0,a)') 'reciprocal grid of ', self%gridpoints, &
            ' points exceeds the limit of ', max_grid_modes, '; reduce --qmax'
         return
      end if

      ! Enumerate the reciprocal grid twice: the first pass only counts the
      ! modes inside the q range, the second fills the arrays sized for them.
      ! This avoids full-grid temporary arrays (which would double the peak
      ! memory of large grids).
      n_keep = 0
      do pass = 1, 2
         if (pass == 2) then
            self%nmodes = n_keep
            allocate (self%hkl(3, n_keep), self%qvec(3, n_keep), self%qlen(n_keep), &
                      self%gidx(n_keep))
            n_keep = 0
         end if
         do p3 = 1, self%modes(3)
            h(3) = p3 - (self%modes(3) + 1)/2
            do p2 = 1, self%modes(2)
               h(2) = p2 - (self%modes(2) + 1)/2
               do p1 = 1, self%modes(1)
                  h(1) = p1 - (self%modes(1) + 1)/2
                  if (all(h == 0)) cycle
                  q = frame%cell%b(:, 1)*real(h(1), rk) + frame%cell%b(:, 2)*real(h(2), rk) &
                      + frame%cell%b(:, 3)*real(h(3), rk)
                  ql = sqrt(sum(q*q))
                  if (ql > self%qmax .or. ql < self%qmin) cycle
                  n_keep = n_keep + 1
                  if (pass == 2) then
                     self%hkl(:, n_keep) = h
                     self%qvec(:, n_keep) = q
                     self%qlen(n_keep) = ql
                     self%gidx(n_keep) = int(p1, lk) + int(self%modes(1), lk) &
                        *(int(p2, lk) - 1 + int(self%modes(2), lk)*(int(p3, lk) - 1))
                  end if
               end do
            end do
         end do
      end do
      if (n_keep == 0) then
         ierr = 1
         message = 'no reciprocal lattice points fall inside the requested q range'
         return
      end if

      ! Shell assignment and per-mode normalization.
      self%shell_dq = (self%qmax - self%qmin)/real(self%nq, rk)
      if (self%shell_dq <= 0.0_rk) then
         ierr = 1
         message = 'q range must be positive'
         return
      end if
      allocate (self%shell(n_keep))
      do i = 1, n_keep
         self%shell(i) = int((self%qlen(i) - self%qmin)/self%shell_dq) + 1
         self%shell(i) = min(max(self%shell(i), 1), self%nq)
      end do
      allocate (self%num(self%nq), self%den(self%nq), self%shell_count(self%nq))
      allocate (self%mode_values(self%nmodes))
      if (self%partials) then
         allocate (self%partial_num(self%ntypes, self%ntypes, self%nq))
         self%partial_num = 0.0_rk
         allocate (self%pair_label(self%ntypes, self%ntypes))
         self%pair_label = ' '
         do t = 1, self%ntypes
            do i = t, self%ntypes
               self%pair_label(t, i) = scheme%pair_label(t, i)
            end do
         end do
      end if
      self%num = 0.0_rk
      self%den = 0.0_rk
      self%shell_count = 0
      self%mode_values = 0.0_rk
      if (self%want_grid) then
         allocate (self%gnum(self%gridpoints), self%gden(self%gridpoints))
         self%gnum = 0.0_rk
         self%gden = 0.0_rk
      end if
      do i = 1, n_keep
         dmode = mode_denominator(self%norm, scheme, self%natoms, type_counts, self%qlen(i))
         self%den(self%shell(i)) = self%den(self%shell(i)) + dmode
         self%shell_count(self%shell(i)) = self%shell_count(self%shell(i)) + 1
         if (self%want_grid) self%gden(self%gidx(i)) = dmode
      end do

      call self%method_setup(frame, scheme, ierr, message)
   end subroutine sf_configure

   !> Per-mode normalization denominator.
   pure real(rk) function mode_denominator(norm, scheme, natoms, type_counts, q) result(d)
      integer, intent(in) :: norm
      type(weight_scheme_t), intent(in) :: scheme
      integer(ik), intent(in) :: natoms
      integer(ik), intent(in) :: type_counts(:)
      real(rk), intent(in) :: q
      real(rk) :: w, sum_w, sum_w2
      integer :: t

      if (norm == norm_natom) then
         d = real(natoms, rk)
         return
      end if
      sum_w = 0.0_rk
      sum_w2 = 0.0_rk
      do t = 1, size(type_counts)
         if (type_counts(t) <= 0) cycle
         w = scheme%amplitude(int(t, ik), q)
         sum_w = sum_w + real(type_counts(t), rk)*w
         sum_w2 = sum_w2 + real(type_counts(t), rk)*w*w
      end do
      if (norm == norm_mean) then
         d = sum_w*sum_w/real(natoms, rk)
      else
         d = sum_w2
      end if
   end function mode_denominator

   !> Add |rho|^2 of the current frame to the shell (and grid) accumulators.
   !!
   !! `amps` is indexed like the flat reciprocal grid, i.e. amps(gidx(i)) is the
   !! amplitude of the i-th kept mode.
   subroutine sf_accumulate_modes(self, amps)
      class(structure_factor_t), intent(inout) :: self
      complex(c_double_complex), intent(in) :: amps(:)
      integer :: im
      integer(lk) :: g

      !$omp parallel do schedule(static)
      do im = 1, int(self%nmodes)
         g = self%gidx(im)
         self%mode_values(im) = real(amps(g)*conjg(amps(g)), rk)
      end do
      !$omp end parallel do
      call self%accumulate_values(self%mode_values)
   end subroutine sf_accumulate_modes

   !> Accumulate the unweighted partial sums Re(rho_a rho_b*) of one frame.
   !!
   !! `rho_species(g, a)` holds the amplitude of species `a` on the flat
   !! reciprocal grid (the same layout as `accumulate_modes` expects).  The
   !! partials are stored unweighted so that the same numbers serve for every
   !! weighting scheme and normalization.
   subroutine sf_accumulate_partials(self, rho_species)
      class(structure_factor_t), intent(inout) :: self
      complex(c_double_complex), intent(in) :: rho_species(:, :)
      integer :: ia, ib, im, ia_max
      integer(ik), allocatable :: local_partial(:, :, :)
      integer(lk) :: g, s
      real(rk) :: value, contribution

      if (.not. self%partials) return
      ia_max = size(rho_species, 2)
      allocate (local_partial(self%ntypes, self%ntypes, self%nq))
      local_partial = 0.0_rk
      do ia = 1, ia_max
         do ib = ia, ia_max
            do im = 1, int(self%nmodes)
               g = self%gidx(im)
               value = real(rho_species(g, ia)*conjg(rho_species(g, ib)), rk)
               s = self%shell(im)
               local_partial(ia, ib, s) = local_partial(ia, ib, s) + value
            end do
         end do
      end do
      self%partial_num = self%partial_num + local_partial
      deallocate (local_partial)
   end subroutine sf_accumulate_partials

   !> Add |rho|^2 values of the kept modes to the shell and grid accumulators.
   !!
   !! The shell reduction runs in parallel with one accumulator per thread; the
   !! optional per-grid-point accumulation is a plain serial pass (it is a small
   !! fraction of the frame cost and avoids grid-sized thread copies).
   subroutine sf_accumulate_values(self, values)
      class(structure_factor_t), intent(inout) :: self
      real(rk), intent(in) :: values(:)
      real(rk), allocatable :: local_num(:, :)
      integer :: im, sh, tid, nthreads
      integer(lk) :: g

      nthreads = 1
      !$ nthreads = omp_get_max_threads()
      allocate (local_num(self%nq, nthreads))
      local_num = 0.0_rk

      !$omp parallel private(im, sh, tid)
      tid = 1
      !$ tid = omp_get_thread_num() + 1
      !$omp do schedule(static)
      do im = 1, int(self%nmodes)
         sh = self%shell(im)
         local_num(sh, tid) = local_num(sh, tid) + values(im)
      end do
      !$omp end do
      !$omp end parallel

      do tid = 1, nthreads
         do sh = 1, self%nq
            self%num(sh) = self%num(sh) + local_num(sh, tid)
         end do
      end do
      deallocate (local_num)

      if (self%want_grid) then
         do im = 1, int(self%nmodes)
            g = self%gidx(im)
            self%gnum(g) = self%gnum(g) + values(im)
         end do
      end if
      self%nframes = self%nframes + 1
   end subroutine sf_accumulate_values

   !> Hook for methods that fill the output arrays rather than the per-frame
   !! accumulators (the Debye method computes S(q) from histograms at the end).
   subroutine sf_prepare_output(self, scheme, ierr, message)
      class(structure_factor_t), intent(inout) :: self
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      ierr = 0
      message = ''
   end subroutine sf_prepare_output

   !> Normalized S(q) of one shell.
   pure real(rk) function sf_shell_value(self, shell) result(value)
      class(structure_factor_t), intent(in) :: self
      integer, intent(in) :: shell
      real(rk) :: frames

      frames = real(max(self%nframes, 1_lk), rk)
      value = 0.0_rk
      if (self%den(shell) > 0.0_rk) value = self%num(shell)/(frames*self%den(shell))
   end function sf_shell_value

   !> Partial structure factors of one shell, in the OVITO convention
   !!
   !!   S_ab(q) = (1/N) < sum_{j in a} sum_{k in b} sin(q r_jk)/(q r_jk) >
   !!
   !! so that S_aa -> x_a and S_ab -> 0 (a /= b) at large q and
   !! `S(q) = sum_ab (2 - delta_ab) S_ab(q)` for unit weights.  With
   !! `faber_ziman` the Faber-Ziman form is returned instead,
   !!
   !!   A_ab(q) = (S_ab(q) - x_a delta_ab)/(x_a x_b) + 1,
   !!
   !! which tends to 1 for every pair.  Values are given for a <= b.
   pure real(rk) function sf_partial_value(self, ia, ib, shell) result(value)
      class(structure_factor_t), intent(in) :: self
      integer, intent(in) :: ia, ib, shell
      real(rk) :: frames, xa, xb

      frames = real(max(self%nframes, 1_lk), rk)
      value = 0.0_rk
      if (.not. allocated(self%partial_num)) return
      ! The grid method sums over the modes of the shell, so the mode count must
      ! be divided out (the total gets that from `den`).  Methods that evaluate
      ! S(q) directly (Debye) leave shell_count at zero.
      value = self%partial_num(ia, ib, shell) &
              /(frames*real(self%natoms, rk)*real(max(sf_shell_modes(self, shell), 1), rk))
      if (.not. self%faber_ziman) return
      xa = real(sf_count_type_of(self, ia), rk)/real(self%natoms, rk)
      xb = real(sf_count_type_of(self, ib), rk)/real(self%natoms, rk)
      if (xa*xb <= 0.0_rk) then
         value = 0.0_rk
         return
      end if
      value = (value - merge(xa, 0.0_rk, ia == ib))/(xa*xb) + 1.0_rk
   end function sf_partial_value

   !> Number of reciprocal lattice modes that contributed to a shell (0 when the
   !! method does not use the reciprocal grid).
   pure integer function sf_shell_modes(self, shell) result(n)
      class(structure_factor_t), intent(in) :: self
      integer, intent(in) :: shell
      n = 0
      if (allocated(self%shell_count)) then
         if (shell >= 1 .and. shell <= size(self%shell_count)) n = int(self%shell_count(shell))
      end if
   end function sf_shell_modes

   !> Number of atoms of a type id (helper for the partial normalization).
   pure integer function sf_count_type_of(self, type_id) result(n)
      class(structure_factor_t), intent(in) :: self
      integer, intent(in) :: type_id
      n = 0
      if (allocated(self%type_counts)) then
         if (type_id >= 1 .and. type_id <= size(self%type_counts)) n = int(self%type_counts(type_id))
      end if
   end function sf_count_type_of

   !> Column label of a partial, "S(Si-O)" or "A(Si-O)" for Faber-Ziman.
   pure function sf_partial_column(self, ia, ib) result(label)
      class(structure_factor_t), intent(in) :: self
      integer, intent(in) :: ia, ib
      character(len=24) :: label
      character(len=18) :: pair
      character(len=2) :: prefix

      pair = ' '
      if (allocated(self%pair_label)) pair = self%pair_label(ia, ib)
      prefix = 'S('
      if (self%faber_ziman) prefix = 'A('
      label = prefix//trim(pair)//')'
   end function sf_partial_column

   !> Normalized S(q) of one kept reciprocal lattice mode.
   pure real(rk) function sf_grid_value(self, mode) result(value)
      class(structure_factor_t), intent(in) :: self
      integer(lk), intent(in) :: mode
      real(rk) :: frames
      integer(lk) :: g

      frames = real(max(self%nframes, 1_lk), rk)
      value = 0.0_rk
      g = self%gidx(mode)
      if (self%gden(g) > 0.0_rk) value = self%gnum(g)/(frames*self%gden(g))
   end function sf_grid_value

   !> Write the 1D shell table and (optionally) the reciprocal grid table.
   subroutine sf_write_results(self, shell_unit, grid_unit, ierr, message)
      class(structure_factor_t), intent(in) :: self
      integer, intent(in) :: shell_unit, grid_unit
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: qc, value, frames
      integer :: s, t, u
      integer(lk) :: i, g

      ierr = 0
      message = ''
      ! num(:) accumulates every frame while den(:) is the single frame
      ! normalization, so the average over frames is num/(nframes*den).
      frames = real(max(self%nframes, 1_lk), rk)
      if (shell_unit /= no_unit) then
         if (self%partials .and. allocated(self%partial_num)) then
            write (shell_unit, '(a)', advance='no') '# q S(q)'
            do t = 1, self%ntypes
               do u = t, self%ntypes
                  if (sf_count_type_of(self, t) == 0) cycle
                  if (sf_count_type_of(self, u) == 0) cycle
                  write (shell_unit, '(a)', advance='no') ' '//self%partial_column(t, u)
               end do
            end do
            write (shell_unit, '(a)') ''
         else
            write (shell_unit, '(a)') '# q S(q)'
         end if
         do s = 1, self%nq
            qc = self%qmin + (real(s, rk) - 0.5_rk)*self%shell_dq
            value = self%shell_value(s)
            if (self%partials .and. allocated(self%partial_num)) then
               write (shell_unit, '(f14.6,2x,es20.12)', advance='no') qc, value
               do t = 1, self%ntypes
                  do u = t, self%ntypes
                     if (sf_count_type_of(self, t) == 0) cycle
                     if (sf_count_type_of(self, u) == 0) cycle
                     write (shell_unit, '(2x,es20.12)', advance='no') &
                        self%partial_value(t, u, s)
                  end do
               end do
               write (shell_unit, '(a)') ''
            else
               write (shell_unit, '(f14.6,2x,es20.12)') qc, value
            end if
         end do
      end if

      if (grid_unit /= no_unit .and. self%want_grid) then
         write (grid_unit, '(a)') '# qx qy qz S(q)'
         do i = 1, self%nmodes
            value = self%grid_value(i)
            write (grid_unit, '(3(f14.8,2x),es20.12)') self%qvec(1, i), self%qvec(2, i), &
               self%qvec(3, i), value
         end do
      end if
   end subroutine sf_write_results

   ! ---------------------------------------------------------------------
   ! NUFFT method
   ! ---------------------------------------------------------------------

   subroutine nufft_setup(self, frame, scheme, ierr, message)
      class(nufft_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(c_int64_t) :: n_modes(3)
      integer(c_int) :: ier

      ierr = 0
      message = ''

      call sf_prepare_species(self, frame, scheme, self%species_of, self%species_type, &
                              self%amp_const, self%amp_table, self%q_dependent, ierr, message)
      if (ierr /= 0) return

      ! Allocate work arrays and build the FINUFFT plan.
      allocate (self%xa(frame%natoms), self%ya(frame%natoms), self%za(frame%natoms))
      allocate (self%strengths(frame%natoms))
      allocate (self%fk(self%gridpoints))
      self%strengths = (0.0_rk, 0.0_rk)
      self%fk = (0.0_rk, 0.0_rk)
      ! The species sum is only needed when the amplitudes depend on q.
      if (self%q_dependent) then
         allocate (self%total(self%gridpoints))
         self%total = (0.0_rk, 0.0_rk)
      end if
      if (self%partials) then
         allocate (self%rho_species(self%gridpoints, size(self%species_type)))
         if (.not. allocated(self%total)) then
            allocate (self%total(self%gridpoints))
            self%total = (0.0_rk, 0.0_rk)
         end if
      end if

      call finufft_default_opts(self%opts)
      self%opts%modeord = mode_cmcl
      self%opts%nthreads = int(max(self%nthreads, 0), c_int)
      n_modes = int(self%modes, c_int64_t)
      ier = finufft_makeplan(type1, 3_c_int, n_modes, 1_c_int, 1_c_int, &
                             real(self%eps, c_double), self%plan, self%opts)
      if (ier /= 0) then
         ierr = 1
         write (message, '(a,i0)') 'FINUFFT failed to create a plan (ier = ', ier
         message = trim(message)//')'
         return
      end if
   end subroutine nufft_setup

   !> Species present in the frame plus the amplitude of each species at every
   !! kept reciprocal lattice mode.  Shared by the CPU and GPU NUFFT methods.
   subroutine sf_prepare_species(self, frame, scheme, species_of, species_type, &
                                 amp_const, amp_table, q_dependent, ierr, message)
      class(structure_factor_t), intent(in) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer(ik), allocatable, intent(out) :: species_of(:)
      integer(ik), allocatable, intent(out) :: species_type(:)
      real(rk), allocatable, intent(out) :: amp_const(:)
      real(rk), allocatable, intent(out) :: amp_table(:, :)
      logical, intent(out) :: q_dependent
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      logical, allocatable :: present(:)
      integer :: i, t, isp, nspecies

      ierr = 0
      message = ''

      ! Species actually present in the trajectory.
      allocate (present(max(self%ntypes, 1)))
      present = .false.
      do i = 1, frame%natoms
         present(int(frame%type_id(i))) = .true.
      end do
      nspecies = count(present)
      allocate (species_type(nspecies))
      isp = 0
      do t = 1, size(present)
         if (.not. present(t)) cycle
         isp = isp + 1
         species_type(isp) = t
      end do
      allocate (species_of(frame%natoms))
      do i = 1, frame%natoms
         species_of(i) = findloc(species_type, frame%type_id(i), dim=1)
      end do

      ! Amplitudes: constant per species unless the weights depend on q.
      q_dependent = scheme%kind == weight_xray
      if (q_dependent) then
         allocate (amp_table(self%nmodes, nspecies))
         do isp = 1, nspecies
            do i = 1, int(self%nmodes)
               amp_table(i, isp) = scheme%amplitude(species_type(isp), self%qlen(i))
            end do
         end do
         allocate (amp_const(0))
      else
         allocate (amp_const(nspecies))
         do isp = 1, nspecies
            amp_const(isp) = scheme%amplitude(species_type(isp), 0.0_rk)
         end do
         allocate (amp_table(0, 0))
      end if
   end subroutine sf_prepare_species

   subroutine nufft_accumulate(self, frame, scheme, ierr, message)
      class(nufft_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: s(3), amp, value
      real(rk) :: dummy(1), pi
      integer(c_int) :: ier
      integer :: i, im, isp, sh
      integer(lk) :: g

      ierr = 0
      message = ''
      dummy = 0.0_rk
      pi = acos(-1.0_rk)

      ! Scaled coordinates in [-pi, pi) for the type-1 transform.
      do i = 1, frame%natoms
         s = frame%cell%fractional(frame%pos(:, i))
         s = s - floor(s)
         self%xa(i) = two_pi*s(1) - pi
         self%ya(i) = two_pi*s(2) - pi
         self%za(i) = two_pi*s(3) - pi
      end do
      ier = finufft_setpts(self%plan, int(frame%natoms, c_int64_t), self%xa, self%ya, &
                           self%za, 0_c_int64_t, dummy, dummy, dummy)
      if (ier /= 0) then
         ierr = 1
         write (message, '(a,i0)') 'FINUFFT setpts failed (ier = ', ier
         return
      end if

      if (self%q_dependent) then
         ! Partial structure factors need the individual species amplitudes,
         ! which also gives the weighted total in the same pass.
         if (self%partials) then
            do isp = 1, size(self%species_type)
               do i = 1, frame%natoms
                  if (self%species_of(i) == isp) then
                     self%strengths(i) = (1.0_rk, 0.0_rk)
                  else
                     self%strengths(i) = (0.0_rk, 0.0_rk)
                  end if
               end do
               ier = finufft_execute(self%plan, self%strengths, self%fk)
               if (ier /= 0) then
                  ierr = 1
                  write (message, '(a,i0)') 'FINUFFT execute failed (ier = ', ier
                  return
               end if
               self%rho_species(:, isp) = self%fk
            end do
            call self%accumulate_partials(self%rho_species)
            !$omp parallel do schedule(static) private(im, isp, g)
            do im = 1, int(self%nmodes)
               g = self%gidx(im)
               self%total(g) = (0.0_rk, 0.0_rk)
               do isp = 1, size(self%species_type)
                  self%total(g) = self%total(g) &
                                  + cmplx(self%amp_table(im, isp), 0.0_rk, c_double_complex) &
                                    *self%rho_species(g, isp)
               end do
            end do
            !$omp end parallel do
            call self%accumulate_modes(self%total)
            return
         end if
         ! X-ray form factors depend on q, so each species needs its own
         ! transform before they are combined with f_alpha(|q|).
         self%total = (0.0_rk, 0.0_rk)
         do isp = 1, size(self%species_type)
            do i = 1, frame%natoms
               if (self%species_of(i) == isp) then
                  self%strengths(i) = (1.0_rk, 0.0_rk)
               else
                  self%strengths(i) = (0.0_rk, 0.0_rk)
               end if
            end do
            ier = finufft_execute(self%plan, self%strengths, self%fk)
            if (ier /= 0) then
               ierr = 1
               write (message, '(a,i0)') 'FINUFFT execute failed (ier = ', ier
               return
            end if
            !$omp parallel do schedule(static)
            do im = 1, int(self%nmodes)
               g = self%gidx(im)
               self%total(g) = self%total(g) &
                               + cmplx(self%amp_table(im, isp), 0.0_rk, c_double_complex) &
                                 *self%fk(g)
            end do
            !$omp end parallel do
         end do
         call self%accumulate_modes(self%total)
      else
         ! Partial structure factors: compute the species amplitudes first,
         ! then the partials and the weighted total from them.
         if (self%partials) then
            do isp = 1, size(self%species_type)
               do i = 1, frame%natoms
                  if (self%species_of(i) == isp) then
                     self%strengths(i) = (1.0_rk, 0.0_rk)
                  else
                     self%strengths(i) = (0.0_rk, 0.0_rk)
                  end if
               end do
               ier = finufft_execute(self%plan, self%strengths, self%fk)
               if (ier /= 0) then
                  ierr = 1
                  write (message, '(a,i0)') 'FINUFFT execute failed (ier = ', ier
                  return
               end if
               self%rho_species(:, isp) = self%fk
            end do
            call self%accumulate_partials(self%rho_species)
            !$omp parallel do schedule(static) private(im, isp, g)
            do im = 1, int(self%nmodes)
               g = self%gidx(im)
               self%total(g) = (0.0_rk, 0.0_rk)
               do isp = 1, size(self%species_type)
                  self%total(g) = self%total(g) &
                                  + cmplx(self%amp_const(isp), 0.0_rk, c_double_complex) &
                                    *self%rho_species(g, isp)
               end do
            end do
            !$omp end parallel do
            call self%accumulate_modes(self%total)
            return
         end if
         ! Unit and neutron weights do not depend on q, so the per-atom
         ! amplitudes go straight into a single transform:
         !   T(q) = sum_alpha f_alpha rho_alpha(q) = sum_j w_j exp(i q.r_j)
         do i = 1, frame%natoms
            self%strengths(i) = cmplx(self%amp_const(self%species_of(i)), 0.0_rk, &
                                      c_double_complex)
         end do
         ier = finufft_execute(self%plan, self%strengths, self%fk)
         if (ier /= 0) then
            ierr = 1
            write (message, '(a,i0)') 'FINUFFT execute failed (ier = ', ier
            return
         end if
         call self%accumulate_modes(self%fk)
      end if
   end subroutine nufft_accumulate

   subroutine nufft_finalize(self)
      type(nufft_structure_factor_t), intent(inout) :: self
      integer(c_int) :: ier
      if (c_associated(self%plan)) then
         ier = finufft_destroy(self%plan)
         self%plan = c_null_ptr
      end if
   end subroutine nufft_finalize

   ! ---------------------------------------------------------------------
   ! Direct method (reference / fallback)
   ! ---------------------------------------------------------------------

   subroutine direct_setup(self, frame, scheme, ierr, message)
      class(direct_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer :: i, t
      integer(lk) :: work

      ierr = 0
      message = ''
      allocate (self%type_of(frame%natoms))
      self%type_of = frame%type_id
      allocate (self%species_first(self%ntypes + 1))
      self%species_first = 1
      do i = 1, frame%natoms
         self%species_first(int(frame%type_id(i)) + 1) = &
            self%species_first(int(frame%type_id(i)) + 1) + 1
      end do
      do t = 1, self%ntypes
         self%species_first(t + 1) = self%species_first(t + 1) + self%species_first(t) - 1
      end do
      allocate (self%atom_of(frame%natoms))
      block
         integer(ik), allocatable :: cursor(:)
         allocate (cursor(self%ntypes))
         cursor = self%species_first(1:self%ntypes)
         do i = 1, frame%natoms
            t = int(frame%type_id(i))
            self%atom_of(cursor(t)) = int(i, ik)
            cursor(t) = cursor(t) + 1
         end do
      end block
      allocate (self%cart(3, frame%natoms))
      allocate (self%intensity(self%nmodes))
      work = int(frame%natoms, lk)*self%nmodes
      if (work > 2000000000_lk) then
         ierr = 1
         write (message, '(a,i0,a)') 'the direct method would need ', work, &
            ' phase evaluations; use the default NUFFT method or a smaller --qmax'
         return
      end if
   end subroutine direct_setup

   subroutine direct_accumulate(self, frame, scheme, ierr, message)
      class(direct_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: s(3), amp, phase, value
      real(rk) :: qx, qy, qz
      complex(c_double_complex) :: acc, total
      integer :: i, im, isp, sh, idx
      integer(lk) :: g

      ierr = 0
      message = ''
      do i = 1, frame%natoms
         s = frame%cell%fractional(frame%pos(:, i))
         s = s - floor(s)
         self%cart(:, i) = matmul(frame%cell%a, s)
      end do

      !$omp parallel do schedule(static) private(im, isp, idx, i, amp, phase, acc, total, qx, qy, qz)
      do im = 1, int(self%nmodes)
         qx = self%qvec(1, im)
         qy = self%qvec(2, im)
         qz = self%qvec(3, im)
         total = (0.0_rk, 0.0_rk)
         do isp = 1, self%ntypes
            amp = scheme%amplitude(int(isp, ik), self%qlen(im))
            acc = (0.0_rk, 0.0_rk)
            do idx = self%species_first(isp), self%species_first(isp + 1) - 1
               i = self%atom_of(idx)
               phase = qx*self%cart(1, i) + qy*self%cart(2, i) + qz*self%cart(3, i)
               acc = acc + cmplx(cos(phase), sin(phase), c_double_complex)
            end do
            total = total + cmplx(amp, 0.0_rk, c_double_complex)*acc
         end do
         self%intensity(im) = real(total*conjg(total), rk)
      end do
      !$omp end parallel do

      call self%accumulate_values(self%intensity)
   end subroutine direct_accumulate

end module sqc_structure_factor
