! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Debye scattering equation method (real space pair histograms).
!!
!! Instead of transforming the periodic density (the NUFFT method), this method
!! counts atom pairs as a function of their distance for every pair of chemical
!! types and averages those histograms over the trajectory.  S(q) is then
!! evaluated directly from the histograms:
!!
!!   S(q) = [ sum_ab f_a(q) f_b(q) sum_k n_ab[k] sin(q r_k)/(q r_k)
!!            + sum_a N_a f_a(q)^2 ] / W(q)
!!
!! where `n_ab[k]` are the ordered pair counts, `r_k = (k-1/2) dr` the bin
!! centers, the explicit sum is the r = 0 self term (an atom with itself, which
!! no histogram can hold, so it is added analytically) and `W(q)` is the same
!! normalization the grid method uses (`--norm`).  The pair counts also give the
!! total and partial pair distribution functions g(r).
!!
!! The same recipe is used by the `debyer` program (pair histograms plus the
!! Debye sum); see B. E. Warren, "X-ray Diffraction".
module sqc_debye
   use sqc_kinds
   use sqc_dump, only: frame_t
   use sqc_weights, only: weight_scheme_t
   use sqc_neighbour_list, only: neighbour_list_t
   use sqc_structure_factor, only: structure_factor_t, mode_denominator, norm_natom, &
                            norm_mean, norm_self
   use omp_lib, only: omp_get_max_threads, omp_get_thread_num
#ifdef SQC_HAVE_HDF5
   use sqc_hdf5, only: hdf5_write_rdf
#endif
   implicit none
   private

   public :: debye_structure_factor_t, debye_default_dr, debye_default_skin

   !> Default r bin width [A] (reliable up to q ~ pi/(2 dr) ~ 150 1/A).
   real(rk), parameter :: debye_default_dr = 0.01_rk
   !> Default Verlet skin [A].
   real(rk), parameter :: debye_default_skin = 1.0_rk

   !> Debye method state: partial pair histograms plus the neighbour list.
   type, extends(structure_factor_t) :: debye_structure_factor_t
      !> Pair cutoff [A]; 0 means "derive from the cell" during setup.
      real(rk) :: rmax = 0.0_rk
      !> Radial bin width [A] and number of bins.
      real(rk) :: dr = debye_default_dr
      integer :: nbins = 0
      !> Verlet skin [A] (0 rebuilds the pair list every frame).
      real(rk) :: skin = debye_default_skin
      !> Periodicity of the box (read from the dump header).
      logical :: pbc(3) = .true.
      !> Optional pair distribution function file ('' = none).
      character(len=:), allocatable :: rdf_path
      !> Apply the cut-off density correction (debyer's add_cutoff_correction).
      logical :: correct_cutoff = .true.
      !> Type id of each atom and number of atoms per type id.
      integer(ik), allocatable :: type_of(:)
      integer(ik), allocatable :: count_type(:)
      !> Ordered pair counts histogram(a, b, bin), averaged over frames by nhist.
      real(rk), allocatable :: histogram(:, :, :)
      !> Box volume, used for the g(r) normalization.
      real(rk) :: volume = 0.0_rk
      !> Cell list used to enumerate the pairs.
      type(neighbour_list_t) :: nlist
      !> Number of frames whose pairs were histogrammed.
      integer(lk) :: nhist = 0
   contains
      procedure :: method_setup => debye_setup
      procedure :: accumulate_frame => debye_accumulate
      procedure :: prepare_output => debye_prepare_output
      procedure :: write_rdf => debye_write_rdf
      final :: debye_finalize
   end type debye_structure_factor_t

contains

   !> Validate the settings, derive defaults from the cell, allocate the state.
   subroutine debye_setup(self, frame, scheme, ierr, message)
      class(debye_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: widths(3), bound, span(3)
      integer :: i, t

      ierr = 0
      message = ''
      self%natoms = frame%natoms
      self%ntypes = max(scheme%mapped_types(), maxval(frame%type_id))
      self%pbc = frame%pbc
      self%volume = frame%cell%volume

      allocate (self%type_of(frame%natoms), self%count_type(self%ntypes))
      self%type_of = frame%type_id
      self%count_type = 0
      do i = 1, frame%natoms
         t = int(frame%type_id(i))
         if (t < 1 .or. t > self%ntypes) then
            ierr = 1
            write (message, '(a,i0,a)') 'atom type id ', t, ' is outside the range given by -m'
            return
         end if
         self%count_type(t) = self%count_type(t) + 1
      end do

      widths = frame%cell%widths()
      if (any(widths <= 0.0_rk)) then
         ierr = 1
         message = 'degenerate simulation cell'
         return
      end if
      bound = huge(1.0_rk)
      do i = 1, 3
         if (self%pbc(i)) bound = min(bound, widths(i))
      end do

      if (self%rmax <= 0.0_rk) then
         ! Default: half of the smallest periodic side (minimum image), or the
         ! whole configuration when the box is not periodic at all.
         if (bound < huge(1.0_rk)) then
            self%rmax = 0.5_rk*bound
         else
            do i = 1, 3
               span(i) = maxval(frame%pos(i, :)) - minval(frame%pos(i, :))
            end do
            self%rmax = sqrt(sum(span*span))
            if (self%rmax <= 0.0_rk) self%rmax = sqrt(sum(matmul(frame%cell%a, &
               [1.0_rk, 1.0_rk, 1.0_rk])**2))
         end if
      else if (bound < huge(1.0_rk)) then
         if (self%rmax >= bound) then
            ierr = 1
            write (message, '(a,f0.3,a,f0.3,a)') 'the pair cutoff (', self%rmax, &
               ' A) must stay below the smallest periodic box side (', bound, ' A)'
            return
         end if
      end if
      if (self%dr <= 0.0_rk) then
         ierr = 1
         message = 'the radial bin width must be positive'
         return
      end if
      self%nbins = max(1, ceiling(self%rmax/self%dr))
      if (real(self%nbins, rk)*real(self%ntypes, rk)**2 > 2.0e8_rk) then
         ierr = 1
         write (message, '(a,i0,a,i0,a)') 'the pair histogram would need ', self%nbins, &
            ' bins x ', self%ntypes, '^2 type pairs; increase --dr or reduce --rmax'
         return
      end if

      allocate (self%histogram(self%ntypes, self%ntypes, self%nbins))
      self%histogram = 0.0_rk
      self%nhist = 0
      ! This method evaluates S(q) directly at the shell centers, so the
      ! reciprocal mode list of the base type is not used: clear its per-shell
      ! mode counts, otherwise the partial normalization would divide by them.
      if (allocated(self%shell_count)) self%shell_count = 0
      self%nmodes = 0
      call self%nlist%configure(frame%cell, frame%natoms, self%rmax, self%skin, self%pbc, &
                                ierr, message)
   end subroutine debye_setup

   !> Count the ordered atom pairs of one frame into the histograms.
   subroutine debye_accumulate(self, frame, scheme, ierr, message)
      class(debye_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: local_hist(:, :, :, :)
      integer :: i, j, ia, ib, nthreads, tid, bin
      integer(lk) :: k, first, last
      real(rk) :: dr(3), d2, inv_dr

      ierr = 0
      message = ''
      do i = 1, frame%natoms
         if (self%type_of(i) /= frame%type_id(i)) then
            ierr = 1
            message = 'the atom types changed between frames'
            return
         end if
      end do
      call self%nlist%begin_frame(frame%pos, ierr, message)
      if (ierr /= 0) return

      nthreads = 1
      !$ nthreads = omp_get_max_threads()
      allocate (local_hist(self%ntypes, self%ntypes, self%nbins, nthreads))
      local_hist = 0.0_rk
      inv_dr = 1.0_rk/self%dr

      !$omp parallel private(tid, i, k, first, last, j, ia, ib, dr, d2, bin)
      tid = 1
      !$ tid = omp_get_thread_num() + 1
      !$omp do schedule(static)
      do i = 1, frame%natoms
         ia = int(self%type_of(i))
         first = self%nlist%pair_start(i)
         last = self%nlist%pair_start(i + 1) - 1
         do k = first, last
            j = int(self%nlist%pair_atom(k))
            dr = self%nlist%pos(:, j) - self%nlist%pos(:, i) &
                 + matmul(self%nlist%bins, real(self%nlist%pair_shift(:, k), rk))
            d2 = sum(dr*dr)
            if (d2 > self%rmax*self%rmax) cycle
            bin = int(sqrt(d2)*inv_dr) + 1
            if (bin < 1 .or. bin > self%nbins) cycle
            ib = int(self%type_of(j))
            local_hist(ia, ib, bin, tid) = local_hist(ia, ib, bin, tid) + 1.0_rk
         end do
      end do
      !$omp end do
      !$omp end parallel

      self%histogram = self%histogram + sum(local_hist, dim=4)
      deallocate (local_hist)
      self%nhist = self%nhist + 1
      self%nframes = self%nframes + 1
   end subroutine debye_accumulate

   !> Evaluate S(q) from the histograms and fill the shared output arrays.
   !!
   !! The shared writers compute `num(s)/(nframes*den(s))`, so the values are
   !! stored as `num = S(q)*W(q)*nframes` with `den = W(q)`.
   subroutine debye_prepare_output(self, scheme, ierr, message)
      class(debye_structure_factor_t), intent(inout) :: self
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: weight(:)
      real(rk) :: q, r_bin, sinc, pair_sum, self_sum, weight_sum, avg, corr, rho0, qc
      real(rk) :: part_sum, part_corr, xa, xb
      integer :: s, ia, ib, k, t

      ierr = 0
      message = ''
      if (self%nhist <= 0) then
         ierr = 1
         message = 'no frames were accumulated'
         return
      end if
      allocate (weight(self%ntypes))
      do s = 1, self%nq
         q = self%qmin + (real(s, rk) - 0.5_rk)*self%shell_dq
         do t = 1, self%ntypes
            weight(t) = scheme%amplitude(int(t, ik), q)
         end do
         pair_sum = 0.0_rk
         do ia = 1, self%ntypes
            if (self%count_type(ia) == 0) cycle
            do ib = 1, self%ntypes
               if (self%count_type(ib) == 0) cycle
               do k = 1, self%nbins
                  if (self%histogram(ia, ib, k) == 0.0_rk) cycle
                  r_bin = (real(k, rk) - 0.5_rk)*self%dr
                  sinc = sin(q*r_bin)/(q*r_bin)
                  pair_sum = pair_sum + weight(ia)*weight(ib) &
                             *(self%histogram(ia, ib, k)/real(self%nhist, rk))*sinc
               end do
            end do
         end do
         ! r = 0 self term: exactly the "self" denominator of the same weights
         self_sum = mode_denominator(norm_self, scheme, self%natoms, self%count_type, q)
         weight_sum = 0.0_rk
         do t = 1, self%ntypes
            weight_sum = weight_sum + real(self%count_type(t), rk)*weight(t)
         end do
         self%den(s) = mode_denominator(self%norm, scheme, self%natoms, self%count_type, q)
         ! Partial structure factors (unweighted), same convention as the
         ! reciprocal method: (1/N) <sum_jk sin(q r)/(q r)> with the r = 0 self
         ! term for equal species.
         if (self%partials) then
            ! The cut-off correction must be shared by the partials in
            ! proportion to x_a x_b so that they still add up to the total
            ! (sum_ab (2 - delta_ab) x_a x_b = 1).
            part_corr = 0.0_rk
            if (self%correct_cutoff .and. any(self%pbc)) then
               qc = q*self%rmax
               part_corr = 4.0_rk*acos(-1.0_rk)*(real(self%natoms, rk)/self%volume)/(q*q) &
                           *(self%rmax*cos(qc) - sin(qc)/q)
            end if
            do ia = 1, self%ntypes
               if (self%count_type(ia) == 0) cycle
               do ib = ia, self%ntypes
                  if (self%count_type(ib) == 0) cycle
                  part_sum = 0.0_rk
                  do k = 1, self%nbins
                     if (self%histogram(ia, ib, k) == 0.0_rk) cycle
                     r_bin = (real(k, rk) - 0.5_rk)*self%dr
                     part_sum = part_sum + (self%histogram(ia, ib, k) &
                                /real(self%nhist, rk))*sin(q*r_bin)/(q*r_bin)
                  end do
                  if (ia == ib) part_sum = part_sum + real(self%count_type(ia), rk)
                  xa = real(self%count_type(ia), rk)/real(self%natoms, rk)
                  xb = real(self%count_type(ib), rk)/real(self%natoms, rk)
                  part_sum = part_sum + xa*xb*part_corr*real(self%natoms, rk)
                  self%partial_num(ia, ib, s) = part_sum*real(self%nframes, rk)
               end do
            end do
         end if
         ! debyer's density correction for the pairs beyond the cut-off
         ! (add_cutoff_correction): restores the low q behaviour and removes the
         ! truncation ripple.  It is expressed in the per-atom normalization
         ! (debyer's S = pattern/N), hence the factor N when we add it to our
         ! numerator.
         corr = 0.0_rk
         if (self%correct_cutoff .and. any(self%pbc)) then
            rho0 = real(self%natoms, rk)/self%volume
            avg = weight_sum/real(self%natoms, rk)
            qc = q*self%rmax
            corr = avg**2*4.0_rk*acos(-1.0_rk)*rho0/(q*q) &
                   *(self%rmax*cos(qc) - sin(qc)/q)
         end if
         if (self%den(s) > 0.0_rk) then
            ! shell_value() returns num(s)/(nframes*den(s)) = S(q)
            self%num(s) = (pair_sum + self_sum + corr*real(self%natoms, rk)) &
                          *real(self%nframes, rk)
         else
            self%num(s) = 0.0_rk
         end if
      end do
      deallocate (weight)
   end subroutine debye_prepare_output

   !> Write the total and every partial g(r) into one file.
   subroutine debye_write_rdf(self, scheme, path, ierr, message)
      class(debye_structure_factor_t), intent(in) :: self
      type(weight_scheme_t), intent(in) :: scheme
      character(len=*), intent(in) :: path
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: g_partial(:, :, :), g_total(:), r(:), weight(:)
      character(len=18), allocatable :: labels(:)
      real(rk) :: shell, rho_b, f_mean, s_of_q
      integer :: k, ia, ib, u, t, npair, p

      ierr = 0
      message = ''
      allocate (g_partial(self%ntypes, self%ntypes, self%nbins), g_total(self%nbins), &
                r(self%nbins), weight(self%ntypes))
      npair = self%ntypes*(self%ntypes + 1)/2
      allocate (labels(npair))
      p = 0
      do ia = 1, self%ntypes
         do ib = ia, self%ntypes
            p = p + 1
            labels(p) = pair_label(scheme, ia, ib)
         end do
      end do

      do k = 1, self%nbins
         r(k) = (real(k, rk) - 0.5_rk)*self%dr
         shell = 4.0_rk*acos(-1.0_rk)*r(k)**2*self%dr
         do ia = 1, self%ntypes
            do ib = 1, self%ntypes
               g_partial(ia, ib, k) = 0.0_rk
               if (self%count_type(ia) == 0 .or. self%count_type(ib) == 0) cycle
               rho_b = real(self%count_type(ib), rk)/self%volume
               if (rho_b*shell > 0.0_rk) g_partial(ia, ib, k) = &
                  (self%histogram(ia, ib, k)/real(self%nhist, rk)) &
                  /(real(self%count_type(ia), rk)*shell*rho_b)
            end do
         end do
      end do

      ! Total g(r), scattering weighted with the q -> 0 amplitudes.
      do t = 1, self%ntypes
         weight(t) = scheme%amplitude(int(t, ik), 0.0_rk)
      end do
      f_mean = 0.0_rk
      do t = 1, self%ntypes
         f_mean = f_mean + real(self%count_type(t), rk)/real(self%natoms, rk)*weight(t)
      end do
      g_total = 0.0_rk
      if (abs(f_mean) > 0.0_rk) then
         do k = 1, self%nbins
            s_of_q = 0.0_rk
            do ia = 1, self%ntypes
               do ib = 1, self%ntypes
                  if (self%count_type(ia) == 0 .or. self%count_type(ib) == 0) cycle
                  s_of_q = s_of_q + (real(self%count_type(ia), rk)/real(self%natoms, rk)) &
                           *(real(self%count_type(ib), rk)/real(self%natoms, rk)) &
                           *weight(ia)*weight(ib)*g_partial(ia, ib, k)
               end do
            end do
            g_total(k) = s_of_q/f_mean**2
         end do
      end if

      if (index(path, '.h5') > 0 .or. index(path, '.hdf5') > 0) then
         call write_rdf_hdf5(self, labels, path, r, g_total, g_partial, ierr, message)
      else
         open (newunit=u, file=trim(path), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            message = 'cannot write the RDF file "'//trim(path)//'"'
            return
         end if
         write (u, '(a)', advance='no') '# r g(r)'
         do p = 1, npair
            write (u, '(a)', advance='no') ' g('//trim(labels(p))//')'
         end do
         write (u, '(a)') ''
         do k = 1, self%nbins
            write (u, '(f12.6,2x,es18.10)', advance='no') r(k), g_total(k)
            do ia = 1, self%ntypes
               do ib = ia, self%ntypes
                  write (u, '(2x,es18.10)', advance='no') g_partial(ia, ib, k)
               end do
            end do
            write (u, '(a)') ''
         end do
         close (u)
      end if
      deallocate (g_partial, g_total, r, weight, labels)
   end subroutine debye_write_rdf

   !> HDF5 flavour of the RDF output.
   subroutine write_rdf_hdf5(self, labels, path, r, g_total, g_partial, ierr, message)
      class(debye_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: labels(:), path
      real(rk), intent(in) :: r(:), g_total(:), g_partial(:, :, :)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
#ifdef SQC_HAVE_HDF5
      call hdf5_write_rdf(path, r, g_total, g_partial, labels, ierr, message)
#else
      ierr = 1
      message = 'this build has no HDF5 support; use a text file name for --rdf'
#endif
   end subroutine write_rdf_hdf5

   !> Label of a type pair ("Si-O" when an element mapping is available).
   pure function pair_label(scheme, ia, ib) result(label)
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(in) :: ia, ib
      character(len=18) :: label
      character(len=2) :: sa, sb

      sa = ' '
      sb = ' '
      if (ia <= scheme%mapped_types()) sa = scheme%symbols(ia)
      if (ib <= scheme%mapped_types()) sb = scheme%symbols(ib)
      if (len_trim(sa) > 0 .and. len_trim(sb) > 0) then
         label = trim(sa)//'-'//trim(sb)
      else
         write (label, '(i0,a,i0)') ia, '-', ib
      end if
   end function pair_label

   subroutine debye_finalize(self)
      type(debye_structure_factor_t), intent(inout) :: self
      call self%nlist%finalize()
      if (allocated(self%histogram)) deallocate (self%histogram)
      if (allocated(self%type_of)) deallocate (self%type_of)
      if (allocated(self%count_type)) deallocate (self%count_type)
   end subroutine debye_finalize

end module sqc_debye
