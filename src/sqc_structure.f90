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
module sqc_structure
   use sqc_kinds
   use sqc_cell, only: cell_t, two_pi
   use sqc_dump, only: frame_t
   use sqc_weights, only: weight_scheme_t, weight_xray
   use sqc_finufft, only: finufft_opts_t, finufft_default_opts, finufft_makeplan, &
                          finufft_setpts, finufft_execute, finufft_destroy, &
                          mode_cmcl, type1
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_double, c_double_complex, &
                                          c_int, c_int64_t, c_associated
   implicit none
   private

   public :: structure_factor_t, nufft_structure_factor_t, direct_structure_factor_t, &
             method_nufft, method_direct, norm_mean, norm_self, norm_natom, no_unit, &
             sf_prepare_species

   !> Sentinel for "do not write this table".  A plain negative test would be
   !! wrong because OPEN(NEWUNIT=) may hand out negative unit numbers.
   integer, parameter :: no_unit = -huge(1)

   !> Available evaluation methods.
   integer, parameter :: method_nufft = 1
   integer, parameter :: method_direct = 2

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
      !> Number of (frame, mode) contributions per shell.
      integer(lk), allocatable :: shell_count(:)
      !> Accumulated numerator and denominator per grid point.
      real(rk), allocatable :: gnum(:), gden(:)
      !> Frame bookkeeping.
      integer(ik) :: natoms = 0
      integer(ik) :: ntypes = 0
      integer(lk) :: nframes = 0
   contains
      procedure :: configure => sf_configure
      procedure :: write_results => sf_write_results
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
   contains
      procedure :: method_setup => nufft_setup
      procedure :: accumulate_frame => nufft_accumulate
      final :: nufft_finalize
   end type nufft_structure_factor_t

   !> Direct O(N * nmodes) evaluation, used as reference and offline fallback.
   type, extends(structure_factor_t) :: direct_structure_factor_t
      !> LAMMPS type id of each atom.
      integer(ik), allocatable :: type_of(:)
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
      real(rk), allocatable :: qvec_all(:, :), qlen_all(:), qlen_keep(:), qvec_keep(:, :)
      integer, allocatable :: hkl_all(:, :), hkl_keep(:, :)
      integer(lk), allocatable :: gidx_all(:), gidx_keep(:)
      real(rk) :: q(3), ql, dmode, length_a
      integer :: i, p1, p2, p3, h(3), n_keep, t, maxtype
      integer(lk) :: nmax, g

      ierr = 0
      message = ''
      self%natoms = frame%natoms
      maxtype = 1
      do i = 1, frame%natoms
         maxtype = max(maxtype, int(frame%type_id(i)))
      end do
      self%ntypes = max(maxtype, size(scheme%symbols))

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

      nmax = self%gridpoints
      allocate (hkl_all(3, nmax), qvec_all(3, nmax), qlen_all(nmax), gidx_all(nmax))
      n_keep = 0
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
               if (ql > self%qmax) cycle
               if (ql < self%qmin) cycle
               n_keep = n_keep + 1
               hkl_all(:, n_keep) = h
               qvec_all(:, n_keep) = q
               qlen_all(n_keep) = ql
               g = int(p1, lk) + int(self%modes(1), lk)*(int(p2, lk) - 1 &
                   + int(self%modes(2), lk)*(int(p3, lk) - 1))
               gidx_all(n_keep) = g
            end do
         end do
      end do
      if (n_keep == 0) then
         ierr = 1
         message = 'no reciprocal lattice points fall inside the requested q range'
         return
      end if
      self%nmodes = n_keep
      allocate (hkl_keep(3, n_keep), qvec_keep(3, n_keep), qlen_keep(n_keep), &
                gidx_keep(n_keep))
      do i = 1, n_keep
         hkl_keep(:, i) = hkl_all(:, i)
         qvec_keep(:, i) = qvec_all(:, i)
         qlen_keep(i) = qlen_all(i)
         gidx_keep(i) = gidx_all(i)
      end do
      deallocate (hkl_all, qvec_all, qlen_all, gidx_all)
      call move_alloc(hkl_keep, self%hkl)
      call move_alloc(qvec_keep, self%qvec)
      call move_alloc(qlen_keep, self%qlen)
      call move_alloc(gidx_keep, self%gidx)

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
      self%num = 0.0_rk
      self%den = 0.0_rk
      self%shell_count = 0
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

   !> Write the 1D shell table and (optionally) the reciprocal grid table.
   subroutine sf_write_results(self, shell_unit, grid_unit, ierr, message)
      class(structure_factor_t), intent(in) :: self
      integer, intent(in) :: shell_unit, grid_unit
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: qc, value, frames
      integer :: s
      integer(lk) :: i, g

      ierr = 0
      message = ''
      ! num(:) accumulates every frame while den(:) is the single frame
      ! normalization, so the average over frames is num/(nframes*den).
      frames = real(max(self%nframes, 1_lk), rk)
      if (shell_unit /= no_unit) then
         write (shell_unit, '(a)') '# q S(q)'
         do s = 1, self%nq
            qc = self%qmin + (real(s, rk) - 0.5_rk)*self%shell_dq
            value = 0.0_rk
            if (self%den(s) > 0.0_rk) value = self%num(s)/(frames*self%den(s))
            write (shell_unit, '(f14.6,2x,es20.12)') qc, value
         end do
      end if

      if (grid_unit /= no_unit .and. self%want_grid) then
         write (grid_unit, '(a)') '# qx qy qz S(q)'
         do i = 1, self%nmodes
            g = self%gidx(i)
            value = 0.0_rk
            if (self%gden(g) > 0.0_rk) value = self%gnum(g)/(frames*self%gden(g))
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
      allocate (self%total(self%gridpoints))
      self%strengths = (0.0_rk, 0.0_rk)
      self%fk = (0.0_rk, 0.0_rk)
      self%total = (0.0_rk, 0.0_rk)

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

      ! One transform per species, summed with its scattering amplitude.
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
         if (self%q_dependent) then
            !$omp parallel do schedule(static)
            do im = 1, int(self%nmodes)
               g = self%gidx(im)
               self%total(g) = self%total(g) &
                               + cmplx(self%amp_table(im, isp), 0.0_rk, c_double_complex) &
                                 *self%fk(g)
            end do
            !$omp end parallel do
         else
            amp = self%amp_const(isp)
            !$omp parallel do schedule(static)
            do im = 1, int(self%nmodes)
               g = self%gidx(im)
               self%total(g) = self%total(g) &
                               + cmplx(amp, 0.0_rk, c_double_complex)*self%fk(g)
            end do
            !$omp end parallel do
         end if
      end do

      ! Accumulate |rho|^2 into the shell and grid statistics.
      do im = 1, int(self%nmodes)
         g = self%gidx(im)
         value = real(self%total(g)*conjg(self%total(g)), rk)
         sh = self%shell(im)
         self%num(sh) = self%num(sh) + value
         if (self%want_grid) self%gnum(g) = self%gnum(g) + value
      end do
      self%nframes = self%nframes + 1
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
      integer(lk) :: work

      ierr = 0
      message = ''
      allocate (self%type_of(frame%natoms))
      self%type_of = frame%type_id
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
      complex(c_double_complex) :: acc, total
      integer :: i, im, isp, sh, atype
      integer(lk) :: g

      ierr = 0
      message = ''
      do i = 1, frame%natoms
         s = frame%cell%fractional(frame%pos(:, i))
         s = s - floor(s)
         self%cart(:, i) = matmul(frame%cell%a, s)
      end do

      !$omp parallel do schedule(static) private(im, isp, i, amp, phase, acc, total, atype)
      do im = 1, int(self%nmodes)
         total = (0.0_rk, 0.0_rk)
         do isp = 1, self%ntypes
            amp = scheme%amplitude(int(isp, ik), self%qlen(im))
            acc = (0.0_rk, 0.0_rk)
            do i = 1, frame%natoms
               if (self%type_of(i) /= isp) cycle
               phase = self%qvec(1, im)*self%cart(1, i) + self%qvec(2, im)*self%cart(2, i) &
                       + self%qvec(3, im)*self%cart(3, i)
               acc = acc + cmplx(cos(phase), sin(phase), c_double_complex)
            end do
            total = total + cmplx(amp, 0.0_rk, c_double_complex)*acc
         end do
         self%intensity(im) = real(total*conjg(total), rk)
      end do
      !$omp end parallel do

      do im = 1, int(self%nmodes)
         g = self%gidx(im)
         value = self%intensity(im)
         sh = self%shell(im)
         self%num(sh) = self%num(sh) + value
         if (self%want_grid) self%gnum(g) = self%gnum(g) + value
      end do
      self%nframes = self%nframes + 1
   end subroutine direct_accumulate

end module sqc_structure
