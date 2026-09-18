!> Cell list for pair enumeration (triclinic cells, per-direction periodicity).
!!
!! Adapted from the nnap2 neighbour list (`~/develop/nnap2/src/neighbour_list.f90`)
!! to the sqcalc conventions: kinds from `sqc_kinds`, geometry from `sqc_cell`,
!! errors through `ierr`/`message` instead of `stop`, and an interface aimed at
!! pair-distance histograms rather than atomic environments.
!!
!! Working principle (unchanged from nnap2):
!!
!! * the cell is divided into `div(1:3)` cell vectors `bins(:,i) = a_i/div(i)`,
!!   chosen so that no cell is narrower than `reach = cutoff + skin`;
!! * every atom is stored in a linked list per cell (`cell_first`/`cell_next`)
!!   and keeps its position relative to its own cell corner (`rel_pos`);
!! * candidate pairs are found by walking the `(2n+1)^3` cell stencil, adding
!!   the precomputed cell shift and keeping pairs closer than `reach`;
!! * the candidate list is kept between frames (Verlet skin) and rebuilt only
!!   when some atom moved by more than `skin/2`, or when `skin = 0`.
!!
!! Positions are unwrapped internally (smallest-image displacement from frame to
!! frame), so the stored image shifts remain valid when atoms cross cell or box
!! boundaries.  A frame-to-frame motion larger than half the box cannot be
!! unwrapped; that forces a rebuild, which is harmless but costs time.
module sqc_neighbour_list
   use sqc_kinds
   use sqc_cell, only: cell_t, inverse3
   implicit none
   private

   public :: neighbour_list_t

   !> Largest cell shift we are willing to store (the values are small because
   !! positions are unwrapped; this only guards against pathological input).
   integer, parameter :: max_shift = 30000

   !> Cell list with a Verlet pair list.
   type :: neighbour_list_t
      !> Geometry and thresholds.
      type(cell_t) :: cell
      real(rk) :: cutoff = 0.0_rk
      real(rk) :: skin = 0.0_rk
      real(rk) :: reach = 0.0_rk
      real(rk) :: reach_sq = 0.0_rk
      logical :: pbc(3) = .true.
      !> Cell decomposition.
      integer :: div(3) = 1
      integer :: nsearch(3) = 0
      integer :: ncells = 0
      integer :: nstencil = 0
      real(rk) :: bins(3, 3) = 0.0_rk
      integer, allocatable :: cell_first(:)
      integer, allocatable :: cell_next(:)
      integer, allocatable :: cell_index(:, :)
      real(rk), allocatable :: rel_pos(:, :)
      integer, allocatable :: stencil(:, :)
      real(rk), allocatable :: stencil_shift(:, :)
      !> Positions: current (unwrapped) and their in-cell form of the previous
      !! frame, plus the reference copy taken when the pair list was built.
      real(rk), allocatable :: pos(:, :)
      real(rk), allocatable :: pos_wrapped(:, :)
      real(rk), allocatable :: pos_build(:, :)
      !> Candidate pairs in CSR form: atom `i` owns pair_start(i) .. pair_start(i+1)-1.
      integer(lk), allocatable :: pair_start(:)
      integer(ik), allocatable :: pair_atom(:)
      integer(ik), allocatable :: pair_shift(:, :)
      integer :: natoms = 0
      integer :: nbuilds = 0
      logical :: started = .false.
      logical :: have_list = .false.
      integer(lk) :: npairs = 0
   contains
      procedure :: configure => nl_configure
      procedure :: begin_frame => nl_begin_frame
      procedure :: finalize => nl_finalize
   end type neighbour_list_t

contains

   !> Set up the cell decomposition.  `natoms` is the (constant) atom count.
   subroutine nl_configure(self, cell, natoms, cutoff, skin, pbc, ierr, message)
      class(neighbour_list_t), intent(inout) :: self
      type(cell_t), intent(in) :: cell
      integer(ik), intent(in) :: natoms
      real(rk), intent(in) :: cutoff
      real(rk), intent(in) :: skin
      logical, intent(in) :: pbc(3)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: widths(3)
      integer :: i, j, k, n
      integer(lk) :: nstencil

      ierr = 0
      message = ''
      if (cutoff <= 0.0_rk) then
         ierr = 1
         message = 'the pair cutoff must be positive'
         return
      end if
      self%cell = cell
      self%natoms = int(natoms)
      self%cutoff = cutoff
      self%skin = max(skin, 0.0_rk)
      self%reach = self%cutoff + self%skin
      self%reach_sq = self%reach**2
      self%pbc = pbc

      widths = cell%widths()
      if (any(widths <= 0.0_rk)) then
         ierr = 1
         message = 'degenerate cell for the neighbour list'
         return
      end if
      do i = 1, 3
         self%div(i) = max(1, int(widths(i)/self%reach))
         self%bins(:, i) = cell%a(:, i)/real(self%div(i), rk)
         self%nsearch(i) = max(1, ceiling(self%reach/widths(i)*real(self%div(i), rk)))
      end do
      self%ncells = self%div(1)*self%div(2)*self%div(3)

      ! Stencil of neighbouring cells with precomputed Cartesian shifts.
      nstencil = int(self%nsearch(1)*2 + 1, lk)*int(self%nsearch(2)*2 + 1, lk) &
                 *int(self%nsearch(3)*2 + 1, lk)
      if (nstencil > 200000_lk) then
         ierr = 1
         write (message, '(a,i0,a)') 'the cell stencil needs ', nstencil, &
            ' cells; reduce the cutoff or use a larger box'
         return
      end if
      self%nstencil = int(nstencil)
      allocate (self%stencil(3, self%nstencil), self%stencil_shift(3, self%nstencil))
      n = 0
      do k = -self%nsearch(3), self%nsearch(3)
         do j = -self%nsearch(2), self%nsearch(2)
            do i = -self%nsearch(1), self%nsearch(1)
               n = n + 1
               self%stencil(:, n) = [i, j, k]
               self%stencil_shift(:, n) = matmul(self%bins, real([i, j, k], rk))
            end do
         end do
      end do

      allocate (self%cell_first(self%ncells), self%cell_next(self%natoms))
      allocate (self%cell_index(3, self%natoms), self%rel_pos(3, self%natoms))
      allocate (self%pos(3, self%natoms), self%pos_wrapped(3, self%natoms))
      allocate (self%pos_build(3, self%natoms))
      allocate (self%pair_start(self%natoms + 1))
      self%cell_first = 0
      self%cell_next = 0
      self%cell_index = 0
      self%rel_pos = 0.0_rk
      self%pos = 0.0_rk
      self%pos_wrapped = 0.0_rk
      self%pos_build = 0.0_rk
      self%pair_start = 1
      self%started = .false.
      self%have_list = .false.
      self%npairs = 0
      self%nbuilds = 0
   end subroutine nl_configure

   !> Release all arrays (keeps the configuration).
   subroutine nl_finalize(self)
      class(neighbour_list_t), intent(inout) :: self
      if (allocated(self%cell_first)) deallocate (self%cell_first)
      if (allocated(self%cell_next)) deallocate (self%cell_next)
      if (allocated(self%cell_index)) deallocate (self%cell_index)
      if (allocated(self%rel_pos)) deallocate (self%rel_pos)
      if (allocated(self%stencil)) deallocate (self%stencil)
      if (allocated(self%stencil_shift)) deallocate (self%stencil_shift)
      if (allocated(self%pos)) deallocate (self%pos)
      if (allocated(self%pos_wrapped)) deallocate (self%pos_wrapped)
      if (allocated(self%pos_build)) deallocate (self%pos_build)
      if (allocated(self%pair_start)) deallocate (self%pair_start)
      if (allocated(self%pair_atom)) deallocate (self%pair_atom)
      if (allocated(self%pair_shift)) deallocate (self%pair_shift)
      self%have_list = .false.
      self%started = .false.
   end subroutine nl_finalize

   !> Bring the list up to date for the given positions.
   !!
   !! `positions` are the Cartesian coordinates of the frame (any wrapping).
   !! The first call initialises the unwrapped reference; later calls continue it
   !! and rebuild the pair list only when needed.
   subroutine nl_begin_frame(self, positions, ierr, message)
      class(neighbour_list_t), intent(inout) :: self
      real(rk), intent(in) :: positions(:, :)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: s(3), r_new(3), dr(3), ds(3), inv(3, 3), disp, max_disp
      integer :: i
      logical :: rebuild

      ierr = 0
      message = ''
      if (size(positions, 2) /= self%natoms) then
         ierr = 1
         message = 'the neighbour list was configured for a different atom count'
         return
      end if

      if (.not. self%started) then
         do i = 1, self%natoms
            s = wrap_fractional(self%cell%fractional(positions(:, i)))
            self%pos(:, i) = matmul(self%cell%a, s)
            self%pos_wrapped(:, i) = self%pos(:, i)
         end do
         self%started = .true.
         rebuild = .true.
      else
         inv = inverse3(self%cell%a)
         do i = 1, self%natoms
            s = wrap_fractional(self%cell%fractional(positions(:, i)))
            r_new = matmul(self%cell%a, s)
            dr = r_new - self%pos_wrapped(:, i)
            ds = matmul(inv, dr)
            ds = ds - anint(ds)
            dr = matmul(self%cell%a, ds)
            self%pos(:, i) = self%pos(:, i) + dr
            self%pos_wrapped(:, i) = r_new
         end do
         rebuild = .false.
         if (self%skin <= 0.0_rk .or. .not. self%have_list) then
            rebuild = .true.
         else
            max_disp = 0.0_rk
            do i = 1, self%natoms
               disp = sum((self%pos(:, i) - self%pos_build(:, i))**2)
               if (disp > max_disp) max_disp = disp
            end do
            rebuild = max_disp > (0.5_rk*self%skin)**2
         end if
      end if

      if (rebuild) call nl_build(self, ierr, message)
   end subroutine nl_begin_frame

   !> (Re)build the cell lists and the Verlet pair list for the current positions.
   subroutine nl_build(self, ierr, message)
      class(neighbour_list_t), intent(inout) :: self
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer :: i, n, cw(3), ci(3), ofs(3), jatom, ncell
      logical :: inside
      integer(lk) :: k, pos_pair
      real(rk) :: dr(3)

      ierr = 0
      message = ''

      ! Cell index (unwrapped) and position relative to the cell corner.
      self%cell_first = 0
      do i = 1, self%natoms
         ci = floor(self%cell%fractional(self%pos(:, i))*real(self%div, rk))
         self%cell_index(:, i) = ci
         self%rel_pos(:, i) = self%pos(:, i) - matmul(self%bins, real(ci, rk))
         cw = ci
         do n = 1, 3
            if (self%pbc(n)) then
               cw(n) = modulo(ci(n), self%div(n))
            else if (ci(n) < 0 .or. ci(n) >= self%div(n)) then
               cw(n) = -1
            end if
         end do
         if (any(cw < 0) .or. any(cw >= self%div)) then
            ierr = 1
            message = 'an atom lies outside a non-periodic box'
            return
         end if
         ncell = linear_cell(self, cw)
         self%cell_next(i) = self%cell_first(ncell)
         self%cell_first(ncell) = i
      end do

      ! First pass: count the candidate pairs per atom.
      self%pair_start = 0
      do i = 1, self%natoms
         ci = self%cell_index(:, i)
         do k = 1, self%nstencil
            ofs = self%stencil(:, k)
            call cell_in_domain(self, ci + ofs, inside, cw)
            if (.not. inside) cycle
            ncell = linear_cell(self, cw)
            jatom = self%cell_first(ncell)
            do while (jatom > 0)
               if (.not. (i == jatom .and. all(ofs == 0))) then
                  dr = self%rel_pos(:, jatom) - self%rel_pos(:, i) + self%stencil_shift(:, k)
                  if (sum(dr*dr) <= self%reach_sq) self%pair_start(i + 1) = &
                     self%pair_start(i + 1) + 1
               end if
               jatom = self%cell_next(jatom)
            end do
         end do
      end do
      ! prefix sum
      self%pair_start(1) = 1
      do i = 1, self%natoms
         self%pair_start(i + 1) = self%pair_start(i + 1) + self%pair_start(i)
      end do
      self%npairs = self%pair_start(self%natoms + 1) - 1
      if (allocated(self%pair_atom)) deallocate (self%pair_atom)
      if (allocated(self%pair_shift)) deallocate (self%pair_shift)
      allocate (self%pair_atom(self%npairs), self%pair_shift(3, self%npairs))
      ! Second pass: fill.
      pos_pair = 0
      do i = 1, self%natoms
         ci = self%cell_index(:, i)
         do k = 1, self%nstencil
            ofs = self%stencil(:, k)
            call cell_in_domain(self, ci + ofs, inside, cw)
            if (.not. inside) cycle
            ncell = linear_cell(self, cw)
            jatom = self%cell_first(ncell)
            do while (jatom > 0)
               if (.not. (i == jatom .and. all(ofs == 0))) then
                  dr = self%rel_pos(:, jatom) - self%rel_pos(:, i) + self%stencil_shift(:, k)
                  if (sum(dr*dr) <= self%reach_sq) then
                     pos_pair = pos_pair + 1
                     self%pair_atom(pos_pair) = int(jatom, ik)
                     self%pair_shift(:, pos_pair) = int(ci + ofs - self%cell_index(:, jatom), ik)
                     if (maxval(abs(self%pair_shift(:, pos_pair))) > max_shift) then
                        ierr = 1
                        message = 'the image shift of a pair grew too large'
                        return
                     end if
                  end if
               end if
               jatom = self%cell_next(jatom)
            end do
         end do
      end do
      self%pos_build = self%pos
      self%have_list = .true.
      self%nbuilds = self%nbuilds + 1
   end subroutine nl_build

   !> Cell index of a wrapped cell, assuming it is inside the domain.
   pure integer function linear_cell(self, c) result(idx)
      class(neighbour_list_t), intent(in) :: self
      integer, intent(in) :: c(3)
      idx = (c(3)*self%div(2) + c(2))*self%div(1) + c(1) + 1
   end function linear_cell

   !> Wrap a target cell back into the domain (periodic directions) or report
   !! that the cell is outside (non-periodic directions).
   pure subroutine cell_in_domain(self, c, inside, wrapped)
      class(neighbour_list_t), intent(in) :: self
      integer, intent(in) :: c(3)
      logical, intent(out) :: inside
      integer, intent(out) :: wrapped(3)
      integer :: i

      inside = .true.
      do i = 1, 3
         if (self%pbc(i)) then
            wrapped(i) = modulo(c(i), self%div(i))
         else
            wrapped(i) = c(i)
            if (c(i) < 0 .or. c(i) >= self%div(i)) inside = .false.
         end if
      end do
   end subroutine cell_in_domain

   !> Fractional coordinates folded into [0, 1).
   pure function wrap_fractional(s) result(wrapped)
      real(rk), intent(in) :: s(3)
      real(rk) :: wrapped(3)
      wrapped = s - floor(s)
   end function wrap_fractional

end module sqc_neighbour_list
