!> Simulation cell geometry: direct and reciprocal lattice vectors.
module sqc_cell
   use sqc_kinds
   implicit none
   private

   public :: cell_t, two_pi, inverse3

   real(rk), parameter :: two_pi = 2.0_rk*acos(-1.0_rk)

   !> A parallelepiped cell defined by three lattice vectors.
   type :: cell_t
      !> Direct lattice vectors in columns: a(:,1) = a1, a(:,2) = a2, a(:,3) = a3.
      real(rk) :: a(3, 3) = 0.0_rk
      !> Coordinates of the cell origin (lower corner of the box).
      real(rk) :: origin(3) = 0.0_rk
      !> Reciprocal lattice vectors, a_i . b_j = 2 pi delta_ij.
      real(rk) :: b(3, 3) = 0.0_rk
      !> Volume of the cell.
      real(rk) :: volume = 0.0_rk
   contains
      procedure :: init_from_bounds => cell_init_from_bounds
      procedure :: set_vectors => cell_set_vectors
      procedure :: fractional => cell_fractional
      procedure :: widths => cell_widths
   end type cell_t

contains

   !> Build the cell from LAMMPS BOX BOUNDS data.
   !!
   !! lo/hi are the (bound) box limits as written in the dump and tilt holds the
   !! xy, xz, yz tilt factors, which are zero for an orthogonal box.
   subroutine cell_init_from_bounds(self, lo, hi, tilt)
      class(cell_t), intent(inout) :: self
      real(rk), intent(in) :: lo(3), hi(3), tilt(3)
      real(rk) :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz

      xy = tilt(1)
      xz = tilt(2)
      yz = tilt(3)

      ! LAMMPS stores bounding limits for triclinic boxes; recover the true
      ! corners of the restricted triclinic cell from them.
      xlo = lo(1) - min(0.0_rk, xy, xz, xy + xz)
      xhi = hi(1) - max(0.0_rk, xy, xz, xy + xz)
      ylo = lo(2) - min(0.0_rk, yz)
      yhi = hi(2) - max(0.0_rk, yz)
      zlo = lo(3)
      zhi = hi(3)

      self%origin = [xlo, ylo, zlo]
      self%a = 0.0_rk
      self%a(:, 1) = [xhi - xlo, 0.0_rk, 0.0_rk]
      self%a(:, 2) = [xy, yhi - ylo, 0.0_rk]
      self%a(:, 3) = [xz, yz, zhi - zlo]
      call cell_setup(self)
   end subroutine cell_init_from_bounds

   !> Set the cell from three lattice vectors (columns) and refresh the
   !! reciprocal basis and volume.
   subroutine cell_set_vectors(self, vectors)
      class(cell_t), intent(inout) :: self
      real(rk), intent(in) :: vectors(3, 3)
      self%a = vectors
      self%origin = 0.0_rk
      call cell_setup(self)
   end subroutine cell_set_vectors

   !> Invert the direct lattice basis and store the reciprocal basis.
   subroutine cell_setup(self)
      class(cell_t), intent(inout) :: self
      real(rk) :: inv(3, 3)

      inv = inverse3(self%a)
      ! b_j = 2 pi * (row j of A^-1) so that a_i . b_j = 2 pi delta_ij.
      self%b(:, 1) = two_pi*inv(1, :)
      self%b(:, 2) = two_pi*inv(2, :)
      self%b(:, 3) = two_pi*inv(3, :)
      self%volume = abs(determinant3(self%a))
   end subroutine cell_setup

   !> Fractional coordinates of a Cartesian point.
   pure function cell_fractional(self, position) result(s)
      class(cell_t), intent(in) :: self
      real(rk), intent(in) :: position(3)
      real(rk) :: s(3)
      real(rk) :: inv(3, 3)

      inv = inverse3(self%a)
      ! r = A s with the lattice vectors in the columns of A, hence s = A^-1 r.
      s = matmul(inv, position - self%origin)
   end function cell_fractional

   !> Perpendicular widths of the cell, i.e. the distance between the pairs of
   !! lattice planes, in the order of the lattice vectors.
   pure function cell_widths(self) result(widths)
      class(cell_t), intent(in) :: self
      real(rk) :: widths(3)
      real(rk) :: n(3)
      integer :: i, j, k

      if (self%volume <= 0.0_rk) then
         widths = 0.0_rk
         return
      end if
      do i = 1, 3
         j = mod(i, 3) + 1
         k = mod(j, 3) + 1
         n = cross3(self%a(:, j), self%a(:, k))
         widths(i) = self%volume/sqrt(sum(n*n))
      end do
   end function cell_widths

   pure function cross3(u, v) result(w)
      real(rk), intent(in) :: u(3), v(3)
      real(rk) :: w(3)
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
   end function cross3

   pure real(rk) function determinant3(m) result(det)
      real(rk), intent(in) :: m(3, 3)
      det = m(1, 1)*(m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2)) &
            - m(1, 2)*(m(2, 1)*m(3, 3) - m(2, 3)*m(3, 1)) &
            + m(1, 3)*(m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))
   end function determinant3

   !> Explicit 3x3 inverse (the cells involved here are always invertible).
   pure function inverse3(m) result(inv)
      real(rk), intent(in) :: m(3, 3)
      real(rk) :: inv(3, 3), det

      det = determinant3(m)
      inv(1, 1) = (m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2))
      inv(2, 1) = -(m(2, 1)*m(3, 3) - m(2, 3)*m(3, 1))
      inv(3, 1) = (m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))
      inv(1, 2) = -(m(1, 2)*m(3, 3) - m(1, 3)*m(3, 2))
      inv(2, 2) = (m(1, 1)*m(3, 3) - m(1, 3)*m(3, 1))
      inv(3, 2) = -(m(1, 1)*m(3, 2) - m(1, 2)*m(3, 1))
      inv(1, 3) = (m(1, 2)*m(2, 3) - m(1, 3)*m(2, 2))
      inv(2, 3) = -(m(1, 1)*m(2, 3) - m(1, 3)*m(2, 1))
      inv(3, 3) = (m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1))
      inv = inv/det
   end function inverse3

end module sqc_cell
