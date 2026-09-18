!> Simulation cell geometry: direct and reciprocal lattice vectors.
module sqc_cell
   use sqc_kinds
   implicit none
   private

   public :: cell_t, two_pi

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
      procedure :: fractional => cell_fractional
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

   !> Invert the direct lattice basis and store the reciprocal basis.
   subroutine cell_setup(self)
      class(cell_t), intent(inout) :: self
      real(rk) :: inv(3, 3)

      inv = inverse3(self%a)
      ! b_j = 2 pi * (row j of A^-1) so that a_i . b_j = 2 pi delta_ij.
      self%b(1, :) = two_pi*inv(1, :)
      self%b(2, :) = two_pi*inv(2, :)
      self%b(3, :) = two_pi*inv(3, :)
      self%volume = abs(determinant3(self%a))
   end subroutine cell_setup

   !> Fractional coordinates of a Cartesian point.
   pure function cell_fractional(self, position) result(s)
      class(cell_t), intent(in) :: self
      real(rk), intent(in) :: position(3)
      real(rk) :: s(3)
      real(rk) :: inv(3, 3)

      inv = inverse3(self%a)
      s = matmul(transpose(inv), position - self%origin)
   end function cell_fractional

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
