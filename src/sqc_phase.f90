! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Separable phase factors for reciprocal-lattice mode sums.
!!
!! A density amplitude sums one phase per (mode, atom),
!!
!!   rho(q_n) = sum_i exp(i q_n . r_i),      q_n = n_1 b_1 + n_2 b_2 + n_3 b_3,
!!
!! and for a lattice mode that phase factors into three independent parts.
!! With the projections of the atom on the reciprocal basis,
!!
!!   phi_j = b_j . r_i,
!!
!! the phase is a sum of three terms, each of which depends on one index
!! component only,
!!
!!   q_n . r_i = n_1 phi_1 + n_2 phi_2 + n_3 phi_3,
!!
!! so the exponential is a product of three factors:
!!
!!   exp(i q_n . r_i) = t_1(n_1) t_2(n_2) t_3(n_3),    t_j(n) = exp(i n phi_j).
!!
!! Building `t_j` by a geometric recurrence costs one sine/cosine pair and
!! `n_max` complex multiplies per axis, once per atom, and every mode then
!! costs two complex multiplies instead of one sine/cosine pair.  This is the
!! identity Ewald-style codes exploit; because the atoms are off the grid
!! here, the factors come from a recurrence rather than from a lookup table.
!!
!! `t_j` is invariant under `phi_j -> phi_j + 2 pi`, so it does not matter
!! whether the coordinates are wrapped, and the recurrence accumulates only
!! about one ulp per step, i.e. O(n_max) ulp in total.
module sqc_phase
   use sqc_kinds
   use, intrinsic :: iso_c_binding, only: c_double_complex
   implicit none
   private

   public :: phase_tables, phase_factor

   !> Smallest mode count for which the separable evaluation is worth it.
   !!
   !! The factor tables are a fixed cost per atom, so a handful of modes is
   !! cheaper to evaluate with one sine/cosine pair per (mode, atom).  Measured
   !! with four modes: 22.2 s separable against 19.7 s direct, same window.
   integer, parameter, public :: phase_min_modes = 16

contains

   !> Fill the three factor tables of one atom.
   !!
   !! `t(j, n + nmax(j))` receives `exp(i n phi_j)` for `n = -nmax(j) .. nmax(j)`.
   !! The negative half is the conjugate of the positive half, so the table
   !! costs `nmax(j)` complex multiplies per axis.  The table must be at least
   !! `2*maxval(nmax)` wide, which lets the three axes share one array even
   !! when their index ranges differ (a triclinic box).
   subroutine phase_tables(t, phi, nmax)
      complex(c_double_complex), intent(inout) :: t(:, 0:)
      real(rk), intent(in) :: phi(3)
      integer, intent(in) :: nmax(3)
      complex(c_double_complex) :: u, prod
      integer :: j, n, off

      do j = 1, 3
         off = max(nmax(j), 0)
         u = cmplx(cos(phi(j)), sin(phi(j)), c_double_complex)
         t(j, off) = (1.0_rk, 0.0_rk)
         prod = (1.0_rk, 0.0_rk)
         do n = 1, off
            prod = prod*u
            t(j, off + n) = prod
            t(j, off - n) = conjg(prod)
         end do
      end do
   end subroutine phase_tables

   !> `exp(i q_n . r)` for the Miller indices `n`, from the tables of
   !! `phase_tables`.
   pure complex(c_double_complex) function phase_factor(t, nmax, n) result(value)
      complex(c_double_complex), intent(in) :: t(:, 0:)
      integer, intent(in) :: nmax(3)
      integer, intent(in) :: n(3)

      value = t(1, n(1) + nmax(1))*t(2, n(2) + nmax(2))*t(3, n(3) + nmax(3))
   end function phase_factor

end module sqc_phase
