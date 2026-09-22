! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Validation of the two-theta and Lorentz-polarization conversions (`--xrd`).
!!
!! The XRD output maps the momentum transfer of each reciprocal lattice point
!! to a two-theta bin; both conversions are pure functions of q and the
!! wavelength, so their edges are cheap to pin here as well as through the
!! end to end pattern the python suite checks.
program test_xrd
   use sqc_kinds
   use sqc_structure_factor, only: two_theta_deg, lorentz_polarization
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   implicit none

   real(rk), parameter :: pi_x = 3.14159265358979323846264338327950288_rk
   real(rk), parameter :: lambda = 1.541838_rk
   real(rk), parameter :: tolerance = 1.0e-10_rk
   real(rk) :: q, tth, worst, lp
   integer :: i
   real(rk), parameter :: fcc_a = 3.52_rk
   real(rk), parameter :: peaks(3, 3) = reshape([1.0_rk, 1.0_rk, 1.0_rk, &
                                                 2.0_rk, 0.0_rk, 0.0_rk, &
                                                 2.0_rk, 2.0_rk, 0.0_rk], [3, 3])
   real(rk) :: expected(3) = [44.585405_rk, 51.955510_rk, 76.552967_rk]

   worst = 0.0_rk

   ! The first three fcc reflections of the LAMMPS Ni example land on the
   ! textbook angles at Cu Kalpha.
   do i = 1, 3
      q = sqrt(sum(peaks(:, i)**2))*2.0_rk*pi_x/fcc_a
      tth = two_theta_deg(q, lambda)
      worst = max(worst, abs(tth - expected(i))/expected(i))
   end do
   call require(worst < 1.0e-7_rk, 'fcc Ni reflections at Cu Kalpha', worst)

   ! Round trip: the q of the 45 degree line has to come back as 45 degrees.
   q = 4.0_rk*pi_x*sin(pi_x/8.0_rk)/lambda
   call require(abs(two_theta_deg(q, lambda) - 45.0_rk) < 1.0e-9_rk, &
                'two_theta_deg inverts q = 4 pi sin(theta)/lambda', &
                abs(two_theta_deg(q, lambda) - 45.0_rk))

   ! Backscattering is out of range: the argument must be clipped, not NaN.
   call require(two_theta_deg(1.0e9_rk, lambda) == 180.0_rk, &
                'two_theta_deg clips q lambda/4 pi > 1', &
                abs(two_theta_deg(1.0e9_rk, lambda) - 180.0_rk))

   ! Lorentz-polarization: (1 + cos^2 2theta)/(sin^2 theta cos theta), which
   ! is 2 sqrt(2) at 90 degrees and diverges towards 0 and 180 degrees.
   call require(abs(lorentz_polarization(90.0_rk) - 2.0_rk*sqrt(2.0_rk)) < 1.0e-12_rk, &
                'LP(90 deg) = 2 sqrt(2)', &
                abs(lorentz_polarization(90.0_rk) - 2.0_rk*sqrt(2.0_rk)))
   call require(lorentz_polarization(1.0_rk) > lorentz_polarization(90.0_rk) .and. &
                lorentz_polarization(179.0_rk) > lorentz_polarization(90.0_rk), &
                'LP grows towards the ends of the range', 0.0_rk)

   write (output_unit, '(a)') 'ok xrd: two-theta and Lorentz-polarization'

contains

   !> Fail loudly with the value that was out of tolerance.
   subroutine require(condition, label, value)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      real(rk), intent(in) :: value

      if (condition) return
      write (error_unit, '(a,a,a,es12.4)') 'FAIL: ', trim(label), ', deviation ', value
      stop 1
   end subroutine require

end program test_xrd
