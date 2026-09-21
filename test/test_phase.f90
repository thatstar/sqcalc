! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Validation of the separable phase factors (`sqc_phase`).
!!
!! The tables are only worth using if they reproduce the phase they replace,
!! so the check is direct: compare `phase_factor` against a sine/cosine
!! evaluation for random atoms and random Miller indices, and pin the three
!! properties the callers rely on - the conjugate symmetry of the negative
!! half, the `n = 0` entry, and the invariance under `phi -> phi + 2 pi` that
!! makes wrapped and unwrapped coordinates equivalent.
program test_phase
   use sqc_kinds
   use sqc_phase
   use, intrinsic :: iso_c_binding, only: c_double_complex
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   implicit none

   integer, parameter :: nmax = 4            ! |n_i| <= 4, the grid range
   integer, parameter :: width = 2*nmax + 1
   real(rk), parameter :: two_pi = 2.0_rk*acos(-1.0_rk)
   integer, parameter :: ntrial = 20000
   complex(c_double_complex) :: t(3, 0:width - 1), value, expected, wrapped
   real(rk) :: phi(3), phase, r(3), worst, dev, worse
   integer :: nm(3), n(3), i, j, nfail

   nfail = 0
   nm = nmax
   worst = 0.0_rk
   worse = 0.0_rk
   do i = 1, ntrial
      call random_number(phi)
      phi = two_pi*(phi - 0.5_rk)            ! projections may be negative
      call random_number(r)
      n = int(r*real(2*nmax + 1, rk)) - nmax ! indices in [-nmax, nmax]

      call phase_tables(t, phi, nm)
      value = phase_factor(t, nm, n)
      phase = real(n(1), rk)*phi(1) + real(n(2), rk)*phi(2) + real(n(3), rk)*phi(3)
      expected = cmplx(cos(phase), sin(phase), c_double_complex)
      worst = max(worst, abs(value - expected))

      ! phi -> phi + 2 pi has to leave every factor unchanged.
      ! The check is limited by the rounding of phi + 2 pi itself, whose
      ! argument is three times larger, so it gets its own threshold.
      call phase_tables(t, phi + two_pi, nm)
      wrapped = phase_factor(t, nm, n)
      worse = max(worse, abs(wrapped - expected))

      ! The negative half of each table is the conjugate of the positive half.
      do j = 1, 3
         if (n(j) /= 0) then
            call phase_tables(t, phi, nm)
            if (abs(t(j, nmax - n(j)) - conjg(t(j, nmax + n(j)))) > 1.0e-14_rk) then
               nfail = nfail + 1
               write (error_unit, '(a)') 'FAIL: the factor table is not conjugate symmetric'
            end if
         end if
      end do
   end do
   if (abs(phase_factor(t, nm, [0, 0, 0]) - (1.0_rk, 0.0_rk)) > 1.0e-14_rk) then
      nfail = nfail + 1
      write (error_unit, '(a)') 'FAIL: the n = 0 factor is not 1'
   end if
   if (worst > 1.0e-14_rk) then
      nfail = nfail + 1
      write (error_unit, '(a,es12.3)') 'FAIL: factor deviates from the direct phase by ', worst
   end if
   if (worse > 1.0e-13_rk) then
      nfail = nfail + 1
      write (error_unit, '(a,es12.3)') 'FAIL: phi + 2 pi changed the factor by ', worse
   end if

   if (nfail > 0) then
      write (error_unit, '(a,i0,a)') 'FAIL: ', nfail, ' phase checks failed'
      stop 1
   end if
   write (output_unit, '(a,i0,a,es10.3)') 'test_phase: ', ntrial, &
      ' random atoms and indices, worst deviation ', worst

end program test_phase
