! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Validation of the embedded Lebedev rules (`--dyn-q shell`).
!!
!! The rules have to be a quadrature of the unit sphere: `lebedev_points`
!! points, no duplicates, weights summing to 4 pi, and an exact result for
!! every monomial whose total degree does not exceed the order of the rule.
!! That last check is what makes the shell average trustworthy; it fails
!! loudly if a row of the symmetry reduced table is mistyped.
program test_lebedev
   use sqc_kinds
   use sqc_lebedev
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   implicit none

   integer, parameter :: orders(3) = [lebedev_low, lebedev_medium, lebedev_high]
   integer, parameter :: counts(3) = [50, 110, 194]
   real(rk), parameter :: pi = 3.14159265358979323846264338327950288_rk
   real(rk) :: points(3, lebedev_max_points), weights(lebedev_max_points)
   real(rk) :: sum_w, quadrature, exact, worst
   integer :: i, j, iorder, n, a, b, c, ierr
   integer(lk) :: nred
   character(len=256) :: message
   logical :: paired

   worst = 0.0_rk
   do iorder = 1, 3
      n = lebedev_points(orders(iorder))
      if (n /= counts(iorder)) then
         write (error_unit, '(a,i0,a,i0)') 'FAIL: order ', orders(iorder), &
            ' has ', n, ' points in the table'
         stop 1
      end if

      ! Poison the buffer so that unwritten slots show up in the checks below.
      points = 0.0_rk
      weights = -1.0_rk
      call lebedev_rule(orders(iorder), points, weights, ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a)') 'FAIL: lebedev_rule: '//trim(message)
         stop 1
      end if

      sum_w = sum(weights(1:n))
      if (abs(sum_w - 4.0_rk*pi) > 1.0e-12_rk) then
         write (error_unit, '(a,i0,a,es14.6)') 'FAIL: order ', orders(iorder), &
            ' weights sum to ', sum_w
         stop 1
      end if
      if (any(weights(n + 1:) /= -1.0_rk)) then
         write (error_unit, '(a,i0,a)') 'FAIL: order ', orders(iorder), &
            ' wrote more than ', n, ' weights'
         stop 1
      end if

      do i = 1, n
         if (abs(sum(points(:, i)**2) - 1.0_rk) > 1.0e-12_rk) then
            write (error_unit, '(a,i0,a,i0)') 'FAIL: order ', orders(iorder), &
               ' point ', i, ' is not on the unit sphere'
            stop 1
         end if
         do j = 1, i - 1
            if (all(points(:, i) == points(:, j))) then
               write (error_unit, '(a,i0,a,i0)') 'FAIL: order ', orders(iorder), &
                  ' repeats point ', i
               stop 1
            end if
         end do
      end do

      ! Exactness of the rule on every monomial of degree <= its order.
      do a = 0, orders(iorder)
         do b = 0, orders(iorder) - a
            do c = 0, orders(iorder) - a - b
               quadrature = 0.0_rk
               do i = 1, n
                  quadrature = quadrature + weights(i) &
                     *points(1, i)**a*points(2, i)**b*points(3, i)**c
               end do
               exact = sphere_moment(a, b, c)
               worst = max(worst, abs(quadrature - exact)/max(abs(exact), 1.0_rk))
               if (abs(quadrature - exact) > 1.0e-12_rk*max(abs(exact), 1.0_rk)) then
                  write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,es12.4,a,es12.4)') &
                     'FAIL: order ', orders(iorder), ' x^', a, ' y^', b, ' z^', c, &
                     ' gives ', quadrature, ' instead of ', exact
                  stop 1
               end if
            end do
         end do
      end do

      ! Every point has its antipode in the rule with the same weight, which is
      ! what lets a shell run keep one vector of each pair.  The reduced rule
      ! has to stay a quadrature for the even functions the shell averages
      ! (S4 and F_s are squared amplitudes, F and S are read through their
      ! real parts), so check the weight sum and every even monomial again.
      do i = 1, n
         paired = .false.
         do j = 1, n
            if (all(points(:, j) == -points(:, i))) then
               paired = .true.
               if (weights(j) /= weights(i)) then
                  write (error_unit, '(a,i0,a,i0,a)') 'FAIL: order ', orders(iorder), &
                     ' points ', i, ' and its antipode carry different weights'
                  stop 1
               end if
               exit
            end if
         end do
         if (.not. paired) then
            write (error_unit, '(a,i0,a,i0,a)') 'FAIL: order ', orders(iorder), &
               ' point ', i, ' has no antipode in the rule'
            stop 1
         end if
      end do

      nred = n
      call lebedev_reduce_pairs(points, weights, nred)
      if (nred /= n/2) then
         write (error_unit, '(a,i0,a,i0,a,i0)') 'FAIL: order ', orders(iorder), &
            ' reduced to ', nred, ' of ', n, ' directions'
         stop 1
      end if
      sum_w = sum(weights(1:nred))
      if (abs(sum_w - 4.0_rk*pi) > 1.0e-12_rk) then
         write (error_unit, '(a,i0,a,es14.6)') 'FAIL: the reduced order ', orders(iorder), &
            ' weights sum to ', sum_w
         stop 1
      end if
      do a = 0, orders(iorder), 2
         do b = 0, orders(iorder) - a, 2
            do c = 0, orders(iorder) - a - b, 2
               quadrature = 0.0_rk
               do i = 1, nred
                  quadrature = quadrature + weights(i) &
                     *points(1, i)**a*points(2, i)**b*points(3, i)**c
               end do
               exact = sphere_moment(a, b, c)
               if (abs(quadrature - exact) > 1.0e-12_rk*max(abs(exact), 1.0_rk)) then
                  write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,es12.4,a,es12.4)') &
                     'FAIL: reduced order ', orders(iorder), ' x^', a, ' y^', b, ' z^', c, &
                     ' gives ', quadrature, ' instead of ', exact
                  stop 1
               end if
            end do
         end do
      end do
   end do

   write (output_unit, '(a,es10.2)') 'PASS: Lebedev rules and their +- halves are exact, worst ', &
      worst

contains

   !> Integral of `x^a y^b z^c` over the unit sphere (0 for odd powers).
   pure real(rk) function sphere_moment(a, b, c) result(value)
      integer, intent(in) :: a, b, c

      if (mod(a, 2) /= 0 .or. mod(b, 2) /= 0 .or. mod(c, 2) /= 0) then
         value = 0.0_rk
         return
      end if
      value = 4.0_rk*pi*double_factorial(a - 1)*double_factorial(b - 1) &
              *double_factorial(c - 1)/double_factorial(a + b + c + 1)
   end function sphere_moment

   !> n!! with 0!! = (-1)!! = 1.
   pure real(rk) function double_factorial(n) result(value)
      integer, intent(in) :: n
      integer :: k

      value = 1.0_rk
      if (n < 1) return
      do k = n, 1, -2
         value = value*real(k, rk)
      end do
   end function double_factorial

end program test_lebedev
