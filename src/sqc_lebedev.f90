! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Lebedev quadrature on the unit sphere, used by `--q shell:Q,ACC`.
!!
!! A Lebedev rule integrates the spherical harmonics exactly up to its order,
!! so averaging a quantity that varies smoothly over a spherical shell in
!! reciprocal space needs only a few dozen directions.  The rules are stored
!! here in the symmetry reduced form the published tables use: one generating
!! point [a, b, c] with 0 <= a <= b <= c per orbit of the octahedral group,
!! plus a single weight, and the full grid is the orbit of that point under the
!! coordinate permutations and sign changes.  The 19 rows below expand to the
!! 50, 110 and 194 point grids of `scipy.integrate.lebedev_rule(11|17|23)`
!! exactly, which is the reference the shell average is tested against.
!!
!! The weights sum to 4 pi, so the shell average of X(q) is
!! `sum_k w_k X(q_k) / sum_k w_k` (see `dyn_shell_average` in sqc_dynamics).
module sqc_lebedev
   use sqc_kinds
   implicit none
   private

   public :: lebedev_rule, lebedev_points, lebedev_order_from_name, &
             lebedev_low, lebedev_medium, lebedev_high, lebedev_max_points, &
             lebedev_reduce_pairs

   !> Accuracy names of `--q shell:Q,ACC` and the order they select.
   integer, parameter :: lebedev_low = 11
   integer, parameter :: lebedev_medium = 17
   integer, parameter :: lebedev_high = 23

   !> Points of the largest rule (a grid array has to hold this many).
   integer, parameter :: lebedev_max_points = 194

   !> Symmetry reduced rules: order, generating point (a <= b <= c) and weight.
   integer, parameter :: n_lebedev_rows = 19
   integer, parameter :: lebedev_row_order(n_lebedev_rows) = [ &
      11, 11, 11, 11, &
      17, 17, 17, 17, 17, 17, &
      23, 23, 23, 23, 23, 23, 23, 23, 23]
   real(rk), parameter :: lebedev_row_point(3, n_lebedev_rows) = reshape([ &
      0.00000000000000000e+00_rk, 0.00000000000000000e+00_rk, 1.00000000000000000e+00_rk, &
      0.00000000000000000e+00_rk, 7.07106781186547573e-01_rk, 7.07106781186547573e-01_rk, &
      3.01511344577763574e-01_rk, 3.01511344577763574e-01_rk, 9.04534033733290888e-01_rk, &
      5.77350269189625731e-01_rk, 5.77350269189625731e-01_rk, 5.77350269189625731e-01_rk, &
      0.00000000000000000e+00_rk, 0.00000000000000000e+00_rk, 1.00000000000000000e+00_rk, &
      0.00000000000000000e+00_rk, 4.78369028812150210e-01_rk, 8.78158910604066145e-01_rk, &
      1.85115635344736212e-01_rk, 1.85115635344736212e-01_rk, 9.65124035086594056e-01_rk, &
      2.15957291845848443e-01_rk, 6.90421048382292235e-01_rk, 6.90421048382292235e-01_rk, &
      3.95689473055941876e-01_rk, 3.95689473055941876e-01_rk, 8.28769981252592269e-01_rk, &
      5.77350269189625731e-01_rk, 5.77350269189625731e-01_rk, 5.77350269189625731e-01_rk, &
      0.00000000000000000e+00_rk, 0.00000000000000000e+00_rk, 1.00000000000000000e+00_rk, &
      0.00000000000000000e+00_rk, 3.45770219761128317e-01_rk, 9.38319218137591560e-01_rk, &
      0.00000000000000000e+00_rk, 7.07106781186547573e-01_rk, 7.07106781186547573e-01_rk, &
      1.29933544765006709e-01_rk, 1.29933544765006709e-01_rk, 9.82972302707253220e-01_rk, &
      1.59041710538353004e-01_rk, 5.25118572443642018e-01_rk, 8.36036015482458872e-01_rk, &
      2.89246562757543901e-01_rk, 2.89246562757543901e-01_rk, 9.12509096867473724e-01_rk, &
      3.14196994182586287e-01_rk, 6.71297344269522589e-01_rk, 6.71297344269522589e-01_rk, &
      4.44693317871743710e-01_rk, 4.44693317871743710e-01_rk, 7.77493219314767114e-01_rk, &
      5.77350269189625731e-01_rk, 5.77350269189625731e-01_rk, 5.77350269189625731e-01_rk], &
      [3, n_lebedev_rows])
   real(rk), parameter :: lebedev_row_weight(n_lebedev_rows) = [ &
      1.59572960182338713e-01_rk, 2.83685262546379879e-01_rk, 2.53505610897311273e-01_rk, &
      2.65071880146638794e-01_rk, &
      4.81074658513965941e-02_rk, 1.21830917385521376e-01_rk, 1.03191734088330406e-01_rk, &
      1.24945096872513303e-01_rk, 1.20580249028527889e-01_rk, 1.23071735281670175e-01_rk, &
      2.23975506210384659e-02_rk, 6.34833699346415564e-02_rk, 7.18407589348473569e-02_rk, &
      5.16072821665131617e-02_rk, 6.94951574710432202e-02_rk, 6.48203268035104641e-02_rk, &
      7.04810541680701286e-02_rk, 6.93509275937110037e-02_rk, 7.00371986012484904e-02_rk]

contains

   !> Number of points of a rule, 0 when the order is not one of the three.
   pure integer function lebedev_points(order) result(n)
      integer, intent(in) :: order

      select case (order)
      case (lebedev_low)
         n = 50
      case (lebedev_medium)
         n = 110
      case (lebedev_high)
         n = 194
      case default
         n = 0
      end select
   end function lebedev_points

   !> Order selected by an accuracy name, -1 when it is not one of them.
   pure integer function lebedev_order_from_name(name) result(order)
      character(len=*), intent(in) :: name

      select case (trim(name))
      case ('low')
         order = lebedev_low
      case ('medium')
         order = lebedev_medium
      case ('high')
         order = lebedev_high
      case default
         order = -1
      end select
   end function lebedev_order_from_name

   !> Expand a rule into `points(3, :)` (unit vectors) and its `weights(:)`.
   subroutine lebedev_rule(order, points, weights, ierr, message)
      integer, intent(in) :: order
      real(rk), intent(out) :: points(:, :)
      real(rk), intent(out) :: weights(:)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      !> The six coordinate permutations that generate an octahedral orbit.
      integer, parameter :: perm(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, &
                                                  2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
      integer :: want, row, ip, is, k, n
      real(rk) :: rep(3), candidate(3)
      logical :: duplicate

      ierr = 0
      message = ''
      want = lebedev_points(order)
      if (want == 0) then
         ierr = 1
         write (message, '(a,i0,a)') 'unsupported Lebedev order ', order, &
            ' (available: 11, 17 and 23)'
         return
      end if
      if (size(points, 1) < 3 .or. size(points, 2) < want .or. size(weights) < want) then
         ierr = 1
         message = 'the Lebedev output arrays are too small for the requested order'
         return
      end if

      n = 0
      do row = 1, n_lebedev_rows
         if (lebedev_row_order(row) /= order) cycle
         rep = lebedev_row_point(:, row)
         do ip = 1, 6
            do is = 0, 7
               do k = 1, 3
                  candidate(k) = rep(perm(k, ip))
                  if (btest(is, k - 1)) candidate(k) = -candidate(k)
               end do
               ! Members of one orbit are generated exactly from the same
               ! numbers, so the duplicates are equal to well within this
               ! tolerance; distinct directions are more than 0.1 apart.
               duplicate = .false.
               do k = 1, n
                  if (maxval(abs(points(:, k) - candidate)) <= 1.0e-12_rk) then
                     duplicate = .true.
                     exit
                  end if
               end do
               if (duplicate) cycle
               n = n + 1
               points(:, n) = candidate
               weights(n) = lebedev_row_weight(row)
            end do
         end do
      end do
      if (n /= want) then
         ierr = 1
         write (message, '(a,i0,a,i0,a,i0,a)') 'the Lebedev rule of order ', order, &
            ' expanded to ', n, ' points instead of ', want, '; the table is corrupt'
      end if
   end subroutine lebedev_rule

   !> Merge every `+-` pair of a rule into one vector of twice the weight.
   !!
   !! A Lebedev rule is built from orbits of the octahedral group, which
   !! contains the inversion, so the antipode of every point is itself a point
   !! of the rule with the very same weight.  Every quantity the shell average
   !! is applied to is even in `q` - `S_4` and `F_s` are squared amplitudes
   !! and `F`, `S` are read through their real parts - so one vector of a pair
   !! can carry both.  Halving the rule therefore leaves the average unchanged
   !! (to rounding) while halving the modes the shell run accumulates.
   !!
   !! The representative is the one whose first non-zero component is
   !! positive, the same convention `modes_canonical` uses for the lattice
   !! `+-` reduction.  The points and weights are compacted in place and `n`
   !! is set to the number of survivors.
   subroutine lebedev_reduce_pairs(points, weights, n)
      real(rk), intent(inout) :: points(:, :)
      real(rk), intent(inout) :: weights(:)
      integer(lk), intent(inout) :: n
      integer(lk) :: i, keep

      keep = 0
      do i = 1, n
         if (.not. lebedev_pair_head(points(:, i))) cycle
         keep = keep + 1
         points(:, keep) = points(:, i)
         weights(keep) = 2.0_rk*weights(i)
      end do
      n = keep
   end subroutine lebedev_reduce_pairs

   !> True for the member of a `+-` pair that `lebedev_reduce_pairs` keeps.
   pure logical function lebedev_pair_head(point) result(head)
      real(rk), intent(in) :: point(3)
      integer :: k

      do k = 1, 3
         if (abs(point(k)) > 0.0_rk) then
            head = point(k) > 0.0_rk
            return
         end if
      end do
      head = .false.
   end function lebedev_pair_head

end module sqc_lebedev
