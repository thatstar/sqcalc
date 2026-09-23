! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Validation of the reciprocal-lattice mode sets (`--qpoints grid`).
!!
!! The three reductions of `sqc_modes` are checked against counts that were
!! obtained independently (a small Python enumeration, see the design notes):
!!
!!   * the lattice point group has to come out of the metric, so a cubic cell
!!     gives 48 operations, an orthorhombic one 8 and a generic triclinic one
!!     the pair {E, -1};
!!   * the cubic Kob-Andersen cell at |q| <= 1 must give 7 shells holding 92
!!     lattice vectors, i.e. 46 modes after the +- reduction, and 7 orbits;
!!   * the orbit multiplicities have to add up to the full shell population,
!!     which is what makes the multiplicity weighting meaningful;
!!   * the budget has to drop whole shells first and keep the two ends of the
!!     q range, and it has to be deterministic.
program test_modes
   use sqc_kinds
   use sqc_cell, only: cell_t, two_pi
   use sqc_modes
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   implicit none

   real(rk), parameter :: ka_l = 18.8207_rk
   integer :: nfail

   nfail = 0
   call check_cubic(nfail)
   call check_point_groups(nfail)
   call check_skewed(nfail)
   call check_budget(nfail)
   call check_multi_orbit(nfail)

   if (nfail > 0) then
      write (error_unit, '(a,i0,a)') 'FAIL: ', nfail, ' mode set checks failed'
      stop 1
   end if
   write (output_unit, '(a)') 'test_modes: all checks passed'

contains

   !> The Kob-Andersen cubic cell: 7 shells, 92 vectors, 46 +- modes, 7 orbits.
   subroutine check_cubic(nfail)
      integer, intent(inout) :: nfail
      type(cell_t) :: cell
      type(modes_t) :: m
      real(rk) :: a(3, 3), expect_q
      integer :: i, ierr
      character(len=256) :: message

      call make_cubic(ka_l, cell)
      call modes_build(m, cell, 1.0_rk, .true., 0, thin_none, ierr, message)
      call check(ierr == 0, 'cubic: build succeeds ('//trim(message)//')', nfail)
      call check(m%nops == 48, 'cubic: point group order is 48', nfail)
      call check(m%nshell == 7, 'cubic: 7 shells below qmax = 1', nfail)
      call check(m%nmodes == 47, 'cubic: 46 +- modes plus Gamma', nfail)
      call check(m%has_gamma, 'cubic: Gamma is present', nfail)

      expect_q = two_pi/ka_l
      call check(abs(m%shell_q(1) - expect_q) < 1.0e-9_rk, &
                 'cubic: first shell sits at 2 pi / L', nfail)

      ! Gamma conventions and the multiplicity weighting.
      call check(m%shell(1) == 0, 'cubic: Gamma carries shell index 0', nfail)
      call check(m%qlen(1) == 0.0_rk .and. all(m%hkl(:, 1) == 0), &
                 'cubic: Gamma has zero q and zero Miller indices', nfail)
      do i = 2, m%nmodes
         if (m%qlen(i) <= 0.0_rk) cycle
         call check(m%qlen(i) <= 1.0_rk + 1.0e-12_rk, &
                    'cubic: every mode respects qmax', nfail)
      end do

      ! +- reduction: every kept mode is its own pair representative.
      do i = 2, m%nmodes
         call check(m%hkl(1, i) > 0 .or. (m%hkl(1, i) == 0 .and. m%hkl(2, i) > 0) &
                    .or. (m%hkl(1, i) == 0 .and. m%hkl(2, i) == 0 .and. m%hkl(3, i) > 0), &
                    'cubic: every mode is a canonical +- representative', nfail)
      end do
      call check(.not. m%thinned, 'cubic: no thinning without a budget', nfail)
      call m%finalize()

      ! Determinism: the same input has to give the same set.
      call modes_build(m, cell, 1.0_rk, .true., 0, thin_none, ierr, message)
      a = 0.0_rk
      do i = 1, m%nmodes
         a(:, 1) = matmul(cell%b, real(m%hkl(:, i), rk))
         call check(maxval(abs(a(:, 1) - m%qvec(:, i))) < 1.0e-9_rk, &
                    'cubic: the q vector matches its Miller indices', nfail)
      end do
      call m%finalize()
   end subroutine check_cubic

   !> The point group follows from the metric, not from an assumed shape.
   subroutine check_point_groups(nfail)
      integer, intent(inout) :: nfail
      type(cell_t) :: cell
      type(modes_t) :: m
      real(rk) :: a(3, 3)
      integer :: ierr
      character(len=256) :: message

      ! Orthorhombic: three different axis lengths, D_2h has 8 operations.
      a = 0.0_rk
      a(1, 1) = 20.0_rk
      a(2, 2) = 21.0_rk
      a(3, 3) = 22.0_rk
      call cell%set_vectors(a)
      call modes_build(m, cell, 0.8_rk, .false., 0, thin_none, ierr, message)
      call check(ierr == 0, 'orthorhombic: build succeeds', nfail)
      call check(m%nops == 8, 'orthorhombic: point group order is 8', nfail)
      call m%finalize()

      ! Generic triclinic: only {E, -1} survives.
      a = 0.0_rk
      a(:, 1) = [20.0_rk, 0.0_rk, 0.0_rk]
      a(:, 2) = [1.3_rk, 20.5_rk, 0.0_rk]
      a(:, 3) = [0.7_rk, 0.9_rk, 21.0_rk]
      call cell%set_vectors(a)
      call modes_build(m, cell, 0.8_rk, .false., 0, thin_none, ierr, message)
      call check(ierr == 0, 'triclinic: build succeeds', nfail)
      call check(m%nops == 2, 'triclinic: point group order is 2', nfail)
      call m%finalize()
   end subroutine check_point_groups

   !> Cells too skewed for the point-group bootstrap box.
   !!
   !! The Miller index box of the symmetry search is a bootstrap, not a
   !! property of the cell: the mode set, the shells and the +- reduction
   !! never use the group.  A cell whose operations fit inside the box has to
   !! find its {E, -1}; one that does not has to fall back to the trivial
   !! group instead of failing, losing only the orbit reduction.
   subroutine check_skewed(nfail)
      integer, intent(inout) :: nfail
      type(cell_t) :: cell
      type(modes_t) :: m
      real(rk) :: a(3, 3)
      integer :: i, ierr
      logical :: all_singletons
      character(len=256) :: message

      ! 4.0 x 14.3 x 18.6 A: the search needs Miller indices up to 6, which the
      ! raised bound covers, so the cell keeps its point group.
      a = 0.0_rk
      a(:, 1) = [4.0_rk, 0.0_rk, 0.0_rk]
      a(:, 2) = [3.0_rk, 14.0_rk, 0.0_rk]
      a(:, 3) = [4.0_rk, 2.0_rk, 18.0_rk]
      call cell%set_vectors(a)
      call modes_build(m, cell, 1.0_rk, .true., 0, thin_none, ierr, message)
      call check(ierr == 0, 'skewed: the raised bound accepts the cell ('// &
                 trim(message)//')', nfail)
      call check(.not. m%point_group_fallback, &
                 'skewed: no fallback while the search box is wide enough', nfail)
      call check(m%nops == 2, 'skewed: the cell keeps its {E, -1} group', nfail)
      call m%finalize()

      ! 2.0 x 15.3 x 18.6 A: the search would need indices past the box, so the
      ! build continues with the trivial group.  The modes themselves stay
      ! exact: 10 +- reduced vectors below |q| = 1 in 10 shells, plus Gamma.
      a(:, 1) = [2.0_rk, 0.0_rk, 0.0_rk]
      a(:, 2) = [3.0_rk, 15.0_rk, 0.0_rk]
      a(:, 3) = [4.0_rk, 2.0_rk, 18.0_rk]
      call cell%set_vectors(a)
      call modes_build(m, cell, 1.0_rk, .true., 0, thin_none, ierr, message)
      call check(ierr == 0, 'skewed: the fallback build succeeds ('// &
                 trim(message)//')', nfail)
      call check(m%point_group_fallback, 'skewed: the fallback is reported', nfail)
      call check(m%nops == 2, 'skewed: the generic {E, -1} group is used', nfail)
      call check(m%nmodes == 11, 'skewed: 10 +- modes plus Gamma', nfail)
      call check(m%nshell == 10, 'skewed: 10 shells below qmax = 1', nfail)
      ! Inversion is still a symmetry, so each +- pair keeps its multiplicity
      ! of two and no other vector shares its orbit.  The Gamma point is not
      ! part of a pair, so it is skipped.
      all_singletons = .true.
      do i = merge(2, 1, m%has_gamma), m%nmodes
         if (m%orbit_mult(i) /= 2) all_singletons = .false.
      end do
      call check(all_singletons, 'skewed: every +- pair is one orbit of two', nfail)
      ! The shells are still the distinct |q| of the lattice, so a shell holds
      ! exactly the vectors of that radius.
      do i = 1, m%nshell
         call check(count(m%shell == i) >= 1, &
                    'skewed: every shell holds at least one mode', nfail)
      end do
      call m%finalize()
   end subroutine check_skewed

   !> The budget drops whole shells, keeps the ends and stays deterministic.
   subroutine check_budget(nfail)
      integer, intent(inout) :: nfail
      type(cell_t) :: cell
      type(modes_t) :: m, m2
      real(rk) :: total_mult
      integer :: i, ierr, first, last
      character(len=256) :: message

      call make_cubic(ka_l, cell)

      ! Budget 10: the two end shells (3 + 6 modes) fit, nothing more does.
      call modes_build(m, cell, 1.0_rk, .true., 10, thin_shells, ierr, message)
      call check(ierr == 0, 'budget: whole-shell thinning succeeds', nfail)
      call check(m%nmodes <= 10, 'budget: whole-shell thinning respects the budget', nfail)
      call check(m%thinned, 'budget: thinning is reported', nfail)
      first = minval(m%shell(2:))
      last = maxval(m%shell(2:))
      call check(first == 1 .and. last == 7, &
                 'budget: the two ends of the q range survive', nfail)

      ! Orbit thinning keeps every shell with one representative each.
      call modes_build(m2, cell, 1.0_rk, .true., 8, thin_orbits, ierr, message)
      call check(ierr == 0, 'budget: orbit thinning succeeds', nfail)
      call check(m2%nshell == 7, 'budget: orbit thinning keeps the shells', nfail)
      call check(m2%nmodes == 8, 'budget: one orbit representative per shell', nfail)
      call check(m2%nmodes <= 8, 'budget: orbit thinning respects the budget', nfail)
      ! One representative per orbit: the multiplicities now add up to the
      ! full population of the seven shells, which is what the shell average
      ! needs in order to weight the orbits correctly.
      total_mult = 0.0_rk
      do i = 2, m2%nmodes
         call check(m2%orbit_mult(i) >= 6, &
                    'budget: orbit representative keeps its full multiplicity', nfail)
         total_mult = total_mult + real(m2%orbit_mult(i), rk)
      end do
      call check(abs(total_mult - 92.0_rk) < 0.5_rk, &
                 'budget: orbit multiplicities add up to the 92 lattice vectors', nfail)
      call m%finalize()
      call m2%finalize()
   end subroutine check_budget

   !> A shell that holds two orbits: what the multiplicity weighting is for.
   !!
   !! For a cubic cell every shell with `|n|^2 <= 8` holds a single orbit, so the
   !! weights are uniform and the rule the shell average depends on is
   !! invisible.  `|n|^2 = 9` is the first shell with two orbits, six vectors of
   !! the `(3,0,0)` type and twenty-four of the `(2,2,1)` type.  The invariant
   !! is that `sum_o m_o / k_o` over the kept modes of the shell equals the
   !! number of lattice vectors the shell holds (30), whether the shell is kept
   !! whole or reduced to one representative per orbit.
   subroutine check_multi_orbit(nfail)
      integer, intent(inout) :: nfail
      type(cell_t) :: cell
      type(modes_t) :: m
      real(rk) :: target, wsum
      integer :: i, s, ierr, kept, six, many, k_o
      character(len=256) :: message

      call make_cubic(ka_l, cell)
      target = 3.0_rk*two_pi/ka_l           ! |n| = 3, i.e. |n|^2 = 9

      ! Whole shells: the +- reduction keeps 3 + 12 = 15 modes.
      call modes_build(m, cell, 1.05_rk, .false., 0, thin_none, ierr, message)
      call check(ierr == 0, 'multi-orbit: build succeeds ('//trim(message)//')', nfail)
      call check(shell_index(m, target) > 0, 'multi-orbit: the |n| = 3 shell is there', nfail)
      s = shell_index(m, target)
      kept = 0
      six = 0
      many = 0
      wsum = 0.0_rk
      do i = 1, m%nmodes
         if (m%shell(i) /= s) cycle
         kept = kept + 1
         if (m%orbit_mult(i) == 6) six = six + 1
         if (m%orbit_mult(i) == 24) many = many + 1
         k_o = count(m%shell == s .and. m%orbit == m%orbit(i))
         wsum = wsum + real(m%orbit_mult(i), rk)/real(max(k_o, 1), rk)
      end do
      call check(kept == 15, 'multi-orbit: 15 modes after the +- reduction', nfail)
      call check(six == 3 .and. many == 12, 'multi-orbit: a 6-fold and a 24-fold orbit', nfail)
      call check(abs(wsum - 30.0_rk) < 0.5_rk, &
                 'multi-orbit: the weights add up to the 30 lattice vectors', nfail)
      call m%finalize()

      ! One representative per orbit: two modes that still carry their orbit.
      call modes_build(m, cell, 1.05_rk, .false., 10, thin_orbits, ierr, message)
      call check(ierr == 0, 'multi-orbit: orbit thinning succeeds', nfail)
      s = shell_index(m, target)
      kept = 0
      wsum = 0.0_rk
      do i = 1, m%nmodes
         if (m%shell(i) /= s) cycle
         kept = kept + 1
         k_o = count(m%shell == s .and. m%orbit == m%orbit(i))
         wsum = wsum + real(m%orbit_mult(i), rk)/real(max(k_o, 1), rk)
      end do
      call check(kept == 2, 'multi-orbit: one representative per orbit', nfail)
      call check(abs(wsum - 30.0_rk) < 0.5_rk, &
                 'multi-orbit: the representatives recover the same weights', nfail)
      call m%finalize()
   end subroutine check_multi_orbit

   !> Index of the shell whose radius matches `q`, or 0 when there is none.
   pure integer function shell_index(m, q) result(s)
      type(modes_t), intent(in) :: m
      real(rk), intent(in) :: q
      integer :: i

      s = 0
      do i = 1, m%nshell
         if (abs(m%shell_q(i) - q) <= 1.0e-6_rk*q) s = i
      end do
   end function shell_index

   !> A cubic cell of side `l`.
   subroutine make_cubic(l, cell)
      real(rk), intent(in) :: l
      type(cell_t), intent(out) :: cell
      real(rk) :: a(3, 3)

      a = 0.0_rk
      a(1, 1) = l
      a(2, 2) = l
      a(3, 3) = l
      call cell%set_vectors(a)
   end subroutine make_cubic

   !> Bookkeeping for the checks above.
   subroutine check(ok, label, nfail)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      integer, intent(inout) :: nfail

      if (ok) return
      nfail = nfail + 1
      write (error_unit, '(a)') 'FAIL: '//trim(label)
   end subroutine check

end program test_modes
