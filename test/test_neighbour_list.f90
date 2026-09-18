!> Brute force validation of the cell list (distances, images, PBC, skin).
program test_neighbour_list
   use sqc_kinds
   use sqc_cell, only: cell_t, inverse3
   use sqc_neighbour_list, only: neighbour_list_t
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   implicit none

   integer :: nfail

   nfail = 0
   call check_case('orthogonal, periodic', [20.0_rk, 20.0_rk, 20.0_rk], [0.0_rk, 0.0_rk, 0.0_rk], &
                   [.true., .true., .true.], 300, 5.0_rk, 1.0_rk, nfail)
   call check_case('triclinic, periodic', [18.0_rk, 17.0_rk, 21.0_rk], [2.5_rk, -1.5_rk, 3.0_rk], &
                   [.true., .true., .true.], 300, 5.0_rk, 1.0_rk, nfail)
   call check_case('slab, 2D periodic', [18.0_rk, 18.0_rk, 14.0_rk], [0.0_rk, 0.0_rk, 0.0_rk], &
                   [.true., .true., .false.], 300, 4.0_rk, 0.0_rk, nfail)
   call check_case('cluster, no periodicity', [24.0_rk, 24.0_rk, 24.0_rk], [0.0_rk, 0.0_rk, 0.0_rk], &
                   [.false., .false., .false.], 200, 6.0_rk, 0.0_rk, nfail)
   call check_case('small box, self images', [8.0_rk, 8.0_rk, 8.0_rk], [0.0_rk, 0.0_rk, 0.0_rk], &
                   [.true., .true., .true.], 120, 5.0_rk, 0.5_rk, nfail)

   call check_skin_reuse(nfail)

   if (nfail > 0) then
      write (error_unit, '(a,i0,a)') 'FAIL: ', nfail, ' neighbour list checks failed'
      stop 1
   end if
   write (output_unit, '(a)') 'PASS: neighbour list matches brute force'

contains

   !> Build one configuration, compare the pair list with a brute force sum.
   subroutine check_case(label, lengths, tilt, pbc, natoms, cutoff, skin, nfail)
      character(len=*), intent(in) :: label
      real(rk), intent(in) :: lengths(3), tilt(3), cutoff, skin
      logical, intent(in) :: pbc(3)
      integer, intent(in) :: natoms
      integer, intent(inout) :: nfail
      type(cell_t) :: cell
      type(neighbour_list_t) :: nl
      real(rk), allocatable :: pos(:, :)
      real(rk), allocatable :: mine_d(:), ref_d(:)
      integer, allocatable :: mine_i(:), mine_j(:), ref_i(:), ref_j(:)
      integer :: ierr, nmine, nref, i
      character(len=512) :: message

      call make_cell(cell, lengths, tilt)
      call random_positions(pos, natoms, cell, pbc)
      call nl%configure(cell, int(natoms, ik), cutoff, skin, pbc, ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a,a,a)') 'FAIL ', trim(label), ': configure: '//trim(message)
         nfail = nfail + 1
         return
      end if
      call nl%begin_frame(pos, ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a,a,a)') 'FAIL ', trim(label), ': begin_frame: '//trim(message)
         nfail = nfail + 1
         return
      end if

      call list_pairs(nl, cutoff, mine_i, mine_j, mine_d, nmine)
      call brute_force(cell, pos, pbc, cutoff, ref_i, ref_j, ref_d, nref)
      call compare_pairs(label, mine_i, mine_j, mine_d, nmine, ref_i, ref_j, ref_d, nref, nfail)
      write (output_unit, '(a,a,a,i0,a,i0)') '  ok   ', trim(label), ': ', nmine, &
         ' pairs vs brute force ', nref
   end subroutine check_case

   !> Reuse across a smooth trajectory must not lose pairs, and skin=0 must
   !! give exactly the same pairs as a rebuild.
   subroutine check_skin_reuse(nfail)
      integer, intent(inout) :: nfail
      type(cell_t) :: cell
      type(neighbour_list_t) :: nl_skin, nl_fresh
      real(rk), allocatable :: pos(:, :)
      integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
      real(rk), allocatable :: d1(:), d2(:)
      integer :: ierr, n1, n2, frame, n
      character(len=512) :: message

      call make_cell(cell, [20.0_rk, 20.0_rk, 20.0_rk], [0.0_rk, 0.0_rk, 0.0_rk])
      call random_positions(pos, 400, cell, [.true., .true., .true.])
      call nl_skin%configure(cell, 400_ik, 5.0_rk, 1.0_rk, [.true., .true., .true.], &
                             ierr, message)
      if (ierr /= 0) return

      do frame = 1, 20
         call rattle(pos, 0.02_rk)
         ! the list with skin is expected to be reused most of the time
         call nl_skin%begin_frame(pos, ierr, message)
         if (ierr /= 0) then
            write (error_unit, '(a)') 'FAIL skin: '//trim(message)
            nfail = nfail + 1
            return
         end if
      end do

      ! compare the final frame against a list built from scratch (skin = 0)
      call nl_fresh%configure(cell, 400_ik, 5.0_rk, 0.0_rk, [.true., .true., .true.], ierr, message)
      call nl_fresh%begin_frame(pos, ierr, message)
      call list_pairs(nl_skin, 5.0_rk, i1, j1, d1, n1)
      call list_pairs(nl_fresh, 5.0_rk, i2, j2, d2, n2)
      call compare_pairs('skin reuse vs rebuild', i1, j1, d1, n1, i2, j2, d2, n2, nfail)
      n = nl_skin%nbuilds
      if (n >= 20) then
         write (error_unit, '(a,i0,a)') 'FAIL skin: ', n, ' rebuilds in 20 rattled frames'
         nfail = nfail + 1
      else
         write (output_unit, '(a,i0,a)') '  ok   skin reuse: ', n, ' rebuilds over 20 frames'
      end if
   end subroutine check_skin_reuse

   !> Pairs of the list that are within the cutoff, as (i, j, d) triples.
   subroutine list_pairs(nl, cutoff, ii, jj, dd, n)
      type(neighbour_list_t), intent(in) :: nl
      real(rk), intent(in) :: cutoff
      integer, allocatable, intent(out) :: ii(:), jj(:)
      real(rk), allocatable, intent(out) :: dd(:)
      integer, intent(out) :: n
      integer(lk) :: k, s, e
      integer :: i, j
      real(rk) :: dr(3), d2

      allocate (ii(nl%npairs), jj(nl%npairs), dd(nl%npairs))
      n = 0
      do i = 1, nl%natoms
         s = nl%pair_start(i)
         e = nl%pair_start(i + 1) - 1
         do k = s, e
            j = int(nl%pair_atom(k))
            dr = nl%pos(:, j) - nl%pos(:, i) &
                 + matmul(nl%bins, real(nl%pair_shift(:, k), rk))
            d2 = sum(dr*dr)
            if (d2 > cutoff*cutoff) cycle
            n = n + 1
            ii(n) = i
            jj(n) = j
            dd(n) = sqrt(d2)
         end do
      end do
   end subroutine list_pairs

   !> All ordered pairs (including images) closer than the cutoff.
   subroutine brute_force(cell, pos, pbc, cutoff, ii, jj, dd, n)
      type(cell_t), intent(in) :: cell
      real(rk), intent(in) :: pos(:, :)
      logical, intent(in) :: pbc(3)
      real(rk), intent(in) :: cutoff
      integer, allocatable, intent(out) :: ii(:), jj(:)
      real(rk), allocatable, intent(out) :: dd(:)
      integer, intent(out) :: n
      real(rk) :: widths(3), dr(3), d2
      integer :: natoms, i, j, a, b, c, na, nb, nc, pass
      integer, allocatable :: t_i(:), t_j(:)
      real(rk), allocatable :: t_d(:)

      natoms = size(pos, 2)
      widths = cell%widths()
      na = merge(ceiling(cutoff/widths(1)) + 1, 0, pbc(1))
      nb = merge(ceiling(cutoff/widths(2)) + 1, 0, pbc(2))
      nc = merge(ceiling(cutoff/widths(3)) + 1, 0, pbc(3))
      ! pass 1 counts, pass 2 fills (avoids guessing a buffer size)
      do pass = 1, 2
         if (pass == 2) then
            allocate (t_i(n), t_j(n), t_d(n))
            n = 0
         else
            n = 0
         end if
         do i = 1, natoms
            do j = 1, natoms
               do a = -na, na
                  if (.not. pbc(1) .and. a /= 0) cycle
                  do b = -nb, nb
                     if (.not. pbc(2) .and. b /= 0) cycle
                     do c = -nc, nc
                        if (.not. pbc(3) .and. c /= 0) cycle
                        if (i == j .and. a == 0 .and. b == 0 .and. c == 0) cycle
                        dr = pos(:, j) + matmul(cell%a, real([a, b, c], rk)) - pos(:, i)
                        d2 = sum(dr*dr)
                        if (d2 > cutoff*cutoff) cycle
                        n = n + 1
                        if (pass == 2) then
                           t_i(n) = i
                           t_j(n) = j
                           t_d(n) = sqrt(d2)
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
      allocate (ii(n), jj(n), dd(n))
      ii = t_i(1:n)
      jj = t_j(1:n)
      dd = t_d(1:n)
      deallocate (t_i, t_j, t_d)
   end subroutine brute_force

   !> Compare two pair lists as multisets of (i, j, distance).
   subroutine compare_pairs(label, i1, j1, d1, n1, i2, j2, d2, n2, nfail)
      character(len=*), intent(in) :: label
      integer, intent(in) :: i1(:), j1(:), n1, i2(:), j2(:), n2
      real(rk), intent(in) :: d1(:), d2(:)
      integer, intent(inout) :: nfail
      integer, allocatable :: a(:), b(:), c(:), x(:), y(:), z(:)
      integer :: k

      if (n1 /= n2) then
         write (error_unit, '(a,a,a,i0,a,i0)') 'FAIL ', trim(label), ': ', n1, &
            ' pairs in the list, ', n2, ' in the reference'
         nfail = nfail + 1
         return
      end if
      allocate (a(n1), b(n1), c(n1), x(n2), y(n2), z(n2))
      a = i1(1:n1)
      b = j1(1:n1)
      c = nint(d1(1:n1)*10000.0_rk)
      x = i2(1:n2)
      y = j2(1:n2)
      z = nint(d2(1:n2)*10000.0_rk)
      call sort_triples(a, b, c, n1)
      call sort_triples(x, y, z, n2)
      do k = 1, n1
         if (a(k) /= x(k) .or. b(k) /= y(k) .or. abs(c(k) - z(k)) > 1) then
            write (error_unit, '(a,a,a,i0,a,i0,a,i3,a,i5,a,i5)') 'FAIL ', trim(label), &
               ': pair ', k, ' differs: ', a(k), '-', b(k), ' d*1e4=', c(k), ' vs ', z(k)
            nfail = nfail + 1
            exit
         end if
      end do
      deallocate (a, b, c, x, y, z)
   end subroutine compare_pairs

   !> Insertion sort of (i, j, d) triples on i, then j, then d.
   subroutine sort_triples(i, j, d, n)
      integer, intent(inout) :: i(:), j(:), d(:)
      integer, intent(in) :: n
      integer :: k, m, ti, tj, td
      do k = 2, n
         ti = i(k)
         tj = j(k)
         td = d(k)
         m = k - 1
         do while (m >= 1)
            if (i(m) > ti .or. (i(m) == ti .and. (j(m) > tj .or. &
                (j(m) == tj .and. d(m) > td)))) then
               i(m + 1) = i(m)
               j(m + 1) = j(m)
               d(m + 1) = d(m)
               m = m - 1
            else
               exit
            end if
         end do
         i(m + 1) = ti
         j(m + 1) = tj
         d(m + 1) = td
      end do
   end subroutine sort_triples

   !> Cell with the given lattice vector lengths and tilt factors.
   subroutine make_cell(cell, lengths, tilt)
      type(cell_t), intent(out) :: cell
      real(rk), intent(in) :: lengths(3), tilt(3)
      real(rk) :: lo(3), hi(3)

      lo = 0.0_rk
      hi = 0.0_rk
      call cell%init_from_bounds(lo, hi, tilt)
      call cell%set_vectors(reshape([ &
         lengths(1), 0.0_rk, 0.0_rk, &
         tilt(1), lengths(2), 0.0_rk, &
         tilt(2), tilt(3), lengths(3)], [3, 3]))
   end subroutine make_cell

   !> Random positions inside the cell (wrapped).
   subroutine random_positions(pos, natoms, cell, pbc)
      real(rk), allocatable, intent(out) :: pos(:, :)
      integer, intent(in) :: natoms
      type(cell_t), intent(in) :: cell
      logical, intent(in) :: pbc(3)
      real(rk) :: s(3)
      integer :: i

      allocate (pos(3, natoms))
      call random_seed()
      do i = 1, natoms
         call random_number(s)
         if (.not. all(pbc)) s = 0.05_rk + 0.9_rk*s
         pos(:, i) = matmul(cell%a, s) + cell%origin
      end do
   end subroutine random_positions

   !> Small random displacement (much smaller than the skin).
   subroutine rattle(pos, amplitude)
      real(rk), intent(inout) :: pos(:, :)
      real(rk), intent(in) :: amplitude
      real(rk) :: r(3)
      integer :: i

      do i = 1, size(pos, 2)
         call random_number(r)
         pos(:, i) = pos(:, i) + amplitude*(r - 0.5_rk)
      end do
   end subroutine rattle

end program test_neighbour_list
