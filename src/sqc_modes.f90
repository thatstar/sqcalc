! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Reciprocal-lattice mode sets for the grid q sampling (`--q grid`).
!!
!! A grid run samples the reciprocal lattice of the (constant) simulation cell
!! instead of an arbitrary line or sphere,
!!
!!   q = n_1 b_1 + n_2 b_2 + n_3 b_3,     0 < |q| <= qmax,
!!
!! with the Gamma point kept separately.  Only a reciprocal-lattice vector
!! carries a density amplitude that does not depend on how the periodic images
!! are chosen, so a grid avoids the box-form-factor contamination that makes
!! the off-lattice Lebedev shell unusable for `S_4`.
!!
!! Three reductions are available, in increasing order of what they cost:
!!
!! * **+-**: `w_i` is real, hence `W(-q) = conjg(W(q))` and `S_4(-q) =
!!   S_4(q)` identically.  Keeping one vector of every pair therefore halves
!!   the modes and loses nothing at all.
!! * **orbits**: the point group of the periodic lattice,
!!
!!     Aut = { M in GL(3,Z) : M^T G M = G },   G_ij = b_i . b_j,
!!
!!   maps every shell onto itself, and vectors of one orbit have the same
!!   expectation for an isotropic system.  One representative per orbit
!!   therefore estimates the shell average without bias.  The orbit
!!   multiplicity is kept, so a shell made of several orbits can be averaged
!!   with the correct weights (for example `|n|^2 = 9` mixes a 6-fold
!!   `(3,0,0)` orbit with a 24-fold `(2,2,1)` orbit).
!! * **shells**: dropping a whole shell leaves the other shells complete, so
!!   their values and their variances do not change; only the number of fit
!!   points and the covered q range change.
!!
!! The mode budget is spent in that order: whole shells first, because they
!! are free in variance, and orbit thinning last, because it is not.  A shell
!! is the smallest unit with an unambiguous isotropic meaning, so a kept shell
!! is never sampled only in part unless the caller asks for orbit thinning.
!!
!! Both the orbit reduction and the shell average assume that the system is
!! isotropic and in equilibrium; the caller has to check that assumption (the
!! dynamic summary carries a cheap three-direction probe for it).
module sqc_modes
   use sqc_kinds
   use sqc_cell, only: cell_t, two_pi
   implicit none
   private

   public :: modes_t, modes_build, modes_point_group, &
             thin_none, thin_shells, thin_orbits

   !> Keep every +- reduced mode (no thinning beyond `+-`).
   integer, parameter :: thin_none = 0
   !> Drop whole shells before reducing the orbits.
   integer, parameter :: thin_shells = 1
   !> Reduce every orbit to one representative before dropping shells.
   integer, parameter :: thin_orbits = 2

   !> Point group operations we can hold (the cubic group O_h has 48).
   integer, parameter :: max_ops = 48
   !> Largest |n_i| explored when the lattice point group is looked for.  The
   !! bootstrap box has to hold every lattice vector as long as the longest
   !! basis vector, which grows with the skew of the cell: 4 covers
   !! |a_i| <= 3 d_min (d_min the smallest interplanar spacing), 8 covers
   !! about |a_i| <= 7 d_min.  Beyond that modes_build falls back to the
   !! trivial group instead of refusing the cell.
   integer, parameter :: max_index = 8
   !> Size of the Miller index box that is scanned.
   integer, parameter :: max_cand = (2*max_index + 1)**3
   !> Relative tolerance that puts two lattice vectors into one shell.
   real(rk), parameter :: shell_tol = 1.0e-6_rk
   !> Relative tolerance of the metric test M^T G M = G.
   real(rk), parameter :: gram_tol = 1.0e-7_rk
   !> modes_point_group error code: the bootstrap box is too small for this
   !! cell, which the caller may answer with the trivial group.
   integer, parameter :: point_group_unbounded = 2
   !> Guard on the enumerated grid, in the spirit of the static path.
   integer(lk), parameter :: max_grid_points = 400000000_lk

   !> Reciprocal-lattice sampling set of one dynamic run.
   type :: modes_t
      !> Number of modes kept, including Gamma when it was requested.
      integer :: nmodes = 0
      !> Number of shells; shell 0 marks the Gamma point.
      integer :: nshell = 0
      !> How many of those shells survived the mode budget.
      integer :: nshell_kept = 0
      !> Order of the lattice point group.
      integer :: nops = 0
      !> True when the point group could not be found within `max_index` and
      !! the trivial group was used instead: the modes are unaffected, but
      !! orbit thinning cannot reduce anything and the caller should say so.
      logical :: point_group_fallback = .false.
      !> Requested mode budget (0 = unlimited).
      integer :: budget = 0
      !> Thinning policy that was applied.
      integer :: thin_kind = thin_none
      !> True when the Gamma point is part of the set.
      logical :: has_gamma = .false.
      !> True when modes were dropped beyond the +- reduction.
      logical :: thinned = .false.
      !> False when the two end shells alone need more modes than the budget:
      !! the set is then kept anyway, and the caller should say so.
      logical :: budget_met = .true.
      !> Upper bound of |q|.
      real(rk) :: qmax = 0.0_rk
      !> Cartesian q vector and length of every kept mode.
      real(rk), allocatable :: qvec(:, :)
      real(rk), allocatable :: qlen(:)
      !> Representative |q| of every shell; entry 0 belongs to the Gamma point.
      real(rk), allocatable :: shell_q(:)
      !> Miller indices of every kept mode.
      integer, allocatable :: hkl(:, :)
      !> Shell index (0 = Gamma) and orbit multiplicity of every mode.
      integer, allocatable :: shell(:), orbit_mult(:)
      !> Orbit number of every mode, so a caller can tell which modes share an
      !! orbit when it has to weight a shell average.
      integer, allocatable :: orbit(:)
      !> Integer point group operations of the cell.  The columns of `ops(:,:,k)`
      !! are the images of b_1, b_2, b_3, and the index space action is n -> M n.
      integer, allocatable :: ops(:, :, :)
   contains
      procedure :: finalize => modes_finalize
   end type modes_t

contains

   !> Release the storage of a mode set.
   subroutine modes_finalize(self)
      class(modes_t), intent(inout) :: self
      if (allocated(self%qvec)) deallocate (self%qvec)
      if (allocated(self%qlen)) deallocate (self%qlen)
      if (allocated(self%shell_q)) deallocate (self%shell_q)
      if (allocated(self%hkl)) deallocate (self%hkl)
      if (allocated(self%shell)) deallocate (self%shell)
      if (allocated(self%orbit)) deallocate (self%orbit)
      if (allocated(self%orbit_mult)) deallocate (self%orbit_mult)
      if (allocated(self%ops)) deallocate (self%ops)
      self%nmodes = 0
      self%nshell = 0
      self%nops = 0
      self%nshell_kept = 0
      self%point_group_fallback = .false.
      self%thinned = .false.
      self%budget_met = .true.
   end subroutine modes_finalize

   !> Integer automorphisms of the reciprocal metric: the lattice point group.
   !!
   !! Every column of an operation is the image of one basis vector, so the
   !! group is found by looking for triples (v1, v2, v3) of lattice vectors
   !! with v_i . v_j = G_ij.  Only integer index vectors are ever combined, so
   !! the operations are integer by construction; the metric test is the only
   !! floating point comparison involved.
   !!
   !! The Miller index box that is scanned is a bootstrap, not a property of
   !! the cell: a very skewed cell needs large indices before the search is
   !! exhaustive.  When that happens the routine returns
   !! `point_group_unbounded` so the caller can continue with the trivial
   !! group, which costs only the orbit reduction.
   subroutine modes_point_group(cell, ops, nops, ierr, message)
      type(cell_t), intent(in) :: cell
      integer, intent(out) :: ops(3, 3, max_ops)
      integer, intent(out) :: nops
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: g(3, 3), diag(3), blen, tol
      real(rk) :: vi(3), vj(3), vk(3)
      ! A vector whose length matches several basis vectors is listed once per
      ! match, so the lists can hold up to three times the scanned box.
      integer :: cand(3, 3*max_cand), clist(3, 3*max_cand), ncand, nlist(3)
      integer :: bound(3), n(3), i, j, d, i1, i2, i3

      ierr = 0
      message = ''
      nops = 0
      ops = 0
      do i = 1, 3
         do j = 1, 3
            g(i, j) = sum(cell%b(:, i)*cell%b(:, j))
         end do
         diag(i) = g(i, i)
      end do
      if (any(diag <= 0.0_rk)) then
         ierr = 1
         message = 'degenerate cell: the reciprocal basis is not usable'
         return
      end if

      ! n_i = q . a_i / (2 pi) bounds the Miller indices of a vector whose
      ! length does not exceed the longest basis vector.
      blen = sqrt(maxval(diag))
      do i = 1, 3
         bound(i) = int(ceiling(blen*sqrt(sum(cell%a(:, i)**2))/two_pi)) + 1
         if (bound(i) > max_index) then
            ierr = point_group_unbounded
            write (message, '(a,i0,a)') 'the cell needs Miller indices up to ', &
               bound(i), ' to find its point group; the cell is too skewed'
            return
         end if
      end do

      ! Candidate columns, collected per basis vector by matching lengths.
      ncand = 0
      nlist = 0
      do i3 = -bound(3), bound(3)
         do i2 = -bound(2), bound(2)
            do i1 = -bound(1), bound(1)
               n = [i1, i2, i3]
               if (all(n == 0)) cycle
               vi = matmul(cell%b, real(n, rk))
               do d = 1, 3
                  if (abs(sum(vi*vi) - diag(d)) > gram_tol*diag(d)) cycle
                  ncand = ncand + 1
                  if (ncand > 3*max_cand) then
                     ierr = 1
                     message = 'internal: too many candidate vectors for the point group'
                     return
                  end if
                  cand(:, ncand) = n
                  nlist(d) = nlist(d) + 1
                  clist(d, nlist(d)) = ncand
               end do
            end do
         end do
      end do

      ! Triple loop over the columns.  The lists are short (six entries for a
      ! cubic cell), so the triple loop stays cheap.
      do i1 = 1, nlist(1)
         vi = matmul(cell%b, real(cand(:, clist(1, i1)), rk))
         do i2 = 1, nlist(2)
            vj = matmul(cell%b, real(cand(:, clist(2, i2)), rk))
            tol = gram_tol*sqrt(diag(1)*diag(2))
            if (abs(sum(vi*vj) - g(1, 2)) > tol) cycle
            do i3 = 1, nlist(3)
               vk = matmul(cell%b, real(cand(:, clist(3, i3)), rk))
               tol = gram_tol*sqrt(diag(1)*diag(3))
               if (abs(sum(vi*vk) - g(1, 3)) > tol) cycle
               tol = gram_tol*sqrt(diag(2)*diag(3))
               if (abs(sum(vj*vk) - g(2, 3)) > tol) cycle
               if (nops >= max_ops) then
                  ierr = 1
                  message = 'internal: more point group operations than expected'
                  return
               end if
               nops = nops + 1
               ops(:, 1, nops) = cand(:, clist(1, i1))
               ops(:, 2, nops) = cand(:, clist(2, i2))
               ops(:, 3, nops) = cand(:, clist(3, i3))
            end do
         end do
      end do
      if (nops < 2) then
         ierr = 1
         message = 'the lattice point group came out too small; check the cell'
         return
      end if
   end subroutine modes_point_group

   !> Build the mode set of one cell.
   !!
   !! `budget` is the largest number of modes the caller wants to pay for
   !! (0 = unlimited).  When the +- reduced set exceeds it, whole shells are
   !! dropped first (`thin_shells`) or every orbit is reduced to a single
   !! representative first (`thin_orbits`).  The retained shells are spread
   !! uniformly over the q range and always include its two ends, so the lever
   !! arm of a fit survives the thinning.  The Gamma point is never dropped.
   !!
   !! The budget is a target, not a hard limit: if even the two end shells do
   !! not fit (together with Gamma), they are kept anyway and `nmodes` comes
   !! out larger than `budget`.  The caller is expected to compare the two and
   !! say so in its summary.
   !!
   !! A cell too skewed for the point-group search is not an error here: the
   !! build continues with the trivial group and sets
   !! `point_group_fallback`, so only the orbit reduction is lost.
   subroutine modes_build(self, cell, qmax, include_gamma, budget, thin_kind, ierr, message)
      class(modes_t), intent(inout) :: self
      type(cell_t), intent(in) :: cell
      real(rk), intent(in) :: qmax
      logical, intent(in) :: include_gamma
      integer, intent(in) :: budget, thin_kind
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: qlen(:), qvec(:, :), shell_rep(:), per_shell(:)
      integer, allocatable :: hkl(:, :), order(:), shellof(:), orbitof(:), &
                              orbsize(:)
      logical, allocatable :: keep(:), shell_keep(:)
      integer :: nfound, nshell, nkept, i, j, k, nops_local, npick
      logical :: use_orbits

      ierr = 0
      message = ''
      if (qmax <= 0.0_rk) then
         ierr = 1
         message = 'the grid qmax must be positive'
         return
      end if
      call self%finalize()
      self%qmax = qmax
      self%budget = max(budget, 0)
      self%thin_kind = thin_kind
      self%has_gamma = include_gamma

      allocate (self%ops(3, 3, max_ops))
      call modes_point_group(cell, self%ops, nops_local, ierr, message)
      if (ierr == point_group_unbounded) then
         ! The search box is a bootstrap for finding the symmetry operations,
         ! not a requirement of the sampling: the mode set, the shells and the
         ! +- reduction do not use the group at all.  Fall back to {E, -1}, the
         ! group of a generic triclinic cell: inversion is always a symmetry of
         ! the reciprocal lattice, so a +- pair keeps its multiplicity of two,
         ! and every other orbit is a single vector.  The caller reports that
         ! orbit thinning is off.
         ierr = 0
         message = ''
         nops_local = 2
         self%ops = 0
         self%ops(1, 1, 1) = 1
         self%ops(2, 2, 1) = 1
         self%ops(3, 3, 1) = 1
         self%ops(1, 1, 2) = -1
         self%ops(2, 2, 2) = -1
         self%ops(3, 3, 2) = -1
         self%point_group_fallback = .true.
      else if (ierr /= 0) then
         return
      end if
      self%nops = nops_local

      call modes_grid_count(cell, qmax, nfound, ierr, message)
      if (ierr /= 0) return
      allocate (hkl(3, nfound), qvec(3, nfound), qlen(nfound))
      call modes_fill(cell, qmax, nfound, hkl, qvec, qlen)
      allocate (order(nfound))
      call modes_order(hkl, qlen, nfound, order)

      allocate (shellof(nfound), shell_rep(nfound))
      call modes_group_shells(qlen, order, nfound, shellof, shell_rep, nshell)

      allocate (orbitof(nfound), orbsize(nfound))
      call modes_orbits(hkl, order, shellof, nfound, nshell, self%ops, &
                        nops_local, orbitof, orbsize)

      ! Per shell: how many modes remain after +-, and how many orbits exist.
      allocate (per_shell(nshell), shell_keep(nshell))
      per_shell = 0.0_rk
      do i = 1, nfound
         if (modes_canonical(hkl(:, i))) per_shell(shellof(i)) = per_shell(shellof(i)) + 1.0_rk
      end do
      shell_keep = .true.
      nkept = int(sum(per_shell)) + merge(1, 0, include_gamma)

      use_orbits = .false.
      if (self%budget > 0 .and. nkept > self%budget) then
         self%thinned = .true.
         if (thin_kind == thin_orbits) then
           use_orbits = .true.
           call modes_shell_orbits(orbitof, shellof, nfound, per_shell)
           call modes_fit_shells(per_shell, nshell, self%budget, shell_keep, npick, include_gamma)
         else
           call modes_fit_shells(per_shell, nshell, self%budget, shell_keep, npick, include_gamma)
           if (npick < 2) then
              ! whole shells do not fit even as a pair: fall back to one
              ! representative per orbit and try again
              use_orbits = .true.
              call modes_shell_orbits(orbitof, shellof, nfound, per_shell)
              call modes_fit_shells(per_shell, nshell, self%budget, shell_keep, npick, include_gamma)
           end if
         end if
         if (npick == 0) then
            ! No selection meets the budget at all.  Keep the two end shells,
            ! the smallest set that still spans the q range, and let the caller
            ! report that the budget was not met.
            shell_keep = .false.
            shell_keep(1) = .true.
            if (nshell > 1) shell_keep(nshell) = .true.
         end if
      end if

      ! Assemble the kept modes: Gamma first, then the shells in q order.
      allocate (keep(nfound))
      keep = .false.
      do i = 1, nfound
         if (.not. shell_keep(shellof(i))) cycle
         if (.not. modes_canonical(hkl(:, i))) cycle
         if (use_orbits .and. .not. modes_is_orbit_head(hkl, orbitof, i)) cycle
         keep(i) = .true.
      end do
      self%nmodes = count(keep) + merge(1, 0, include_gamma)
      self%nshell = nshell
      self%nshell_kept = count(shell_keep)
      self%budget_met = self%budget <= 0 .or. self%nmodes <= self%budget

      allocate (self%qvec(3, self%nmodes), self%qlen(self%nmodes), &
                self%hkl(3, self%nmodes), self%shell(self%nmodes), &
                self%orbit(self%nmodes), self%orbit_mult(self%nmodes), &
                self%shell_q(0:nshell))
      self%shell_q = 0.0_rk
      if (nshell > 0) self%shell_q(1:nshell) = shell_rep(1:nshell)
      k = 0
      if (include_gamma) then
         k = 1
         self%qvec(:, k) = 0.0_rk
         self%qlen(k) = 0.0_rk
         self%hkl(:, k) = 0
         self%shell(k) = 0
         self%orbit(k) = 0
         self%orbit_mult(k) = 1
      end if
      do j = 1, nfound
         i = order(j)
         if (.not. keep(i)) cycle
         k = k + 1
         self%hkl(:, k) = hkl(:, i)
         self%qvec(:, k) = qvec(:, i)
         self%qlen(k) = qlen(i)
         self%shell(k) = shellof(i)
         self%orbit(k) = orbitof(i)
         self%orbit_mult(k) = orbsize(i)
      end do
      deallocate (qlen, qvec, shell_rep, per_shell, hkl, order, shellof, &
                  orbitof, orbsize, keep, shell_keep)
   end subroutine modes_build

   !> How many reciprocal-lattice vectors lie inside the q range.
   !!
   !! The Grid is enumerated in Miller index space, so the only limit needed
   !! is the index box `2*ceil(qmax*|a_i|/2 pi) + 1` per direction.
   subroutine modes_grid_count(cell, qmax, nfound, ierr, message)
      type(cell_t), intent(in) :: cell
      real(rk), intent(in) :: qmax
      integer, intent(out) :: nfound
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer :: bound(3), i
      integer(lk) :: gridpoints

      ierr = 0
      message = ''
      do i = 1, 3
         bound(i) = int(ceiling(qmax*sqrt(sum(cell%a(:, i)**2))/two_pi)) + 1
      end do
      gridpoints = int(2*bound(1) + 1, lk)*int(2*bound(2) + 1, lk) &
                   *int(2*bound(3) + 1, lk)
      if (gridpoints > max_grid_points) then
         ierr = 1
         write (message, '(a,i0,a,i0,a)') 'the grid needs ', gridpoints, &
            ' points, above the limit of ', max_grid_points, '; reduce qmax'
         return
      end if
      nfound = modes_count(cell, qmax, bound)
      if (nfound == 0) then
         ierr = 1
         message = 'no reciprocal lattice points fall inside qmax'
         return
      end if
   end subroutine modes_grid_count

   !> How many reciprocal-lattice vectors lie inside the q range.
   pure integer function modes_count(cell, qmax, bound) result(n)
      type(cell_t), intent(in) :: cell
      real(rk), intent(in) :: qmax
      integer, intent(in) :: bound(3)
      integer :: i1, i2, i3
      real(rk) :: q(3)

      n = 0
      do i3 = -bound(3), bound(3)
         do i2 = -bound(2), bound(2)
            do i1 = -bound(1), bound(1)
               if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
               q = matmul(cell%b, real([i1, i2, i3], rk))
               if (sum(q*q) <= qmax*qmax) n = n + 1
            end do
         end do
      end do
   end function modes_count

   !> Write the Miller indices and q vector of every vector inside the range.
   subroutine modes_fill(cell, qmax, nfound, hkl, qvec, qlen)
      type(cell_t), intent(in) :: cell
      real(rk), intent(in) :: qmax
      integer, intent(in) :: nfound
      integer, intent(inout) :: hkl(:, :)
      real(rk), intent(inout) :: qvec(:, :), qlen(:)
      integer :: bound(3), i1, i2, i3, n, i
      real(rk) :: q(3), ql

      do i = 1, 3
         bound(i) = int(ceiling(qmax*sqrt(sum(cell%a(:, i)**2))/two_pi)) + 1
      end do
      n = 0
      do i3 = -bound(3), bound(3)
         do i2 = -bound(2), bound(2)
            do i1 = -bound(1), bound(1)
               if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
               q = matmul(cell%b, real([i1, i2, i3], rk))
               ql = sqrt(sum(q*q))
               if (ql > qmax) cycle
               n = n + 1
               hkl(:, n) = [i1, i2, i3]
               qvec(:, n) = q
               qlen(n) = ql
            end do
         end do
      end do
   end subroutine modes_fill

   !> True when mode `a` sorts before mode `b`: by |q|, then by Miller index.
   pure logical function modes_before(hkl, qlen, a, b) result(before)
      integer, intent(in) :: hkl(:, :), a, b
      real(rk), intent(in) :: qlen(:)
      integer :: k
      real(rk) :: scale

      scale = max(qlen(a), qlen(b))
      if (abs(qlen(a) - qlen(b)) > shell_tol*scale) then
         before = qlen(a) < qlen(b)
         return
      end if
      do k = 1, 3
         if (hkl(k, a) /= hkl(k, b)) then
            before = hkl(k, a) < hkl(k, b)
            return
         end if
      end do
      before = .false.
   end function modes_before

   !> Order the modes by |q| with a Shell sort (Fortran has no intrinsic sort).
   subroutine modes_order(hkl, qlen, n, order)
      integer, intent(in) :: hkl(:, :), n
      real(rk), intent(in) :: qlen(:)
      integer, intent(out) :: order(:)
      integer :: i, j, gap, tmp

      do i = 1, n
         order(i) = i
      end do
      gap = n/2
      do while (gap > 0)
         do i = gap + 1, n
            tmp = order(i)
            j = i
            do while (j > gap)
               if (.not. modes_before(hkl, qlen, tmp, order(j - gap))) exit
               order(j) = order(j - gap)
               j = j - gap
            end do
            order(j) = tmp
         end do
         gap = gap/2
      end do
   end subroutine modes_order

   !> Group the ordered modes into shells of equal |q|.
   subroutine modes_group_shells(qlen, order, n, shellof, shell_rep, nshell)
      real(rk), intent(in) :: qlen(:)
      integer, intent(in) :: order(:), n
      integer, intent(out) :: shellof(:)
      real(rk), intent(out) :: shell_rep(:)
      integer, intent(out) :: nshell
      integer :: i

      nshell = 0
      shell_rep = 0.0_rk
      shellof = 0
      do i = 1, n
         if (nshell == 0) then
            nshell = 1
            shell_rep(1) = qlen(order(i))
         else if (abs(qlen(order(i)) - shell_rep(nshell)) > shell_tol*shell_rep(nshell)) then
            nshell = nshell + 1
            shell_rep(nshell) = qlen(order(i))
         end if
         shellof(order(i)) = nshell
      end do
   end subroutine modes_group_shells

   !> True for the member of a +- pair that the code keeps.
   pure logical function modes_canonical(n) result(keep)
      integer, intent(in) :: n(3)
      integer :: k

      do k = 1, 3
         if (n(k) /= 0) then
            keep = n(k) > 0
            return
         end if
      end do
      keep = .false.
   end function modes_canonical

   !> True when `i` is the first canonical member of its orbit.
   pure logical function modes_is_orbit_head(hkl, orbitof, i) result(head)
      integer, intent(in) :: hkl(:, :), orbitof(:), i
      integer :: j

      do j = 1, i - 1
         if (orbitof(j) /= orbitof(i)) cycle
         if (.not. modes_canonical(hkl(:, j))) cycle
         head = .false.
         return
      end do
      head = .true.
   end function modes_is_orbit_head

   !> Partition every shell into orbits of the lattice point group.
   !!
   !! `orbitof` is an orbit number per mode (numbered in order of first
   !! appearance) and `orbsize` its multiplicity, i.e. the number of lattice
   !! vectors it contains.  The Gamma point is not part of the enumeration.
   subroutine modes_orbits(hkl, order, shellof, n, nshell, ops, nops, orbitof, orbsize)
      integer, intent(in) :: hkl(:, :), order(:), shellof(:), n, nshell, nops
      integer, intent(in) :: ops(3, 3, max_ops)
      integer, intent(out) :: orbitof(:), orbsize(:)
      integer :: i, j, op, s, nimg(3), norb

      orbitof = 0
      orbsize = 0
      norb = 0
      do s = 1, nshell
         do i = 1, n
            if (shellof(order(i)) /= s) cycle
            if (orbitof(order(i)) /= 0) cycle
            norb = norb + 1
            orbitof(order(i)) = norb
            do op = 1, nops
               nimg = matmul(ops(:, :, op), hkl(:, order(i)))
               do j = 1, n
                  if (shellof(order(j)) /= s) cycle
                  if (all(hkl(:, order(j)) == nimg)) then
                     orbitof(order(j)) = norb
                     exit
                  end if
               end do
            end do
         end do
      end do
      do i = 1, n
         if (orbitof(i) <= 0) cycle
         orbsize(i) = count(orbitof(1:n) == orbitof(i))
      end do
   end subroutine modes_orbits

   !> Replace the per-shell counts with the number of orbits per shell.
   subroutine modes_shell_orbits(orbitof, shellof, n, per_shell)
      integer, intent(in) :: orbitof(:), shellof(:), n
      real(rk), intent(inout) :: per_shell(:)
      integer :: i

      per_shell = 0.0_rk
      do i = 1, n
         if (orbitof(i) <= 0) cycle
         if (any(orbitof(1:i - 1) == orbitof(i))) cycle
         per_shell(shellof(i)) = per_shell(shellof(i)) + 1.0_rk
      end do
   end subroutine modes_shell_orbits

   !> Pick the largest set of shells, spread uniformly over the list and
   !! including its two ends, whose mode count fits the budget.
   !!
   !! `npick` comes back as 0 when no non-empty selection fits, which lets the
   !! caller walk down its ladder of reductions (whole shells, then one
   !! representative per orbit, then the two end shells as the floor) instead
   !! of giving up on the first rung.
   subroutine modes_fit_shells(per_shell, nshell, budget, shell_keep, npick, include_gamma)
      real(rk), intent(in) :: per_shell(:)
      integer, intent(in) :: nshell, budget
      logical, intent(in) :: include_gamma
      logical, intent(out) :: shell_keep(:)
      integer, intent(out) :: npick
      integer, allocatable :: pick(:)
      integer :: k, j, nout
      real(rk) :: total
      integer :: gamma

      gamma = merge(1, 0, include_gamma)
      allocate (pick(nshell))
      shell_keep = .false.
      npick = 0
      do k = nshell, 1, -1
         call modes_uniform_pick(nshell, k, pick, nout)
         total = real(gamma, rk)
         do j = 1, nout
            total = total + per_shell(pick(j))
         end do
         if (budget > 0 .and. total > real(budget, rk)) cycle
         npick = nout
         exit
      end do
      if (npick == 0) then
         ! nothing fits: the caller tries the next reduction of the ladder
      else
         do j = 1, npick
            shell_keep(pick(j)) = .true.
         end do
      end if
      deallocate (pick)
   end subroutine modes_fit_shells

   !> Uniformly spaced shell indices from 1 to `nshell`, both ends included.
   subroutine modes_uniform_pick(nshell, k, pick, nout)
      integer, intent(in) :: nshell, k
      integer, intent(out) :: pick(:)
      integer, intent(out) :: nout
      integer :: j, idx, prev

      nout = 0
      prev = 0
      do j = 1, k
         if (k <= 1) then
            idx = 1
         else
            idx = int(nint(1.0_rk + real(j - 1, rk)*real(nshell - 1, rk)/real(k - 1, rk)))
         end if
         if (idx == prev) cycle
         nout = nout + 1
         pick(nout) = idx
         prev = idx
      end do
   end subroutine modes_uniform_pick

end module sqc_modes
