! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Dynamic structure factor S(q,w) along a line in reciprocal space.
!!
!! The method follows the recipe used by dynasor and MDANSE: evaluate the
!! density amplitude rho_a(q,t) of every species at a few q points, correlate
!! it in time and Fourier transform the correlation.  Nothing here uses the
!! reciprocal grid:
!!
!!   * the q line is arbitrary (off lattice q included), so the amplitudes are
!!     summed directly, O(N * nq) per frame - for the ~100 q points of a line
!!     scan that is far cheaper than transforming the whole grid;
!!   * the correlation is accumulated in a ring buffer, so the memory does not
!!     depend on the length of the trajectory.
!!
!! Definitions (W(q) is the same normalization the static table uses):
!!
!!   C_ab(q,l) = <rho_a(q,t+l) rho_b*(q,t)>            l = 0..maxframes
!!   F(q,l)    = sum_ab w_a(q) w_b(q) C_ab(q,l) / W(q)
!!   S(q,w_n)  = (dt/2pi) sum_l F(q,l) exp(i w_n l dt)      (even extension)
!!
!! The sum over time origins carries its own count per lag (`C_cnt`), which is
!! the "nanmean" bookkeeping: a trajectory whose length is not a multiple of
!! the window needs no padding and its ragged tail cannot bias a lag.  With
!! the one-sided folding the zeroth moment is exact,
!!
!!   sum_n S(q,w_n) dw = F(q,0) = S(q) = the OUTPUT table of the same run,
!!
!! which is the primary consistency test of the method.
module sqc_dynamics
   use sqc_kinds
   use sqc_output, only: format_text, format_hdf5
   use sqc_cell, only: cell_t
   use sqc_dump, only: frame_t
   use sqc_weights, only: weight_scheme_t
   use sqc_structure_factor, only: structure_factor_t, sf_shared_setup, sf_alloc_partials, &
                                   sf_prepare_species, mode_denominator, norm_self, &
                                   norm_natom, no_unit
   use sqc_lebedev, only: lebedev_rule, lebedev_points, lebedev_reduce_pairs
   use sqc_modes, only: modes_t, modes_build, thin_none, thin_shells, thin_orbits
   use sqc_phase, only: phase_tables, phase_factor, phase_min_modes
   use, intrinsic :: iso_c_binding, only: c_double_complex
   use, intrinsic :: iso_fortran_env, only: output_unit
   use omp_lib, only: omp_get_max_threads, omp_get_thread_num
#ifdef SQC_HAVE_HDF5
   use sqc_hdf5, only: hdf5_write_dynamics, hdf5_write_chi4, hdf5_write_fqt_self, &
                       hdf5_write_msd
#endif
   implicit none
   private

   public :: dynamics_structure_factor_t, &
             dyn_q_none, dyn_q_line, dyn_q_shell, dyn_q_grid, dyn_q_single, &
             dyn_q_powder

   !> q sampling modes of a dynamic run (`--qpoints`).
   integer, parameter :: dyn_q_none = 0
   integer, parameter :: dyn_q_line = 1
   integer, parameter :: dyn_q_shell = 2
   integer, parameter :: dyn_q_grid = 3
   integer, parameter :: dyn_q_single = 4
   !> A window of lattice shells around a requested |q|, averaged into one row.
   integer, parameter :: dyn_q_powder = 5

   !> Which quantity a writer is asked for.
   integer, parameter :: dyn_sqw = 1
   integer, parameter :: dyn_fqt = 2
   integer, parameter :: dyn_s4 = 3

   !> Direct summation over a q line plus a multi-origin time correlation.
   type, extends(structure_factor_t) :: dynamics_structure_factor_t
      !> Time between two dumped frames [same unit as `--dt`].
      real(rk) :: frame_dt = 0.0_rk
      !> Correlation window in frames (the largest lag kept).
      integer :: maxframes = 0
      !> Frames between two consecutive time origins.
      integer :: lag_stride = 1
   !> The q line: `nintervals` intervals from `s0` to `s1` [1/A] along
   !! `direction` (unnormalized), through the Gamma point.
      integer :: nintervals = 0
      real(rk) :: s0 = 0.0_rk, s1 = 0.0_rk
      real(rk) :: direction(3) = 0.0_rk
      !> q sampling (dyn_q_none, dyn_q_line or dyn_q_shell) and the shell
      !! radius with the order of its Lebedev rule.
      integer :: q_mode = dyn_q_none
      real(rk) :: shell_q = 0.0_rk
      integer :: shell_order = 0
      !> Reciprocal-lattice grid sampling: qmax, mode budget, thinning policy,
      !! and what the build actually produced (for the run summary).
      real(rk) :: grid_qmax = 0.0_rk
      integer :: grid_budget = 0
      integer :: grid_thin = thin_none
      integer :: grid_nops = 0
      integer :: grid_nshell = 0
      !> True when the skew of the cell exceeded the point-group search box,
      !! so the modes are exact but the orbit reduction is unavailable.
      logical :: grid_group_fallback = .false.
      integer :: grid_nshell_kept = 0
      logical :: grid_thinned = .false.
      logical :: grid_budget_met = .true.
      !> Lattice shell and orbit multiplicity of every mode (grid sampling).
      integer, allocatable :: grid_shell(:), grid_orbit(:), orbit_mult(:)
      !> Miller indices of every mode and the range they span; allocated for a
      !! lattice sampling, which is what allows the separable phase evaluation.
      integer, allocatable :: mode_index(:, :)
      integer :: index_max(3) = 0
      !> Phase gradients: the argument of the factor table of axis `j` is
      !! `phi_j = recip(:, j) . r`, with `r` the position (rho) or the
      !! displacement (F_s).  A reciprocal-lattice mode fills the columns
      !! with the reciprocal basis `b_j`; a q line puts its step `dq*u` in
      !! the first column and leaves the other two at zero.  One phase-factor
      !! table per thread follows, plus the per-thread accumulators.
      real(rk) :: recip(3, 3) = 0.0_rk
      !> Mode-independent start of a q line: the first mode's phase is
      !! `phase_base_vec . r`, so the table of axis 1 has to be built from
      !! `exp(i phase_base_vec.r)` rather than from 1.  Unset, and therefore
      !! exactly 1, for a lattice sampling.
      logical :: phase_base_on = .false.
      real(rk) :: phase_base_vec(3) = 0.0_rk
      integer :: phase_threads = 1
      complex(c_double_complex), allocatable :: phase_pool(:, :, :)
      !> One workspace per thread: every worker owns a whole lag at a time, so
      !! the mode sums of a lag never share a partial result with another lag.
      complex(c_double_complex), allocatable :: acc_thread(:, :)
      complex(c_double_complex), allocatable :: fs_thread(:, :, :), rho_thread(:, :, :)
      !> Single lattice vector sampling (`--qpoints single`): the Miller indices
      !! and the q vector they build.
      integer :: single_index(3) = 0
      real(rk) :: single_q(3) = 0.0_rk
      !> Modes of the grid before the rows were collapsed onto the shells, and
      !! whether the per-mode rows were kept instead.
      integer :: grid_nmodes = 0
      !> Lattice-shell window (`dyn_q_powder`): the requested |q|, the target
      !! number of lattice vectors in the window, and the half width when the
      !! caller fixed one instead.
      real(rk) :: powder_q = 0.0_rk
      integer :: powder_modes = 0
      real(rk) :: powder_dq = 0.0_rk
      logical :: powder_dq_given = .false.
      !> What the window came out to be, for the header and the summary: the
      !! multiplicity weighted mean |q| of its vectors, the half width that
      !! was used, how many lattice vectors and shells it holds, and whether
      !! it had to be widened to the nearest shell.
      real(rk) :: powder_qmean = 0.0_rk
      real(rk) :: powder_dq_used = 0.0_rk
      integer :: powder_vectors = 0
      integer :: powder_nshell = 0
      logical :: powder_nearest = .false.
      logical :: keep_modes = .false.
      !> Isotropy probe of a grid run: at the lag of the chi4 peak, the S4
      !! values of the modes that share the smallest lattice shell.  Their
      !! expectation is equal for an equilibrium isotropic system, so their
      !! spread is a cheap smoke test for the orbit reduction and for the
      !! shell average.
      logical :: probe_valid = .false.
      real(rk) :: probe_tau = 0.0_rk
      real(rk) :: probe_mean = 0.0_rk
      real(rk) :: probe_spread = 0.0_rk
      integer :: probe_modes = 0
      !> Quadrature weight of every q mode: 1 on a line, the Lebedev weights on
      !! a shell (the shell average is `sum_k w_k X_k / sum_k w_k`).
      real(rk), allocatable :: mode_weight(:)
      !> Optional output files (empty = not written) and their format.
      character(len=:), allocatable :: sqw_path, fqt_path, fqt_self_path
      integer :: sqw_format = format_text
      integer :: fqt_format = format_text
      integer :: fqt_self_format = format_text
      !> Self intermediate scattering function F_s(q,t).
      logical :: fqt_self_enabled = .false.
      !> Four-point structure factor and average overlap / chi4.
      logical :: s4_enabled = .false.
      logical :: chi4_enabled = .false.
      real(rk) :: s4_cutoff = 0.0_rk
      real(rk) :: s4_cutoff2 = 0.0_rk
      character(len=:), allocatable :: s4_path, chi4_path
      integer :: s4_format = format_text
      integer :: chi4_format = format_text
      !> Mean squared displacement MSD(tau) and its per-species split.
      logical :: msd_enabled = .false.
      character(len=:), allocatable :: msd_path
      integer :: msd_format = format_text
      !> True when F(q,t)/S(q,w) buffers and transforms are needed.
      logical :: coherent_enabled = .true.
      !> Position buffer limit [GB, 10^9 bytes] for --s4/--chi4.
      real(rk) :: buffer_limit_gb = 2.0_rk
      !> Frame stride of the S4/chi4 trajectory (1 = every dump frame).
      integer :: stride = 1
      !> Number of stride steps in the effective S4 window.
      integer :: nsteps = 0
      !> Effective S4 window, nsteps*stride (<= maxframes).
      integer :: effective_maxframes = 0
      !> Ring buffer of the last frames' positions for the overlap function,
      !! `(3, atom, lag)`.  The lag is the last index on purpose: the first
      !! index varies fastest, so the scan of one lag walks the atoms
      !! contiguously and streams instead of fetching a cache line per atom.
      real(rk), allocatable :: pos_buffer(:, :, :)
      !> Atoms inside the overlap cutoff, one list per thread: a worker builds
      !! the list of the lag it owns and no other thread can touch it.
      integer(ik), allocatable :: over_seg(:, :)
      !> S4 accumulators: sum |W|^2 and sum W per (mode, lag).
      real(rk), allocatable :: s4_asum(:, :)
      complex(c_double_complex), allocatable :: s4_bsum(:, :)
      !> Self-function accumulators: sum exp(i q.dr) per (mode, lag, species).
      complex(c_double_complex), allocatable :: fs_sum(:, :, :)
      !> chi4 accumulators over the scalar overlap W0.
      real(rk), allocatable :: chi4_asum(:), chi4_bsum(:)
      !> MSD accumulator: sum |r_i(t0+tau) - r_i(t0)|^2 per (species, lag).
      real(rk), allocatable :: msd_sum(:, :)
      !> Time origins per S4 lag and the S4 time axis.
      integer(lk), allocatable :: sample_cnt(:, :)
      real(rk), allocatable :: sample_tau(:)
      !> Results: S4(q,t), average overlap Q(t) and chi4(t).
      real(rk), allocatable :: s4(:, :)
      real(rk), allocatable :: overlap(:), chi4(:)
      !> Results: total and per-species mean squared displacement.
      real(rk), allocatable :: msd_total(:), msd_partial(:, :)
      !> Results: total and per-species F_s(q,t).
      real(rk), allocatable :: fqt_self_total(:, :)
      real(rk), allocatable :: fqt_self_partial(:, :, :)
      !> Species present and their grouping by contiguous atom ranges.
      integer(ik), allocatable :: species_of(:), species_type(:)
      integer(ik), allocatable :: atom_of(:), species_first(:)
      integer :: nspecies = 0
      !> Amplitude w_a(q) of each species at each q (weights applied at output).
      real(rk), allocatable :: amp(:, :)
      logical :: q_dependent = .false.
      real(rk), allocatable :: amp_const(:), amp_table(:, :)
      !> Ring buffer rho(nmodes, nspecies, 0:maxframes).
      complex(c_double_complex), allocatable :: rho(:, :, :)
      !> Correlation sums and counts per (mode, lag, species, species).
      complex(c_double_complex), allocatable :: csum(:, :, :, :)
      integer(lk), allocatable :: ccnt(:, :)
      !> Zero lag sums over the whole trajectory (the static OUTPUT table).
      complex(c_double_complex), allocatable :: c0sum(:, :, :)
      !> Results: the spectra and the intermediate scattering function.
      real(rk), allocatable :: omega(:), sqw(:, :), sqw_partial(:, :, :)
      real(rk), allocatable :: tau(:), ftau(:, :), ftau_partial(:, :, :)
      !> Frequency grid spacing, pi/(maxframes*frame_dt).
      real(rk) :: dw = 0.0_rk
   contains
      procedure :: configure => dyn_configure
      procedure :: method_setup => dyn_setup
      procedure :: accumulate_frame => dyn_accumulate
      procedure :: prepare_output => dyn_prepare_output
      procedure :: write_sqw => dyn_write_sqw
      procedure :: write_fqt => dyn_write_fqt
      procedure :: write_fqt_self => dyn_write_fqt_self
      procedure :: write_s4 => dyn_write_s4
      procedure :: write_chi4 => dyn_write_chi4
      procedure :: write_msd => dyn_write_msd
   end type dynamics_structure_factor_t

contains

   !> Build the q line and the shared accumulators from the first frame.
   subroutine dyn_configure(self, frame, scheme, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: u(3), norm_u, s_i
      !! Lattice-shell window bookkeeping (dyn_q_powder).
      real(rk) :: qmax_use, dq_use, pi
      integer :: budget_use, thin_use
      real(rk), allocatable :: kept_q(:, :), kept_w(:)
      integer(lk) :: npts
      integer :: i
      type(modes_t) :: grid

      ierr = 0
      message = ''
      pi = acos(-1.0_rk)
      call sf_shared_setup(self, frame, scheme, ierr, message)
      if (ierr /= 0) return

      select case (self%q_mode)
      case (dyn_q_line)
         if (self%nintervals < 1) then
            ierr = 1
            message = '--qpoints line needs at least one interval on the q line'
            return
         end if
         norm_u = sqrt(sum(self%direction**2))
         if (norm_u <= 0.0_rk) then
            ierr = 1
            message = '--qpoints line needs a non-zero direction, e.g. 1,1,0'
            return
         end if
         u = self%direction/norm_u

         ! q_i = s_i u, anchored at the Gamma point.
         self%nmodes = int(self%nintervals + 1, lk)
         allocate (self%qvec(3, self%nmodes), self%qlen(self%nmodes), &
                   self%shell(self%nmodes), self%gidx(self%nmodes))
         do i = 1, int(self%nmodes)
            s_i = self%s0 + real(i - 1, rk)*(self%s1 - self%s0)/real(self%nintervals, rk)
            self%qvec(:, i) = s_i*u
            self%qlen(i) = s_i
            self%shell(i) = i
            self%gidx(i) = int(i, lk)
         end do
         ! The shared writer prints qmin + (shell-1/2)*shell_dq, so these two
         ! make it print the q line itself.
         self%shell_dq = (self%s1 - self%s0)/real(self%nintervals, rk)
         self%qmin = self%s0 - 0.5_rk*self%shell_dq
         self%qmax = self%s1 + 0.5_rk*self%shell_dq
         ! Separable phases of the line.  The scale steps by a constant, so
         ! q_i . r = s0*(u.r) + i*dq*(u.r): a mode-independent base times the
         ! i-th power of one factor, i.e. one table on the first axis with the
         ! other two held at index 0.  The gradient columns hold dq*u, and the
         ! base carries the s0*u part, both from the Cartesian direction
         ! because a line is not a lattice vector.
         allocate (self%mode_index(3, self%nmodes))
         do i = 1, int(self%nmodes)
            self%mode_index(:, i) = [i - 1, 0, 0]
         end do
         self%index_max = [self%nintervals, 0, 0]
         self%recip = 0.0_rk
         self%recip(:, 1) = self%shell_dq*u
         self%phase_base_on = abs(self%s0) > 0.0_rk
         self%phase_base_vec = self%s0*u
      case (dyn_q_shell)
         ! Every direction of the shell |q| = Q, on a Lebedev rule.  The modes
         ! are averaged into a single output row once they are assembled.
         npts = int(lebedev_points(self%shell_order), lk)
         if (npts < 1) then
            ierr = 1
            write (message, '(a,i0,a)') 'unsupported Lebedev order ', self%shell_order, &
               ' for --qpoints shell (use low, medium or high)'
            return
         end if
         self%nmodes = npts
         allocate (self%qvec(3, npts), self%mode_weight(npts))
         call lebedev_rule(self%shell_order, self%qvec, self%mode_weight, ierr, message)
         if (ierr /= 0) return
         ! The rule is symmetric under q -> -q with equal weights and every
         ! quantity averaged here is even in q, so one vector of each pair
         ! carries the pair with twice the weight; the average is unchanged
         ! and the shell costs half the modes.
         call lebedev_reduce_pairs(self%qvec, self%mode_weight, self%nmodes)
         ! The rule filled arrays sized for the full rule, so build the ones
         ! the run keeps at the reduced length: the shell average sums the
         ! weight array, and entries past nmodes would otherwise be summed as
         ! well.
         allocate (kept_q(3, self%nmodes), kept_w(self%nmodes))
         kept_q = self%qvec(:, 1:self%nmodes)
         kept_w = self%mode_weight(1:self%nmodes)
         call move_alloc(kept_q, self%qvec)
         call move_alloc(kept_w, self%mode_weight)
         allocate (self%qlen(self%nmodes), self%shell(self%nmodes), self%gidx(self%nmodes))
         self%qvec = self%shell_q*self%qvec
         self%qlen = self%shell_q
         do i = 1, int(self%nmodes)
            self%shell(i) = i
            self%gidx(i) = int(i, lk)
         end do
         ! The shell collapses to one row at |q| = Q, so it has no width.
         self%shell_dq = 0.0_rk
         self%qmin = self%shell_q
         self%qmax = self%shell_q
      case (dyn_q_grid, dyn_q_powder)
         ! Reciprocal-lattice sampling.  Only a lattice vector carries a
         ! density amplitude that is independent of how the periodic images
         ! are chosen, so this is the sampling the OZ fit of S4(q,t) needs;
         ! the modes are one output row each, exactly like the q line.  The
         ! `powder` variant averages a window of shells around the requested
         ! |q| into a single row instead (see below).
         qmax_use = self%grid_qmax
         dq_use = 0.0_rk
         budget_use = self%grid_budget
         thin_use = self%grid_thin
         if (self%q_mode == dyn_q_powder) then
            ! The window half width follows from the target number of lattice
            ! vectors: the reciprocal-lattice density is V/(2 pi)^3 and a
            ! window of half width dQ at |q| = Q holds about
            ! V Q^2 dQ / pi^2 of them.  The cap keeps the q resolution from
            ! collapsing in a small box; a caller who fixed dq gets it as is.
            dq_use = self%powder_dq
            if (.not. self%powder_dq_given) then
               dq_use = min(self%powder_q/20.0_rk, &
                            pi**2*real(max(self%powder_modes, 1), rk)/ &
                            (frame%cell%volume*self%powder_q**2))
            end if
            qmax_use = self%powder_q + dq_use
            ! Thinning would bias the average the window is meant to take, and
            ! the window width is the knob that controls the cost.
            budget_use = 0
            thin_use = thin_none
         end if
         call modes_build(grid, frame%cell, qmax_use, self%q_mode == dyn_q_grid, &
                          budget_use, thin_use, ierr, message)
         if (ierr /= 0) then
            if (self%q_mode == dyn_q_powder .and. &
                index(message, 'no reciprocal lattice points') > 0) then
               write (message, '(a,f0.4,a,f0.4,a)') 'no reciprocal-lattice vector within |q| = ', &
                  self%powder_q - dq_use, ' .. ', self%powder_q + dq_use, &
                  ' 1/A; the box cannot resolve that shell'
            end if
            return
         end if
         if (grid%nmodes < 1) then
            ierr = 1
            message = 'the reciprocal-lattice grid came out empty'
            return
         end if
         call dyn_adopt_modes(self, grid, frame%cell)
         if (self%q_mode == dyn_q_powder) then
            call dyn_powder_select(self, grid, frame%cell, dq_use, ierr, message)
            if (ierr /= 0) return
         end if
         call grid%finalize()
      case (dyn_q_single)
         ! One reciprocal-lattice vector, built from its Miller indices so that
         ! it sits exactly on the lattice.  The row is labelled by |q|, while
         ! the indices, the vector and its length all go into the table header.
         self%nmodes = 1
         allocate (self%qvec(3, 1), self%qlen(1), self%shell(1), self%gidx(1))
         self%single_q = matmul(frame%cell%b, real(self%single_index, rk))
         self%qlen(1) = sqrt(sum(self%single_q**2))
         ! qvec is what the density amplitude uses, so it keeps the true vector
         ! here; the rows are relabelled by |q| once the accumulation is done.
         self%qvec(:, 1) = self%single_q
         self%shell(1) = 1
         self%gidx(1) = 1
         self%shell_dq = 0.0_rk
         self%qmin = self%qlen(1)
         self%qmax = self%qlen(1)
      case default
         ! No q points at all: the scalar overlap Q(t)/chi4(t) is all that is
         ! left, and it needs neither the density amplitudes nor a table.
         self%nmodes = 0
         self%shell_dq = 0.0_rk
      end select

      self%nq = int(self%nmodes)
      allocate (self%num(self%nq), self%den(self%nq), self%shell_count(self%nq))
      allocate (self%mode_values(self%nmodes))
      self%num = 0.0_rk
      self%den = 0.0_rk
      self%mode_values = 0.0_rk
      self%shell_count = 1
      do i = 1, int(self%nmodes)
         self%den(i) = mode_denominator(self%norm, scheme, self%natoms, self%type_counts, &
                                        self%qlen(i))
      end do
      if (self%nmodes > 0) call sf_alloc_partials(self, scheme)

      call self%method_setup(frame, scheme, ierr, message)
   end subroutine dyn_configure

   !> Copy a reciprocal-lattice mode set into the run.
   !!
   !! The modes, their shells and orbits, the Miller indices and the reciprocal
   !! basis all come from the grid builder, so every lattice sampling (grid,
   !! powder, single) shares this one adoption step.  The caller finalizes the
   !! builder afterwards; a sampling that has to re-enumerate calls this again
   !! with the new set.
   subroutine dyn_adopt_modes(self, grid, cell)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(modes_t), intent(in) :: grid
      type(cell_t), intent(in) :: cell
      integer :: i

      self%nmodes = int(grid%nmodes, lk)
      self%qvec = grid%qvec
      self%qlen = grid%qlen
      self%shell_radii = grid%qlen
      self%orbit_mult = grid%orbit_mult
      self%grid_shell = grid%shell
      self%grid_orbit = grid%orbit
      self%shell = [(i, i = 1, int(self%nmodes))]
      self%gidx = [(int(i, lk), i = 1, int(self%nmodes))]
      ! The lattice shells are not evenly spaced, so the rows carry their own
      ! |q| instead of the qmin + (s-1/2) dq formula.
      self%label_by_shell_radius = .true.
      self%shell_dq = 0.0_rk
      self%qmin = minval(self%qlen)
      self%qmax = maxval(self%qlen)
      self%grid_nops = grid%nops
      self%grid_nshell = grid%nshell
      self%grid_group_fallback = grid%point_group_fallback
      self%grid_nshell_kept = grid%nshell_kept
      self%grid_thinned = grid%thinned
      self%grid_budget_met = grid%budget_met
      self%grid_nmodes = int(self%nmodes)
      ! The separable phase evaluation needs the Miller indices and the
      ! reciprocal basis, both of which the grid builder already has.
      self%mode_index = grid%hkl
      do i = 1, 3
         self%index_max(i) = maxval(abs(grid%hkl(i, :)))
      end do
      self%recip = cell%b
   end subroutine dyn_adopt_modes

   !> Weight the modes of a powder window and fill its report fields.
   !!
   !! `grid` holds the lattice vectors within `dq` of the requested |q|, which
   !! is the band the run pays for.  Only whole shells can be averaged, so a
   !! window that falls between two of them has to fall back on the nearest
   !! one.  That search re-enumerates the lattice up to `2 Q`: the shell above
   !! the request is not part of the window build, and without it the nearest
   !! shell below the request would win even when a closer one sits just above.
   subroutine dyn_powder_select(self, grid, cell, dq, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(modes_t), intent(inout) :: grid
      type(cell_t), intent(in) :: cell
      real(rk), intent(inout) :: dq
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: stol, dnear, wsum, qsum
      integer :: i, s, sbest, nkeep_o, nshell_used
      logical :: averaged

      ierr = 0
      message = ''
      stol = 1.0e-6_rk*max(self%powder_q, 1.0_rk)
      sbest = 0
      dnear = huge(1.0_rk)
      do s = 1, grid%nshell
         if (abs(grid%shell_q(s) - self%powder_q) < dnear) then
            dnear = abs(grid%shell_q(s) - self%powder_q)
            sbest = s
         end if
      end do
      self%powder_nearest = dnear > dq + stol
      if (self%powder_nearest) then
         ! The window held no shell of its own: widen the search and average
         ! exactly the nearest shell, rather than a band that would also take
         ! in whatever else lies between it and the request.
         call modes_build(grid, cell, 2.0_rk*self%powder_q, .false., 0, &
                          thin_none, ierr, message)
         ! The first build already found a vector within `Q + dq` and this bound
         ! is wider, so the enumeration cannot come back empty here: whatever
         ! failure arrives (the grid size limit, say) is its own and keeps its
         ! own message.
         if (ierr /= 0) return
         call dyn_adopt_modes(self, grid, cell)
         dnear = huge(1.0_rk)
         sbest = 0
         do s = 1, grid%nshell
            if (abs(grid%shell_q(s) - self%powder_q) < dnear) then
               dnear = abs(grid%shell_q(s) - self%powder_q)
               sbest = s
            end if
         end do
         dq = dnear + stol
      end if
      if (sbest < 1) then
         ierr = 1
         message = 'the lattice shell window came out empty'
         return
      end if

      ! Keep every mode but weight only the ones the row averages: the
      ! single-row collapse of dyn_shell_average then averages exactly them,
      ! and the sum rules stay exact because numerator and denominator carry
      ! the same weights.  The weight of a mode is its share of the orbit it
      ! stands for (an orbit lies in one shell, so its shares add up to the
      ! multiplicity of the vector), which makes the row an average over the
      ! lattice vectors of the window rather than over the modes.
      allocate (self%mode_weight(self%nmodes))
      self%mode_weight = 0.0_rk
      do i = 1, int(self%nmodes)
         s = self%grid_shell(i)
         if (s < 1) cycle
         if (self%powder_nearest) then
            averaged = s == sbest
         else
            averaged = abs(grid%shell_q(s) - self%powder_q) <= dq + stol
         end if
         if (.not. averaged) cycle
         nkeep_o = count(self%grid_orbit == self%grid_orbit(i))
         self%mode_weight(i) = real(self%orbit_mult(i), rk)/real(max(nkeep_o, 1), rk)
      end do
      wsum = sum(self%mode_weight)
      if (wsum <= 0.0_rk) then
         ierr = 1
         message = 'the lattice shell window came out empty'
         return
      end if
      qsum = 0.0_rk
      do i = 1, int(self%nmodes)
         qsum = qsum + self%mode_weight(i)*self%qlen(i)
      end do
      self%powder_qmean = qsum/wsum
      self%powder_dq_used = dq
      self%powder_vectors = nint(wsum)
      nshell_used = 0
      do s = 1, grid%nshell
         if (any(self%mode_weight > 0.0_rk .and. self%grid_shell == s)) &
            nshell_used = nshell_used + 1
      end do
      self%powder_nshell = nshell_used
      ! The row carries the weighted mean |q|.  The writers label a row as
      ! qmin + (s - 1/2) dq when the shells are evenly spaced (what the Lebedev
      ! shell does), so the mean goes there; the window itself is reported
      ! through the powder fields.
      self%shell_q = self%powder_qmean
      self%shell_dq = 0.0_rk
      self%qmin = self%powder_qmean
      self%qmax = self%powder_qmean
      self%label_by_shell_radius = .false.
   end subroutine dyn_powder_select

   !> Species bookkeeping, weights per q and the correlation buffers.
   subroutine dyn_setup(self, frame, scheme, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(ik), allocatable :: cursor(:)
      integer(lk) :: limit_bytes, need_bytes
      integer :: i, t, isp, astat
      character(len=256) :: amsg

      ierr = 0
      message = ''
      ! Validate the overlap request before the coherent ring buffers are
      ! allocated, so an oversized --maxframes fails with the S4 message
      ! instead of exhausting memory in the unrelated F(q,t) accumulators.
      if (self%s4_enabled .or. self%chi4_enabled .or. self%fqt_self_enabled .or. &
          self%msd_enabled) then
         if (self%s4_enabled .or. self%chi4_enabled) then
            if (self%s4_cutoff <= 0.0_rk) then
               ierr = 1
               message = '--s4-cutoff must be a positive number'
               return
            end if
            self%s4_cutoff2 = self%s4_cutoff**2
         end if
         if (self%stride < 1) then
            ierr = 1
            message = '--stride must be a positive integer'
            return
         end if
         self%nsteps = self%maxframes/self%stride
         if (self%nsteps < 1) then
            ierr = 1
            message = '--stride cannot exceed --maxframes'
            return
         end if
         self%effective_maxframes = self%nsteps*self%stride
         limit_bytes = int(self%buffer_limit_gb*1.0e9_rk, lk)
         need_bytes = 24_lk*int(self%natoms, lk)*int(self%nsteps + 1, lk)
         if (need_bytes > limit_bytes) then
            ierr = 1
            write (message, '(a,f0.4,a,i0,a,f0.4,a,i0,a)') &
               'the S4/chi4/F_s/MSD position buffer needs ', real(need_bytes, rk)/1.0e9_rk, &
               ' GB (', need_bytes, ' bytes), above the ', self%buffer_limit_gb, &
               ' GB (', limit_bytes, ' bytes) limit; raise --buffer-limit, '// &
               'increase --stride or reduce --maxframes'
            return
         end if
      end if
      call sf_prepare_species(self, frame, scheme, self%species_of, self%species_type, &
                              self%amp_const, self%amp_table, self%q_dependent, ierr, message)
      if (ierr /= 0) return
      self%nspecies = size(self%species_type)

      ! Contiguous atom ranges per species, as in the direct method.
      allocate (self%species_first(self%nspecies + 1))
      self%species_first = 1
      do i = 1, frame%natoms
         t = int(self%species_of(i))
         self%species_first(t + 1) = self%species_first(t + 1) + 1
      end do
      do t = 1, self%nspecies
         self%species_first(t + 1) = self%species_first(t + 1) + self%species_first(t) - 1
      end do
      allocate (self%atom_of(frame%natoms))
      allocate (cursor(self%nspecies))
      cursor = self%species_first(1:self%nspecies)
      do i = 1, frame%natoms
         t = int(self%species_of(i))
         self%atom_of(cursor(t)) = int(i, ik)
         cursor(t) = cursor(t) + 1
      end do
      deallocate (cursor)

      ! Weight of every species at every q, applied when the results are
      ! assembled (never inside the per frame sum).
      ! A run without --qpoints has no c0sum and no rho at all.
      if (self%nmodes > 0) then
         allocate (self%amp(self%nmodes, self%nspecies))
         if (self%q_dependent) then
            self%amp = self%amp_table
         else
            do isp = 1, self%nspecies
               self%amp(:, isp) = self%amp_const(isp)
            end do
         end if

         if (self%coherent_enabled) then
            allocate (self%rho(self%nmodes, self%nspecies, 0:self%maxframes))
            self%rho = (0.0_rk, 0.0_rk)
            allocate (self%csum(self%nmodes, 0:self%maxframes, self%nspecies, self%nspecies))
            self%csum = (0.0_rk, 0.0_rk)
         else
            ! The static S(q) table only needs the current frame; the ring
            ! buffer and the coherent multi-origin correlation belong to
            ! F(q,t)/S(q,w).
            allocate (self%rho(self%nmodes, self%nspecies, 0:0))
            self%rho = (0.0_rk, 0.0_rk)
         end if
         allocate (self%c0sum(self%nmodes, self%nspecies, self%nspecies))
         self%c0sum = (0.0_rk, 0.0_rk)
         allocate (self%ccnt(self%nmodes, 0:self%maxframes))
         self%ccnt = 0
      end if

      ! --- four point structure factor / overlap accumulators ---------------
      ! S4(q,t), chi4(t), F_s(q,t) and MSD(t) all need the positions of the
      ! frames that are still inside the correlation window.  The list of
      ! overlapping atoms is built once per lag and reused by every q mode.
      if (self%s4_enabled .or. self%chi4_enabled .or. self%fqt_self_enabled .or. &
          self%msd_enabled) then
         allocate (self%pos_buffer(3, int(self%natoms), 0:self%nsteps), &
                   stat=astat, errmsg=amsg)
         if (astat /= 0) then
            ierr = 1
            message = 'cannot allocate the position buffer for S4/chi4/F_s/MSD: '//trim(amsg)
            return
         end if
         self%pos_buffer = 0.0_rk
         ! The origin count is the same at every q mode, so one row is enough
         ! for the writers; keep a row even without q points.
         allocate (self%sample_cnt(max(int(self%nmodes), 1), 0:self%nsteps))
         self%sample_cnt = 0
      end if
      ! Every worker owns a whole lag at a time, so the overlap list it builds
      ! belongs to it alone and needs no merge afterwards.  The thread count is
      ! the one the driver pinned before configure.
      self%phase_threads = max(omp_get_max_threads(), 1)
      if (self%s4_enabled .or. self%chi4_enabled) then
         allocate (self%over_seg(int(self%natoms), self%phase_threads))
         allocate (self%chi4_asum(0:self%nsteps), self%chi4_bsum(0:self%nsteps))
         self%chi4_asum = 0.0_rk
         self%chi4_bsum = 0.0_rk
      end if
      if (self%msd_enabled) then
         allocate (self%msd_sum(self%nspecies, 0:self%nsteps))
         self%msd_sum = 0.0_rk
      end if
      if (self%s4_enabled) then
         allocate (self%acc_thread(self%nmodes, self%phase_threads))
         self%acc_thread = (0.0_rk, 0.0_rk)
         allocate (self%s4_asum(self%nmodes, 0:self%nsteps))
         allocate (self%s4_bsum(self%nmodes, 0:self%nsteps))
         self%s4_asum = 0.0_rk
         self%s4_bsum = (0.0_rk, 0.0_rk)
      end if
      if (self%fqt_self_enabled) then
         allocate (self%fs_thread(self%nmodes, self%nspecies, self%phase_threads))
         self%fs_thread = (0.0_rk, 0.0_rk)
         allocate (self%fs_sum(self%nmodes, 0:self%nsteps, self%nspecies))
         self%fs_sum = (0.0_rk, 0.0_rk)
      end if
      ! Separable phase evaluation: one factor table per thread, filled and
      ! used by the worker that owns the lag or the frame being summed.
      if (allocated(self%mode_index)) then
         allocate (self%phase_pool(3, 0:2*maxval(self%index_max), self%phase_threads))
         self%phase_pool = (0.0_rk, 0.0_rk)
         if (self%coherent_enabled) then
            allocate (self%rho_thread(self%nmodes, self%nspecies, self%phase_threads))
            self%rho_thread = (0.0_rk, 0.0_rk)
         end if
      end if
   end subroutine dyn_setup

   !> Add one frame: the density amplitudes and their time correlations.
   subroutine dyn_accumulate(self, frame, scheme, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      complex(c_double_complex) :: acc
      real(rk) :: qx, qy, qz, phase, dx, dy, dz, rr(3)
      integer :: t_frame, slot, sample_slot, sample, l, j, oslot, im, isp, jsp, i, k, nover
      integer :: tid

      ierr = 0
      message = ''
      if (.not. frame%unwrapped) then
         ierr = 1
         message = 'the dynamic method needs unwrapped coordinates'
         return
      end if
      if (frame%natoms /= self%natoms) then
         ierr = 1
         message = 'the number of atoms changed between frames'
         return
      end if

      t_frame = int(self%nframes)
      if (self%coherent_enabled) then
         slot = mod(t_frame, self%maxframes + 1)
      else
         slot = 0
      end if

      ! Without q points there are no amplitudes and nothing to correlate; the
      ! overlap accumulators below are the whole calculation.
      if (self%nmodes > 0) then
         ! --- rho_a(q,t) by direct summation over the q line ----------------
         if (allocated(self%rho_thread) .and. self%nmodes >= phase_min_modes) then
            ! Separable phases, as in the S4 branch: the atom loop is the
            ! parallel one, and every mode costs two complex multiplies.
            self%rho_thread = (0.0_rk, 0.0_rk)
            !$omp parallel private(i, im, isp, tid)
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            !$omp do schedule(static)
            do i = 1, int(self%natoms)
               call dyn_atom_phase(self, frame%pos(:, i), tid)
               isp = int(self%species_of(i))
               do im = 1, int(self%nmodes)
                  self%rho_thread(im, isp, tid) = self%rho_thread(im, isp, tid) &
                     + phase_factor(self%phase_pool(:, :, tid), self%index_max, &
                                    self%mode_index(:, im))
               end do
            end do
            !$omp end do
            !$omp end parallel
            do tid = 2, self%phase_threads
               do isp = 1, self%nspecies
                  do im = 1, int(self%nmodes)
                     self%rho_thread(im, isp, 1) = self%rho_thread(im, isp, 1) &
                        + self%rho_thread(im, isp, tid)
                  end do
               end do
            end do
            do isp = 1, self%nspecies
               do im = 1, int(self%nmodes)
                  self%rho(im, isp, slot) = self%rho_thread(im, isp, 1)
               end do
            end do
         else
         !$omp parallel do schedule(static) private(im, isp, i, acc, phase, qx, qy, qz)
         do im = 1, int(self%nmodes)
            qx = self%qvec(1, im)
            qy = self%qvec(2, im)
            qz = self%qvec(3, im)
            do isp = 1, self%nspecies
               acc = (0.0_rk, 0.0_rk)
               do i = self%species_first(isp), self%species_first(isp + 1) - 1
                  phase = qx*frame%pos(1, self%atom_of(i)) + qy*frame%pos(2, self%atom_of(i)) &
                          + qz*frame%pos(3, self%atom_of(i))
                  acc = acc + cmplx(cos(phase), sin(phase), c_double_complex)
               end do
               self%rho(im, isp, slot) = acc
            end do
         end do
         !$omp end parallel do
         end if

         ! --- zero lag sums over the whole trajectory -----------------------
         do isp = 1, self%nspecies
            do jsp = 1, self%nspecies
               do im = 1, int(self%nmodes)
                  self%c0sum(im, isp, jsp) = self%c0sum(im, isp, jsp) &
                     + self%rho(im, isp, slot)*conjg(self%rho(im, jsp, slot))
               end do
            end do
         end do

         ! --- multi-origin correlation of the lags --------------------------
         if (self%coherent_enabled) then
            !$omp parallel do schedule(static) private(l, oslot, im, isp, jsp)
            do l = 0, min(t_frame, self%maxframes)
               if (mod(t_frame - l, self%lag_stride) /= 0) cycle
               oslot = mod(t_frame - l, self%maxframes + 1)
               do im = 1, int(self%nmodes)
                  do isp = 1, self%nspecies
                     do jsp = 1, self%nspecies
                        self%csum(im, l, isp, jsp) = self%csum(im, l, isp, jsp) &
                           + self%rho(im, isp, slot)*conjg(self%rho(im, jsp, oslot))
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         end if
         do l = 0, min(t_frame, self%maxframes)
            if (mod(t_frame - l, self%lag_stride) /= 0) cycle
            self%ccnt(:, l) = self%ccnt(:, l) + 1
         end do
      end if

      ! --- self F_s(q,t) and the four-point structure factor ---------------
      ! The S4/chi4/F_s trajectory is subsampled by stride: only every
      ! stride-th dump frame contributes, and lag j corresponds to
      ! j*stride original frames.  The coherent F(q,t)/S(q,w) path above
      ! still uses every dump frame.  F_s always includes every atom; the
      ! overlap cutoff is only used by S4/chi4.
      if ((self%s4_enabled .or. self%chi4_enabled .or. self%fqt_self_enabled .or. &
           self%msd_enabled) .and. mod(t_frame, self%stride) == 0) then
         sample = t_frame/self%stride
         sample_slot = mod(sample, self%nsteps + 1)
         self%pos_buffer(:, :, sample_slot) = frame%pos
         ! One lag per thread.  Every lag is an independent sum into its own
         ! row of the accumulators, so two workers never share a partial
         ! result and none of them has to wait for another: the schedule
         ! cannot change an answer, a lag is summed in the same order it would
         ! be alone, and the region is entered once per frame instead of once
         ! per lag.  `dynamic` spreads the uneven work - a lag that still
         ! overlaps every atom against one that overlaps few - without
         ! touching the result.  The scan of the overlap list rides along with
         ! the lag that owns it: it is a streaming pass over memory that came
         ! out no faster when it was split across threads on its own, so it is
         ! not split separately.
         !
         ! The lags are only spread when they carry a mode sum to pay for it.
         ! A chi4- or MSD-only lag is nothing but a streaming scan, and
         ! spreading that over threads cost a third of its time on the
         ! machines this was timed on, so those cases run on one thread.
         !$omp parallel do if(self%s4_enabled .or. self%fqt_self_enabled) &
         !$omp& schedule(dynamic, 1) &
         !$omp& private(oslot, nover, k, i, im, isp, tid, dx, dy, dz, rr, &
         !$omp&         qx, qy, qz, phase, acc)
         do j = 0, min(sample, self%nsteps)
            if (mod(t_frame - j*self%stride, self%lag_stride) /= 0) cycle
            oslot = mod(sample - j, self%nsteps + 1)
            tid = 1
            !$ tid = omp_get_thread_num() + 1

            ! Displacement of every atom.  S4/chi4 additionally build the
            ! overlap list, while F_s needs the displacement phase of every
            ! atom; the cutoff never enters F_s.
            nover = 0
            if (self%s4_enabled .or. self%chi4_enabled) then
               do i = 1, int(self%natoms)
                  dx = frame%pos(1, i) - self%pos_buffer(1, i, oslot)
                  dy = frame%pos(2, i) - self%pos_buffer(2, i, oslot)
                  dz = frame%pos(3, i) - self%pos_buffer(3, i, oslot)
                  if (dx*dx + dy*dy + dz*dz <= self%s4_cutoff2) then
                     nover = nover + 1
                     self%over_seg(nover, tid) = int(i, ik)
                  end if
               end do
            end if

            if (self%fqt_self_enabled) then
               if (allocated(self%phase_pool) .and. self%nmodes >= phase_min_modes) then
                  ! Separable phases of the *displacement*: the same identity
                  ! with phi_j = b_j . dr_i, so F_s costs two complex multiplies
                  ! per (mode, atom) as well.
                  self%fs_thread(:, :, tid) = (0.0_rk, 0.0_rk)
                  do i = 1, int(self%natoms)
                     rr(1) = frame%pos(1, i) - self%pos_buffer(1, i, oslot)
                     rr(2) = frame%pos(2, i) - self%pos_buffer(2, i, oslot)
                     rr(3) = frame%pos(3, i) - self%pos_buffer(3, i, oslot)
                     call dyn_atom_phase(self, rr, tid)
                     isp = int(self%species_of(i))
                     do im = 1, int(self%nmodes)
                        self%fs_thread(im, isp, tid) = self%fs_thread(im, isp, tid) &
                           + phase_factor(self%phase_pool(:, :, tid), self%index_max, &
                                          self%mode_index(:, im))
                     end do
                  end do
                  do isp = 1, self%nspecies
                     do im = 1, int(self%nmodes)
                        self%fs_sum(im, j, isp) = self%fs_sum(im, j, isp) &
                           + self%fs_thread(im, isp, tid)
                     end do
                  end do
               else
                  do im = 1, int(self%nmodes)
                     qx = self%qvec(1, im)
                     qy = self%qvec(2, im)
                     qz = self%qvec(3, im)
                     do isp = 1, self%nspecies
                        acc = (0.0_rk, 0.0_rk)
                        do k = self%species_first(isp), self%species_first(isp + 1) - 1
                           i = int(self%atom_of(k))
                           phase = qx*(frame%pos(1, i) - self%pos_buffer(1, i, oslot)) &
                                   + qy*(frame%pos(2, i) - self%pos_buffer(2, i, oslot)) &
                                   + qz*(frame%pos(3, i) - self%pos_buffer(3, i, oslot))
                           acc = acc + cmplx(cos(phase), sin(phase), c_double_complex)
                        end do
                        self%fs_thread(im, isp, tid) = acc
                     end do
                  end do
                  do isp = 1, self%nspecies
                     do im = 1, int(self%nmodes)
                        self%fs_sum(im, j, isp) = self%fs_sum(im, j, isp) &
                           + self%fs_thread(im, isp, tid)
                     end do
                  end do
               end if
            end if

            if (self%s4_enabled) then
               if (nover > 0 .and. allocated(self%mode_index) .and. &
                   self%nmodes >= phase_min_modes) then
                  ! Separable phases: one factor table per atom replaces the
                  ! sine and cosine of every (mode, atom) pair, and every mode
                  ! costs two complex multiplies.
                  self%acc_thread(:, tid) = (0.0_rk, 0.0_rk)
                  do k = 1, nover
                     i = int(self%over_seg(k, tid))
                     call dyn_atom_phase(self, self%pos_buffer(:, i, oslot), tid)
                     do im = 1, int(self%nmodes)
                        self%acc_thread(im, tid) = self%acc_thread(im, tid) &
                           + phase_factor(self%phase_pool(:, :, tid), self%index_max, &
                                          self%mode_index(:, im))
                     end do
                  end do
               else if (nover > 0) then
                  do im = 1, int(self%nmodes)
                     qx = self%qvec(1, im)
                     qy = self%qvec(2, im)
                     qz = self%qvec(3, im)
                     acc = (0.0_rk, 0.0_rk)
                     do k = 1, nover
                        i = int(self%over_seg(k, tid))
                        phase = qx*self%pos_buffer(1, i, oslot) &
                                + qy*self%pos_buffer(2, i, oslot) &
                                + qz*self%pos_buffer(3, i, oslot)
                        acc = acc + cmplx(cos(phase), sin(phase), c_double_complex)
                     end do
                     self%acc_thread(im, tid) = acc
                  end do
               else
                  self%acc_thread(:, tid) = (0.0_rk, 0.0_rk)
               end if
               do im = 1, int(self%nmodes)
                  self%s4_bsum(im, j) = self%s4_bsum(im, j) + self%acc_thread(im, tid)
                  self%s4_asum(im, j) = self%s4_asum(im, j) &
                     + real(self%acc_thread(im, tid), rk)**2 &
                     + aimag(self%acc_thread(im, tid))**2
               end do
            end if

            if (self%s4_enabled .or. self%chi4_enabled) then
               self%chi4_bsum(j) = self%chi4_bsum(j) + real(nover, rk)
               self%chi4_asum(j) = self%chi4_asum(j) + real(nover, rk)**2
            end if

            ! MSD(tau): the squared displacement of every atom, split by
            ! species.  It never uses the cutoff and, like chi4, is a plain
            ! streaming pass, so the lag that owns it does the whole sum.
            if (self%msd_enabled) then
               do i = 1, int(self%natoms)
                  dx = frame%pos(1, i) - self%pos_buffer(1, i, oslot)
                  dy = frame%pos(2, i) - self%pos_buffer(2, i, oslot)
                  dz = frame%pos(3, i) - self%pos_buffer(3, i, oslot)
                  isp = int(self%species_of(i))
                  self%msd_sum(isp, j) = self%msd_sum(isp, j) + dx*dx + dy*dy + dz*dz
               end do
            end if
            self%sample_cnt(:, j) = self%sample_cnt(:, j) + 1
         end do
         !$omp end parallel do
      end if

      self%nframes = self%nframes + 1
   end subroutine dyn_accumulate

   !> Fill the phase-factor tables of one atom's contribution.
   !!
   !! `r` is the position of the atom for the density amplitude and its
   !! displacement for `F_s`; both enter only as the three projections
   !! `phi_j = recip(:, j) . r` that the separable evaluation needs.  A q line
   !! additionally carries the mode-independent start `s0` of its scale, which
   !! becomes the base of the first table - no lattice sampling has one, so
   !! the base-free build is the common case and stays a branch away.
   subroutine dyn_atom_phase(self, r, tid)
      class(dynamics_structure_factor_t), intent(inout) :: self
      real(rk), intent(in) :: r(3)
      integer, intent(in) :: tid
      complex(c_double_complex) :: base(3)
      real(rk) :: phi(3), d
      integer :: im

      do im = 1, 3
         phi(im) = self%recip(1, im)*r(1) + self%recip(2, im)*r(2) &
                 + self%recip(3, im)*r(3)
      end do
      if (self%phase_base_on) then
         d = self%phase_base_vec(1)*r(1) + self%phase_base_vec(2)*r(2) &
           + self%phase_base_vec(3)*r(3)
         base = (1.0_rk, 0.0_rk)
         base(1) = cmplx(cos(d), sin(d), c_double_complex)
         call phase_tables(self%phase_pool(:, :, tid), phi, self%index_max, base)
      else
         call phase_tables(self%phase_pool(:, :, tid), phi, self%index_max)
      end if
   end subroutine dyn_atom_phase

   !> Assemble the static table, F(q,tau) and the spectra.
   subroutine dyn_prepare_output(self, scheme, ierr, message)
      class(dynamics_structure_factor_t), intent(inout) :: self
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: f(:)
      real(rk) :: frames, value
      integer :: im, l, j, isp, jsp, t, u, p, n_pairs, n

      ierr = 0
      message = ''
      if (self%nframes <= self%maxframes) then
         write (message, '(a,i0,a,i0,a)') 'the trajectory has ', self%nframes, &
            ' frames but --maxframes is ', self%maxframes, &
            '; at least maxframes+1 frames are needed'
         ierr = 1
         return
      end if
      frames = real(self%nframes, rk)
      allocate (self%tau(0:self%maxframes))
      do l = 0, self%maxframes
         self%tau(l) = real(l, rk)*self%frame_dt
      end do
      if (self%coherent_enabled) then
         self%dw = acos(-1.0_rk)/(real(self%maxframes, rk)*self%frame_dt)
         allocate (self%omega(0:self%maxframes))
         do l = 0, self%maxframes
            self%omega(l) = real(l, rk)*self%dw
         end do
         allocate (self%sqw(self%nmodes, 0:self%maxframes))
         allocate (self%ftau(self%nmodes, 0:self%maxframes))
         self%sqw = 0.0_rk
         self%ftau = 0.0_rk
      end if

      ! Static table: |rho(q)|^2 averaged over every frame, with the same
      ! normalization as the other methods.
      do im = 1, int(self%nmodes)
         value = 0.0_rk
         do isp = 1, self%nspecies
            do jsp = 1, self%nspecies
               value = value + self%amp(im, isp)*self%amp(im, jsp) &
                       *real(self%c0sum(im, isp, jsp), rk)
            end do
         end do
         self%num(im) = value
      end do
      ! A run without q points has no partial columns (they are allocated with
      ! the shells in dyn_configure).
      if (self%partials .and. allocated(self%partial_num)) then
         self%partial_num = 0.0_rk
         do isp = 1, self%nspecies
            do jsp = 1, self%nspecies
               t = int(self%species_type(isp))
               u = int(self%species_type(jsp))
               do im = 1, int(self%nmodes)
                  self%partial_num(t, u, im) = real(self%c0sum(im, isp, jsp), rk)
               end do
            end do
         end do
      end if

      if (self%coherent_enabled) then
         ! F(q,tau), and the partial intermediate functions.
         n_pairs = 0
         do t = 1, self%ntypes
            do u = t, self%ntypes
               if (.not. pair_present(self, t, u)) cycle
               n_pairs = n_pairs + 1
            end do
         end do
         if (self%partials) then
            allocate (self%ftau_partial(self%nmodes, 0:self%maxframes, max(n_pairs, 1)))
            allocate (self%sqw_partial(self%nmodes, 0:self%maxframes, max(n_pairs, 1)))
            self%ftau_partial = 0.0_rk
            self%sqw_partial = 0.0_rk
         else
            allocate (self%ftau_partial(0, 0, 0), self%sqw_partial(0, 0, 0))
         end if

         allocate (f(0:self%maxframes))
         do im = 1, int(self%nmodes)
            ! total
            do l = 0, self%maxframes
               if (l == 0) then
                  value = 0.0_rk
                  do isp = 1, self%nspecies
                     do jsp = 1, self%nspecies
                        value = value + self%amp(im, isp)*self%amp(im, jsp) &
                                *real(self%c0sum(im, isp, jsp), rk)/frames
                     end do
                  end do
               else
                  value = 0.0_rk
                  if (self%ccnt(im, l) > 0) then
                     do isp = 1, self%nspecies
                        do jsp = 1, self%nspecies
                           value = value + self%amp(im, isp)*self%amp(im, jsp) &
                                   *real(self%csum(im, l, isp, jsp), rk)/real(self%ccnt(im, l), rk)
                        end do
                     end do
                  end if
               end if
               if (self%den(im) > 0.0_rk) then
                  self%ftau(im, l) = value/self%den(im)
               else
                  self%ftau(im, l) = 0.0_rk
               end if
            end do
            f = self%ftau(im, :)
            call dyn_spectrum(f, self%frame_dt, self%sqw(im, :))

            ! partials, in the OVITO (1/N) convention of the static columns
            if (self%partials) then
               p = 0
               do t = 1, self%ntypes
                  do u = t, self%ntypes
                     if (.not. pair_present(self, t, u)) cycle
                     p = p + 1
                     isp = species_index(self, t)
                     jsp = species_index(self, u)
                     do l = 0, self%maxframes
                        if (l == 0) then
                           value = real(self%c0sum(im, isp, jsp), rk)/frames
                        else
                           value = 0.0_rk
                           if (self%ccnt(im, l) > 0) then
                              value = real(self%csum(im, l, isp, jsp), rk)/real(self%ccnt(im, l), rk)
                              ! The cross partial has to be symmetrized, otherwise
                              ! the OVITO sum rule S = S_aa + 2 S_ab + S_bb would
                              ! hold only at zero lag (C_ab and C_ba differ by the
                              ! finite-trajectory noise).
                              if (isp /= jsp) value = 0.5_rk*(value &
                                 + real(self%csum(im, l, jsp, isp), rk)/real(self%ccnt(im, l), rk))
                           end if
                        end if
                        self%ftau_partial(im, l, p) = value/real(self%natoms, rk)
                     end do
                     f = self%ftau_partial(im, :, p)
                     call dyn_spectrum(f, self%frame_dt, self%sqw_partial(im, :, p))
                  end do
               end do
            end if
         end do
         deallocate (f)
      end if

      ! --- four point structure factor and the scalar overlap ----------------
      ! The connected covariance is evaluated with the unbiased (n-1)
      ! denominator; a lag with a single time origin has no fluctuation and is
      ! written as zero (the HDF5 counts show why).
      if (self%s4_enabled) then
         allocate (self%s4(self%nmodes, 0:self%nsteps))
         self%s4 = 0.0_rk
         do im = 1, int(self%nmodes)
            do j = 0, self%nsteps
               n = int(self%sample_cnt(im, j))
               if (n > 1) then
                  value = (self%s4_asum(im, j) &
                           - abs(self%s4_bsum(im, j))**2/real(n, rk))/real(n - 1, rk)
                  self%s4(im, j) = value/real(self%natoms, rk)
               end if
            end do
         end do
      end if
      if (self%chi4_enabled) then
         allocate (self%overlap(0:self%nsteps), self%chi4(0:self%nsteps))
         self%overlap = 0.0_rk
         self%chi4 = 0.0_rk
         do j = 0, self%nsteps
            n = int(self%sample_cnt(1, j))
            if (n > 0) then
               self%overlap(j) = (self%chi4_bsum(j)/real(n, rk))/real(self%natoms, rk)
            end if
            if (n > 1) then
               value = (self%chi4_asum(j) &
                        - self%chi4_bsum(j)**2/real(n, rk))/real(n - 1, rk)
               self%chi4(j) = value/real(self%natoms, rk)
            end if
         end do
      end if
      if (self%fqt_self_enabled) then
         allocate (self%fqt_self_partial(self%nmodes, 0:self%nsteps, self%nspecies))
         allocate (self%fqt_self_total(self%nmodes, 0:self%nsteps))
         self%fqt_self_partial = 0.0_rk
         self%fqt_self_total = 0.0_rk
         do im = 1, int(self%nmodes)
            do j = 0, self%nsteps
               n = int(self%sample_cnt(im, j))
               if (n > 0) then
                  do isp = 1, self%nspecies
                     value = real(self%fs_sum(im, j, isp), rk)/real(n, rk) &
                             /real(self%natoms, rk)
                     self%fqt_self_partial(im, j, isp) = value
                     self%fqt_self_total(im, j) = self%fqt_self_total(im, j) + value
                  end do
               end if
            end do
         end do
      end if
      ! MSD(tau): the origin count is the same at every q mode, so the first
      ! row of sample_cnt carries it.  The per-species columns use the same
      ! 1/N normalization as F_s, so they add up to the total.
      if (self%msd_enabled) then
         allocate (self%msd_partial(self%nspecies, 0:self%nsteps))
         allocate (self%msd_total(0:self%nsteps))
         self%msd_partial = 0.0_rk
         self%msd_total = 0.0_rk
         do j = 0, self%nsteps
            n = int(self%sample_cnt(1, j))
            if (n > 0) then
               do isp = 1, self%nspecies
                  value = self%msd_sum(isp, j)/real(n, rk)/real(self%natoms, rk)
                  self%msd_partial(isp, j) = value
                  self%msd_total(j) = self%msd_total(j) + value
               end do
            end if
         end do
      end if
      if (self%s4_enabled .or. self%chi4_enabled .or. self%fqt_self_enabled .or. &
          self%msd_enabled) then
         allocate (self%sample_tau(0:self%nsteps))
         do j = 0, self%nsteps
           self%sample_tau(j) = real(j*self%stride, rk)*self%frame_dt
        end do
      end if

      ! A shell has one |q| and many directions, and a lattice window averages
      ! a set of shells: both collapse into the single row every writer
      ! expects, weighted by the mode weights the setup left behind.
      if (self%q_mode == dyn_q_shell .or. self%q_mode == dyn_q_powder) &
         call dyn_shell_average(self)
      if (self%q_mode == dyn_q_grid .and. self%s4_enabled) call dyn_isotropy_probe(self)
      ! Every downstream use of a grid wants the isotropic average: the OZ fit
      ! of S4(q,t), the powder average of F(q,t), the isotropic S(q,w).  A run
      ! that needs the individual lattice vectors asks for them with
      ! --keep-modes.
      if (self%q_mode == dyn_q_grid .and. .not. self%keep_modes) call dyn_collapse_shells(self)
      ! A single lattice vector is accumulated at its true q, which is what the
      ! density amplitude needs; its rows are labelled by |q| instead, with the
      ! indices and the vector left to the table header.
      if (self%q_mode == dyn_q_single) self%qvec(:, 1) = [0.0_rk, 0.0_rk, self%qlen(1)]
   end subroutine dyn_prepare_output

   !> Collapse the grid tables onto one row per lattice shell.
   !!
   !! The rows are weighted by the orbit multiplicities divided by how many
   !! members of the orbit survived the thinning.  That is what makes a shell
   !! that holds several lattice orbits come out right: for `|n|^2 = 9` the
   !! 6-fold `(3,0,0)` orbit has to count six times against the 24-fold
   !! `(2,2,1)` one.  After the `+-` reduction every orbit keeps half of its
   !! members, so the weights are uniform and the row is the plain mean; with
   !! `--thin orbits` one representative per orbit survives and it carries
   !! the full multiplicity of its orbit.  Either way the result estimates the
   !! same isotropic shell average.
   !!
   !! This works on the assembled tables, so it costs no extra accumulation.
   !! The Gamma point is a row of its own, with unit weight.
   subroutine dyn_collapse_shells(self)
      class(dynamics_structure_factor_t), intent(inout) :: self
      real(rk), allocatable :: weight(:), wsum(:), row_q(:)
      real(rk), allocatable :: r2a(:, :), r3a(:, :, :)
      real(rk), allocatable :: r1a(:), r1b(:), v1(:)
      real(rk), allocatable :: pnum(:, :, :)
      integer, allocatable :: row_of(:)
      integer(lk), allocatable :: cnt2(:, :)
      integer :: nrow, im, r, prev, nkeep_o

      if (.not. allocated(self%grid_shell)) return
      if (self%nmodes <= 0) return

      ! Row of every mode: the shells are contiguous and sorted by |q|, and the
      ! Gamma point comes first with shell index 0.
      allocate (row_of(self%nmodes), row_q(self%nmodes))
      nrow = 0
      prev = -1
      do im = 1, int(self%nmodes)
         if (self%grid_shell(im) /= prev) then
            nrow = nrow + 1
            prev = self%grid_shell(im)
            row_q(nrow) = self%qlen(im)
         end if
         row_of(im) = nrow
      end do
      if (nrow >= int(self%nmodes)) then
         deallocate (row_of, row_q)
         return
      end if

      ! Weight of every mode: how many lattice vectors its orbit stands for.
      allocate (weight(self%nmodes), wsum(nrow))
      wsum = 0.0_rk
      do im = 1, int(self%nmodes)
         if (self%grid_shell(im) == 0) then
            weight(im) = 1.0_rk
         else
            nkeep_o = count(self%grid_orbit == self%grid_orbit(im))
            weight(im) = real(self%orbit_mult(im), rk)/real(max(nkeep_o, 1), rk)
         end if
         wsum(row_of(im)) = wsum(row_of(im)) + weight(im)
      end do

      ! Collapse every assembled table, then swap it back into the type.
      if (allocated(self%s4)) then
         call collapse_lag_table(self%s4, row_of, weight, wsum, r2a)
         call move_alloc(r2a, self%s4)
      end if
      if (allocated(self%ftau)) then
         call collapse_lag_table(self%ftau, row_of, weight, wsum, r2a)
         call move_alloc(r2a, self%ftau)
      end if
      if (allocated(self%sqw)) then
         call collapse_lag_table(self%sqw, row_of, weight, wsum, r2a)
         call move_alloc(r2a, self%sqw)
      end if
      if (allocated(self%fqt_self_total)) then
         call collapse_lag_table(self%fqt_self_total, row_of, weight, wsum, r2a)
         call move_alloc(r2a, self%fqt_self_total)
      end if
      if (allocated(self%ftau_partial)) then
         call collapse_lag_pair(self%ftau_partial, row_of, weight, wsum, r3a)
         call move_alloc(r3a, self%ftau_partial)
      end if
      if (allocated(self%sqw_partial)) then
         call collapse_lag_pair(self%sqw_partial, row_of, weight, wsum, r3a)
         call move_alloc(r3a, self%sqw_partial)
      end if
      if (allocated(self%fqt_self_partial)) then
         call collapse_lag_pair(self%fqt_self_partial, row_of, weight, wsum, r3a)
         call move_alloc(r3a, self%fqt_self_partial)
      end if

      ! The static table: num is accumulated per mode and den is the single
      ! frame normalization, so both take the same weights and their ratio is
      ! the weighted shell average.
      if (allocated(self%num) .and. allocated(self%den)) then
         allocate (r1a(nrow), r1b(nrow))
         r1a = 0.0_rk
         r1b = 0.0_rk
         do im = 1, int(self%nmodes)
            r = row_of(im)
            r1a(r) = r1a(r) + weight(im)*self%num(im)
            r1b(r) = r1b(r) + weight(im)*self%den(im)
         end do
         call move_alloc(r1a, self%num)
         call move_alloc(r1b, self%den)
      end if
      if (allocated(self%shell_count)) then
         allocate (v1(nrow))
         v1 = 0.0_rk
         do im = 1, int(self%nmodes)
            ! The partial columns divide by this count, so it has to be the
            ! number of lattice vectors the row averages: the weighted sum of
            ! the modes, not their number.  A +- reduced mode stands for two or
            ! more vectors (its `weight` is the orbit share), and the partial
            ! numerator carries the same weights.
            v1(row_of(im)) = v1(row_of(im)) + weight(im)
         end do
         deallocate (self%shell_count)
         allocate (self%shell_count(nrow))
         self%shell_count = int(v1, lk)
         deallocate (v1)
      end if
      if (self%partials .and. allocated(self%partial_num)) then
         allocate (pnum(self%ntypes, self%ntypes, nrow))
         pnum = 0.0_rk
         do im = 1, int(self%nmodes)
            pnum(:, :, row_of(im)) = pnum(:, :, row_of(im)) + weight(im)*self%partial_num(:, :, im)
         end do
         call move_alloc(pnum, self%partial_num)
      end if

      ! Row labels and bookkeeping.
      deallocate (self%qlen)
      allocate (self%qlen(nrow))
      self%qlen = row_q(1:nrow)
      deallocate (row_q)
      do r = 1, nrow
         self%qvec(:, r) = [0.0_rk, 0.0_rk, self%qlen(r)]
         self%shell(r) = r
         self%gidx(r) = int(r, lk)
      end do
      if (allocated(self%shell_radii)) then
         deallocate (self%shell_radii)
         allocate (self%shell_radii(nrow))
         self%shell_radii = self%qlen
      end if
      if (allocated(self%grid_shell)) deallocate (self%grid_shell)
      if (allocated(self%grid_orbit)) deallocate (self%grid_orbit)
      if (allocated(self%orbit_mult)) deallocate (self%orbit_mult)

      ! The origin counts are the same for every mode of a run, but the HDF5
      ! writer groups them per row, so they follow the rows down as well.
      if (allocated(self%sample_cnt)) then
         allocate (cnt2(nrow, 0:size(self%sample_cnt, 2) - 1))
         do r = 1, nrow
            cnt2(r, :) = self%sample_cnt(1, :)
         end do
         call move_alloc(cnt2, self%sample_cnt)
      end if
      if (allocated(self%ccnt)) then
         allocate (cnt2(nrow, 0:size(self%ccnt, 2) - 1))
         do r = 1, nrow
            cnt2(r, :) = self%ccnt(1, :)
         end do
         call move_alloc(cnt2, self%ccnt)
      end if
      self%nmodes = nrow
      self%nq = nrow
      deallocate (row_of, weight, wsum)
   end subroutine dyn_collapse_shells

   !> Weighted row sum of a table with one lag axis.
   subroutine collapse_lag_table(src, row_of, weight, wsum, dst)
      real(rk), intent(in) :: src(:, 0:)
      integer, intent(in) :: row_of(:)
      real(rk), intent(in) :: weight(:), wsum(:)
      real(rk), allocatable, intent(out) :: dst(:, :)
      integer :: im, r, k

      allocate (dst(size(wsum), 0:ubound(src, 2)))
      dst = 0.0_rk
      do im = 1, size(row_of)
         r = row_of(im)
         do k = 0, ubound(src, 2)
            dst(r, k) = dst(r, k) + weight(im)*src(im, k)
         end do
      end do
      do r = 1, size(wsum)
         if (wsum(r) > 0.0_rk) dst(r, :) = dst(r, :)/wsum(r)
      end do
   end subroutine collapse_lag_table

   !> Weighted row sum of a table with a lag axis and a pair/species axis.
   subroutine collapse_lag_pair(src, row_of, weight, wsum, dst)
      real(rk), intent(in) :: src(:, 0:, :)
      integer, intent(in) :: row_of(:)
      real(rk), intent(in) :: weight(:), wsum(:)
      real(rk), allocatable, intent(out) :: dst(:, :, :)
      integer :: im, r, k, p

      allocate (dst(size(wsum), 0:ubound(src, 2), size(src, 3)))
      dst = 0.0_rk
      do im = 1, size(row_of)
         r = row_of(im)
         do p = 1, size(src, 3)
            do k = 0, ubound(src, 2)
               dst(r, k, p) = dst(r, k, p) + weight(im)*src(im, k, p)
            end do
         end do
      end do
      do r = 1, size(wsum)
         if (wsum(r) > 0.0_rk) dst(r, :, :) = dst(r, :, :)/wsum(r)
      end do
   end subroutine collapse_lag_pair

   !> Spread of S4 among the symmetry related modes of the first shell.
   !!
   !! The modes of one lattice shell are related by the point group of the
   !! periodic box, so an equilibrium isotropic system has to give them the
   !! same S4 in expectation.  This is the assumption behind both the orbit
   !! reduction and the shell average, and it is cheap to check: the probe
   !! reports their spread at the lag of the chi4 peak.
   !!
   !! The spread is not an error bar on its own.  It has to be compared with
   !! the run-to-run scatter of a single mode (a second run with another
   !! `--lag`, or the two halves of the trajectory), which is the recipe the
   !! skill reference gives.  A spread that stays far above that scatter, and
   !! that does not shrink with the system size, means the trajectory is not
   !! equilibrated or the system is not isotropic, and then neither the orbit
   !! reduction nor an isotropic correlation length may be trusted.
   subroutine dyn_isotropy_probe(self)
      class(dynamics_structure_factor_t), intent(inout) :: self
      real(rk) :: value, best, vmin, vmax, vsum
      integer :: im, j, jpeak, n, count

      self%probe_valid = .false.
      if (.not. allocated(self%grid_shell)) return
      if (.not. allocated(self%s4)) return

      ! The peak of the scalar overlap fluctuation sets the time scale the
      ! four-point analysis is quoted at.
      jpeak = 0
      best = -huge(1.0_rk)
      do j = 0, self%nsteps
         n = int(self%sample_cnt(1, j))
         if (n <= 1) cycle
         value = (self%chi4_asum(j) - self%chi4_bsum(j)**2/real(n, rk))/real(n - 1, rk)
         if (value > best) then
            best = value
            jpeak = j
         end if
      end do

      vmin = huge(1.0_rk)
      vmax = -huge(1.0_rk)
      vsum = 0.0_rk
      count = 0
      do im = 1, int(self%nmodes)
         if (self%grid_shell(im) /= 1) cycle
         value = self%s4(im, jpeak)
         vmin = min(vmin, value)
         vmax = max(vmax, value)
         vsum = vsum + value
         count = count + 1
      end do
      if (count < 2 .or. vsum <= 0.0_rk) return
      self%probe_valid = .true.
      self%probe_modes = count
      self%probe_tau = self%sample_tau(jpeak)
      self%probe_mean = vsum/real(count, rk)
      self%probe_spread = (vmax - vmin)/self%probe_mean
   end subroutine dyn_isotropy_probe

   !> Collapse the per-mode shell results into one quadrature averaged row.
   !!
   !! Every mode of the shell has the same |q|, so each output is the average
   !! `sum_k w_k X(q_k) / sum_k w_k` over the Lebedev directions and the shell
   !! can share the single-row layout of a one point q line.
   !!
   !! The same collapse serves the lattice-shell window (`dyn_q_powder`), where
   !! the modes are the reciprocal-lattice vectors of a band of shells and
   !! `w_k` is the orbit multiplicity of each of them: the modes the window
   !! excludes carry a zero weight, so the sum is the multiplicity weighted
   !! average over the lattice vectors inside the window.  The row is labelled
   !! with `shell_q`, which the setup sets to the weighted mean |q|.
   !!
   !! The static table is stored as the unnormalized `num(q)` over a per-mode
   !! `den(q)`, but the average is taken over the *normalized* S(q) of every
   !! mode: that is the definition of the shell average, and it is the quantity
   !! F(q,0) and the zeroth moment of S(q,w) refer to.  It also stays correct
   !! should a weight scheme ever make `den` vary within one shell - for the
   !! unit, neutron and X-ray weights `den` depends only on |q|, so here it is
   !! the same for every mode.  The other arrays (F(q,tau), S4 and F_s) are
   !! already normalized per mode, and the spectra are transformed from the
   !! averaged F(q,tau) so that `sum_n S(q,w_n) dw = F(q,0) = S(q)` still holds
   !! exactly.
   subroutine dyn_shell_average(self)
      class(dynamics_structure_factor_t), intent(inout) :: self
      real(rk), allocatable :: keep1(:), keep2(:, :), keep3(:, :, :)
      integer(lk), allocatable :: kj1(:), kj2(:, :)
      integer(ik), allocatable :: ki1(:)
      real(rk) :: wsum, value
      integer :: im, l, j, t, u, p, isp

      wsum = sum(self%mode_weight)

      ! Static table and its partial columns.
      ! `num` is the sum over frames, so num(q)/den(q) is `frames` times S(q);
      ! the writer divides by the frame count, which is therefore not repeated
      ! here.
      value = 0.0_rk
      do im = 1, int(self%nmodes)
         if (self%den(im) > 0.0_rk) value = value &
            + self%mode_weight(im)*self%num(im)/self%den(im)
      end do
      self%num(1) = value/wsum
      self%den(1) = 1.0_rk
      self%shell_count(1) = 1
      if (self%partials .and. allocated(self%partial_num)) then
         do t = 1, self%ntypes
            do u = t, self%ntypes
               value = 0.0_rk
               do im = 1, int(self%nmodes)
                  value = value + self%mode_weight(im)*self%partial_num(t, u, im)
               end do
               self%partial_num(t, u, 1) = value/wsum
            end do
         end do
      end if

      ! F(q,tau) and its partials, followed by the spectra of the average.
      if (allocated(self%ftau)) then
         do l = 0, self%maxframes
            value = 0.0_rk
            do im = 1, int(self%nmodes)
               value = value + self%mode_weight(im)*self%ftau(im, l)
            end do
            self%ftau(1, l) = value/wsum
         end do
         call dyn_spectrum(self%ftau(1, :), self%frame_dt, self%sqw(1, :))
      end if
      if (allocated(self%ftau_partial)) then
         do p = 1, size(self%ftau_partial, 3)
            do l = 0, self%maxframes
               value = 0.0_rk
               do im = 1, int(self%nmodes)
                  value = value + self%mode_weight(im)*self%ftau_partial(im, l, p)
               end do
               self%ftau_partial(1, l, p) = value/wsum
            end do
            call dyn_spectrum(self%ftau_partial(1, :, p), self%frame_dt, &
                              self%sqw_partial(1, :, p))
         end do
      end if

      ! S4(q,tau) and F_s(q,tau) need no normalization by hand.
      if (allocated(self%s4)) then
         do j = 0, self%nsteps
            value = 0.0_rk
            do im = 1, int(self%nmodes)
               value = value + self%mode_weight(im)*self%s4(im, j)
            end do
            self%s4(1, j) = value/wsum
         end do
      end if
      if (allocated(self%fqt_self_total)) then
         do j = 0, self%nsteps
            value = 0.0_rk
            do im = 1, int(self%nmodes)
               value = value + self%mode_weight(im)*self%fqt_self_total(im, j)
            end do
            self%fqt_self_total(1, j) = value/wsum
            do isp = 1, self%nspecies
               value = 0.0_rk
               do im = 1, int(self%nmodes)
                  value = value + self%mode_weight(im)*self%fqt_self_partial(im, j, isp)
               end do
               self%fqt_self_partial(1, j, isp) = value/wsum
            end do
         end do
      end if

      ! One row from here on.  The writers take the size of their tables from
      ! nmodes (the HDF5 ones reshape the result arrays with it), so every
      ! result has to shrink with the mode axis, not just be overwritten.
      self%nmodes = 1
      self%nq = 1
      if (allocated(self%sample_cnt)) then
         call move_alloc(self%sample_cnt, kj2)
         allocate (self%sample_cnt(1, 0:self%nsteps))
         self%sample_cnt(1, :) = kj2(1, :)
      end if
      if (allocated(self%ccnt)) then
         call move_alloc(self%ccnt, kj2)
         allocate (self%ccnt(1, 0:self%maxframes))
         self%ccnt(1, :) = kj2(1, :)
      end if
      call move_alloc(self%num, keep1)
      allocate (self%num(1))
      self%num(1) = keep1(1)
      call move_alloc(self%den, keep1)
      allocate (self%den(1))
      self%den(1) = keep1(1)
      call move_alloc(self%qlen, keep1)
      allocate (self%qlen(1))
      self%qlen(1) = self%shell_q
      call move_alloc(self%mode_values, keep1)
      allocate (self%mode_values(1))
      self%mode_values(1) = keep1(1)
      call move_alloc(self%mode_weight, keep1)
      allocate (self%mode_weight(1))
      self%mode_weight(1) = keep1(1)
      call move_alloc(self%shell_count, kj1)
      allocate (self%shell_count(1))
      self%shell_count(1) = 1
      call move_alloc(self%gidx, kj1)
      allocate (self%gidx(1))
      self%gidx(1) = 1
      call move_alloc(self%shell, ki1)
      allocate (self%shell(1))
      self%shell(1) = 1
      call move_alloc(self%qvec, keep2)
      allocate (self%qvec(3, 1))
      self%qvec(:, 1) = [0.0_rk, 0.0_rk, self%shell_q]
      if (allocated(self%ftau)) then
         call move_alloc(self%ftau, keep2)
         allocate (self%ftau(1, 0:self%maxframes))
         self%ftau(1, :) = keep2(1, :)
      end if
      if (allocated(self%sqw)) then
         call move_alloc(self%sqw, keep2)
         allocate (self%sqw(1, 0:self%maxframes))
         self%sqw(1, :) = keep2(1, :)
      end if
      if (allocated(self%s4)) then
         call move_alloc(self%s4, keep2)
         allocate (self%s4(1, 0:self%nsteps))
         self%s4(1, :) = keep2(1, :)
      end if
      if (allocated(self%fqt_self_total)) then
         call move_alloc(self%fqt_self_total, keep2)
         allocate (self%fqt_self_total(1, 0:self%nsteps))
         self%fqt_self_total(1, :) = keep2(1, :)
      end if
      if (allocated(self%partial_num)) then
         call move_alloc(self%partial_num, keep3)
         allocate (self%partial_num(self%ntypes, self%ntypes, 1))
         self%partial_num(:, :, 1) = keep3(:, :, 1)
      end if
      if (allocated(self%ftau_partial)) then
         p = size(self%ftau_partial, 3)
         call move_alloc(self%ftau_partial, keep3)
         allocate (self%ftau_partial(1, 0:self%maxframes, p))
         ! The array is allocated with a zero pair axis when --no-partials
         ! asks for no columns.
         if (p > 0) self%ftau_partial(1, :, :) = keep3(1, :, :)
      end if
      if (allocated(self%sqw_partial)) then
         p = size(self%sqw_partial, 3)
         call move_alloc(self%sqw_partial, keep3)
         allocate (self%sqw_partial(1, 0:self%maxframes, p))
         if (p > 0) self%sqw_partial(1, :, :) = keep3(1, :, :)
      end if
      if (allocated(self%fqt_self_partial)) then
         call move_alloc(self%fqt_self_partial, keep3)
         allocate (self%fqt_self_partial(1, 0:self%nsteps, self%nspecies))
         self%fqt_self_partial(1, :, :) = keep3(1, :, :)
      end if
      ! The per-mode bookkeeping of a lattice sampling does not describe a
      ! single row: drop it rather than leave arrays whose length the type no
      ! longer agrees on.  Both are unallocated on the Lebedev path, so this
      ! is a no-op there.
      if (allocated(self%shell_radii)) deallocate (self%shell_radii)
      if (allocated(self%grid_shell)) deallocate (self%grid_shell)
      if (allocated(self%grid_orbit)) deallocate (self%grid_orbit)
      if (allocated(self%orbit_mult)) deallocate (self%orbit_mult)
   end subroutine dyn_shell_average

   !> One-sided spectrum of a real, even correlation function.
   !!
   !! `f(l)` is F(q, l*dt) for l = 0..L; the even extension of length 2L is
   !! transformed with the trapezoid rule and the interior points are doubled
   !! (S(q,-w) = S(q,w)), so that sum_n spec(n)*dw = f(0) holds exactly.
   pure subroutine dyn_spectrum(f, dt, spec)
      real(rk), intent(in) :: f(0:)
      real(rk), intent(in) :: dt
      real(rk), intent(out) :: spec(0:)
      real(rk) :: acc, pi
      integer :: l_max, n, l

      pi = acos(-1.0_rk)
      l_max = ubound(f, 1)
      do n = 0, l_max
         acc = f(0)
         do l = 1, l_max - 1
            acc = acc + 2.0_rk*f(l)*cos(pi*real(n*l, rk)/real(l_max, rk))
         end do
         if (l_max > 0) then
            if (mod(n, 2) == 0) then
               acc = acc + f(l_max)
            else
               acc = acc - f(l_max)
            end if
         end if
         spec(n) = acc*dt/(2.0_rk*pi)
      end do
      do n = 1, l_max - 1
         spec(n) = 2.0_rk*spec(n)
      end do
   end subroutine dyn_spectrum

   !> Number of atoms of a LAMMPS type id (0 when the type is absent).
   pure integer function count_pair(self, type_id) result(n)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer, intent(in) :: type_id
      n = 0
      if (allocated(self%type_counts)) then
         if (type_id >= 1 .and. type_id <= size(self%type_counts)) n = int(self%type_counts(type_id))
      end if
   end function count_pair

   !> True when a pair of type ids has the partial column the base writer uses.
   pure logical function pair_present(self, ia, ib) result(present)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer, intent(in) :: ia, ib
      present = count_pair(self, ia) > 0 .and. count_pair(self, ib) > 0
   end function pair_present

   !> Number of type pairs that get a column (`npair` of the base writer).
   pure integer function npairs_out(self) result(n)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer :: t, u
      n = 0
      do t = 1, self%ntypes
         do u = t, self%ntypes
            if (.not. pair_present(self, t, u)) cycle
            n = n + 1
         end do
      end do
   end function npairs_out

   !> Species index of a LAMMPS type id (0 when the type is absent).
   pure integer function species_index(self, type_id) result(idx)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer, intent(in) :: type_id
      idx = 0
      if (allocated(self%species_type)) then
         idx = findloc(self%species_type, int(type_id, ik), dim=1)
      end if
   end function species_index

   !> Short label of one species, e.g. F_s(Si) or F_s(2).
   pure subroutine species_label(scheme, type_id, label, prefix)
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(in) :: type_id
      character(len=*), intent(out) :: label
      character(len=*), intent(in), optional :: prefix
      character(len=8) :: pfx
      character(len=8) :: sym

      pfx = 'F_s'
      if (present(prefix)) pfx = trim(prefix)
      sym = ' '
      if (allocated(scheme%species)) then
         if (type_id >= 1 .and. type_id <= size(scheme%species)) sym = scheme%species(type_id)
      end if
      if (len_trim(sym) > 0) then
         label = trim(pfx)//'('//trim(sym)//')'
      else
         write (label, '(a,i0,a)') trim(pfx)//'(', type_id, ')'
      end if
   end subroutine species_label

   !> The S(q,w) table (text or HDF5).
   subroutine dyn_write_sqw(self, path, format, scheme, input, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input
      integer, intent(in) :: format
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message

      call dyn_write_table(self, path, format, scheme, input, dyn_sqw, 'omega', &
                           self%omega, self%sqw, self%sqw_partial, ierr, message)
   end subroutine dyn_write_sqw

   !> The F(q,tau) table (text or HDF5).
   subroutine dyn_write_fqt(self, path, format, scheme, input, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input
      integer, intent(in) :: format
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message

      call dyn_write_table(self, path, format, scheme, input, dyn_fqt, 'tau', &
                           self%tau, self%ftau, self%ftau_partial, ierr, message)
   end subroutine dyn_write_fqt

   !> The self intermediate scattering function F_s(q,tau).
   subroutine dyn_write_fqt_self(self, path, format, scheme, input, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input
      integer, intent(in) :: format
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: q4(:, :)
      character(len=18), allocatable :: labels(:)
      character(len=18) :: slabel
      integer :: unit, im, j, isp

      ierr = 0
      message = ''
      if (format == format_hdf5) then
#ifdef SQC_HAVE_HDF5
         allocate (q4(self%nmodes, 4), labels(self%nspecies))
         do im = 1, int(self%nmodes)
            q4(im, 1:3) = self%qvec(:, im)
            q4(im, 4) = self%qlen(im)
         end do
         do isp = 1, self%nspecies
            call species_label(scheme, int(self%species_type(isp)), labels(isp))
         end do
         call hdf5_write_fqt_self(path, q4, self%sample_tau, self%fqt_self_total, &
                                  self%fqt_self_partial, labels, self%sample_cnt, self%nframes, &
                                  self%frame_dt, self%maxframes, self%lag_stride, self%stride, &
                                  self%effective_maxframes, ierr, message)
         deallocate (q4, labels)
#else
         ierr = 1
         message = 'this build has no HDF5 support; use --format text'
#endif
         return
      end if

      if (trim(path) == '-') then
         unit = output_unit
      else
         open (newunit=unit, file=trim(path), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            message = 'cannot write the F_s output "'//trim(path)//'"'
            return
         end if
      end if

      write (unit, '(a)') '# sqcalc 0.1.0 self intermediate scattering function F_s(q,t)'
      write (unit, '(a)') '# input '//trim(input)//'  weight unit (self)  norm unit'
      write (unit, '(a,f12.6,a,i0,a,i0,a,i0)') '# frame_dt ', self%frame_dt, &
         '  maxframes ', self%maxframes, '  lag ', self%lag_stride, '  nframes ', self%nframes
      call dyn_write_sampling(self, unit)
      write (unit, '(a,i0,a,f0.6,a)') '# nlag ', self%nsteps + 1, &
         '  tau 0 .. ', real(self%effective_maxframes, rk)*self%frame_dt, ' (time unit of dt)'
      write (unit, '(a,i0,a,i0,a,i0,a)') '# stride ', self%stride, &
         ' dump frames; effective maxframes ', self%effective_maxframes, ' (requested ', &
         self%maxframes, ')'
      write (unit, '(a)', advance='no') '# qx qy qz tau F_s(q,t)'
      do isp = 1, self%nspecies
         call species_label(scheme, int(self%species_type(isp)), slabel)
         write (unit, '(a)', advance='no') ' '//trim(slabel)
      end do
      write (unit, '(a)') ''
      do im = 1, int(self%nmodes)
         do j = 0, self%nsteps
            write (unit, '(3(f14.8,2x),f16.8,2x,es20.12)', advance='no') &
               self%qvec(1, im), self%qvec(2, im), self%qvec(3, im), self%sample_tau(j), &
               self%fqt_self_total(im, j)
            do isp = 1, self%nspecies
               write (unit, '(2x,es20.12)', advance='no') self%fqt_self_partial(im, j, isp)
            end do
            write (unit, '(a)') ''
         end do
      end do
      if (unit /= output_unit) close (unit)
   end subroutine dyn_write_fqt_self

   !> The total S4(q,tau) table (text or HDF5).
   subroutine dyn_write_s4(self, path, format, scheme, input, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input
      integer, intent(in) :: format
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: empty_part(:, :, :)

      allocate (empty_part(0, 0, 0))
      call dyn_write_table(self, path, format, scheme, input, dyn_s4, 'tau', &
                           self%sample_tau, self%s4, empty_part, ierr, message)
      deallocate (empty_part)
   end subroutine dyn_write_s4

   !> The average overlap Q(t) and the dynamic susceptibility chi4(t).
   subroutine dyn_write_chi4(self, path, format, input, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input
      integer, intent(in) :: format
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer :: unit, l

      ierr = 0
      message = ''
      if (format == format_hdf5) then
#ifdef SQC_HAVE_HDF5
         call hdf5_write_chi4(path, self%sample_tau, self%overlap, self%chi4, self%sample_cnt(1, :), &
                              self%nframes, self%frame_dt, self%maxframes, self%lag_stride, &
                              self%s4_cutoff, self%stride, self%effective_maxframes, ierr, message)
#else
         ierr = 1
         message = 'this build has no HDF5 support; use --format text'
#endif
         return
      end if

      if (trim(path) == '-') then
         unit = output_unit
      else
         open (newunit=unit, file=trim(path), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            message = 'cannot write the chi4 output "'//trim(path)//'"'
            return
         end if
      end if

      write (unit, '(a)') '# sqcalc 0.1.0 average overlap Q(t) and four-point '// &
         'susceptibility chi4(t)'
      write (unit, '(a)') '# input '//trim(input)
      write (unit, '(a,f12.6,a,f12.6,a,i0,a,i0,a,i0,a,i0)') '# overlap ', self%s4_cutoff, &
         '  frame_dt ', self%frame_dt, '  stride ', self%stride, '  maxframes ', &
         self%effective_maxframes, '  lag ', self%lag_stride, '  nframes ', self%nframes
      write (unit, '(a)') '# the overlap uses unit weights; Q and chi4 are the q = 0 limit'
      write (unit, '(a)') '# tau Q(t) chi4(t)'
      do l = 0, self%nsteps
         write (unit, '(f16.8,2x,es20.12,2x,es20.12)') self%sample_tau(l), self%overlap(l), self%chi4(l)
      end do
      if (unit /= output_unit) close (unit)
   end subroutine dyn_write_chi4

   !> Mean squared displacement MSD(tau) and its per-species split.
   subroutine dyn_write_msd(self, path, format, scheme, input, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input
      integer, intent(in) :: format
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=18), allocatable :: labels(:)
      character(len=18) :: slabel
      integer :: unit, j, isp

      ierr = 0
      message = ''
      if (format == format_hdf5) then
#ifdef SQC_HAVE_HDF5
         allocate (labels(self%nspecies))
         do isp = 1, self%nspecies
            call species_label(scheme, int(self%species_type(isp)), labels(isp), 'MSD')
         end do
         call hdf5_write_msd(path, self%sample_tau, self%msd_total, self%msd_partial, &
                             labels, self%sample_cnt(1, :), self%partials, self%nframes, &
                             self%frame_dt, self%maxframes, self%lag_stride, self%stride, &
                             self%effective_maxframes, ierr, message)
         deallocate (labels)
#else
         ierr = 1
         message = 'this build has no HDF5 support; use --format text'
#endif
         return
      end if

      if (trim(path) == '-') then
         unit = output_unit
      else
         open (newunit=unit, file=trim(path), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            message = 'cannot write the MSD output "'//trim(path)//'"'
            return
         end if
      end if

      write (unit, '(a)') '# sqcalc 0.1.0 mean squared displacement MSD(t)'
      write (unit, '(a)') '# input '//trim(input)//'  weight unit (self)  norm unit'
      write (unit, '(a,f12.6,a,i0,a,i0,a,i0)') '# frame_dt ', self%frame_dt, &
         '  maxframes ', self%maxframes, '  lag ', self%lag_stride, '  nframes ', self%nframes
      write (unit, '(a,i0,a,f0.6,a)') '# nlag ', self%nsteps + 1, &
         '  tau 0 .. ', real(self%effective_maxframes, rk)*self%frame_dt, ' (time unit of dt)'
      write (unit, '(a,i0,a,i0,a,i0,a)') '# stride ', self%stride, &
         ' dump frames; effective maxframes ', self%effective_maxframes, ' (requested ', &
         self%maxframes, ')'
      write (unit, '(a)') '# MSD(t) = N^-1 <sum_i |r_i(t0+t) - r_i(t0)|^2>_t0 '// &
         '(unit weights); one species column, the total is their sum'
      write (unit, '(a)', advance='no') '# tau MSD(t)'
      if (self%partials) then
         do isp = 1, self%nspecies
            call species_label(scheme, int(self%species_type(isp)), slabel, 'MSD')
            write (unit, '(a)', advance='no') ' '//trim(slabel)
         end do
      end if
      write (unit, '(a)') ''
      do j = 0, self%nsteps
         write (unit, '(f16.8,2x,es20.12)', advance='no') self%sample_tau(j), self%msd_total(j)
         if (self%partials) then
            do isp = 1, self%nspecies
               write (unit, '(2x,es20.12)', advance='no') self%msd_partial(isp, j)
            end do
         end if
         write (unit, '(a)') ''
      end do
      if (unit /= output_unit) close (unit)
   end subroutine dyn_write_msd

   !> Text or HDF5 writer shared by the spectra and F(q,tau).
   subroutine dyn_write_table(self, path, format, scheme, input, kind, axis_name, &
                              axis, spec, part, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      character(len=*), intent(in) :: path, input, axis_name
      integer, intent(in) :: format, kind
      type(weight_scheme_t), intent(in) :: scheme
      real(rk), intent(in) :: axis(0:), spec(:, 0:), part(:, 0:, :)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer :: unit, im, n, p, t, u

      ierr = 0
      message = ''
      if (format == format_hdf5) then
#ifdef SQC_HAVE_HDF5
         call dyn_write_hdf5(self, scheme, path, kind, axis_name, axis, spec, part, ierr, message)
#else
         ierr = 1
         message = 'this build has no HDF5 support; write the dynamic tables as text'
#endif
         return
      end if

      if (trim(path) == '-') then
         unit = output_unit
      else
         open (newunit=unit, file=trim(path), status='replace', action='write', iostat=ierr)
         if (ierr /= 0) then
            message = 'cannot write the dynamic output "'//trim(path)//'"'
            return
         end if
      end if

      call dyn_header(self, scheme, input, kind, axis_name, unit)
      do im = 1, int(self%nmodes)
         do n = 0, ubound(axis, 1)
            write (unit, '(3(f14.8,2x),f16.8,2x,es20.12)', advance='no') &
               self%qvec(1, im), self%qvec(2, im), self%qvec(3, im), axis(n), spec(im, n)
            if (self%partials .and. kind /= dyn_s4) then
               p = 0
               do t = 1, self%ntypes
                  do u = t, self%ntypes
                     if (.not. pair_present(self, t, u)) cycle
                     p = p + 1
                     write (unit, '(2x,es20.12)', advance='no') part(im, n, p)
                  end do
               end do
            end if
            write (unit, '(a)') ''
         end do
      end do
      if (unit /= output_unit) close (unit)
   end subroutine dyn_write_table

   !> Metadata preamble of the text tables.
   subroutine dyn_header(self, scheme, input, kind, axis_name, unit)
      class(dynamics_structure_factor_t), intent(in) :: self
      type(weight_scheme_t), intent(in) :: scheme
      character(len=*), intent(in) :: input, axis_name
      integer, intent(in) :: kind, unit
      integer :: t, u
      character(len=24) :: label

      select case (kind)
      case (dyn_sqw)
         write (unit, '(a)') '# sqcalc 0.1.0 dynamic structure factor S(q,w)'
      case (dyn_fqt)
         write (unit, '(a)') '# sqcalc 0.1.0 intermediate scattering function F(q,t)'
      case default
         write (unit, '(a)') '# sqcalc 0.1.0 four-point structure factor S4(q,t)'
      end select
      if (kind == dyn_s4) then
         write (unit, '(a)') '# input '//trim(input)//'  weight unit (overlap)  norm unit'
      else
         write (unit, '(a)') '# input '//trim(input)//'  weight '//trim(scheme%label())// &
            '  norm '//trim(norm_name(self%norm))
      end if
      if (self%q_dependent .and. kind /= dyn_s4) then
         write (unit, '(a)') '# note the amplitudes depend on q (X-ray form factors)'
      end if
      write (unit, '(a,f12.6,a,i0,a,i0,a,i0)') '# frame_dt ', self%frame_dt, &
         '  maxframes ', self%maxframes, '  lag ', self%lag_stride, '  nframes ', self%nframes
      call dyn_write_sampling(self, unit)
      select case (kind)
      case (dyn_sqw)
         write (unit, '(a,i0,a,f0.6,a,f0.6,a)') '# nomega ', self%maxframes + 1, &
            '  omega 0 .. ', real(self%maxframes, rk)*self%dw, ' dw ', self%dw, &
            ' (1/time unit of dt, one-sided)'
         write (unit, '(a)') '# S(q,-w) = S(q,w);  sum_n S(q,w_n)*dw = F(q,0) '// &
            '= the OUTPUT table of this run'
      case (dyn_fqt)
         write (unit, '(a,i0,a,f0.6,a)') '# nlag ', self%maxframes + 1, &
            '  tau 0 .. ', real(self%maxframes, rk)*self%frame_dt, ' (time unit of dt)'
         write (unit, '(a)') '# F(q,0) = the OUTPUT table of this run; the spectra are '// &
            'its Fourier transform'
      case default
         write (unit, '(a,i0,a,f0.6,a)') '# nlag ', self%nsteps + 1, &
            '  tau 0 .. ', real(self%effective_maxframes, rk)*self%frame_dt, ' (time unit of dt)'
         write (unit, '(a,i0,a,i0,a,i0,a)') '# stride ', self%stride, &
            ' dump frames; effective maxframes ', self%effective_maxframes, ' (requested ', &
            self%maxframes, ')'
         write (unit, '(a,f12.6,a)') '# overlap cutoff a = ', self%s4_cutoff, &
            ' (same length unit as the dump)'
         write (unit, '(a)') '# S4 is the unbiased connected overlap covariance; '// &
            'lim_q->0 S4 = chi4(t)'
      end select
      write (unit, '(a)', advance='no') '# qx qy qz '//trim(axis_name)//' '
      select case (kind)
      case (dyn_sqw)
         write (unit, '(a)', advance='no') 'S(q,w)'
      case (dyn_fqt)
         write (unit, '(a)', advance='no') 'F(q,t)'
      case default
         write (unit, '(a)', advance='no') 'S4(q,t)'
      end select
      if (self%partials .and. kind /= dyn_s4) then
         do t = 1, self%ntypes
            do u = t, self%ntypes
               if (.not. pair_present(self, t, u)) cycle
               label = self%partial_column(t, u)
               write (unit, '(a)', advance='no') ' '//trim(label)
            end do
         end do
      end if
      write (unit, '(a)') ''
   end subroutine dyn_header

   !> The q sampling preamble shared by the dynamic text tables.
   subroutine dyn_write_sampling(self, unit)
      class(dynamics_structure_factor_t), intent(in) :: self
      integer, intent(in) :: unit
      real(rk) :: unorm

      if (self%q_mode == dyn_q_shell) then
         write (unit, '(a,f12.6,a,i0,a,i0,a,i0,a)') '# q shell |q| = ', self%shell_q, &
            ' 1/A: Lebedev average over every direction (order ', self%shell_order, ', ', &
            lebedev_points(self%shell_order), ' points, ', self%nmodes, &
            ' after the +- reduction)'
         return
      end if
      if (self%q_mode == dyn_q_powder) then
         write (unit, '(a,f0.6,a,f0.6,a,i0,a,i0,a)') '# q powder shell: window ', &
            self%powder_q - self%powder_dq_used, ' ~ ', self%powder_q + self%powder_dq_used, &
            ' 1/A holds ', self%powder_vectors, ' lattice vectors in ', &
            self%powder_nshell, ' shells'
         write (unit, '(a,f0.6,a,f0.6,a,f0.6,a,f0.6,a)') '# requested |q| ', self%powder_q, &
            ', mean |q| ', self%powder_qmean, ', dq ', self%powder_dq_used, ', offset ', &
            self%powder_qmean - self%powder_q, ' 1/A'
         if (self%powder_nearest) then
            write (unit, '(a)') '# the window held no shell of its own; the nearest one was used'
         end if
         return
      end if
      if (self%q_mode == dyn_q_grid) then
         if (self%keep_modes) then
            write (unit, '(a,f10.6,a,i0,a,i0,a,i0,a)') '# q grid |q| <= ', self%grid_qmax, &
               ' 1/A: one row per reciprocal lattice vector, ', self%nmodes, ' modes in ', &
               self%grid_nshell_kept, ' shells, lattice point group of order ', self%grid_nops
         else
            write (unit, '(a,f10.6,a,i0,a,i0,a,i0,a)') '# q grid |q| <= ', self%grid_qmax, &
               ' 1/A: isotropic average per lattice shell, ', self%nmodes, ' shells from ', &
               self%grid_nmodes, ' lattice vectors, point group of order ', self%grid_nops
            write (unit, '(a)') '# each row is the shell average over the reciprocal lattice '// &
               'vectors of that |q|, weighted by their orbit multiplicities'
         end if
         if (self%grid_thinned) then
            write (unit, '(a)') '# note: the mode budget thinned the grid, see the run summary'
         end if
         return
      end if
      if (self%q_mode == dyn_q_single) then
         write (unit, '(a,3(i0,1x),a,3(f10.6,1x),a,f10.6,a)') '# q single n = (', &
            self%single_index, ') -> q = (', self%single_q, ') 1/A, |q| = ', self%qlen(1), ' 1/A'
         return
      end if
      unorm = sqrt(sum(self%direction**2))
      write (unit, '(a,i0,a,f10.6,a,f10.6,a)') '# q line ', self%nintervals, &
         ' intervals: |q| ', self%s0, ' .. ', self%s1, ' 1/A (through Gamma)'
      write (unit, '(a,3(f12.8,1x))') '# direction ', self%direction/unorm
   end subroutine dyn_write_sampling

#ifdef SQC_HAVE_HDF5
   !> HDF5 flavour: one group per quantity, mirroring the shell/rdf writers.
   subroutine dyn_write_hdf5(self, scheme, path, kind, axis_name, axis, spec, part, ierr, message)
      class(dynamics_structure_factor_t), intent(in) :: self
      type(weight_scheme_t), intent(in) :: scheme
      character(len=*), intent(in) :: path, axis_name
      integer, intent(in) :: kind
      real(rk), intent(in) :: axis(0:), spec(:, 0:), part(:, 0:, :)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk), allocatable :: q4(:, :)
      character(len=18), allocatable :: labels(:)
      character(len=:), allocatable :: group
      character(len=32) :: wlabel, nlabel
      character(len=8) :: sampling
      real(rk) :: overlap
      integer(lk), allocatable :: count_use(:, :)
      integer :: im, p, t, u, npair, maxframes_use, stride_use, effective_maxframes_use

      ierr = 0
      message = ''
      overlap = -1.0_rk
      select case (self%q_mode)
      case (dyn_q_grid)
         sampling = 'grid'
      case (dyn_q_powder)
         sampling = 'powder'
      case (dyn_q_single)
         sampling = 'single'
      case (dyn_q_shell)
         sampling = 'shell'
      case default
         sampling = 'line'
      end select
      maxframes_use = self%maxframes
      stride_use = -1
      effective_maxframes_use = -1
      select case (kind)
      case (dyn_fqt)
         group = 'fqt'
         wlabel = trim(scheme%label())
         nlabel = trim(norm_name(self%norm))
      case (dyn_s4)
         group = 's4'
         wlabel = 'unit'
         nlabel = 'unit'
         overlap = self%s4_cutoff
         stride_use = self%stride
         effective_maxframes_use = self%effective_maxframes
      case default
         group = 'sqw'
         wlabel = trim(scheme%label())
         nlabel = trim(norm_name(self%norm))
      end select
      allocate (q4(self%nmodes, 4))
      do im = 1, int(self%nmodes)
         q4(im, 1:3) = self%qvec(:, im)
         q4(im, 4) = self%qlen(im)
      end do
      npair = npairs_out(self)
      allocate (labels(max(npair, 1)))
      labels = ' '
      p = 0
      do t = 1, self%ntypes
         do u = t, self%ntypes
            if (.not. pair_present(self, t, u)) cycle
            p = p + 1
            labels(p) = scheme%pair_label(t, u)
         end do
      end do
      if (kind == dyn_s4) then
         allocate (count_use(size(self%sample_cnt, 1), 0:size(self%sample_cnt, 2) - 1))
         count_use = self%sample_cnt
      else
         allocate (count_use(size(self%ccnt, 1), 0:size(self%ccnt, 2) - 1))
         count_use = self%ccnt
      end if
      call hdf5_write_dynamics(path, group, axis_name, q4, axis, spec, part, labels, count_use, &
                               self%nframes, self%frame_dt, maxframes_use, self%lag_stride, &
                               trim(wlabel), trim(nlabel), overlap, stride_use, &
                               effective_maxframes_use, trim(sampling), self%grid_budget, &
                               self%grid_thinned, self%single_index, self%powder_q, &
                               self%powder_dq_used, self%powder_vectors, ierr, message)
      deallocate (q4, labels, count_use)
   end subroutine dyn_write_hdf5
#endif

   !> Text label of a normalization code.
   pure function norm_name(norm) result(label)
      integer, intent(in) :: norm
      character(len=8) :: label
      select case (norm)
      case (norm_self)
         label = 'self'
      case (norm_natom)
         label = 'n'
      case default
         label = 'mean'
      end select
   end function norm_name

end module sqc_dynamics
