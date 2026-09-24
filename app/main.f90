! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> sqcalc: total structure factors from LAMMPS dump trajectories.
program sqcalc
   use sqc_kinds
   use sqc_options
   use sqc_output, only: format_hdf5
   use sqc_dump
   use sqc_weights
   use sqc_cell, only: cell_t
   use sqc_elements, only: element_table
   use sqc_structure_factor
   use sqc_dynamics, only: dynamics_structure_factor_t, &
                           dyn_q_none, dyn_q_line, dyn_q_shell, dyn_q_grid, dyn_q_single, &
                           dyn_q_powder
   use sqc_debye, only: debye_structure_factor_t
   use sqc_finufft, only: finufft_opts_is_consistent
#ifdef SQC_ENABLE_CUDA
   use sqc_gpu, only: cufinufft_structure_factor_t
#endif
#ifdef SQC_HAVE_HDF5
   use sqc_hdf5, only: hdf5_write_results, hdf5_write_xrd
#endif
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none

   type(options_t) :: opts
   class(frame_reader_t), allocatable :: reader
   class(structure_factor_t), allocatable :: method
   type(frame_t) :: frame
   character(len=512) :: message, accum
   integer :: ierr, shell_unit, grid_unit, xrd_unit, tick, tick_rate, progress_step
   integer(ik) :: natoms, step
   real(rk) :: ref_a(3, 3), elapsed, wall0, wall1, eps_used
   logical :: first_frame

   call parse_options(opts, ierr, message)
   if (opts%show_help) then
      call print_usage(output_unit, opts%command)
      stop
   end if
   if (opts%show_version) then
      write (output_unit, '(a)') program_version
      stop
   end if
   if (ierr /= 0) then
      write (error_unit, '(a)') 'sqcalc: '//trim(message)
      call print_usage(error_unit, opts%command)
      stop 2
   end if

   ! --- thread setup ----------------------------------------------------
   if (opts%threads <= 0) opts%threads = omp_get_max_threads()
   call omp_set_num_threads(opts%threads)

   if (opts%method == method_nufft .and. .not. finufft_opts_is_consistent()) then
      write (error_unit, '(a)') 'sqcalc: the linked FINUFFT library uses an incompatible '// &
         'finufft_opts layout; rebuild against the vendored FINUFFT version'
      stop 3
   end if

   ! --- open the trajectory and inspect the first frame -----------------
   allocate (lammps_dump_reader_t :: reader)
   select type (reader)
   type is (lammps_dump_reader_t)
      call reader%open(opts%input, ierr)
   end select
   if (ierr /= 0) then
      write (error_unit, '(a)') 'sqcalc: cannot open dump file "'//opts%input//'"'
      stop 4
   end if

   call reader%next_frame(frame, ierr)
   if (ierr /= 0) then
      write (error_unit, '(a,i0,a)') 'sqcalc: no complete frame found in "'//opts%input// &
         '" (reader error ', ierr, ')'
      stop 5
   end if
   natoms = frame%natoms
   ref_a = frame%cell%a
   call validate_mapping(opts%scheme, frame%type_id)

   if (.not. opts%quiet) call report_setup(opts, frame)

   ! --- create the requested method -------------------------------------
   ! The subcommand picks the method family; --method picks the variant inside
   ! the static family (the dyn subcommand has no --method).
   select case (opts%command)
   case (cmd_dyn)
      allocate (dynamics_structure_factor_t :: method)
      select type (method)
      type is (dynamics_structure_factor_t)
         method%q_mode = opts%q_mode
         method%nintervals = opts%q_intervals
         method%s0 = opts%q_s0
         method%s1 = opts%q_s1
         method%direction = opts%q_dir
         method%shell_q = opts%shell_q
         method%shell_order = opts%shell_order
         method%grid_qmax = opts%grid_qmax
         method%grid_budget = opts%modes
         method%grid_thin = opts%thin
         method%powder_q = opts%powder_q
         method%powder_modes = opts%powder_modes
         method%powder_dq = opts%powder_dq
         method%powder_dq_given = opts%powder_dq_given
         method%keep_modes = opts%keep_modes
         method%single_index = opts%single_index
         method%maxframes = opts%maxframes
         method%lag_stride = opts%lag_stride
         method%sqw_format = opts%sqw_format
         method%fqt_format = opts%fqt_format
         method%fqt_self_format = opts%fqt_self_format
         if (allocated(opts%sqw_output)) method%sqw_path = opts%sqw_output
         if (allocated(opts%fqt_output)) method%fqt_path = opts%fqt_output
         if (allocated(opts%fqt_self_output)) then
            method%fqt_self_path = opts%fqt_self_output
            method%fqt_self_enabled = .true.
         end if
         method%coherent_enabled = allocated(opts%sqw_output) .or. allocated(opts%fqt_output)
         method%s4_cutoff = opts%s4_cutoff
         method%buffer_limit_gb = opts%buffer_limit_gb
         method%stride = opts%stride
         method%s4_format = opts%s4_format
         method%chi4_format = opts%chi4_format
         method%msd_format = opts%msd_format
         method%ngp_format = opts%ngp_format
         if (allocated(opts%s4_output)) then
            method%s4_path = opts%s4_output
            method%s4_enabled = .true.
         end if
         if (allocated(opts%chi4_output)) then
            method%chi4_path = opts%chi4_output
            method%chi4_enabled = .true.
         end if
         if (allocated(opts%msd_output)) then
            method%msd_path = opts%msd_output
            method%msd_enabled = .true.
         end if
         if (allocated(opts%ngp_output)) then
            method%ngp_path = opts%ngp_output
            method%ngp_enabled = .true.
         end if
      end select
   case (cmd_static)
      select case (opts%method)
      case (method_direct)
         allocate (direct_structure_factor_t :: method)
      case (method_debye)
         allocate (debye_structure_factor_t :: method)
         select type (method)
         type is (debye_structure_factor_t)
            method%rmax = opts%rmax
            method%dr = opts%dr
            method%skin = opts%skin
            method%correct_cutoff = .not. opts%no_cutoff_correction
            ! --format overrides the suffix rule of g(r) and the pair entropy.
            if (opts%output_format_given) method%output_format = opts%output_format
            if (allocated(opts%rdf_output)) method%rdf_path = opts%rdf_output
            if (allocated(opts%pair_entropy_output)) method%pair_entropy_path = &
               opts%pair_entropy_output
            if (allocated(opts%s2_accum_output)) method%s2_accum_path = opts%s2_accum_output
         end select
      case default
         if (opts%device == device_gpu) then
#ifdef SQC_ENABLE_CUDA
            allocate (cufinufft_structure_factor_t :: method)
            select type (method)
            type is (cufinufft_structure_factor_t)
               method%device = opts%gpu_id
               method%single_precision = opts%precision == precision_single
            end select
#else
            write (error_unit, '(a)') 'sqcalc: this build has no GPU support; '// &
               'configure with -DSQC_ENABLE_CUDA=ON'
            stop 16
#endif
         else
            allocate (nufft_structure_factor_t :: method)
         end if
      end select
   end select
   method%norm = opts%norm
   method%nq = opts%nq
   method%qmin = opts%qmin
   method%qmax = opts%qmax
   eps_used = opts%eps
   if (opts%device == device_gpu .and. opts%precision == precision_single .and. &
       .not. opts%eps_given) then
      ! float32 cannot reach the double precision default; 1e-5 is a good and
      ! fast choice for S(q) (cufinufft itself would clamp anything smaller).
      eps_used = 1.0e-5_rk
   end if
   method%eps = eps_used
   method%nthreads = opts%threads
   method%want_grid = opts%want_grid
   method%partials = opts%partials
   method%faber_ziman = opts%faber_ziman
   method%xrd_enabled = allocated(opts%xrd_output)
   method%xrd_lambda = opts%xrd_lambda
   method%xrd_2theta_min = opts%xrd_2theta_min
   method%xrd_2theta_max = opts%xrd_2theta_max
   method%xrd_step = opts%xrd_step
   method%xrd_lp = opts%lp

   call method%configure(frame, opts%scheme, ierr, message)
   if (ierr /= 0) then
      write (error_unit, '(a)') 'sqcalc: '//trim(message)
      stop 6
   end if
   if (.not. opts%quiet) call report_grid(method)

   ! --- accumulate over the trajectory ----------------------------------
   call system_clock(tick, tick_rate)
   wall0 = real(tick, rk)/real(tick_rate, rk)
   progress_step = max(1, 200)
   first_frame = .true.
   do
      if (.not. first_frame) then
         call reader%next_frame(frame, ierr)
         if (ierr == 1) exit
         if (ierr == 10) then
            write (error_unit, '(a)') 'sqcalc: the atom order (id) changes between frames; '// &
               'sort the dump (for example dump_modify ... sort id) or keep the order fixed'
            stop 24
         end if
         if (ierr == 11) then
            write (error_unit, '(a)') 'sqcalc: the number of atoms changes between frames'
            stop 25
         end if
         if (ierr /= 0) then
            write (error_unit, '(a,i0,a,i0,a)') 'sqcalc: malformed dump near line ', &
               reader_line(reader), ' (reader error ', ierr, ')'
            stop 7
         end if
         if (frame%natoms /= natoms) then
            write (error_unit, '(a)') 'sqcalc: number of atoms changes between frames'
            stop 8
         end if
         if (maxval(abs(frame%cell%a - ref_a)) > 1.0e-6_rk*max(1.0_rk, maxval(abs(ref_a)))) then
            write (error_unit, '(a)') 'sqcalc: the simulation box changes between frames'
            stop 9
         end if
      end if
      first_frame = .false.

      call method%accumulate_frame(frame, opts%scheme, ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a)') 'sqcalc: '//trim(message)
         stop 10
      end if
      if (.not. opts%quiet .and. mod(method%nframes, int(progress_step, lk)) == 0) then
         write (error_unit, '(a,i0,a)') '  processed ', method%nframes, ' frames'
      end if
   end do
   call reader%close()

   ! --- the dynamic method needs a constant time between frames ----------
   if (opts%command == cmd_dyn) then
      step = reader_step_stride(reader)
      if (.not. reader_uniform_steps(reader) .or. step <= 0) then
         write (error_unit, '(a)') 'sqcalc: the dump must sample the trajectory at a '// &
            'constant timestep interval for the dyn subcommand'
         stop 21
      end if
      select type (method)
      type is (dynamics_structure_factor_t)
         method%frame_dt = opts%dt*real(step, rk)
         if (.not. opts%quiet) then
            write (error_unit, '(a,f0.6,a,f0.6,a,i0,a)') '  frame dt   : ', method%frame_dt, &
               ' = dt ', opts%dt, ' x ', step, ' steps'
            write (error_unit, '(a,f0.4,a)') '  omega max  : ', &
               acos(-1.0_rk)/method%frame_dt, ' 1/time (Nyquist, from the dump interval)'
         end if
      end select
      if (reader_ambiguous(reader) > 0 .and. .not. opts%quiet) then
         write (error_unit, '(a,i0,a)') '  note       : ', reader_ambiguous(reader), &
            ' atom displacements were too large to unwrap reliably; dump more often'
      end if
   end if

   call method%prepare_output(opts%scheme, ierr, message)
   if (ierr /= 0) then
      write (error_unit, '(a)') 'sqcalc: '//trim(message)
      stop 19
   end if
   call system_clock(tick, tick_rate)
   wall1 = real(tick, rk)/real(tick_rate, rk)
   elapsed = wall1 - wall0

   ! --- write the results ------------------------------------------------
#ifndef SQC_HAVE_HDF5
   if (opts%want_grid .and. opts%grid_format == format_hdf5) then
      write (error_unit, '(a)') 'sqcalc: this build has no HDF5 support; '// &
         'configure with -DSQC_ENABLE_HDF5=ON or use --format text'
      stop 17
   end if
   if (allocated(opts%xrd_output) .and. opts%xrd_format == format_hdf5) then
      write (error_unit, '(a)') 'sqcalc: this build has no HDF5 support; '// &
         'configure with -DSQC_ENABLE_HDF5=ON or name the XRD file .txt'
      stop 17
   end if
#endif
   if (opts%method == method_debye .and. allocated(opts%rdf_output)) then
      select type (method)
      type is (debye_structure_factor_t)
         call method%write_rdf(opts%scheme, opts%rdf_output, ierr, message)
      end select
      if (ierr /= 0) then
         write (error_unit, '(a)') 'sqcalc: '//trim(message)
         stop 20
      end if
   end if
   if (opts%method == method_debye .and. allocated(opts%pair_entropy_output)) then
      select type (method)
      type is (debye_structure_factor_t)
         accum = ''
         if (allocated(method%s2_accum_path)) accum = method%s2_accum_path
         call method%write_pair_entropy(opts%scheme, method%pair_entropy_path, trim(accum), &
                                        opts%input, ierr, message)
      end select
      if (ierr /= 0) then
         write (error_unit, '(a)') 'sqcalc: '//trim(message)
         stop 30
      end if
   end if
   ! A run without --q has no S(q) table to write.
   shell_unit = no_unit
   if (allocated(opts%output)) then
      call open_output(opts%output, shell_unit, ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a)') 'sqcalc: '//trim(message)
         stop 11
      end if
      ! A shell average is the single row at |q| = Q; label it before the
      ! shared table header names the columns.
      if (opts%command == cmd_dyn .and. opts%q_mode == dyn_q_shell) then
         select type (method)
         type is (dynamics_structure_factor_t)
            write (shell_unit, '(a,f12.6,a)') '# shell average over |q| = ', &
               method%shell_q, ' 1/A (Lebedev quadrature, one row)'
         end select
      end if
      if (opts%command == cmd_dyn .and. opts%q_mode == dyn_q_powder) then
         select type (method)
         type is (dynamics_structure_factor_t)
            write (shell_unit, '(a,f0.6,a,f0.6,a,i0,a,i0,a,f0.6,a)') &
               '# powder shell average: window ', method%powder_q - method%powder_dq_used, &
               ' ~ ', method%powder_q + method%powder_dq_used, ' 1/A, ', &
               method%powder_vectors, ' lattice vectors in ', method%powder_nshell, &
               ' shells, mean |q| = ', method%powder_qmean, ' 1/A'
         end select
      end if
      if (opts%command == cmd_dyn .and. opts%q_mode == dyn_q_grid) then
         select type (method)
         type is (dynamics_structure_factor_t)
            if (method%keep_modes) then
               write (shell_unit, '(a,i0,a,f0.4,a)') '# reciprocal lattice grid: ', &
                  method%nmodes, ' lattice vectors with |q| <= ', method%grid_qmax, &
                  ' 1/A (one row per vector)'
            else
               write (shell_unit, '(a,i0,a,f0.4,a)') '# reciprocal lattice grid: ', &
                  method%nmodes, ' rows with |q| <= ', method%grid_qmax, &
                  ' 1/A (one row per lattice shell, Gamma included)'
            end if
         end select
      end if
      if (opts%command == cmd_dyn .and. opts%q_mode == dyn_q_single) then
         select type (method)
         type is (dynamics_structure_factor_t)
            write (shell_unit, '(a,3(i0,1x),a,f0.4,a)') '# single q from n = (', &
               method%single_index, ') at |q| = ', method%qlen(1), ' 1/A (one row)'
         end select
      end if
   end if
   grid_unit = no_unit
   if (opts%want_grid) then
      if (opts%grid_format /= format_hdf5) then
         call open_output(opts%grid_output, grid_unit, ierr, message)
         if (ierr /= 0) then
            write (error_unit, '(a)') 'sqcalc: '//trim(message)
            stop 12
         end if
      end if
   end if
   call method%write_results(shell_unit, grid_unit, ierr, message)
   if (ierr /= 0) then
      write (error_unit, '(a)') 'sqcalc: '//trim(message)
      stop 13
   end if
   ! --- powder XRD pattern -----------------------------------------------
   xrd_unit = no_unit
   if (allocated(opts%xrd_output) .and. opts%xrd_format /= format_hdf5) then
      call open_output(opts%xrd_output, xrd_unit, ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a)') 'sqcalc: '//trim(message)
         stop 21
      end if
   end if
   call method%write_xrd(xrd_unit, ierr, message)
   if (ierr /= 0) then
      write (error_unit, '(a)') 'sqcalc: '//trim(message)
      stop 21
   end if
#ifdef SQC_HAVE_HDF5
   if (allocated(opts%xrd_output) .and. opts%xrd_format == format_hdf5) then
      call hdf5_write_xrd(opts%xrd_output, method, ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a)') 'sqcalc: '//trim(message)
         stop 21
      end if
   end if
#endif
   if (opts%command == cmd_dyn) then
      select type (method)
      type is (dynamics_structure_factor_t)
         if (allocated(method%sqw_path)) then
            call method%write_sqw(method%sqw_path, method%sqw_format, opts%scheme, opts%input, &
                                  ierr, message)
            if (ierr /= 0) then
               write (error_unit, '(a)') 'sqcalc: '//trim(message)
               stop 22
            end if
         end if
         if (allocated(method%fqt_path)) then
            call method%write_fqt(method%fqt_path, method%fqt_format, opts%scheme, opts%input, &
                                  ierr, message)
            if (ierr /= 0) then
               write (error_unit, '(a)') 'sqcalc: '//trim(message)
               stop 23
            end if
         end if
         if (allocated(method%fqt_self_path)) then
            call method%write_fqt_self(method%fqt_self_path, method%fqt_self_format, &
                                       opts%scheme, opts%input, ierr, message)
            if (ierr /= 0) then
               write (error_unit, '(a)') 'sqcalc: '//trim(message)
               stop 28
            end if
         end if
         if (allocated(method%s4_path)) then
            call method%write_s4(method%s4_path, method%s4_format, opts%scheme, opts%input, &
                                 ierr, message)
            if (ierr /= 0) then
               write (error_unit, '(a)') 'sqcalc: '//trim(message)
               stop 26
            end if
         end if
         if (allocated(method%chi4_path)) then
            call method%write_chi4(method%chi4_path, method%chi4_format, opts%input, &
                                   ierr, message)
            if (ierr /= 0) then
               write (error_unit, '(a)') 'sqcalc: '//trim(message)
               stop 27
            end if
         end if
         if (allocated(method%msd_path)) then
            call method%write_msd(method%msd_path, method%msd_format, opts%scheme, &
                                  opts%input, ierr, message)
            if (ierr /= 0) then
               write (error_unit, '(a)') 'sqcalc: '//trim(message)
               stop 31
            end if
         end if
         if (allocated(method%ngp_path)) then
            call method%write_ngp(method%ngp_path, method%ngp_format, opts%scheme, &
                                  opts%input, ierr, message)
            if (ierr /= 0) then
               write (error_unit, '(a)') 'sqcalc: '//trim(message)
               stop 32
            end if
         end if
      end select
   end if
   if (opts%want_grid .and. opts%grid_format == format_hdf5) then
#ifdef SQC_HAVE_HDF5
      call hdf5_write_results(opts%grid_output, method, opts%scheme, natoms, ref_a, eps_used, &
                              device_name(opts), precision_name(opts), ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a)') 'sqcalc: '//trim(message)
         stop 18
      end if
#endif
   end if
   if (shell_unit /= no_unit .and. shell_unit /= output_unit) close (shell_unit)
   if (grid_unit /= no_unit .and. grid_unit /= output_unit) close (grid_unit)
   if (xrd_unit /= no_unit .and. xrd_unit /= output_unit) close (xrd_unit)

   if (.not. opts%quiet) then
      if (opts%command == cmd_dyn) then
         select type (method)
         type is (dynamics_structure_factor_t)
            if (method%probe_valid) then
               write (error_unit, '(a,i0,a,f0.2,a,f0.4)') '  isotropy   : ', &
                  method%probe_modes, ' modes of the first shell spread by ', &
                  100.0_rk*method%probe_spread, ' % at tau = ', method%probe_tau
               write (error_unit, '(a)') '  note       : compare that spread with the run-to-run '// &
                  'scatter of a single mode; if it stays far above it the system is not '// &
                  'isotropic or not equilibrated, and the orbit reduction does not apply'
            end if
         end select
      end if
      write (error_unit, '(a,i0,a,a,a,f0.2,a,i0,a)') 'averaged ', method%nframes, ' frames (', &
         trim(method_description(opts)), ', ', elapsed, ' s wall, threads = ', &
         opts%threads, ')'
   end if

contains

   !> Element symbols of every type id must be known for weighted schemes.
   subroutine validate_mapping(scheme, type_id)
      type(weight_scheme_t), intent(in) :: scheme
      integer(ik), intent(in) :: type_id(:)
      logical, allocatable :: seen(:)
      integer :: i, t, idx

      if (scheme%kind == weight_unit) return
      allocate (seen(maxval(type_id)))
      seen = .false.
      do i = 1, size(type_id)
         seen(int(type_id(i))) = .true.
      end do
      do t = 1, size(seen)
         if (.not. seen(t)) cycle
         if (t > scheme%mapped_types()) then
            write (error_unit, '(a,i0,a)') 'sqcalc: atom type ', t, &
               ' appears in the dump but is missing from the element mapping'
            stop 14
         end if
         idx = scheme%element(t)
         if (idx == 0) then
            write (error_unit, '(a,i0,a)') 'sqcalc: atom type ', t, &
               ' has no species in the mapping'
            stop 14
         end if
         if (scheme%kind == weight_neutron .and. .not. element_table%has_neutron(idx)) then
            write (error_unit, '(a,a,a)') 'sqcalc: no neutron scattering length for ', &
               trim(scheme%species(t)), '; use --weight xray or unit'
            stop 15
         end if
         if (scheme%kind == weight_xray .and. scheme%species_index(t) == 0) then
            write (error_unit, '(a,a,a)') 'sqcalc: no X-ray form factor for ', &
               trim(scheme%species(t)), '; use --weight neutron or unit'
            stop 15
         end if
      end do
   end subroutine validate_mapping

   subroutine report_setup(opt, frame)
      type(options_t), intent(in) :: opt
      type(frame_t), intent(in) :: frame
      type(cell_t) :: cell
      integer :: i
      character(len=8) :: species

      cell = frame%cell
      write (error_unit, '(a)') 'sqcalc: '//program_version
      write (error_unit, '(a,a)') '  input      : ', trim(opt%input)
      write (error_unit, '(a,i0)') '  atoms      : ', frame%natoms
      write (error_unit, '(a,3(f0.4,1x),a,f0.2)') '  cell       : ', &
         sqrt(sum(cell%a(:, 1)**2)), sqrt(sum(cell%a(:, 2)**2)), sqrt(sum(cell%a(:, 3)**2)), &
         ' volume ', cell%volume
      write (error_unit, '(a,a)') '  weights    : ', trim(opt%scheme%label())
      if (opt%scheme%has_mapping()) then
         write (error_unit, '(a)', advance='no') '  mapping    : '
         do i = 1, opt%scheme%mapped_types()
            species = opt%scheme%species(i)
            if (len_trim(species) == 0) cycle
            write (error_unit, '(a,i0,a,a,a)', advance='no') ' ', i, ':', trim(species), ','
         end do
         write (error_unit, '(a)') ''
      end if
      write (error_unit, '(a,i0)') '  threads    : ', opt%threads
      if (opt%device == device_gpu) write (error_unit, '(a,i0)') '  gpu device : ', opt%gpu_id
      if (opt%method == method_debye) then
         if (opt%rmax > 0.0_rk) then
            write (error_unit, '(a,f0.3,a)') '  pair cutoff: ', opt%rmax, ' A'
         else
            write (error_unit, '(a)') '  pair cutoff: half the smallest periodic side'
         end if
         write (error_unit, '(a,f0.4,a,f0.3,a)') '  r bins     : ', opt%dr, ' A (skin ', &
            opt%skin, ' A)'
         if (opt%dr > 0.5_rk*acos(-1.0_rk)/opt%qmax) then
            write (error_unit, '(a,f0.2,a)') '  note       : the radial bin width limits the '// &
               'reliable range to q ~ ', 0.5_rk*acos(-1.0_rk)/opt%dr, ' 1/A'
         end if
      end if
      if (opt%device == device_gpu) then
         if (opt%precision == precision_single) then
            write (error_unit, '(a)') '  precision  : single (float32, cufinufftf)'
         else
            write (error_unit, '(a)') '  precision  : double (float64, cufinufft)'
         end if
      end if
   end subroutine report_setup

   subroutine report_grid(m)
      class(structure_factor_t), intent(in) :: m
      select type (m)
      type is (dynamics_structure_factor_t)
         select case (m%q_mode)
         case (dyn_q_shell)
            write (error_unit, '(a,f0.4,a,i0,a,i0,a)') '  q shell    : |q| = ', m%shell_q, &
               ' 1/A, Lebedev order ', m%shell_order, ' (', m%nmodes, &
               ' directions, +- pairs merged)'
         case (dyn_q_powder)
            write (error_unit, '(a,f0.4,a,f0.4,a,i0,a,i0,a,f0.4,a)') '  q powder   : |q| = ', &
               m%powder_q, ' 1/A +- ', m%powder_dq_used, ', ', m%powder_vectors, &
               ' lattice vectors in ', m%powder_nshell, ' shells, mean |q| = ', &
               m%powder_qmean, ' 1/A'
            if (m%powder_nearest) then
               write (error_unit, '(a)') '  note       : the window held no shell of its '// &
                  'own; the nearest one was used'
            end if
         case (dyn_q_grid)
            if (m%grid_nshell_kept < m%grid_nshell) then
               write (error_unit, '(a,f0.4,a,i0,a,i0,a,i0,a,i0,a)') '  q grid     : |q| <= ', &
                  m%grid_qmax, ' 1/A, ', m%nmodes, ' lattice modes in ', m%grid_nshell_kept, &
                  ' of ', m%grid_nshell, ' shells (point group order ', m%grid_nops, ')'
            else
               write (error_unit, '(a,f0.4,a,i0,a,i0,a,i0,a)') '  q grid     : |q| <= ', &
                  m%grid_qmax, ' 1/A, ', m%nmodes, ' lattice modes in ', m%grid_nshell, &
                  ' shells (point group order ', m%grid_nops, ')'
            end if
            if (m%grid_thinned .and. .not. m%grid_budget_met) then
               write (error_unit, '(a,i0,a,i0,a)') '  note       : --modes ', m%grid_budget, &
                  ' could not be met; the two end shells alone need ', m%nmodes, ' modes'
            end if
            if (m%grid_group_fallback) then
               write (error_unit, '(a)') '  note       : the lattice point group needs '// &
                  'larger Miller indices than the search box; the modes are exact, '// &
                  'but orbit thinning is off'
            end if
         case (dyn_q_single)
            write (error_unit, '(a,3(i0,1x),a,f0.4,a)') '  q single   : n = (', &
               m%single_index, ') -> |q| = ', m%qlen(1), ' 1/A'
         case (dyn_q_none)
            write (error_unit, '(a)') '  q sampling : none (only the overlap, MSD and NGP outputs)'
         case default
            write (error_unit, '(a,i0,a,f0.4,a,f0.4,a)') '  q line     : ', &
               m%nintervals + 1, ' points from ', m%s0, ' to ', m%s1, ' 1/A'
         end select
         write (error_unit, '(a,i0,a,i0)') '  window     : ', m%maxframes, &
            ' frames, origin lag ', m%lag_stride
         if (m%s4_enabled .or. m%chi4_enabled .or. m%fqt_self_enabled .or. &
             m%msd_enabled .or. m%ngp_enabled) then
            if (m%s4_enabled .or. m%chi4_enabled) then
               write (error_unit, '(a,f0.4,a)') '  overlap    : cutoff ', m%s4_cutoff, &
                  ' (S4/chi4 use unit weights)'
            end if
            write (error_unit, '(a,i0,a,i0,a,i0,a)') '  stride     : ', m%stride, &
               ' dump frames, effective maxframes ', m%effective_maxframes, ' (requested ', &
               m%maxframes, ')'
            if (m%effective_maxframes /= m%maxframes) then
               write (error_unit, '(a)') '  note       : the S4/chi4/F_s/MSD/NGP window was reduced to a '// &
                  'multiple of --stride; the coherent F(q,t)/S(q,w) window is unchanged'
            end if
            write (error_unit, '(a,f0.4,a,f0.4,a)') '  buffer     : ', &
               24.0_rk*real(m%natoms, rk)*real(m%nsteps + 1, rk)/1.0e9_rk, &
               ' GB (limit ', m%buffer_limit_gb, ' GB)'
         end if
         if (m%fqt_self_enabled) then
            write (error_unit, '(a)') '  F_s        : self intermediate scattering function '// &
               '(unit weights)'
         end if
         if (m%msd_enabled) then
            write (error_unit, '(a)') '  MSD        : mean squared displacement (unit weights)'
         end if
         if (m%ngp_enabled) then
            write (error_unit, '(a)') '  NGP        : non-Gaussian parameter alpha2 (unit weights)'
         end if
      class default
         write (error_unit, '(a,3(i0,1x))') '  grid modes : ', m%modes
         if (m%expected_bytes > 0.0_rk) then
            write (error_unit, '(a,i0,a,f0.1,a)') '  grid points: ', m%gridpoints, &
               ' (~', m%expected_bytes/1.0e6_rk, ' MB peak, estimated)'
         else
            write (error_unit, '(a,i0)') '  grid points: ', m%gridpoints
         end if
         write (error_unit, '(a,i0,a,f0.4,a,f0.4,a)') '  q range    : ', m%nmodes, &
            ' modes in [', m%qmin, ', ', m%qmax, '] 1/A'
         if (m%xrd_enabled) then
            write (error_unit, '(a,f0.6,a,f0.4,a,f0.4,a,f0.5,a)') '  xrd        : lambda ', &
               m%xrd_lambda, ' A, 2theta ', m%xrd_2theta_min, '-', m%xrd_2theta_max, &
               ' deg, step ', m%xrd_step, ' deg'
            write (error_unit, '(a,i0,a,i0,a)') '  xrd bins   : ', m%xrd_bins, &
               ' bins, ', sf_xrd_empty_bins(m), ' without a reciprocal lattice point'
         end if
      end select
   end subroutine report_grid

   !> Open the requested output target; "-" means standard output.
   subroutine open_output(path, unit, ierr, msg)
      character(len=*), intent(in) :: path
      integer, intent(out) :: unit
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: msg

      ierr = 0
      msg = ''
      if (trim(path) == '-') then
         unit = output_unit
         return
      end if
      open (newunit=unit, file=trim(path), status='replace', action='write', iostat=ierr)
      if (ierr /= 0) then
         unit = -1
         msg = 'cannot write output file "'//trim(path)//'"'
      end if
   end subroutine open_output

   function method_description(opt) result(text)
      type(options_t), intent(in) :: opt
      character(len=:), allocatable :: text
      if (opt%command == cmd_dyn) then
         select case (opt%q_mode)
         case (dyn_q_line)
            text = 'dynamic structure factor (direct summation on a q line)'
         case (dyn_q_shell)
            text = 'dynamic structure factor (Lebedev average on a q shell)'
         case (dyn_q_grid)
            text = 'dynamic structure factor on the reciprocal grid'
         case (dyn_q_single)
            text = 'dynamic structure factor at one lattice vector'
         case (dyn_q_powder)
            text = 'dynamic structure factor averaged over a window of lattice shells'
         case default
            text = 'overlap dynamics (average overlap and chi4)'
         end select
      else
         select case (opt%method)
         case (method_direct)
            text = 'direct summation'
         case (method_debye)
            text = 'Debye pair histograms'
         case default
            if (opt%device == device_gpu) then
               if (opt%precision == precision_single) then
                  text = 'NUFFT on GPU, float32 (cufinufftf/cuFFT)'
               else
                  text = 'NUFFT on GPU, float64 (cufinufft/cuFFT)'
               end if
            else
               text = 'NUFFT on CPU (FINUFFT)'
            end if
         end select
      end if
   end function method_description

   !> Name of the device used for the transform (written to HDF5 metadata).
   function device_name(opt) result(text)
      type(options_t), intent(in) :: opt
      character(len=:), allocatable :: text
      if (opt%device == device_gpu) then
         text = 'gpu'
      else
         text = 'cpu'
      end if
   end function device_name

   !> Name of the transform precision (written to HDF5 metadata).
   function precision_name(opt) result(text)
      type(options_t), intent(in) :: opt
      character(len=:), allocatable :: text
      if (opt%device == device_gpu .and. opt%precision == precision_single) then
         text = 'single'
      else
         text = 'double'
      end if
   end function precision_name

   !> Current line number of the concrete reader (for error messages).
   integer function reader_line(r) result(line)
      class(frame_reader_t), intent(in) :: r
      select type (r)
      type is (lammps_dump_reader_t)
         line = r%lineno
      class default
         line = 0
      end select
   end function reader_line

   !> Timestep increment between consecutive frames (0 when not determined).
   integer(ik) function reader_step_stride(r) result(step)
      class(frame_reader_t), intent(in) :: r
      step = 0
      select type (r)
      type is (lammps_dump_reader_t)
         step = r%step_stride
      end select
   end function reader_step_stride

   !> True while every frame was written after the same number of steps.
   logical function reader_uniform_steps(r) result(uniform)
      class(frame_reader_t), intent(in) :: r
      uniform = .false.
      select type (r)
      type is (lammps_dump_reader_t)
         uniform = r%uniform_timestep
      end select
   end function reader_uniform_steps

   !> Number of atom displacements that were too large to unwrap reliably.
   integer function reader_ambiguous(r) result(n)
      class(frame_reader_t), intent(in) :: r
      n = 0
      select type (r)
      type is (lammps_dump_reader_t)
         n = r%n_ambiguous
      end select
   end function reader_ambiguous

end program sqcalc
