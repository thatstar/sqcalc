! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> sqcalc: total structure factors from LAMMPS dump trajectories.
program sqcalc
   use sqc_kinds
   use sqc_options
   use sqc_dump
   use sqc_weights
   use sqc_cell, only: cell_t
   use sqc_elements, only: element_table
   use sqc_structure_factor
   use sqc_debye, only: debye_structure_factor_t
   use sqc_finufft, only: finufft_opts_is_consistent
#ifdef SQC_ENABLE_CUDA
   use sqc_gpu, only: cufinufft_structure_factor_t
#endif
#ifdef SQC_HAVE_HDF5
   use sqc_hdf5, only: hdf5_write_results
#endif
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none

   type(options_t) :: opts
   class(frame_reader_t), allocatable :: reader
   class(structure_factor_t), allocatable :: method
   type(frame_t) :: frame
   character(len=512) :: message
   integer :: ierr, shell_unit, grid_unit, tick, tick_rate, progress_step
   integer(ik) :: natoms
   real(rk) :: ref_a(3, 3), elapsed, wall0, wall1, eps_used
   logical :: first_frame

   call parse_options(opts, ierr, message)
   if (opts%show_help) then
      call print_usage(output_unit)
      stop
   end if
   if (opts%show_version) then
      write (output_unit, '(a)') program_version
      stop
   end if
   if (ierr /= 0) then
      write (error_unit, '(a)') 'sqcalc: '//trim(message)
      call print_usage(error_unit)
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
         if (allocated(opts%rdf_output)) method%rdf_path = opts%rdf_output
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
   if (opts%want_grid .and. opts%grid_format == grid_format_hdf5) then
      write (error_unit, '(a)') 'sqcalc: this build has no HDF5 support; '// &
         'configure with -DSQC_ENABLE_HDF5=ON or use --grid-format text'
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
   call open_output(opts%output, shell_unit, ierr, message)
   if (ierr /= 0) then
      write (error_unit, '(a)') 'sqcalc: '//trim(message)
      stop 11
   end if
   grid_unit = no_unit
   if (opts%want_grid) then
      if (opts%grid_format /= grid_format_hdf5) then
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
   if (opts%want_grid .and. opts%grid_format == grid_format_hdf5) then
#ifdef SQC_HAVE_HDF5
      call hdf5_write_results(opts%grid_output, method, opts%scheme, natoms, ref_a, eps_used, &
                              device_name(opts), precision_name(opts), ierr, message)
      if (ierr /= 0) then
         write (error_unit, '(a)') 'sqcalc: '//trim(message)
         stop 18
      end if
#endif
   end if
   if (shell_unit /= output_unit) close (shell_unit)
   if (grid_unit /= no_unit .and. grid_unit /= output_unit) close (grid_unit)

   if (.not. opts%quiet) then
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
               ' has no element symbol in the mapping'
            stop 14
         end if
         if (scheme%kind == weight_neutron .and. .not. element_table%has_neutron(idx)) then
            write (error_unit, '(a,a,a)') 'sqcalc: no neutron scattering length for ', &
               trim(element_table%symbol_of(idx)), '; use --weight xray or unit'
            stop 15
         end if
         if (scheme%kind == weight_xray .and. .not. element_table%has_xray(idx)) then
            write (error_unit, '(a,a,a)') 'sqcalc: no X-ray form factor for ', &
               trim(element_table%symbol_of(idx)), '; use --weight neutron or unit'
            stop 15
         end if
      end do
   end subroutine validate_mapping

   subroutine report_setup(opt, frame)
      type(options_t), intent(in) :: opt
      type(frame_t), intent(in) :: frame
      type(cell_t) :: cell
      integer :: i
      character(len=2) :: symbol

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
            symbol = opt%scheme%symbols(i)
            if (len_trim(symbol) == 0) cycle
            write (error_unit, '(a,i0,a,a,a)', advance='no') ' ', i, ':', trim(symbol), ','
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
      write (error_unit, '(a,3(i0,1x))') '  grid modes : ', m%modes
      write (error_unit, '(a,i0)') '  grid points: ', m%gridpoints
      write (error_unit, '(a,i0,a,f0.4,a,f0.4,a)') '  q range    : ', m%nmodes, &
         ' modes in [', m%qmin, ', ', m%qmax, '] 1/A'
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

end program sqcalc
