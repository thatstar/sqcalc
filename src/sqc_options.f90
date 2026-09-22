! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Command line handling for sqcalc.
module sqc_options
   use sqc_kinds
   use sqc_weights, only: weight_scheme_t, weight_unit, weight_neutron, weight_xray, &
                          scheme_from_name
   use sqc_structure_factor, only: method_nufft, method_direct, method_debye, norm_mean, &
                            norm_self, norm_natom, method_dynamic
   use sqc_debye, only: debye_default_dr, debye_default_skin
   use sqc_dynamics, only: dyn_q_none, dyn_q_line, dyn_q_shell, dyn_q_grid, dyn_q_single
   use sqc_modes, only: thin_none, thin_shells, thin_orbits
   use sqc_lebedev, only: lebedev_points, lebedev_order_from_name, &
                          lebedev_low, lebedev_medium, lebedev_high
   implicit none
   private

   public :: options_t, parse_options, print_usage, program_version, device_cpu, device_gpu, &
             precision_double, precision_single, grid_format_text, grid_format_hdf5

   !> Where the NUFFT transforms run.
   integer, parameter :: device_cpu = 0
   integer, parameter :: device_gpu = 1

   !> Precision of the GPU transform.
   integer, parameter :: precision_double = 0
   integer, parameter :: precision_single = 1

   !> Format of the reciprocal grid table.
   integer, parameter :: grid_format_text = 0
   integer, parameter :: grid_format_hdf5 = 1

   character(len=*), parameter :: program_version = 'sqcalc 0.1.0'

   !> pi and the degree to radian conversion (the XRD options are in degrees).
   real(rk), parameter :: pi = 3.14159265358979323846_rk
   real(rk), parameter :: deg2rad = pi/180.0_rk

   !> Everything configurable from the command line.
   type :: options_t
      character(len=:), allocatable :: input
      character(len=:), allocatable :: output
      character(len=:), allocatable :: grid_output
      integer :: threads = 0
      integer :: method = method_nufft
      logical :: method_given = .false.
      integer :: norm = norm_mean
      integer :: device = device_cpu
      integer :: gpu_id = 0
      integer :: precision = precision_double
      !> Debye method settings (--rmax, --dr, --skin, --rdf).
      real(rk) :: rmax = 0.0_rk
      real(rk) :: dr = debye_default_dr
      real(rk) :: skin = debye_default_skin
      character(len=:), allocatable :: rdf_output
      character(len=:), allocatable :: pair_entropy_output
      character(len=:), allocatable :: s2_accum_output
      logical :: rmax_given = .false.
      logical :: dr_given = .false.
      logical :: skin_given = .false.
      !> Disable the Debye cut-off density correction.
      logical :: no_cutoff_correction = .false.
      !> Write partial structure factor columns (default: when -m is given).
      logical :: partials = .true.
      logical :: partials_given = .false.
      !> Use the Faber-Ziman form for the partials.
      logical :: faber_ziman = .false.
      !> True when the user gave --eps explicitly.
      logical :: eps_given = .false.
      integer :: nq = 500
      real(rk) :: qmin = 0.0_rk
      real(rk) :: qmax = 20.0_rk
      !> True when the user set the q sampling himself (the XRD output owns it).
      logical :: qmin_given = .false.
      logical :: qmax_given = .false.
      logical :: nq_given = .false.
      real(rk) :: eps = 1.0e-9_rk
      logical :: want_grid = .false.
      integer :: grid_format = grid_format_text
      logical :: grid_format_given = .false.
      !> Powder XRD output (--xrd FILE) and the settings it needs: the
      !! wavelength [dump length unit], the two-theta range [deg], the bin
      !! width [deg, 0 = derive from the box] and the LP/weight switches.
      character(len=:), allocatable :: xrd_output
      real(rk) :: xrd_lambda = 0.0_rk
      logical :: xrd_lambda_given = .false.
      real(rk) :: xrd_2theta_min = 1.0_rk
      real(rk) :: xrd_2theta_max = 179.0_rk
      logical :: xrd_range_given = .false.
      real(rk) :: xrd_step = 0.0_rk
      logical :: xrd_step_given = .false.
      integer :: xrd_format = grid_format_text
      logical :: lp = .true.
      logical :: lp_given = .false.
      !> True when -w was given; --xrd implies the x-ray weights otherwise.
      logical :: weight_given = .false.
      logical :: quiet = .false.
      logical :: show_help = .false.
      logical :: show_version = .false.
      !> Dynamic run (--dyn) and the q sampling that feeds it (--dyn-q).
      logical :: dynamic = .false.
      integer :: dyn_q_mode = dyn_q_none
      logical :: dyn_q_given = .false.
      !> line: NINT intervals, scale from S0 to S1 along (DX, DY, DZ).
      integer :: dyn_intervals = 0
      real(rk) :: dyn_s0 = 0.0_rk
      real(rk) :: dyn_s1 = 0.0_rk
      real(rk) :: dyn_dir(3) = 0.0_rk
      !> shell: |q| radius and the order of the Lebedev rule it selects.
      real(rk) :: dyn_shell_q = 0.0_rk
      integer :: dyn_shell_order = 0
      !> grid: the upper bound of |q|, the mode budget and the thinning policy.
      real(rk) :: dyn_grid_qmax = 0.0_rk
      !> single: the Miller indices of the one lattice vector to sample.
      integer :: dyn_single(3) = 0
      integer :: dyn_modes = 0
      integer :: dyn_thin = thin_none
      logical :: dyn_modes_given = .false.
      logical :: dyn_thin_given = .false.
      logical :: dyn_keep_modes = .false.
      !> MD time step [time units] and the correlation window settings.
      real(rk) :: dt = 0.0_rk
      integer :: maxframes = 0
      integer :: lag_stride = 1
      logical :: lag_given = .false.
      !> Optional dynamic outputs.
      character(len=:), allocatable :: sqw_output
      character(len=:), allocatable :: fqt_output
      character(len=:), allocatable :: fqt_self_output
      !> Unified format override; otherwise each file infers from its suffix.
      integer :: dyn_format = grid_format_text
      logical :: dyn_format_given = .false.
      integer :: sqw_format = grid_format_text
      integer :: fqt_format = grid_format_text
      integer :: fqt_self_format = grid_format_text
      !> Four-point structure factor and average overlap / chi4 outputs.
      character(len=:), allocatable :: s4_output, chi4_output
      !> Mean squared displacement output (--msd).
      character(len=:), allocatable :: msd_output
      real(rk) :: s4_cutoff = 0.0_rk
      logical :: s4_cutoff_given = .false.
      real(rk) :: buffer_limit_gb = 2.0_rk
      logical :: buffer_limit_given = .false.
      integer :: stride = 1
      logical :: stride_given = .false.
      integer :: s4_format = grid_format_text
      integer :: chi4_format = grid_format_text
      integer :: msd_format = grid_format_text
      type(weight_scheme_t) :: scheme
   end type options_t

contains

   !> Parse the argument vector into an options_t.
   subroutine parse_options(self, ierr, message)
      type(options_t), intent(inout) :: self
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=256) :: arg, name, value, value2, positional(4)
      integer :: i, nargs, npos, eq, kind
      logical :: has_inline, needs_value

      ierr = 0
      message = ''
      npos = 0
      nargs = command_argument_count()
      i = 1
      do while (i <= nargs)
         call get_command_argument(i, arg)
         arg = trim(adjustl(arg))
         if (len_trim(arg) == 0) then
            i = i + 1
            cycle
         end if

         ! A bare "-" is the stdout marker and counts as positional.
         if (arg(1:1) == '-' .and. len_trim(arg) > 1) then
            eq = index(arg, '=')
            if (eq > 0) then
               name = arg(:eq - 1)
               value = arg(eq + 1:)
               has_inline = .true.
            else
               name = arg
               value = ''
               has_inline = .false.
            end if

            needs_value = .true.
            select case (trim(name))
            case ('-h', '--help')
               self%show_help = .true.
               needs_value = .false.
            case ('-v', '--version')
               self%show_version = .true.
               needs_value = .false.
            case ('-q', '--quiet')
               self%quiet = .true.
               needs_value = .false.
            case ('--dyn')
               if (has_inline) then
                  ierr = 1
                  message = '--dyn takes no value; pick the q sampling with --dyn-q, e.g. '// &
                     '--dyn-q line:100,0.5,20,1,1,0'
                  return
               end if
               self%dynamic = .true.
               needs_value = .false.
            case ('--no-cutoff-correction')
               self%no_cutoff_correction = .true.
               needs_value = .false.
            case ('--lp')
               self%lp = .true.
               self%lp_given = .true.
               needs_value = .false.
            case ('--no-lp')
               self%lp = .false.
               self%lp_given = .true.
               needs_value = .false.
            case ('-fz', '--faber-ziman')
               self%faber_ziman = .true.
               self%partials = .true.
               needs_value = .false.
            case ('--partials')
               self%partials = .true.
               self%partials_given = .true.
               needs_value = .false.
            case ('--no-partials')
               self%partials = .false.
               self%partials_given = .true.
               needs_value = .false.
            case ('--dyn-keep-modes')
               self%dyn_keep_modes = .true.
               needs_value = .false.
            case ('-i', '--input', '--mapping', '-m', '-w', '--weight', '-t', '--threads', &
                  '--qmin', '--qmax', '--nq', '--eps', '--method', '--norm', '--grid', &
                  '--device', '--gpu-id', '--precision', '--grid-format', &
                  '--xrd', '--xrd-lambda', '--xrd-range', '--xrd-step', &
                  '--rmax', '--dr', '--skin', '--rdf', &
                  '--pair-entropy', '--s2-accum', &
                  '--dyn-q', '--dt', '--maxframes', '--lag', '--sqw', '--fqt', '--dyn-format', &
                  '--fqt-self', '--s4', '--chi4', '--msd', '--s4-cutoff', '--buffer-limit', &
                  '--stride', &
                  '--dyn-modes', '--dyn-thin')
               if (.not. has_inline) then
                  if (i + 1 > nargs) then
                     ierr = 1
                     message = 'missing value for option '//trim(name)
                     return
                  end if
                  i = i + 1
                  call get_command_argument(i, value)
               end if
               ! --xrd-range takes the two limits as separate arguments.
               if (trim(name) == '--xrd-range') then
                  if (has_inline) then
                     ierr = 1
                     message = '--xrd-range takes two values: MIN2TH MAX2TH'
                     return
                  end if
                  if (i + 1 > nargs) then
                     ierr = 1
                     message = 'missing value for option '//trim(name)
                     return
                  end if
                  i = i + 1
                  call get_command_argument(i, value2)
               end if
            case default
               ierr = 1
               message = 'unknown option "'//trim(name)//'"'
               return
            end select

            if (needs_value) then
               select case (trim(name))
               case ('-i', '--input')
                  self%input = trim(value)
               case ('-m', '--mapping')
                  call self%scheme%set_mapping(trim(value), ierr, message)
                  if (ierr /= 0) return
               case ('-w', '--weight')
                  kind = scheme_from_name(value)
                  if (kind < 0) then
                     ierr = 1
                     message = 'unknown weight scheme "'//trim(value)//'" (use unit, neutron or xray)'
                     return
                  end if
                  self%scheme%kind = kind
                  self%weight_given = .true.
               case ('-t', '--threads')
                  read (value, *, iostat=ierr) self%threads
                  if (ierr /= 0 .or. self%threads < 1) then
                     ierr = 1
                     message = 'thread count must be a positive integer'
                     return
                  end if
               case ('--qmin')
                  read (value, *, iostat=ierr) self%qmin
                  if (ierr /= 0 .or. self%qmin < 0.0_rk) then
                     ierr = 1
                     message = '--qmin must be a non-negative number'
                     return
                  end if
                  self%qmin_given = .true.
               case ('--qmax')
                  read (value, *, iostat=ierr) self%qmax
                  if (ierr /= 0 .or. self%qmax <= 0.0_rk) then
                     ierr = 1
                     message = '--qmax must be a positive number'
                     return
                  end if
                  self%qmax_given = .true.
               case ('--nq')
                  read (value, *, iostat=ierr) self%nq
                  if (ierr /= 0 .or. self%nq < 1) then
                     ierr = 1
                     message = '--nq must be a positive integer'
                     return
                  end if
                  self%nq_given = .true.
               case ('--eps')
                  read (value, *, iostat=ierr) self%eps
                  if (ierr /= 0 .or. self%eps <= 0.0_rk) then
                     ierr = 1
                     message = '--eps must be a positive number'
                     return
                  end if
                  self%eps_given = .true.
               case ('--precision')
                  select case (trim(value))
                  case ('double', 'float64', 'f64')
                     self%precision = precision_double
                  case ('single', 'float', 'float32', 'f32')
                     self%precision = precision_single
                  case default
                     ierr = 1
                     message = 'unknown precision "'//trim(value)//'" (use single or double)'
                     return
                  end select
               case ('--method')
                  select case (trim(value))
                  case ('nufft', 'finufft')
                     self%method = method_nufft
                  case ('direct')
                     self%method = method_direct
                  case ('debye')
                     self%method = method_debye
                  case default
                     ierr = 1
                     message = 'unknown method "'//trim(value)//'" (use nufft, direct or debye)'
                     return
                  end select
                  self%method_given = .true.
               case ('--rmax')
                  read (value, *, iostat=ierr) self%rmax
                  if (ierr /= 0 .or. self%rmax <= 0.0_rk) then
                     ierr = 1
                     message = '--rmax must be a positive number'
                     return
                  end if
                  self%rmax_given = .true.
               case ('--dr')
                  read (value, *, iostat=ierr) self%dr
                  if (ierr /= 0 .or. self%dr <= 0.0_rk) then
                     ierr = 1
                     message = '--dr must be a positive number'
                     return
                  end if
                  self%dr_given = .true.
               case ('--skin')
                  read (value, *, iostat=ierr) self%skin
                  if (ierr /= 0 .or. self%skin < 0.0_rk) then
                     ierr = 1
                     message = '--skin must be zero or positive'
                     return
                  end if
                  self%skin_given = .true.
               case ('--rdf')
                  self%rdf_output = trim(value)
               case ('--pair-entropy')
                  self%pair_entropy_output = trim(value)
               case ('--s2-accum')
                  self%s2_accum_output = trim(value)
               case ('--dyn-q')
                  call parse_dyn_q(self, trim(value), ierr, message)
                  if (ierr /= 0) return
               case ('--dyn-modes')
                  read (value, *, iostat=ierr) self%dyn_modes
                  if (ierr /= 0 .or. self%dyn_modes < 0) then
                     ierr = 1
                     message = '--dyn-modes must be a non-negative integer'
                     return
                  end if
                  self%dyn_modes_given = .true.
               case ('--dyn-thin')
                  select case (trim(value))
                  case ('shells')
                     self%dyn_thin = thin_shells
                  case ('orbits')
                     self%dyn_thin = thin_orbits
                  case default
                     ierr = 1
                     message = '--dyn-thin wants "shells" or "orbits"'
                     return
                  end select
                  self%dyn_thin_given = .true.
               case ('--dt')
                  read (value, *, iostat=ierr) self%dt
                  if (ierr /= 0 .or. self%dt <= 0.0_rk) then
                     ierr = 1
                     message = '--dt must be a positive number (the MD time step)'
                     return
                  end if
               case ('--maxframes')
                  read (value, *, iostat=ierr) self%maxframes
                  if (ierr /= 0 .or. self%maxframes < 1) then
                     ierr = 1
                     message = '--maxframes must be a positive integer'
                     return
                  end if
               case ('--lag')
                  read (value, *, iostat=ierr) self%lag_stride
                  if (ierr /= 0 .or. self%lag_stride < 1) then
                     ierr = 1
                     message = '--lag must be a positive integer'
                     return
                  end if
                  self%lag_given = .true.
               case ('--sqw')
                  self%sqw_output = trim(value)
               case ('--fqt')
                  self%fqt_output = trim(value)
               case ('--fqt-self')
                  self%fqt_self_output = trim(value)
               case ('--dyn-format')
                  select case (trim(value))
                  case ('text', 'txt', 'ascii')
                     self%dyn_format = grid_format_text
                  case ('hdf5', 'h5', 'hdf')
                     self%dyn_format = grid_format_hdf5
                  case default
                     ierr = 1
                     message = 'unknown dynamic format "'//trim(value)//'" (use text or hdf5)'
                     return
                  end select
                  self%dyn_format_given = .true.
               case ('--s4')
                  self%s4_output = trim(value)
               case ('--chi4')
                  self%chi4_output = trim(value)
               case ('--msd')
                  self%msd_output = trim(value)
               case ('--s4-cutoff')
                  read (value, *, iostat=ierr) self%s4_cutoff
                  if (ierr /= 0 .or. self%s4_cutoff <= 0.0_rk) then
                     ierr = 1
                     message = '--s4-cutoff must be a positive number'
                     return
                  end if
                  self%s4_cutoff_given = .true.
               case ('--buffer-limit')
                  read (value, *, iostat=ierr) self%buffer_limit_gb
                  if (ierr /= 0 .or. self%buffer_limit_gb <= 0.0_rk) then
                     ierr = 1
                     message = '--buffer-limit must be a positive number of GB'
                     return
                  end if
                  self%buffer_limit_given = .true.
               case ('--stride')
                  read (value, *, iostat=ierr) self%stride
                  if (ierr /= 0 .or. self%stride < 1) then
                     ierr = 1
                     message = '--stride must be a positive integer'
                     return
                  end if
                  self%stride_given = .true.
               case ('--device')
                  select case (trim(value))
                  case ('cpu')
                     self%device = device_cpu
                  case ('gpu', 'cuda')
                     self%device = device_gpu
                  case default
                     ierr = 1
                     message = 'unknown device "'//trim(value)//'" (use cpu or gpu)'
                     return
                  end select
               case ('--gpu-id')
                  read (value, *, iostat=ierr) self%gpu_id
                  if (ierr /= 0 .or. self%gpu_id < 0) then
                     ierr = 1
                     message = 'GPU id must be a non-negative integer'
                     return
                  end if
               case ('--norm')
                  select case (trim(value))
                  case ('mean', 'fz')
                     self%norm = norm_mean
                  case ('self')
                     self%norm = norm_self
                  case ('n', 'natom')
                     self%norm = norm_natom
                  case default
                     ierr = 1
                     message = 'unknown normalization "'//trim(value)//'" (use mean, self or n)'
                     return
                  end select
               case ('--grid')
                  self%want_grid = .true.
                  self%grid_output = trim(value)
               case ('--grid-format')
                  select case (trim(value))
                  case ('text', 'txt', 'ascii')
                     self%grid_format = grid_format_text
                  case ('hdf5', 'h5', 'hdf')
                     self%grid_format = grid_format_hdf5
                  case default
                     ierr = 1
                     message = 'unknown grid format "'//trim(value)//'" (use text or hdf5)'
                     return
                  end select
                  self%grid_format_given = .true.
               case ('--xrd')
                  self%xrd_output = trim(value)
               case ('--xrd-lambda')
                  read (value, *, iostat=ierr) self%xrd_lambda
                  if (ierr /= 0 .or. self%xrd_lambda <= 0.0_rk) then
                     ierr = 1
                     message = '--xrd-lambda must be a positive wavelength'
                     return
                  end if
                  self%xrd_lambda_given = .true.
               case ('--xrd-range')
                  read (value, *, iostat=ierr) self%xrd_2theta_min
                  if (ierr == 0) read (value2, *, iostat=ierr) self%xrd_2theta_max
                  if (ierr /= 0) then
                     ierr = 1
                     message = '--xrd-range takes two numbers: MIN2TH MAX2TH in degrees'
                     return
                  end if
                  self%xrd_range_given = .true.
               case ('--xrd-step')
                  read (value, *, iostat=ierr) self%xrd_step
                  if (ierr /= 0 .or. self%xrd_step <= 0.0_rk) then
                     ierr = 1
                     message = '--xrd-step must be a positive number of degrees'
                     return
                  end if
                  self%xrd_step_given = .true.
               end select
            end if
         else
            npos = npos + 1
            if (npos > size(positional)) then
               ierr = 1
               message = 'too many positional arguments'
               return
            end if
            positional(npos) = trim(arg)
         end if
         i = i + 1
      end do

      if (self%show_help .or. self%show_version) return
      ! A stray value after --dyn is the old command line, which carried the q
      ! line itself; name the replacement instead of treating it as OUTPUT.
      if (self%dynamic .and. self%dyn_q_mode == dyn_q_none) then
         do i = 1, npos
            if (index(positional(i), ',') > 0) then
               ierr = 1
               message = '--dyn takes no value; pick the q sampling with --dyn-q, e.g. '// &
                  '--dyn-q line:'//trim(positional(i))
               return
            end if
         end do
      end if
      if (.not. allocated(self%input)) then
         ierr = 1
         message = 'missing input dump file (-i)'
         return
      end if
      if (self%qmax <= self%qmin) then
         ierr = 1
         message = 'qmax must be larger than qmin'
         return
      end if
      ! --- powder XRD output ------------------------------------------------
      if (allocated(self%xrd_output)) then
         if (self%dynamic) then
            ierr = 1
            message = '--xrd is a static output; it cannot be combined with --dyn'
            return
         end if
         if (self%method == method_debye) then
            ierr = 1
            message = '--method debye sums the orientation-averaged Debye intensity, which is a '// &
               'different quantity from the reciprocal lattice sum --xrd reports (it carries '// &
               'no multiplicity); use the default nufft method, or --method direct as its reference'
            return
         end if
         if (.not. self%xrd_lambda_given) then
            ierr = 1
            message = '--xrd needs --xrd-lambda LAMBDA (the incident wavelength)'
            return
         end if
         if (self%qmin_given .or. self%qmax_given) then
            ierr = 1
            message = '--xrd takes the q range from --xrd-lambda and --xrd-range; '// &
               'do not pass --qmin or --qmax with it'
            return
         end if
         if (self%xrd_2theta_min <= 0.0_rk .or. self%xrd_2theta_max >= 180.0_rk .or. &
             self%xrd_2theta_max <= self%xrd_2theta_min) then
            ierr = 1
            message = '--xrd-range needs 0 < MIN2TH < MAX2TH < 180 degrees'
            return
         end if
         ! A powder x-ray pattern unless the user picked a weighting scheme;
         ! the neutron and unit schemes give legitimate patterns too.
         if (.not. self%weight_given) self%scheme%kind = weight_xray
         ! q = 4 pi sin(theta)/lambda at the two ends of the requested range.
         self%qmin = 4.0_rk*pi*sin(0.5_rk*deg2rad*self%xrd_2theta_min)/self%xrd_lambda
         self%qmax = 4.0_rk*pi*sin(0.5_rk*deg2rad*self%xrd_2theta_max)/self%xrd_lambda
      else if (self%xrd_lambda_given .or. self%xrd_range_given .or. self%xrd_step_given .or. &
               self%lp_given) then
         ierr = 1
         message = '--xrd-lambda, --xrd-range, --xrd-step and --lp/--no-lp '// &
            'belong to --xrd FILE'
         return
      end if
      if (self%scheme%kind /= weight_unit .and. .not. self%scheme%has_mapping()) then
         ierr = 1
         message = 'weighted schemes need an element mapping, e.g. -m 1:Si,2:O'
         return
      end if
      if (self%device == device_gpu .and. self%method == method_direct) then
         ierr = 1
         message = 'the direct method runs on the CPU; use --device cpu or --method nufft'
         return
      end if
      if (self%precision == precision_single .and. self%device == device_cpu) then
         ierr = 1
         message = 'single precision is only available on the GPU path (--device gpu)'
         return
      end if
      if (self%dynamic) then
         if (self%method_given) then
            ierr = 1
            message = '--dyn selects the dynamic method; do not pass --method as well'
            return
         end if
         self%method = method_dynamic
         if (self%dt <= 0.0_rk) then
            ierr = 1
            message = '--dyn needs --dt DT, the time step of the trajectory'
            return
         end if
         if (self%maxframes < 1) then
            ierr = 1
            message = '--dyn needs --maxframes L, the correlation window in frames'
            return
         end if
         if (self%want_grid) then
            ierr = 1
            message = 'the dynamic method samples a q line or shell; --grid is not available'
            return
         end if
         if (self%device == device_gpu) then
            ierr = 1
            message = 'the dynamic method runs on the CPU; use --device cpu'
            return
         end if
         if (allocated(self%s4_output) .or. allocated(self%chi4_output) .or. &
             allocated(self%fqt_self_output) .or. allocated(self%msd_output)) then
            if ((allocated(self%s4_output) .or. allocated(self%chi4_output)) .and. &
                .not. self%s4_cutoff_given) then
               ierr = 1
               message = '--s4 and --chi4 need --s4-cutoff A'
               return
            end if
            if (self%s4_cutoff_given .and. .not. (allocated(self%s4_output) .or. &
                allocated(self%chi4_output))) then
               ierr = 1
               message = '--s4-cutoff belongs to --s4 or --chi4'
               return
            end if
            if (self%stride > self%maxframes) then
               ierr = 1
               message = '--stride cannot exceed --maxframes'
               return
            end if
         else if (self%s4_cutoff_given .or. self%buffer_limit_given .or. self%stride_given) then
            ierr = 1
            message = '--s4-cutoff, --buffer-limit and --stride need --s4, --chi4, '// &
               '--fqt-self or --msd'
            return
         end if
         if ((self%dyn_modes_given .or. self%dyn_thin_given) .and. &
             self%dyn_q_mode /= dyn_q_grid) then
            ierr = 1
            message = '--dyn-modes and --dyn-thin need --dyn-q grid:QMAX'
            return
         end if
         if (self%dyn_keep_modes .and. self%dyn_q_mode /= dyn_q_grid) then
            ierr = 1
            message = '--dyn-keep-modes needs --dyn-q grid:QMAX'
            return
         end if
      else if (self%dyn_q_given) then
         ierr = 1
         message = '--dyn-q belongs to --dyn; add --dyn to the command line'
         return
      else if (allocated(self%sqw_output) .or. allocated(self%fqt_output) .or. &
               allocated(self%fqt_self_output) .or. allocated(self%s4_output) .or. &
               allocated(self%chi4_output) .or. allocated(self%msd_output)) then
         ierr = 1
         message = '--sqw, --fqt, --fqt-self, --s4, --chi4 and --msd belong to --dyn'
         return
      else if (self%dt > 0.0_rk .or. self%maxframes > 0 .or. self%lag_given .or. &
               self%dyn_format_given .or. self%s4_cutoff_given .or. &
               self%buffer_limit_given .or. self%stride_given .or. &
               self%dyn_modes_given .or. self%dyn_thin_given .or. self%dyn_keep_modes) then
         ierr = 1
         message = '--dt, --maxframes, --lag, --dyn-format, --s4-cutoff, --buffer-limit '// &
            'and --stride, --dyn-modes and --dyn-thin belong to --dyn'
         return
      end if
      if (self%method == method_debye) then
         if (self%want_grid) then
            ierr = 1
            message = 'the Debye method evaluates S(q) directly; --grid is not available'
            return
         end if
         if (self%device == device_gpu) then
            ierr = 1
            message = 'the Debye method runs on the CPU; use --device cpu'
            return
         end if
         if (allocated(self%s2_accum_output) .and. .not. allocated(self%pair_entropy_output)) then
            ierr = 1
            message = '--s2-accum needs --pair-entropy'
            return
         end if
      else
         if (self%rmax_given .or. self%dr_given .or. self%skin_given .or. &
             allocated(self%rdf_output) .or. allocated(self%pair_entropy_output) .or. &
             allocated(self%s2_accum_output)) then
            ierr = 1
            message = '--rmax, --dr, --skin, --rdf, --pair-entropy and --s2-accum '// &
               'belong to --method debye'
            return
         end if
      end if

      ! The positional S(q) table is judged once the flags are consistent: a
      ! run without q points has nothing to tabulate.
      if (self%dynamic .and. self%dyn_q_mode == dyn_q_none) then
         if (allocated(self%sqw_output) .or. allocated(self%fqt_output) .or. &
             allocated(self%fqt_self_output) .or. allocated(self%s4_output)) then
            ierr = 1
            message = '--sqw, --fqt, --fqt-self and --s4 need a q sampling; add '// &
               '--dyn-q line:... or --dyn-q shell:...'
            return
         end if
         if (.not. allocated(self%chi4_output) .and. .not. allocated(self%msd_output)) then
            ierr = 1
            message = '--dyn-q - computes only --chi4 and --msd; add --chi4 or --msd FILE, '// &
               'or a q sampling'
            return
         end if
      end if
      if (npos > 1) then
         ierr = 1
         message = 'only one output argument is allowed'
         return
      end if
      if (npos == 0) then
         ! Only a run that writes some other table can do without the S(q)
         ! table: a dynamic run without a q sampling, or an XRD-only run.
         if (.not. (self%dynamic .and. self%dyn_q_mode == dyn_q_none) .and. &
             .not. allocated(self%xrd_output)) then
            ierr = 1
            message = 'missing output argument (use - for stdout)'
            return
         end if
      else
         if (self%dynamic .and. self%dyn_q_mode == dyn_q_none) then
            ierr = 1
            message = '--dyn-q - computes no S(q) table; drop the OUTPUT argument'
            return
         end if
         self%output = trim(positional(1))
      end if
      if (.not. self%want_grid) self%grid_output = ''
      ! Default the grid format from the file name.
      if (self%want_grid .and. .not. self%grid_format_given) then
         if (ends_with(self%grid_output, '.h5') .or. ends_with(self%grid_output, '.hdf5')) then
            self%grid_format = grid_format_hdf5
         end if
      end if
      ! The XRD table infers its format from the file name as well.
      if (allocated(self%xrd_output)) then
         if (ends_with(self%xrd_output, '.h5') .or. ends_with(self%xrd_output, '.hdf5')) then
            self%xrd_format = grid_format_hdf5
         end if
      end if
      if (self%dyn_format_given) then
         self%sqw_format = self%dyn_format
         self%fqt_format = self%dyn_format
         self%fqt_self_format = self%dyn_format
         self%s4_format = self%dyn_format
         self%chi4_format = self%dyn_format
         self%msd_format = self%dyn_format
      else
         if (allocated(self%sqw_output)) then
            if (ends_with(self%sqw_output, '.h5') .or. ends_with(self%sqw_output, '.hdf5')) &
               self%sqw_format = grid_format_hdf5
         end if
         if (allocated(self%fqt_output)) then
            if (ends_with(self%fqt_output, '.h5') .or. ends_with(self%fqt_output, '.hdf5')) &
               self%fqt_format = grid_format_hdf5
         end if
         if (allocated(self%fqt_self_output)) then
            if (ends_with(self%fqt_self_output, '.h5') .or. &
                ends_with(self%fqt_self_output, '.hdf5')) self%fqt_self_format = grid_format_hdf5
         end if
         if (allocated(self%s4_output)) then
            if (ends_with(self%s4_output, '.h5') .or. ends_with(self%s4_output, '.hdf5')) &
               self%s4_format = grid_format_hdf5
         end if
         if (allocated(self%chi4_output)) then
            if (ends_with(self%chi4_output, '.h5') .or. ends_with(self%chi4_output, '.hdf5')) &
               self%chi4_format = grid_format_hdf5
         end if
         if (allocated(self%msd_output)) then
            if (ends_with(self%msd_output, '.h5') .or. ends_with(self%msd_output, '.hdf5')) &
               self%msd_format = grid_format_hdf5
         end if
      end if
   end subroutine parse_options

   !> Parse the q sampling of a dynamic run (`--dyn-q SPEC`).
   !!
   !! `-` keeps the time axis without sampling q at all, so only the scalar
   !! overlap Q(t) and chi4(t) can be computed.  `line:NINT,S0,S1,DX,DY,DZ` is
   !! the density amplitude along a line in reciprocal space through Gamma, and
   !! `shell:Q,ACC` averages every direction of the shell |q| = Q on a Lebedev
   !! grid of the requested accuracy.
   subroutine parse_dyn_q(self, spec, ierr, message)
      type(options_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message

      ierr = 0
      message = ''
      if (trim(spec) == '-') then
         self%dyn_q_mode = dyn_q_none
      else if (starts_with(spec, 'line:')) then
         call parse_dyn_line(self, spec(6:), ierr, message)
         if (ierr /= 0) return
      else if (starts_with(spec, 'shell:')) then
         call parse_dyn_shell(self, spec(7:), ierr, message)
         if (ierr /= 0) return
      else if (starts_with(spec, 'grid:')) then
         call parse_dyn_grid(self, spec(6:), ierr, message)
         if (ierr /= 0) return
      else if (starts_with(spec, 'single:')) then
         call parse_dyn_single(self, spec(8:), ierr, message)
         if (ierr /= 0) return
      else
         ierr = 1
         message = '--dyn-q wants "-", "line:NINT,S0,S1,DX,DY,DZ" or '// &
            '"shell:Q,low|medium|high", "grid:QMAX" or "single:N1,N2,N3"'
         return
      end if
      self%dyn_q_given = .true.
   end subroutine parse_dyn_q

   !> Parse "NINT,S0,S1,DX,DY,DZ" of `--dyn-q line`.
   !!
   !! NINT is the number of intervals of the q line (so NINT+1 q points), S0/S1
   !! the range of its scale in 1/A, and DX,DY,DZ the (unnormalized) direction.
   subroutine parse_dyn_line(self, spec, ierr, message)
      type(options_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=48) :: fields(6)
      real(rk) :: values(6)
      integer :: j, n

      ierr = 0
      message = ''
      fields = ' '
      values = 0.0_rk
      call split_fields(spec, fields, n)
      if (n /= 6) then
         ierr = 1
         message = '--dyn-q line wants NINT,S0,S1,DX,DY,DZ, e.g. '// &
            '--dyn-q line:100,0.5,20,1,1,0'
         return
      end if
      do j = 1, 6
         read (fields(j), *, iostat=ierr) values(j)
         if (ierr /= 0) then
            ierr = 1
            message = 'cannot read "'//trim(fields(j))//'" as a number in --dyn-q line'
            return
         end if
      end do
      self%dyn_intervals = nint(values(1))
      self%dyn_s0 = values(2)
      self%dyn_s1 = values(3)
      self%dyn_dir = values(4:6)
      if (self%dyn_intervals < 1) then
         ierr = 1
         message = '--dyn-q line needs at least one interval on the q line'
         return
      end if
      if (self%dyn_s1 <= self%dyn_s0 .or. self%dyn_s0 < 0.0_rk) then
         ierr = 1
         message = '--dyn-q line needs 0 <= S0 < S1 for the scale of the q line'
         return
      end if
      if (sum(self%dyn_dir**2) <= 0.0_rk) then
         ierr = 1
         message = '--dyn-q line needs a non-zero direction, e.g. 1,1,0'
         return
      end if
      self%dyn_q_mode = dyn_q_line
   end subroutine parse_dyn_line

   !> Parse "Q,ACC" of `--dyn-q shell`.
   !!
   !! Q is the radius of the shell in 1/A and ACC one of the accuracy names of
   !! the Lebedev rules: low, medium or high.
   subroutine parse_dyn_shell(self, spec, ierr, message)
      type(options_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=48) :: fields(2)
      real(rk) :: radius
      integer :: order, n

      ierr = 0
      message = ''
      fields = ' '
      call split_fields(spec, fields, n)
      if (n /= 2) then
         ierr = 1
         message = '--dyn-q shell wants Q,low|medium|high, e.g. '// &
            '--dyn-q shell:2.5,medium'
         return
      end if
      read (fields(1), *, iostat=ierr) radius
      if (ierr /= 0) then
         ierr = 1
         message = 'cannot read "'//trim(fields(1))//'" as the |q| of --dyn-q shell'
         return
      end if
      if (radius <= 0.0_rk) then
         ierr = 1
         message = '--dyn-q shell needs a positive |q| radius'
         return
      end if
      order = lebedev_order_from_name(trim(fields(2)))
      if (order < 0) then
         ierr = 1
         message = 'unknown --dyn-q shell accuracy "'//trim(fields(2))// &
            '" (use low, medium or high)'
         return
      end if
      self%dyn_shell_q = radius
      self%dyn_shell_order = order
      self%dyn_q_mode = dyn_q_shell
   end subroutine parse_dyn_shell

   !> Parse "QMAX" of `--dyn-q grid`.
   !!
   !! QMAX is the upper bound of |q| in 1/A.  Every reciprocal-lattice vector
   !! with 0 < |q| <= QMAX is sampled, together with the Gamma point.  Only
   !! lattice vectors carry a density amplitude that is independent of how the
   !! periodic images are chosen, so a grid is free of the box-form-factor
   !! contamination of the off-lattice Lebedev shell; its mode count follows
   !! from the cell and can be capped with `--dyn-modes`.
   subroutine parse_dyn_grid(self, spec, ierr, message)
      type(options_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: qmax

      ierr = 0
      message = ''
      read (spec, *, iostat=ierr) qmax
      if (ierr /= 0) then
         ierr = 1
         message = 'cannot read "'//trim(spec)//'" as the qmax of --dyn-q grid'
         return
      end if
      if (qmax <= 0.0_rk) then
         ierr = 1
         message = '--dyn-q grid needs a positive qmax'
         return
      end if
      self%dyn_grid_qmax = qmax
      self%dyn_q_mode = dyn_q_grid
   end subroutine parse_dyn_grid

   !> Parse "N1,N2,N3" of `--dyn-q single`.
   !!
   !! The three integers are the Miller indices of one reciprocal-lattice
   !! vector of the dump box,
   !!
   !!   q = n1 b1 + n2 b2 + n3 b3.
   !!
   !! Building q from the integers keeps it exactly on the lattice, which is
   !! what S4 needs: a hand written |q| would sit a few ulps away, the density
   !! amplitude would stop being independent of the periodic images, and the
   !! box form factor would come back.
   subroutine parse_dyn_single(self, spec, ierr, message)
      type(options_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=48) :: fields(4)
      integer :: n, i, values(3)

      ierr = 0
      message = ''
      fields = ' '
      call split_fields(spec, fields, n)
      if (n /= 3) then
         ierr = 1
         message = '--dyn-q single wants N1,N2,N3, e.g. --dyn-q single:1,0,0'
         return
      end if
      do i = 1, 3
         if (.not. is_integer_field(fields(i))) then
            ierr = 1
            message = 'cannot read "'//trim(fields(i))//'" as an integer of --dyn-q single'
            return
         end if
         read (fields(i), *, iostat=ierr) values(i)
         if (ierr /= 0) then
            ierr = 1
            message = 'the index "'//trim(fields(i))//'" of --dyn-q single is out of range'
            return
         end if
      end do
      self%dyn_single = values
      self%dyn_q_mode = dyn_q_single
   end subroutine parse_dyn_single

   !> True when `text` is an optionally signed integer and nothing else.
   !!
   !! List directed input would read "1.5" as 1 and silently ignore the rest,
   !! so the field is checked character by character first.
   pure logical function is_integer_field(text) result(ok)
      character(len=*), intent(in) :: text
      integer :: i, first, last
      logical :: seen_digit

      ok = .false.
      last = len_trim(text)
      first = 0
      do i = 1, last
         if (text(i:i) /= ' ') then
            first = i
            exit
         end if
      end do
      if (first == 0) return
      seen_digit = .false.
      do i = first, last
         select case (text(i:i))
         case ('+', '-')
            if (i /= first) return
         case ('0':'9')
            seen_digit = .true.
         case default
            return
         end select
      end do
      ok = seen_digit
   end function is_integer_field

   !> Split a comma separated option value into at most `size(fields)` tokens.
   !!
   !! `n` is the number of tokens, so a caller that expects a fixed count also
   !! rejects a value with too many fields (`n` becomes `size(fields) + 1`).
   subroutine split_fields(text, fields, n)
      character(len=*), intent(in) :: text
      character(len=*), intent(out) :: fields(:)
      integer, intent(out) :: n
      integer :: i, start, last
      logical :: at_end

      fields = ' '
      n = 0
      start = 1
      last = len_trim(text)
      do i = 1, last + 1
         ! Fortran does not short-circuit .or., so the end of the string has to
         ! be tested separately before text(i:i) is evaluated.
         at_end = i > last
         if (.not. at_end) at_end = text(i:i) == ','
         if (at_end) then
            if (i > start) then
               n = n + 1
               if (n > size(fields)) exit
               fields(n) = text(start:i - 1)
            end if
            start = i + 1
         end if
      end do
   end subroutine split_fields

   !> Case sensitive prefix test for the `--dyn-q` keywords.
   pure logical function starts_with(text, prefix) result(found)
      character(len=*), intent(in) :: text, prefix
      integer :: n

      n = len(prefix)
      found = .false.
      if (len(text) < n) return
      found = text(1:n) == prefix
   end function starts_with

   !> Case sensitive file suffix test used for the grid format default.
   pure logical function ends_with(text, suffix) result(found)
      character(len=*), intent(in) :: text, suffix
      integer :: n, m
      n = len_trim(text)
      m = len(suffix)
      found = n >= m .and. text(n - m + 1:n) == suffix
   end function ends_with

   subroutine print_usage(unit)
      integer, intent(in) :: unit
      write (unit, '(a)') program_version
      write (unit, '(a)') ''
      write (unit, '(a)') 'Total structure factor S(q) from LAMMPS dump files.'
      write (unit, '(a)') ''
      write (unit, '(a)') 'usage: sqcalc -i DUMP [options] OUTPUT'
      write (unit, '(a)') ''
      write (unit, '(a)') 'OUTPUT is the shell averaged S(q) table, use - for stdout.'
      write (unit, '(a)') ''
      write (unit, '(a)') 'options:'
      write (unit, '(a)') '  -i, --input FILE    LAMMPS dump trajectory (required)'
      write (unit, '(a)') '  -m, --mapping LIST  LAMMPS type id to element or species, e.g.'
      write (unit, '(a)') '                      1:Si,2:O or 1:Si4+,2:O2- (ions, see the skill docs)'
      write (unit, '(a)') '  -w, --weight SCHEME unit (default), neutron or xray'
      write (unit, '(a)') '  -t, --threads N     OpenMP threads (default: all available)'
      write (unit, '(a)') '      --qmin VALUE    smallest |q| in the output [1/A] (default 0)'
      write (unit, '(a)') '      --qmax VALUE    largest |q| in the output [1/A] (default 20)'
      write (unit, '(a)') '      --nq N          number of q shells (default 500)'
      write (unit, '(a)') '      --grid FILE     also write S(q) on every reciprocal lattice point'
      write (unit, '(a)') '      --grid-format NAME  text (default) or hdf5 (.h5/.hdf5 implies hdf5)'
      write (unit, '(a)') '      --xrd FILE      also write a powder XRD pattern (.h5 = HDF5):'
      write (unit, '(a)') '                      I(2theta) = sum over the reciprocal lattice points'
      write (unit, '(a)') '                      in each bin of |rho(q)|^2 LP(2theta), per atom'
      write (unit, '(a)') '      --xrd-lambda VALUE  incident wavelength [dump length unit]'
      write (unit, '(a)') '      --xrd-range MIN MAX two-theta range [deg] (default 1 179)'
      write (unit, '(a)') '      --xrd-step VALUE    two-theta bin width [deg] (default: from'
      write (unit, '(a)') '                      the box, at most as coarse as the 2theta ='
      write (unit, '(a)') '                      120 deg spacing, so that no bin is empty)'
      write (unit, '(a)') '      --lp, --no-lp   apply (default) or drop the Lorentz-polarization'
      write (unit, '(a)') '                      factor of the XRD pattern'
      write (unit, '(a)') '      --method NAME   nufft (default) or direct'
      write (unit, '(a)') '                      debye: real space pair histograms'
      write (unit, '(a)') '      --rmax VALUE    Debye pair cutoff [A] (default: half the'
      write (unit, '(a)') '                      smallest periodic box side)'
      write (unit, '(a)') '      --dr VALUE      Debye radial bin width [A] (default 0.01)'
      write (unit, '(a)') '      --skin VALUE    Verlet skin for the pair list [A] (default 1.0)'
      write (unit, '(a)') '      --rdf FILE      total and partial g(r) in one file (.h5 = HDF5)'
      write (unit, '(a)') '      --pair-entropy FILE  total and partial pair entropy S2/kB'
      write (unit, '(a)') '      --s2-accum FILE  S2(r) accumulation curve for tail extrapolation'
      write (unit, '(a)') '      --no-cutoff-correction  disable the Debye cut-off density correction'
      write (unit, '(a)') '      --dyn           keep the time axis: dynamic structure factor and'
      write (unit, '(a)') '                      the four-point structure factor and overlap'
      write (unit, '(a)') '      --dyn-q SPEC    q sampling of --dyn; default "-" (no q points):'
      write (unit, '(a)') '                      -   no q points: only --chi4 and --msd are available'
      write (unit, '(a)') '                      line:NINT,S0,S1,DX,DY,DZ  q line through Gamma,'
      write (unit, '(a)') '                      NINT intervals, scale S0..S1 [1/A], direction'
      write (unit, '(a)') '                      shell:Q,low|medium|high  Lebedev average on |q| = Q'
      write (unit, '(a,i0,a,i0,a,i0,a)') '                      (low, medium and high are the ', &
         lebedev_points(lebedev_low), ', ', lebedev_points(lebedev_medium), ' and ', &
         lebedev_points(lebedev_high), ' point rules)'
      write (unit, '(a)') '                      grid:QMAX  every reciprocal lattice vector |q| <= QMAX'
      write (unit, '(a)') '                      single:N1,N2,N3  one lattice vector of the box'
      write (unit, '(a)') '      --dt VALUE      MD time step of the trajectory (with --dyn)'
      write (unit, '(a)') '      --dyn-modes N   mode budget of grid:QMAX, 0 = unlimited (default)'
      write (unit, '(a)') '      --dyn-thin NAME order of the grid thinning: shells (default) or orbits'
      write (unit, '(a)') '      --dyn-keep-modes  write the grid rows per lattice vector instead of'
      write (unit, '(a)') '                      the |q| shell average (for diagnostics)'
      write (unit, '(a)') '      --maxframes N   correlation window in frames (with --dyn)'
      write (unit, '(a)') '      --lag N         frames between consecutive time origins (default 1)'
      write (unit, '(a)') '      --sqw FILE      S(q,w) spectra, one row per (q,w) (.h5 = HDF5)'
      write (unit, '(a)') '      --fqt FILE      coherent F(q,t) intermediate scattering function'
      write (unit, '(a)') '      --fqt-self FILE self F_s(q,t) intermediate scattering function'
      write (unit, '(a)') '      --dyn-format NAME  text (default) or hdf5 for all dynamic outputs'
      write (unit, '(a)') '      --s4 FILE       S4(q,t) four-point structure factor (.h5 = HDF5)'
      write (unit, '(a)') '      --chi4 FILE     Q(t) and chi4(t) average overlap / susceptibility'
      write (unit, '(a)') '      --msd FILE      MSD(t) mean squared displacement (.h5 = HDF5)'
      write (unit, '(a)') '      --s4-cutoff A   overlap cutoff for --s4 and --chi4 [dump length unit]'
      write (unit, '(a)') '      --buffer-limit GB  position buffer limit for S4/chi4/F_s/MSD (default 2.0)'
      write (unit, '(a)') '      --stride N      use every N-th frame for S4/chi4/F_s/MSD (default 1)'
      write (unit, '(a)') '  -fz, --faber-ziman  partials in the Faber-Ziman normalization'
      write (unit, '(a)') '      --partials      write the partial structure factor columns (default)'
      write (unit, '(a)') '      --no-partials   do not write partial structure factor columns'
      write (unit, '(a)') '      --device NAME   cpu (default) or gpu (cufinufft + cuFFT)'
      write (unit, '(a)') '      --gpu-id N      CUDA device to use (default 0)'
      write (unit, '(a)') '      --precision NAME  double (default) or single (float32 GPU)'
      write (unit, '(a)') '      --norm NAME     mean (default), self or n'
      write (unit, '(a)') '      --eps VALUE     NUFFT tolerance (default 1e-9)'
      write (unit, '(a)') '  -q, --quiet         do not write progress information to stderr'
      write (unit, '(a)') '  -h, --help          show this help'
      write (unit, '(a)') '  -v, --version       show the program version'
   end subroutine print_usage

end module sqc_options
