! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Command line handling for sqcalc.
module sqc_options
   use sqc_kinds
   use sqc_output, only: format_text, format_hdf5
   use sqc_weights, only: weight_scheme_t, weight_unit, weight_neutron, weight_xray, &
                          scheme_from_name
   use sqc_structure_factor, only: method_nufft, method_direct, method_debye, norm_mean, &
                            norm_self, norm_natom
   use sqc_debye, only: debye_default_dr, debye_default_skin
   use sqc_dynamics, only: dyn_q_none, dyn_q_line, dyn_q_shell, dyn_q_grid, dyn_q_single, &
                           dyn_q_powder
   use sqc_modes, only: thin_none, thin_shells, thin_orbits
   use sqc_lebedev, only: lebedev_points, lebedev_order_from_name, &
                          lebedev_low, lebedev_medium, lebedev_high
   implicit none
   private

   public :: options_t, parse_options, print_usage, program_version, device_cpu, device_gpu, &
             precision_double, precision_single, cmd_none, cmd_static, cmd_dyn

   !> Which subcommand the run belongs to; cmd_none is a bare -h/--version.
   integer, parameter :: cmd_none = -1
   integer, parameter :: cmd_static = 0
   integer, parameter :: cmd_dyn = 1

   !> Which subcommands accept an option (bit set), see flag_scope.
   integer, parameter :: scope_static = 1
   integer, parameter :: scope_dyn = 2
   integer, parameter :: scope_both = 3

   !> Where the NUFFT transforms run.
   integer, parameter :: device_cpu = 0
   integer, parameter :: device_gpu = 1

   !> Precision of the GPU transform.
   integer, parameter :: precision_double = 0
   integer, parameter :: precision_single = 1

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
      integer :: grid_format = format_text
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
      integer :: xrd_format = format_text
      logical :: lp = .true.
      logical :: lp_given = .false.
      !> True when -w was given; --xrd implies the x-ray weights otherwise.
      logical :: weight_given = .false.
      logical :: quiet = .false.
      logical :: show_help = .false.
      logical :: show_version = .false.
      !> Subcommand: cmd_static, cmd_dyn, or cmd_none for a bare -h/--version.
      integer :: command = cmd_none
      !> The q sampling of the dyn subcommand.
      integer :: q_mode = dyn_q_none
      !> line: NINT intervals, scale from S0 to S1 along (DX, DY, DZ).
      integer :: q_intervals = 0
      real(rk) :: q_s0 = 0.0_rk
      real(rk) :: q_s1 = 0.0_rk
      real(rk) :: q_dir(3) = 0.0_rk
      !> shell: |q| radius and the order of the Lebedev rule it selects.
      real(rk) :: shell_q = 0.0_rk
      integer :: shell_order = 0
      !> grid: the upper bound of |q|, the mode budget and the thinning policy.
      real(rk) :: grid_qmax = 0.0_rk
      !> powder: the requested |q|, the target number of lattice vectors in
      !! the window, and the half width when one was given.
      real(rk) :: powder_q = 0.0_rk
      integer :: powder_modes = 50
      real(rk) :: powder_dq = 0.0_rk
      logical :: powder_dq_given = .false.
      !> single: the Miller indices of the one lattice vector to sample.
      integer :: single_index(3) = 0
      integer :: modes = 0
      integer :: thin = thin_none
      logical :: modes_given = .false.
      logical :: thin_given = .false.
      logical :: keep_modes = .false.
      !> MD time step [time units] and the correlation window settings.
      real(rk) :: dt = 0.0_rk
      integer :: maxframes = 0
      integer :: lag_stride = 1
      logical :: lag_given = .false.
      !> Optional dynamic outputs.
      character(len=:), allocatable :: sq_output
      character(len=:), allocatable :: sqw_output
      character(len=:), allocatable :: fqt_output
      character(len=:), allocatable :: fqt_self_output
      !> --format override, for every table the run writes; otherwise each
      !! file infers the format from its suffix.
      integer :: output_format = format_text
      logical :: output_format_given = .false.
      integer :: sqw_format = format_text
      integer :: fqt_format = format_text
      integer :: fqt_self_format = format_text
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
      integer :: s4_format = format_text
      integer :: chi4_format = format_text
      integer :: msd_format = format_text
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

      ! --- the subcommand is the first argument ----------------------------
      if (nargs < 1) then
         ierr = 1
         message = 'expected a subcommand: "sqcalc static ..." or "sqcalc dyn ..."'
         return
      end if
      call get_command_argument(1, arg)
      arg = trim(adjustl(arg))
      select case (arg)
      case ('static')
         self%command = cmd_static
      case ('dyn')
         self%command = cmd_dyn
      case ('-h', '--help')
         self%show_help = .true.
         return
      case ('-v', '--version')
         self%show_version = .true.
         return
      case default
         ierr = 1
         message = 'expected a subcommand as the first argument, got "'//trim(arg)// &
            '"; use "sqcalc static ..." or "sqcalc dyn ..." (see sqcalc -h)'
         return
      end select

      i = 2
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

            ! A flag of the other subcommand is refused here, so the message
            ! can name the subcommand that does take it.
            if (iand(flag_scope(trim(name)), &
                     merge(scope_static, scope_dyn, self%command == cmd_static)) == 0) then
               ierr = 1
               if (self%command == cmd_static) then
                  message = trim(name)//' is not a static option; it belongs to '// &
                     '"sqcalc dyn"'
               else
                  message = trim(name)//' is not a dyn option; it belongs to '// &
                     '"sqcalc static"'
               end if
               return
            end if

            needs_value = .true.
            select case (trim(name))
            case ('-h', '--help')
               self%show_help = .true.
               needs_value = .false.
            case ('-v', '--version')
               self%show_version = .true.
               needs_value = .false.
            case ('--quiet')
               self%quiet = .true.
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
            case ('--keep-modes')
               self%keep_modes = .true.
               needs_value = .false.
            case ('-i', '--input', '--mapping', '-m', '-w', '--weight', '-t', '--threads', &
                  '--qmin', '--qmax', '--nq', '--eps', '--method', '--norm', '--grid', &
                  '--device', '--gpu-id', '--precision', &
                  '--xrd', '--xrd-lambda', '--xrd-range', '--xrd-step', &
                  '--rmax', '--dr', '--skin', '--rdf', &
                  '--pair-entropy', '--s2-accum', &
                  '-q', '--qpoints', '--dt', '--maxframes', '--lag', '--sqw', '--fqt', '--format', &
                  '--fqt-self', '--s4', '--chi4', '--msd', '--s4-cutoff', '--buffer-limit', &
                  '--stride', &
                  '--modes', '--thin', '--sq')
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
               case ('-q', '--qpoints')
                  call parse_q_sampling(self, trim(value), ierr, message)
                  if (ierr /= 0) return
               case ('--modes')
                  read (value, *, iostat=ierr) self%modes
                  if (ierr /= 0 .or. self%modes < 0) then
                     ierr = 1
                     message = '--modes must be a non-negative integer'
                     return
                  end if
                  self%modes_given = .true.
               case ('--thin')
                  select case (trim(value))
                  case ('shells')
                     self%thin = thin_shells
                  case ('orbits')
                     self%thin = thin_orbits
                  case default
                     ierr = 1
                     message = '--thin wants "shells" or "orbits"'
                     return
                  end select
                  self%thin_given = .true.
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
               case ('--sq')
                  self%sq_output = trim(value)
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
               case ('--format')
                  select case (trim(value))
                  case ('text', 'txt', 'ascii')
                     self%output_format = format_text
                  case ('hdf5', 'h5', 'hdf')
                     self%output_format = format_hdf5
                  case default
                     ierr = 1
                     message = 'unknown format "'//trim(value)//'" (use text or hdf5)'
                     return
                  end select
                  self%output_format_given = .true.
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
      if (self%command == cmd_dyn) then
         if (self%dt <= 0.0_rk) then
            ierr = 1
            message = 'sqcalc dyn needs --dt DT, the time step of the trajectory'
            return
         end if
         if (self%maxframes < 1) then
            ierr = 1
            message = 'sqcalc dyn needs --maxframes L, the correlation window in frames'
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
         if ((self%modes_given .or. self%thin_given) .and. &
             self%q_mode /= dyn_q_grid) then
            ierr = 1
            message = '--modes and --thin need --qpoints grid:QMAX'
            return
         end if
         if (self%keep_modes .and. self%q_mode /= dyn_q_grid) then
            ierr = 1
            message = '--keep-modes needs --qpoints grid:QMAX'
            return
         end if
      else
         ! --- the static subcommand -----------------------------------------
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
            if (allocated(self%s2_accum_output) .and. .not. &
                allocated(self%pair_entropy_output)) then
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
      end if

      ! The shell averaged S(q) table is the positional argument of static and
      ! --sq of dyn; a dyn run without q points tabulates nothing.
      if (self%command == cmd_dyn) then
         if (self%q_mode == dyn_q_none) then
            if (allocated(self%sqw_output) .or. allocated(self%fqt_output) .or. &
                allocated(self%fqt_self_output) .or. allocated(self%s4_output)) then
               ierr = 1
               message = '--sqw, --fqt, --fqt-self and --s4 need a q sampling; add '// &
                  '--qpoints line:... or --qpoints shell:...'
               return
            end if
            if (allocated(self%sq_output)) then
               ierr = 1
               message = '--sq needs a q sampling; add --qpoints line:... or --qpoints shell:...'
               return
            end if
            if (.not. allocated(self%chi4_output) .and. .not. allocated(self%msd_output)) then
               ierr = 1
               message = 'a run without a q sampling computes only --chi4 and --msd; '// &
                  'add --chi4 or --msd FILE, or a q sampling'
               return
            end if
         end if
         if (allocated(self%sq_output)) self%output = self%sq_output
         if (npos > 0) then
            ierr = 1
            message = 'the dyn subcommand takes no positional argument; write the S(q) '// &
               'table with --sq FILE'
            return
         end if
      else
         if (npos > 1) then
            ierr = 1
            message = 'only one output argument is allowed'
            return
         end if
         if (npos == 0) then
            ! Only a run that writes some other table can do without the S(q)
            ! table: an XRD-only run.
            if (.not. allocated(self%xrd_output)) then
               ierr = 1
               message = 'missing output argument (use - for stdout)'
               return
            end if
         else
            self%output = trim(positional(1))
         end if
      end if
      if (.not. self%want_grid) self%grid_output = ''
      ! --format names the format of every table the run writes; without it
      ! each file infers the format from its suffix.
      if (self%output_format_given) then
         self%grid_format = self%output_format
         self%xrd_format = self%output_format
         self%sqw_format = self%output_format
         self%fqt_format = self%output_format
         self%fqt_self_format = self%output_format
         self%s4_format = self%output_format
         self%chi4_format = self%output_format
         self%msd_format = self%output_format
      else
         if (self%want_grid) then
            if (ends_with(self%grid_output, '.h5') .or. ends_with(self%grid_output, '.hdf5')) &
               self%grid_format = format_hdf5
         end if
         if (allocated(self%xrd_output)) then
            if (ends_with(self%xrd_output, '.h5') .or. ends_with(self%xrd_output, '.hdf5')) &
               self%xrd_format = format_hdf5
         end if
         if (allocated(self%sqw_output)) then
            if (ends_with(self%sqw_output, '.h5') .or. ends_with(self%sqw_output, '.hdf5')) &
               self%sqw_format = format_hdf5
         end if
         if (allocated(self%fqt_output)) then
            if (ends_with(self%fqt_output, '.h5') .or. ends_with(self%fqt_output, '.hdf5')) &
               self%fqt_format = format_hdf5
         end if
         if (allocated(self%fqt_self_output)) then
            if (ends_with(self%fqt_self_output, '.h5') .or. &
                ends_with(self%fqt_self_output, '.hdf5')) self%fqt_self_format = format_hdf5
         end if
         if (allocated(self%s4_output)) then
            if (ends_with(self%s4_output, '.h5') .or. ends_with(self%s4_output, '.hdf5')) &
               self%s4_format = format_hdf5
         end if
         if (allocated(self%chi4_output)) then
            if (ends_with(self%chi4_output, '.h5') .or. ends_with(self%chi4_output, '.hdf5')) &
               self%chi4_format = format_hdf5
         end if
         if (allocated(self%msd_output)) then
            if (ends_with(self%msd_output, '.h5') .or. ends_with(self%msd_output, '.hdf5')) &
               self%msd_format = format_hdf5
         end if
      end if
   end subroutine parse_options

   !> Parse the q sampling of the dyn subcommand (`--qpoints SPEC`).
   !!
   !! `line:NINT,S0,S1,DX,DY,DZ` is the density amplitude along a line in
   !! reciprocal space through Gamma, and `shell:Q,ACC` averages every
   !! direction of the shell |q| = Q on a Lebedev grid of the requested
   !! accuracy.  Leaving `--qpoints` out samples no q at all, which leaves only the
   !! scalar overlap Q(t) and chi4(t).
   subroutine parse_q_sampling(self, spec, ierr, message)
      type(options_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message

      ierr = 0
      message = ''
      if (starts_with(spec, 'line:')) then
         call parse_q_line(self, spec(6:), ierr, message)
         if (ierr /= 0) return
      else if (starts_with(spec, 'shell:')) then
         call parse_q_shell(self, spec(7:), ierr, message)
         if (ierr /= 0) return
      else if (starts_with(spec, 'grid:')) then
         call parse_q_grid(self, spec(6:), ierr, message)
         if (ierr /= 0) return
      else if (starts_with(spec, 'powder:')) then
         call parse_q_powder(self, spec(8:), ierr, message)
         if (ierr /= 0) return
      else if (starts_with(spec, 'single:')) then
         call parse_q_single(self, spec(8:), ierr, message)
         if (ierr /= 0) return
      else
         ierr = 1
         message = '--qpoints wants "line:NINT,S0,S1,DX,DY,DZ", '// &
            '"shell:Q,low|medium|high", "grid:QMAX", "powder:Q[,M|,dq=VALUE]" '// &
            'or "single:N1,N2,N3"'
         return
      end if
   end subroutine parse_q_sampling

   !> Parse "NINT,S0,S1,DX,DY,DZ" of `--qpoints line`.
   !!
   !! NINT is the number of intervals of the q line (so NINT+1 q points), S0/S1
   !! the range of its scale in 1/A, and DX,DY,DZ the (unnormalized) direction.
   subroutine parse_q_line(self, spec, ierr, message)
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
         message = '--qpoints line wants NINT,S0,S1,DX,DY,DZ, e.g. '// &
            '--qpoints line:100,0.5,20,1,1,0'
         return
      end if
      do j = 1, 6
         read (fields(j), *, iostat=ierr) values(j)
         if (ierr /= 0) then
            ierr = 1
            message = 'cannot read "'//trim(fields(j))//'" as a number in --qpoints line'
            return
         end if
      end do
      self%q_intervals = nint(values(1))
      self%q_s0 = values(2)
      self%q_s1 = values(3)
      self%q_dir = values(4:6)
      if (self%q_intervals < 1) then
         ierr = 1
         message = '--qpoints line needs at least one interval on the q line'
         return
      end if
      if (self%q_s1 <= self%q_s0 .or. self%q_s0 < 0.0_rk) then
         ierr = 1
         message = '--qpoints line needs 0 <= S0 < S1 for the scale of the q line'
         return
      end if
      if (sum(self%q_dir**2) <= 0.0_rk) then
         ierr = 1
         message = '--qpoints line needs a non-zero direction, e.g. 1,1,0'
         return
      end if
      self%q_mode = dyn_q_line
   end subroutine parse_q_line

   !> Parse "Q,ACC" of `--qpoints shell`.
   !!
   !! Q is the radius of the shell in 1/A and ACC one of the accuracy names of
   !! the Lebedev rules: low, medium or high.
   subroutine parse_q_shell(self, spec, ierr, message)
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
         message = '--qpoints shell wants Q,low|medium|high, e.g. '// &
            '--qpoints shell:2.5,medium'
         return
      end if
      read (fields(1), *, iostat=ierr) radius
      if (ierr /= 0) then
         ierr = 1
         message = 'cannot read "'//trim(fields(1))//'" as the |q| of --qpoints shell'
         return
      end if
      if (radius <= 0.0_rk) then
         ierr = 1
         message = '--qpoints shell needs a positive |q| radius'
         return
      end if
      order = lebedev_order_from_name(trim(fields(2)))
      if (order < 0) then
         ierr = 1
         message = 'unknown --qpoints shell accuracy "'//trim(fields(2))// &
            '" (use low, medium or high)'
         return
      end if
      self%shell_q = radius
      self%shell_order = order
      self%q_mode = dyn_q_shell
   end subroutine parse_q_shell

   !> Parse "QMAX" of `--qpoints grid`.
   !!
   !! QMAX is the upper bound of |q| in 1/A.  Every reciprocal-lattice vector
   !! with 0 < |q| <= QMAX is sampled, together with the Gamma point.  Only
   !! lattice vectors carry a density amplitude that is independent of how the
   !! periodic images are chosen, so a grid is free of the box-form-factor
   !! contamination of the off-lattice Lebedev shell; its mode count follows
   !! from the cell and can be capped with `--modes`.
   subroutine parse_q_grid(self, spec, ierr, message)
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
         message = 'cannot read "'//trim(spec)//'" as the qmax of --qpoints grid'
         return
      end if
      if (qmax <= 0.0_rk) then
         ierr = 1
         message = '--qpoints grid needs a positive qmax'
         return
      end if
      self%grid_qmax = qmax
      self%q_mode = dyn_q_grid
   end subroutine parse_q_grid

   !> Parse "Q[,M][,dq=VALUE]" of `--qpoints powder`.
   !!
   !! Q is the radius of the shell in 1/A.  Every reciprocal-lattice vector
   !! with |q| inside Q +- dq is averaged into one row with its multiplicity,
   !! which is the sampling the coherent quantities need: only lattice vectors
   !! carry amplitudes that are independent of the periodic-image convention.
   !! M is the target number of lattice vectors in the window (50 by default)
   !! and `dq=VALUE` fixes the half width instead, in 1/A.
   subroutine parse_q_powder(self, spec, ierr, message)
      type(options_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=48) :: fields(3)
      real(rk) :: radius, dq
      integer :: n, j, nmodes
      logical :: modes_given

      ierr = 0
      message = ''
      fields = ' '
      modes_given = .false.
      call split_fields(spec, fields, n)
      if (n < 1 .or. n > 3) then
         ierr = 1
         message = '--qpoints powder wants Q, or Q,M, or Q,dq=VALUE, e.g. '// &
            '--qpoints powder:2.5,80'
         return
      end if
      read (fields(1), *, iostat=ierr) radius
      if (ierr /= 0) then
         ierr = 1
         message = 'cannot read "'//trim(fields(1))//'" as the |q| of --qpoints powder'
         return
      end if
      if (radius <= 0.0_rk) then
         ierr = 1
         message = '--qpoints powder needs a positive |q| radius'
         return
      end if
      do j = 2, n
         if (starts_with(trim(fields(j)), 'dq=')) then
            read (fields(j)(4:), *, iostat=ierr) dq
            if (ierr /= 0 .or. dq <= 0.0_rk) then
               ierr = 1
               message = '--qpoints powder needs a positive dq= window half width'
               return
            end if
            self%powder_dq = dq
            self%powder_dq_given = .true.
         else
            read (fields(j), *, iostat=ierr) nmodes
            if (ierr /= 0 .or. nmodes < 1) then
               ierr = 1
               message = '--qpoints powder wants a positive number of lattice '// &
                  'vectors, e.g. powder:2.5,80'
               return
            end if
            self%powder_modes = nmodes
            modes_given = .true.
         end if
      end do
      ! M and dq= describe the same knob from two sides, so a line that carries
      ! both would silently drop one of them.
      if (modes_given .and. self%powder_dq_given) then
         ierr = 1
         message = '--qpoints powder takes either a number of lattice vectors '// &
            'or dq=, not both (M sets the width, dq= fixes it)'
         return
      end if
      self%powder_q = radius
      self%q_mode = dyn_q_powder
   end subroutine parse_q_powder

   !> Parse "N1,N2,N3" of `--qpoints single`.
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
   subroutine parse_q_single(self, spec, ierr, message)
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
         message = '--qpoints single wants N1,N2,N3, e.g. --qpoints single:1,0,0'
         return
      end if
      do i = 1, 3
         if (.not. is_integer_field(fields(i))) then
            ierr = 1
            message = 'cannot read "'//trim(fields(i))//'" as an integer of --qpoints single'
            return
         end if
         read (fields(i), *, iostat=ierr) values(i)
         if (ierr /= 0) then
            ierr = 1
            message = 'the index "'//trim(fields(i))//'" of --qpoints single is out of range'
            return
         end if
      end do
      self%single_index = values
      self%q_mode = dyn_q_single
   end subroutine parse_q_single

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

   !> Which subcommands accept an option.
   !!
   !! Unknown options are reported as accepted by both, so that the caller
   !! can still name them in its "unknown option" message.
   pure integer function flag_scope(name) result(scope)
      character(len=*), intent(in) :: name

      select case (trim(name))
      case ('--method', '--qmin', '--qmax', '--nq', '--eps', &
            '--device', '--gpu-id', '--precision', '--grid', &
            '--xrd', '--xrd-lambda', '--xrd-range', '--xrd-step', &
            '--lp', '--no-lp', '-fz', '--faber-ziman', &
            '--rmax', '--dr', '--skin', '--rdf', '--pair-entropy', &
            '--s2-accum', '--no-cutoff-correction')
         scope = scope_static
      case ('-q', '--qpoints', '--dt', '--maxframes', '--lag', '--modes', '--thin', &
            '--keep-modes', '--sq', '--sqw', '--fqt', '--fqt-self', &
            '--s4', '--chi4', '--msd', '--s4-cutoff', '--buffer-limit', '--stride')
         scope = scope_dyn
      case default
         scope = scope_both
      end select
   end function flag_scope

   !> Case sensitive prefix test for the `--qpoints` keywords.
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

   !> Print the option summary.
   !!
   !! Without COMMAND the subcommand list and the global options are printed;
   !! with cmd_static or cmd_dyn the global options come first and the options
   !! of that subcommand follow, so nothing is listed that the subcommand
   !! would refuse.
   subroutine print_usage(unit, command)
      integer, intent(in) :: unit
      integer, intent(in), optional :: command
      integer :: cmd

      cmd = cmd_none
      if (present(command)) cmd = command

      write (unit, '(a)') program_version
      write (unit, '(a)') ''
      write (unit, '(a)') 'Total structure factor S(q) from LAMMPS dump files.'
      write (unit, '(a)') ''
      if (cmd == cmd_none) then
         write (unit, '(a)') 'usage: sqcalc static [global options] [options] OUTPUT'
         write (unit, '(a)') '       sqcalc dyn    [global options] [options]'
         write (unit, '(a)') ''
         write (unit, '(a)') 'subcommands:'
         write (unit, '(a)') '  static           the quantities averaged over the frames: the'
         write (unit, '(a)') '                   shell averaged S(q) table, g(r) and the'
         write (unit, '(a)') '                   powder XRD pattern'
         write (unit, '(a)') '  dyn              the time axis is kept: S(q,w), F(q,t),'
         write (unit, '(a)') '                   S4(q,t), the overlap Q(t)/chi4(t) and the MSD'
         write (unit, '(a)') ''
         write (unit, '(a)') 'Run "sqcalc static -h" or "sqcalc dyn -h" for the options of'
         write (unit, '(a)') 'one subcommand.'
         write (unit, '(a)') ''
      else if (cmd == cmd_static) then
         write (unit, '(a)') 'usage: sqcalc static [global options] [options] OUTPUT'
         write (unit, '(a)') ''
         write (unit, '(a)') 'OUTPUT is the shell averaged S(q) table, use - for stdout.'
         write (unit, '(a)') ''
      else
         write (unit, '(a)') 'usage: sqcalc dyn [global options] [options]'
         write (unit, '(a)') ''
         write (unit, '(a)') 'dyn takes no positional argument; every table is named by'
         write (unit, '(a)') 'an option.'
         write (unit, '(a)') ''
      end if

      write (unit, '(a)') 'global options:'
      write (unit, '(a)') '  -i, --input FILE    LAMMPS dump trajectory (required)'
      write (unit, '(a)') '  -m, --mapping LIST  LAMMPS type id to element or species, e.g.'
      write (unit, '(a)') '                      1:Si,2:O or 1:Si4+,2:O2- (ions, see the skill docs)'
      write (unit, '(a)') '  -w, --weight SCHEME unit (default), neutron or xray'
      write (unit, '(a)') '  -t, --threads N     OpenMP threads (default: all available)'
      write (unit, '(a)') '      --norm NAME     mean (default), self or n'
      write (unit, '(a)') '      --partials      write the partial structure factor columns (default)'
      write (unit, '(a)') '      --no-partials   do not write partial structure factor columns'
      write (unit, '(a)') '      --format NAME   text (default) or hdf5 for the files that offer'
      write (unit, '(a)') '                      both: the grid, the XRD pattern, g(r), the pair'
      write (unit, '(a)') '                      entropy and the dynamic tables (the S(q) table'
      write (unit, '(a)') '                      itself is always text); without it a .h5/.hdf5'
      write (unit, '(a)') '                      name means hdf5'
      write (unit, '(a)') '      --quiet         do not write progress information to stderr'
      write (unit, '(a)') '  -h, --help          show this help'
      write (unit, '(a)') '  -v, --version       show the program version'

      if (cmd == cmd_none) return

      if (cmd == cmd_static) then
         write (unit, '(a)') ''
         write (unit, '(a)') 'static options:'
         write (unit, '(a)') '      --method NAME   nufft (default) or direct'
         write (unit, '(a)') '                      debye: real space pair histograms'
         write (unit, '(a)') '      --qmin VALUE    smallest |q| in the output [1/A] (default 0)'
         write (unit, '(a)') '      --qmax VALUE    largest |q| in the output [1/A] (default 20)'
         write (unit, '(a)') '      --nq N          number of q shells (default 500)'
         write (unit, '(a)') '      --eps VALUE     NUFFT tolerance (default 1e-9)'
         write (unit, '(a)') '      --device NAME   cpu (default) or gpu (cufinufft + cuFFT)'
         write (unit, '(a)') '      --gpu-id N      CUDA device to use (default 0)'
         write (unit, '(a)') '      --precision NAME  double (default) or single (float32 GPU)'
         write (unit, '(a)') '      --grid FILE     also write S(q) on every reciprocal lattice point'
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
         write (unit, '(a)') '      --rmax VALUE    Debye pair cutoff [A] (default: half the'
         write (unit, '(a)') '                      smallest periodic box side)'
         write (unit, '(a)') '      --dr VALUE      Debye radial bin width [A] (default 0.01)'
         write (unit, '(a)') '      --skin VALUE    Verlet skin for the pair list [A] (default 1.0)'
         write (unit, '(a)') '      --rdf FILE      total and partial g(r) in one file (.h5 = HDF5)'
         write (unit, '(a)') '      --pair-entropy FILE  total and partial pair entropy S2/kB'
         write (unit, '(a)') '      --s2-accum FILE  S2(r) accumulation curve for tail extrapolation'
         write (unit, '(a)') '      --no-cutoff-correction  disable the Debye cut-off density correction'
         write (unit, '(a)') '  -fz, --faber-ziman  partials in the Faber-Ziman normalization'
      else
         write (unit, '(a)') ''
         write (unit, '(a)') 'dyn options:'
         write (unit, '(a)') '  -q, --qpoints SPEC  q sampling, without it no q is sampled:'
         write (unit, '(a)') '                      line:NINT,S0,S1,DX,DY,DZ  q line through Gamma,'
         write (unit, '(a)') '                      NINT intervals, scale S0..S1 [1/A], direction'
         write (unit, '(a)') '                      shell:Q,low|medium|high  Lebedev average on |q| = Q'
         write (unit, '(a,i0,a,i0,a,i0,a)') '                      (low, medium and high are the ', &
            lebedev_points(lebedev_low), ', ', lebedev_points(lebedev_medium), ' and ', &
            lebedev_points(lebedev_high), ' point rules)'
         write (unit, '(a)') '                      grid:QMAX  every reciprocal lattice vector |q| <= QMAX'
         write (unit, '(a)') '                      powder:Q[,M|,dq=VALUE]  the lattice vectors with'
         write (unit, '(a)') '                      |q| = Q +- dq, averaged into one row (M, default'
         write (unit, '(a)') '                      50, sets dq; dq=VALUE fixes the half width)'
         write (unit, '(a)') '                      single:N1,N2,N3  one lattice vector of the box'
         write (unit, '(a)') '      --sq FILE       the shell averaged S(q) table of a --qpoints run'
         write (unit, '(a)') '      --dt VALUE      MD time step of the trajectory (required)'
         write (unit, '(a)') '      --maxframes N   correlation window in frames (required)'
         write (unit, '(a)') '      --lag N         frames between consecutive time origins (default 1)'
         write (unit, '(a)') '      --modes N       mode budget of --qpoints grid:QMAX, 0 = unlimited (default)'
         write (unit, '(a)') '      --thin NAME     order of the grid thinning: shells (default) or orbits'
         write (unit, '(a)') '      --keep-modes    write the grid rows per lattice vector instead of'
         write (unit, '(a)') '                      the |q| shell average (for diagnostics)'
         write (unit, '(a)') '      --sqw FILE      S(q,w) spectra, one row per (q,w) (.h5 = HDF5)'
         write (unit, '(a)') '      --fqt FILE      coherent F(q,t) intermediate scattering function'
         write (unit, '(a)') '      --fqt-self FILE self F_s(q,t) intermediate scattering function'
         write (unit, '(a)') '      --s4 FILE       S4(q,t) four-point structure factor (.h5 = HDF5)'
         write (unit, '(a)') '      --chi4 FILE     Q(t) and chi4(t) average overlap / susceptibility'
         write (unit, '(a)') '      --msd FILE      MSD(t) mean squared displacement (.h5 = HDF5)'
         write (unit, '(a)') '      --s4-cutoff A   overlap cutoff for --s4 and --chi4 [dump length unit]'
         write (unit, '(a)') '      --buffer-limit GB  position buffer limit for S4/chi4/F_s/MSD (default 2.0)'
         write (unit, '(a)') '      --stride N      use every N-th frame for S4/chi4/F_s/MSD (default 1)'
      end if
   end subroutine print_usage

end module sqc_options
