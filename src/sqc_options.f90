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
      real(rk) :: eps = 1.0e-9_rk
      logical :: want_grid = .false.
      integer :: grid_format = grid_format_text
      logical :: grid_format_given = .false.
      logical :: quiet = .false.
      logical :: show_help = .false.
      logical :: show_version = .false.
      !> Dynamic structure factor (--dyn and its parameters).
      logical :: dynamic = .false.
      !> q line: NINT intervals, scale from S0 to S1 along (DX, DY, DZ).
      integer :: dyn_intervals = 0
      real(rk) :: dyn_s0 = 0.0_rk
      real(rk) :: dyn_s1 = 0.0_rk
      real(rk) :: dyn_dir(3) = 0.0_rk
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
      real(rk) :: s4_cutoff = 0.0_rk
      logical :: s4_cutoff_given = .false.
      real(rk) :: buffer_limit_gb = 2.0_rk
      logical :: buffer_limit_given = .false.
      integer :: stride = 1
      logical :: stride_given = .false.
      integer :: s4_format = grid_format_text
      integer :: chi4_format = grid_format_text
      type(weight_scheme_t) :: scheme
   end type options_t

contains

   !> Parse the argument vector into an options_t.
   subroutine parse_options(self, ierr, message)
      type(options_t), intent(inout) :: self
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=256) :: arg, name, value, positional(4)
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
            case ('--no-cutoff-correction')
               self%no_cutoff_correction = .true.
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
            case ('-i', '--input', '--mapping', '-m', '-w', '--weight', '-t', '--threads', &
                  '--qmin', '--qmax', '--nq', '--eps', '--method', '--norm', '--grid', &
                  '--device', '--gpu-id', '--precision', '--grid-format', &
                  '--rmax', '--dr', '--skin', '--rdf', &
                  '--pair-entropy', '--s2-accum', &
                  '--dyn', '--dt', '--maxframes', '--lag', '--sqw', '--fqt', '--dyn-format', &
                  '--fqt-self', '--s4', '--chi4', '--s4-cutoff', '--buffer-limit', '--stride')
               if (.not. has_inline) then
                  if (i + 1 > nargs) then
                     ierr = 1
                     message = 'missing value for option '//trim(name)
                     return
                  end if
                  i = i + 1
                  call get_command_argument(i, value)
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
               case ('--qmax')
                  read (value, *, iostat=ierr) self%qmax
                  if (ierr /= 0 .or. self%qmax <= 0.0_rk) then
                     ierr = 1
                     message = '--qmax must be a positive number'
                     return
                  end if
               case ('--nq')
                  read (value, *, iostat=ierr) self%nq
                  if (ierr /= 0 .or. self%nq < 1) then
                     ierr = 1
                     message = '--nq must be a positive integer'
                     return
                  end if
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
               case ('--dyn')
                  call parse_dyn(self, trim(value), ierr, message)
                  if (ierr /= 0) return
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
      if (npos < 1) then
         ierr = 1
         message = 'missing output argument (use - for stdout)'
         return
      end if
      if (npos > 1) then
         ierr = 1
         message = 'only one output argument is allowed'
         return
      end if
      self%output = trim(positional(1))
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
            message = 'the dynamic method samples a q line; --grid is not available'
            return
         end if
         if (self%device == device_gpu) then
            ierr = 1
            message = 'the dynamic method runs on the CPU; use --device cpu'
            return
         end if
         if (allocated(self%s4_output) .or. allocated(self%chi4_output) .or. &
             allocated(self%fqt_self_output)) then
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
            message = '--s4-cutoff, --buffer-limit and --stride need --s4, --chi4 or --fqt-self'
            return
         end if
      else if (allocated(self%sqw_output) .or. allocated(self%fqt_output) .or. &
               allocated(self%fqt_self_output) .or. allocated(self%s4_output) .or. &
               allocated(self%chi4_output)) then
         ierr = 1
         message = '--sqw, --fqt, --fqt-self, --s4 and --chi4 belong to --dyn'
         return
      else if (self%dt > 0.0_rk .or. self%maxframes > 0 .or. self%lag_given .or. &
               self%dyn_format_given .or. self%s4_cutoff_given .or. &
               self%buffer_limit_given .or. self%stride_given) then
         ierr = 1
         message = '--dt, --maxframes, --lag, --dyn-format, --s4-cutoff, --buffer-limit '// &
            'and --stride belong to --dyn'
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
      if (.not. self%want_grid) self%grid_output = ''
      ! Default the grid format from the file name.
      if (self%want_grid .and. .not. self%grid_format_given) then
         if (ends_with(self%grid_output, '.h5') .or. ends_with(self%grid_output, '.hdf5')) then
            self%grid_format = grid_format_hdf5
         end if
      end if
      if (self%dyn_format_given) then
         self%sqw_format = self%dyn_format
         self%fqt_format = self%dyn_format
         self%fqt_self_format = self%dyn_format
         self%s4_format = self%dyn_format
         self%chi4_format = self%dyn_format
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
      end if
   end subroutine parse_options

   !> Parse the packed `--dyn` specification "NINT,S0,S1,DX,DY,DZ".
   !!
   !! NINT is the number of intervals of the q line (so NINT+1 q points), S0/S1
   !! the range of its scale in 1/A, and DX,DY,DZ the (unnormalized) direction.
   subroutine parse_dyn(self, spec, ierr, message)
      type(options_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=48) :: fields(6)
      real(rk) :: values(6)
      integer :: i, j, n, start, last
      logical :: at_end

      ierr = 0
      message = ''
      fields = ' '
      values = 0.0_rk
      n = 0
      start = 1
      last = len_trim(spec)
      do i = 1, last + 1
         ! Fortran does not short-circuit .or., so the end of the string has to
         ! be tested separately before spec(i:i) is evaluated.
         at_end = i > last
         if (.not. at_end) at_end = spec(i:i) == ','
         if (at_end) then
            if (i > start) then
               n = n + 1
               if (n > 6) exit
               fields(n) = spec(start:i - 1)
            end if
            start = i + 1
         end if
      end do
      if (n /= 6) then
         ierr = 1
         message = '--dyn wants NINT,S0,S1,DX,DY,DZ, e.g. --dyn 100,0.5,20,1,1,0'
         return
      end if
      do j = 1, 6
         read (fields(j), *, iostat=ierr) values(j)
         if (ierr /= 0) then
            ierr = 1
            message = 'cannot read "'//trim(fields(j))//'" as a number in --dyn'
            return
         end if
      end do
      self%dyn_intervals = nint(values(1))
      self%dyn_s0 = values(2)
      self%dyn_s1 = values(3)
      self%dyn_dir = values(4:6)
      if (self%dyn_intervals < 1) then
         ierr = 1
         message = '--dyn needs at least one interval on the q line'
         return
      end if
      if (self%dyn_s1 <= self%dyn_s0 .or. self%dyn_s0 < 0.0_rk) then
         ierr = 1
         message = '--dyn needs 0 <= S0 < S1 for the scale of the q line'
         return
      end if
      if (sum(self%dyn_dir**2) <= 0.0_rk) then
         ierr = 1
         message = '--dyn needs a non-zero direction, e.g. 1,1,0'
         return
      end if
      self%dynamic = .true.
   end subroutine parse_dyn

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
      write (unit, '(a)') '  -m, --mapping LIST  LAMMPS type id to element symbol, e.g. 1:Si,2:O'
      write (unit, '(a)') '  -w, --weight SCHEME unit (default), neutron or xray'
      write (unit, '(a)') '  -t, --threads N     OpenMP threads (default: all available)'
      write (unit, '(a)') '      --qmin VALUE    smallest |q| in the output [1/A] (default 0)'
      write (unit, '(a)') '      --qmax VALUE    largest |q| in the output [1/A] (default 20)'
      write (unit, '(a)') '      --nq N          number of q shells (default 500)'
      write (unit, '(a)') '      --grid FILE     also write S(q) on every reciprocal lattice point'
      write (unit, '(a)') '      --grid-format NAME  text (default) or hdf5 (.h5/.hdf5 implies hdf5)'
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
      write (unit, '(a)') '      --dyn SPEC      dynamic structure factor S(q,w) along a q line:'
      write (unit, '(a)') '                      NINT,S0,S1,DX,DY,DZ, e.g. 100,0.5,20,1,1,0'
      write (unit, '(a)') '                      NINT intervals, scale S0..S1 [1/A], direction'
      write (unit, '(a)') '      --dt VALUE      MD time step of the trajectory (with --dyn)'
      write (unit, '(a)') '      --maxframes N   correlation window in frames (with --dyn)'
      write (unit, '(a)') '      --lag N         frames between consecutive time origins (default 1)'
      write (unit, '(a)') '      --sqw FILE      S(q,w) spectra, one row per (q,w) (.h5 = HDF5)'
      write (unit, '(a)') '      --fqt FILE      coherent F(q,t) intermediate scattering function'
      write (unit, '(a)') '      --fqt-self FILE self F_s(q,t) intermediate scattering function'
      write (unit, '(a)') '      --dyn-format NAME  text (default) or hdf5 for all dynamic outputs'
      write (unit, '(a)') '      --s4 FILE       S4(q,t) four-point structure factor (.h5 = HDF5)'
      write (unit, '(a)') '      --chi4 FILE     Q(t) and chi4(t) average overlap / susceptibility'
      write (unit, '(a)') '      --s4-cutoff A   overlap cutoff for --s4 and --chi4 [dump length unit]'
      write (unit, '(a)') '      --buffer-limit GB  position buffer limit for S4/chi4/F_s (default 2.0)'
      write (unit, '(a)') '      --stride N      use every N-th frame for S4/chi4/F_s (default 1)'
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
