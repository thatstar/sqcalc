!> Command line handling for sqcalc.
module sqc_options
   use sqc_kinds
   use sqc_weights, only: weight_scheme_t, weight_unit, weight_neutron, weight_xray, &
                          scheme_from_name
   use sqc_structure_factor, only: method_nufft, method_direct, method_debye, norm_mean, &
                            norm_self, norm_natom
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
      integer :: norm = norm_mean
      integer :: device = device_cpu
      integer :: gpu_id = 0
      integer :: precision = precision_double
      !> Debye method settings (--rmax, --dr, --skin, --rdf).
      real(rk) :: rmax = 0.0_rk
      real(rk) :: dr = debye_default_dr
      real(rk) :: skin = debye_default_skin
      character(len=:), allocatable :: rdf_output
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
                  '--rmax', '--dr', '--skin', '--rdf')
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
      else
         if (self%rmax_given .or. self%dr_given .or. self%skin_given .or. &
             allocated(self%rdf_output)) then
            ierr = 1
            message = '--rmax, --dr, --skin and --rdf belong to --method debye'
            return
         end if
      end if
      ! Partial columns need to know which type is which element.
      if (self%partials .and. .not. self%scheme%has_mapping()) self%partials = .false.
      if (.not. self%want_grid) self%grid_output = ''
      ! Default the grid format from the file name.
      if (self%want_grid .and. .not. self%grid_format_given) then
         if (ends_with(self%grid_output, '.h5') .or. ends_with(self%grid_output, '.hdf5')) then
            self%grid_format = grid_format_hdf5
         end if
      end if
   end subroutine parse_options

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
      write (unit, '(a)') '      --rmax VALUE    Debye pair cutoff [1/A] (default: half the'
      write (unit, '(a)') '                      smallest periodic box side)'
      write (unit, '(a)') '      --dr VALUE      Debye radial bin width [1/A] (default 0.01)'
      write (unit, '(a)') '      --skin VALUE    Verlet skin for the pair list [1/A] (default 1.0)'
      write (unit, '(a)') '      --rdf FILE      total and partial g(r) in one file (.h5 = HDF5)'
      write (unit, '(a)') '      --no-cutoff-correction  disable the Debye cut-off density correction'
      write (unit, '(a)') '  -fz, --faber-ziman  partials in the Faber-Ziman normalization'
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
