!> Per-atom scattering weights: unit, neutron or X-ray.
module sqc_weights
   use sqc_kinds
   use sqc_elements, only: element_table, normalise_symbol
   implicit none
   private

   public :: weight_scheme_t, weight_unit, weight_neutron, weight_xray, &
             scheme_from_name

   !> Kinds of weighting available on the command line.
   integer, parameter :: weight_unit = 0
   integer, parameter :: weight_neutron = 1
   integer, parameter :: weight_xray = 2

   !> Weighting scheme plus the type-id to element mapping it needs.
   type :: weight_scheme_t
      !> One of weight_unit, weight_neutron, weight_xray.
      integer :: kind = weight_unit
      !> Number of LAMMPS atom types covered by the mapping.
      integer(ik) :: ntypes = 0
      !> Canonical element symbol of each type id (blank when unmapped).
      character(len=2), allocatable :: symbols(:)
      !> Index of each type id in the element table (0 when unmapped).
      integer, allocatable :: element(:)
   contains
      procedure :: grow => scheme_grow
      procedure :: set_mapping => scheme_set_mapping
      procedure :: has_mapping => scheme_has_mapping
      procedure :: mapped_types => scheme_mapped_types
      procedure :: amplitude => scheme_amplitude
      procedure :: amplitude_of_element => scheme_amplitude_of_element
      procedure :: label => scheme_label
   end type weight_scheme_t

contains

   !> Parse a mapping specification such as "1:Si,2:O".
   subroutine scheme_set_mapping(self, spec, ierr, message)
      class(weight_scheme_t), intent(inout) :: self
      character(len=*), intent(in) :: spec
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=64), allocatable :: items(:)
      character(len=:), allocatable :: item, left, right
      integer :: nitems, i, colon, type_id, idx

      ierr = 0
      message = ''
      call split_char(spec, ',', items, nitems)
      if (nitems == 0) then
         ierr = 1
         message = 'empty element mapping'
         return
      end if

      do i = 1, nitems
         item = trim(adjustl(items(i)))
         if (len_trim(item) == 0) cycle
         colon = index(item, ':')
         if (colon <= 1 .or. colon == len(item)) then
            ierr = 1
            message = 'cannot parse mapping entry "'//trim(item)//'" (expected TYPE:SYMBOL)'
            return
         end if
         left = trim(adjustl(item(:colon - 1)))
         right = trim(adjustl(item(colon + 1:)))
         read (left, *, iostat=ierr) type_id
         if (ierr /= 0 .or. type_id < 1) then
            ierr = 1
            message = 'invalid atom type id "'//left//'" in the element mapping'
            return
         end if
         idx = element_table%index_of(right)
         if (idx == 0) then
            ierr = 1
            message = 'unknown element symbol "'//right//'" in the element mapping'
            return
         end if
         call self%grow(type_id)
         self%symbols(type_id) = element_table%symbol_of(idx)
         self%element(type_id) = idx
      end do
   end subroutine scheme_set_mapping

   !> Resize the per-type arrays, preserving what is already stored.
   subroutine scheme_grow(self, ntypes)
      class(weight_scheme_t), intent(inout) :: self
      integer, intent(in) :: ntypes
      character(len=2), allocatable :: symbols(:)
      integer, allocatable :: element(:)
      integer :: n

      if (ntypes <= self%ntypes) return
      n = max(ntypes, max(self%ntypes*2, 4))
      allocate(symbols(n), element(n))
      symbols = ' '
      element = 0
      if (self%ntypes > 0) then
         symbols(1:self%ntypes) = self%symbols(1:self%ntypes)
         element(1:self%ntypes) = self%element(1:self%ntypes)
      end if
      call move_alloc(symbols, self%symbols)
      call move_alloc(element, self%element)
      self%ntypes = n
   end subroutine scheme_grow

   pure logical function scheme_has_mapping(self) result(mapped)
      class(weight_scheme_t), intent(in) :: self
      mapped = self%ntypes > 0
   end function scheme_has_mapping

   !> Number of mapped LAMMPS type ids (0 when no mapping was given).
   !!
   !! The per-type arrays are only allocated once a mapping is supplied, so all
   !! consumers must ask here instead of looking at size(symbols) directly.
   pure integer function scheme_mapped_types(self) result(n)
      class(weight_scheme_t), intent(in) :: self
      if (allocated(self%symbols)) then
         n = size(self%symbols)
      else
         n = 0
      end if
   end function scheme_mapped_types

   !> Scattering amplitude of one LAMMPS type id at momentum transfer q.
   pure real(rk) function scheme_amplitude(self, type_id, q) result(weight)
      class(weight_scheme_t), intent(in) :: self
      integer(ik), intent(in) :: type_id
      real(rk), intent(in) :: q
      integer :: idx

      select case (self%kind)
      case (weight_neutron, weight_xray)
         idx = 0
         if (type_id >= 1 .and. type_id <= self%ntypes) idx = self%element(type_id)
         weight = self%amplitude_of_element(idx, q)
      case default
         weight = 1.0_rk
      end select
   end function scheme_amplitude

   !> Scattering amplitude of an element table index at momentum transfer q.
   pure real(rk) function scheme_amplitude_of_element(self, idx, q) result(weight)
      class(weight_scheme_t), intent(in) :: self
      integer, intent(in) :: idx
      real(rk), intent(in) :: q

      select case (self%kind)
      case (weight_neutron)
         weight = element_table%neutron_b(idx)
      case (weight_xray)
         weight = element_table%xray_f(idx, q)
      case default
         weight = 1.0_rk
      end select
   end function scheme_amplitude_of_element

   !> Human readable description used in diagnostics.
   pure function scheme_label(self) result(label)
      class(weight_scheme_t), intent(in) :: self
      character(len=:), allocatable :: label
      select case (self%kind)
      case (weight_neutron)
         label = 'neutron bound coherent scattering lengths'
      case (weight_xray)
         label = 'X-ray atomic form factors f(q)'
      case default
         label = 'unit weights'
      end select
   end function scheme_label

   !> Map a command line scheme name to its integer code.
   integer function scheme_from_name(name) result(kind)
      character(len=*), intent(in) :: name
      character(len=16) :: text
      integer :: i, code

      text = name
      do i = 1, len_trim(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) text(i:i) = achar(code + 32)
      end do
      select case (trim(text))
      case ('unit', 'none', '1')
         kind = weight_unit
      case ('neutron', 'n')
         kind = weight_neutron
      case ('xray', 'x-ray', 'x_ray', 'x')
         kind = weight_xray
      case default
         kind = -1
      end select
   end function scheme_from_name

   !> Split a string on a single character delimiter.
   subroutine split_char(text, delimiter, items, nitems)
      character(len=*), intent(in) :: text, delimiter
      character(len=64), allocatable, intent(out) :: items(:)
      integer, intent(out) :: nitems
      integer :: pos, start, n

      n = len_trim(text)
      allocate (items(max(n, 1)))
      nitems = 0
      pos = 1
      do while (pos <= n)
         start = pos
         do while (pos <= n .and. text(pos:pos) /= delimiter)
            pos = pos + 1
         end do
         nitems = nitems + 1
         items(nitems) = text(start:pos - 1)
         pos = pos + 1
      end do
      if (n == 0) nitems = 0
   end subroutine split_char

end module sqc_weights
