! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Per-atom scattering weights: unit, neutron or X-ray.
module sqc_weights
   use sqc_kinds
   use sqc_elements, only: element_table
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
      !> Canonical species label of each type id (blank when unmapped).
      character(len=8), allocatable :: species(:)
      !> X-ray species index of each type id (0 for an element with no row).
      integer, allocatable :: species_index(:)
      !> Element table index of each type id (0 when unmapped).
      integer, allocatable :: element(:)
      !> Largest type id that was actually mapped (0 when none).
      integer(ik) :: n_mapped = 0
   contains
      procedure :: grow => scheme_grow
      procedure :: set_mapping => scheme_set_mapping
      procedure :: has_mapping => scheme_has_mapping
      procedure :: mapped_types => scheme_mapped_types
      procedure :: amplitude => scheme_amplitude
      procedure :: amplitude_of_species => scheme_amplitude_of_species
      procedure :: label => scheme_label
      procedure :: pair_label => scheme_pair_label
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
      character(len=8) :: label
      integer :: nitems, i, colon, type_id, idx, element

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
         call element_table%find_species(right, idx, element, label, ierr, message)
         if (ierr /= 0) return
         call self%grow(type_id)
         self%species(type_id) = label
         self%species_index(type_id) = idx
         self%element(type_id) = element
         self%n_mapped = max(self%n_mapped, int(type_id, ik))
      end do
   end subroutine scheme_set_mapping

   !> Resize the per-type arrays, preserving what is already stored.
   subroutine scheme_grow(self, ntypes)
      class(weight_scheme_t), intent(inout) :: self
      integer, intent(in) :: ntypes
      character(len=8), allocatable :: species(:)
      integer, allocatable :: species_index(:), element(:)
      integer :: n

      if (ntypes <= self%ntypes) return
      n = max(ntypes, max(self%ntypes*2, 4))
      allocate(species(n), species_index(n), element(n))
      species = ' '
      species_index = 0
      element = 0
      if (self%ntypes > 0) then
         species(1:self%ntypes) = self%species(1:self%ntypes)
         species_index(1:self%ntypes) = self%species_index(1:self%ntypes)
         element(1:self%ntypes) = self%element(1:self%ntypes)
      end if
      call move_alloc(species, self%species)
      call move_alloc(species_index, self%species_index)
      call move_alloc(element, self%element)
      self%ntypes = n
   end subroutine scheme_grow

   pure logical function scheme_has_mapping(self) result(mapped)
      class(weight_scheme_t), intent(in) :: self
      mapped = self%ntypes > 0
   end function scheme_has_mapping

   !> Number of mapped LAMMPS type ids (0 when no mapping was given).
   !!
   !! The per-type arrays are grown in blocks, so their size is the capacity and
   !! not the number of mapped types; consumers must ask here instead of looking
   !! at size(symbols) directly.
   pure integer function scheme_mapped_types(self) result(n)
      class(weight_scheme_t), intent(in) :: self
      n = int(self%n_mapped)
   end function scheme_mapped_types

   !> Scattering amplitude of one LAMMPS type id at momentum transfer q.
   pure real(rk) function scheme_amplitude(self, type_id, q) result(weight)
      class(weight_scheme_t), intent(in) :: self
      integer(ik), intent(in) :: type_id
      real(rk), intent(in) :: q
      integer :: idx

      select case (self%kind)
      case (weight_neutron)
         idx = 0
         if (type_id >= 1 .and. type_id <= self%ntypes) idx = self%element(type_id)
         weight = element_table%neutron_b(idx)
      case (weight_xray)
         idx = 0
         if (type_id >= 1 .and. type_id <= self%ntypes) idx = self%species_index(type_id)
         weight = element_table%species_xray_f(idx, q)
      case default
         weight = 1.0_rk
      end select
   end function scheme_amplitude

   !> Scattering amplitude of a species index at momentum transfer q.
   !!
   !! Neutron lengths are nuclear and do not depend on the electronic state, so
   !! the change of an ion label applies to the X-ray form factor alone.
   pure real(rk) function scheme_amplitude_of_species(self, idx, q) result(weight)
      class(weight_scheme_t), intent(in) :: self
      integer, intent(in) :: idx
      real(rk), intent(in) :: q

      select case (self%kind)
      case (weight_neutron)
         weight = element_table%neutron_b(element_table%species_element_of(idx))
      case (weight_xray)
         weight = element_table%species_xray_f(idx, q)
      case default
         weight = 1.0_rk
      end select
   end function scheme_amplitude_of_species

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

   !> Label of a type pair, "Si-O" when the types are mapped and "1-2" otherwise.
   pure function scheme_pair_label(self, ia, ib) result(label)
      class(weight_scheme_t), intent(in) :: self
      integer, intent(in) :: ia, ib
      character(len=18) :: label
      character(len=8) :: sa, sb

      sa = ' '
      sb = ' '
      if (allocated(self%species)) then
         if (ia >= 1 .and. ia <= size(self%species)) sa = self%species(ia)
         if (ib >= 1 .and. ib <= size(self%species)) sb = self%species(ib)
      end if
      if (len_trim(sa) > 0 .and. len_trim(sb) > 0) then
         label = trim(sa)//'-'//trim(sb)
      else
         write (label, '(i0,a,i0)') ia, '-', ib
      end if
   end function scheme_pair_label

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
         ! note: Fortran does not short-circuit .and., so the comparison must
         ! not be evaluated once pos has run past the end of the string
         do
            if (pos > n) exit
            if (text(pos:pos) == delimiter) exit
            pos = pos + 1
         end do
         nitems = nitems + 1
         items(nitems) = text(start:pos - 1)
         pos = pos + 1
      end do
      if (n == 0) nitems = 0
   end subroutine split_char

end module sqc_weights
