! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Access to the periodic table, the X-ray species and the neutron lengths.
!!
!! The element table has one row per element and carries the mass and the
!! bound coherent neutron scattering length.  X-ray form factors live in a
!! second table indexed by *species*: the neutral atom of an element is the
!! species with charge 0 and no valence flag, while labels such as "O2-",
!! "Si4+" or "Cval" select an ion or a valence-state parameterization of the
!! same element.  Both tables come from the same IT92 source; see
!! `src/sqc_element_data.f90`.
module sqc_elements
   use sqc_kinds
   use sqc_element_data, only: n_elements, n_species, n_species_alias, element_z, &
                               element_symbol, element_name, element_mass, &
                               element_neutron_b, element_has_neutron, &
                               element_neutral_species, species_label, species_element, &
                               species_charge, species_valence, species_xray_a, &
                               species_xray_b, species_xray_c, species_alias_from, &
                               species_alias_to
   implicit none
   private

   public :: element_table_t, element_table, normalise_symbol

   !> Stateless lookup object for the bundled element and species data.
   type :: element_table_t
   contains
      procedure, nopass :: size => table_size
      procedure, nopass :: index_of => table_index_of
      procedure, nopass :: symbol_of => table_symbol_of
      procedure, nopass :: name_of => table_name_of
      procedure, nopass :: mass_of => table_mass_of
      procedure, nopass :: z_of => table_z_of
      procedure, nopass :: has_neutron => table_has_neutron
      procedure, nopass :: neutron_b => table_neutron_b
      procedure, nopass :: species_count => table_species_count
      procedure, nopass :: neutral_species => table_neutral_species
      procedure, nopass :: find_species => table_find_species
      procedure, nopass :: species_label_of => table_species_label_of
      procedure, nopass :: species_element_of => table_species_element_of
      procedure, nopass :: species_charge_of => table_species_charge_of
      procedure, nopass :: species_is_valence => table_species_is_valence
      procedure, nopass :: species_xray_f => table_species_xray_f
   end type element_table_t

   !> Shared instance used by the rest of the program.
   type(element_table_t), parameter :: element_table = element_table_t()

contains

   pure integer function table_size() result(n)
      n = n_elements
   end function table_size

   pure integer function table_species_count() result(n)
      n = n_species
   end function table_species_count

   !> Canonical spelling of an element label: "SI" and "si" both become "Si".
   pure function normalise_symbol(text) result(symbol)
      character(len=*), intent(in) :: text
      character(len=2) :: symbol
      integer :: i, n

      symbol = ' '
      n = min(len_trim(text), 2)
      if (n < 1) return
      do i = 1, n
         if (i == 1) then
            symbol(i:i) = upper_char(text(i:i))
         else
            symbol(i:i) = lower_char(text(i:i))
         end if
      end do
   end function normalise_symbol

   pure function upper_char(ch) result(out)
      character(len=1), intent(in) :: ch
      character(len=1) :: out
      integer :: code
      code = iachar(ch)
      if (code >= iachar('a') .and. code <= iachar('z')) then
         out = achar(code - 32)
      else
         out = ch
      end if
   end function upper_char

   pure function lower_char(ch) result(out)
      character(len=1), intent(in) :: ch
      character(len=1) :: out
      integer :: code
      code = iachar(ch)
      if (code >= iachar('A') .and. code <= iachar('Z')) then
         out = achar(code + 32)
      else
         out = ch
      end if
   end function lower_char

   !> Case insensitive comparison of two labels, ignoring trailing blanks.
   pure logical function same_text(first, second) result(equal)
      character(len=*), intent(in) :: first, second
      integer :: i, n

      equal = .false.
      n = len_trim(first)
      if (len_trim(second) /= n) return
      do i = 1, n
         if (lower_char(first(i:i)) /= lower_char(second(i:i))) return
      end do
      equal = .true.
   end function same_text

   !> Lower case copy of a label, in a buffer large enough for the suffixes.
   pure function lower_text(text) result(out)
      character(len=*), intent(in) :: text
      character(len=16) :: out
      integer :: i, n

      out = ' '
      n = min(len_trim(text), len(out))
      do i = 1, n
         out(i:i) = lower_char(text(i:i))
      end do
   end function lower_text

   !> Index into the element table, or 0 when the symbol is unknown.
   pure integer function table_index_of(symbol) result(idx)
      character(len=*), intent(in) :: symbol
      character(len=2) :: wanted
      integer :: i

      wanted = normalise_symbol(symbol)
      idx = 0
      if (len_trim(wanted) == 0) return
      do i = 1, n_elements
         if (element_symbol(i) == wanted) then
            idx = i
            return
         end if
      end do
   end function table_index_of

   pure function table_symbol_of(idx) result(symbol)
      integer, intent(in) :: idx
      character(len=2) :: symbol
      symbol = '  '
      if (idx >= 1 .and. idx <= n_elements) symbol = element_symbol(idx)
   end function table_symbol_of

   pure function table_name_of(idx) result(name)
      integer, intent(in) :: idx
      character(len=16) :: name
      name = ''
      if (idx >= 1 .and. idx <= n_elements) name = element_name(idx)
   end function table_name_of

   pure real(rk) function table_mass_of(idx) result(mass)
      integer, intent(in) :: idx
      mass = 0.0_rk
      if (idx >= 1 .and. idx <= n_elements) mass = element_mass(idx)
   end function table_mass_of

   pure integer function table_z_of(idx) result(z)
      integer, intent(in) :: idx
      z = 0
      if (idx >= 1 .and. idx <= n_elements) z = element_z(idx)
   end function table_z_of

   pure logical function table_has_neutron(idx) result(found)
      integer, intent(in) :: idx
      found = idx >= 1 .and. idx <= n_elements .and. element_has_neutron(idx)
   end function table_has_neutron

   !> Bound coherent neutron scattering length [fm].
   pure real(rk) function table_neutron_b(idx) result(b)
      integer, intent(in) :: idx
      b = 0.0_rk
      if (idx >= 1 .and. idx <= n_elements) b = element_neutron_b(idx)
   end function table_neutron_b

   !> Neutral species of an element, or 0 when it has no X-ray row.
   pure integer function table_neutral_species(element) result(idx)
      integer, intent(in) :: element
      idx = 0
      if (element >= 1 .and. element <= n_elements) idx = element_neutral_species(element)
   end function table_neutral_species

   !> Canonical label of a species index.
   pure function table_species_label_of(idx) result(label)
      integer, intent(in) :: idx
      character(len=8) :: label
      label = ''
      if (idx >= 1 .and. idx <= n_species) label = species_label(idx)
   end function table_species_label_of

   !> Element table index of a species index.
   pure integer function table_species_element_of(idx) result(element)
      integer, intent(in) :: idx
      element = 0
      if (idx >= 1 .and. idx <= n_species) element = species_element(idx)
   end function table_species_element_of

   !> Formal charge of a species index (0 for neutral and valence states).
   pure integer function table_species_charge_of(idx) result(charge)
      integer, intent(in) :: idx
      charge = 0
      if (idx >= 1 .and. idx <= n_species) charge = species_charge(idx)
   end function table_species_charge_of

   !> True for a valence-state parameterization such as "Cval".
   pure logical function table_species_is_valence(idx) result(valence)
      integer, intent(in) :: idx
      valence = .false.
      if (idx >= 1 .and. idx <= n_species) valence = species_valence(idx)
   end function table_species_is_valence

   !> X-ray atomic form factor in electron units at momentum transfer q [1/A].
   !!
   !! f(q) = sum_i a_i exp(-b_i (q/4pi)^2) + c   (IT92 / Cromer-Mann)
   pure real(rk) function table_species_xray_f(idx, q) result(f)
      integer, intent(in) :: idx
      real(rk), intent(in) :: q
      real(rk) :: stol2
      integer :: i

      f = 0.0_rk
      if (idx < 1 .or. idx > n_species) return
      stol2 = (q/(4.0_rk*acos(-1.0_rk)))**2
      f = species_xray_c(idx)
      do i = 1, 4
         f = f + species_xray_a(i, idx)*exp(-species_xray_b(i, idx)*stol2)
      end do
   end function table_species_xray_f

   !> Resolve a species label such as "Si", "O2-" or "Cval".
   !!
   !! On return `ierr` is 0 for a resolved label, 1 for an unknown element and
   !! 2 for a known element that has no such charge state; `message` then
   !! names the element and lists the states that are tabulated.  A bare
   !! element symbol of an element without an X-ray row resolves with
   !! `sidx = 0`, which the caller reports when X-ray weights are requested.
   subroutine table_find_species(text, sidx, element, label, ierr, message)
      character(len=*), intent(in) :: text
      integer, intent(out) :: sidx
      integer, intent(out) :: element
      character(len=*), intent(out) :: label
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      character(len=len(text)) :: wanted
      character(len=2) :: symbol
      character(len=8) :: canonical
      character(len=256) :: states
      character(len=len(text)) :: tail
      integer :: i, charge, idx
      logical :: valence, ok

      sidx = 0
      element = 0
      label = ''
      ierr = 0
      message = ''
      wanted = trim(adjustl(text))
      do i = 1, n_species_alias
         if (same_text(wanted, species_alias_from(i))) wanted = species_alias_to(i)
      end do

      call split_species_name(wanted, symbol, charge, valence, tail, ok)
      if (.not. ok) then
         ierr = 1
         message = 'unknown element symbol "'//trim(adjustl(text))//'" in the element mapping'
         return
      end if
      element = table_index_of(symbol)
      if (charge == 0 .and. .not. valence) then
         sidx = element_neutral_species(element)
         label = symbol
         return
      end if

      canonical = species_name(symbol, charge, valence)
      do idx = 1, n_species
         if (same_text(canonical, species_label(idx))) then
            sidx = idx
            label = species_label(idx)
            return
         end if
      end do

      ierr = 2
      call available_states(element, states)
      message = 'element "'//trim(symbol)//'" has no charge state "'//trim(adjustl(text))// &
                '"; available: '//trim(states)
      label = symbol
   end subroutine table_find_species

   !> Split "O2-", "Si4+" or "Cval" into element, charge and valence flag.
   !!
   !! A neutral label leaves both the charge and the valence flag false.  The
   !! element is the longest of the one and two letter symbols that leaves a
   !! suffix the grammar accepts, so "Cval" is carbon, not a "Cv" element.
   pure subroutine split_species_name(text, symbol, charge, valence, tail, ok)
      character(len=*), intent(in) :: text
      character(len=2), intent(out) :: symbol
      integer, intent(out) :: charge
      logical, intent(out) :: valence
      character(len=*), intent(out) :: tail
      logical, intent(out) :: ok
      character(len=2) :: two
      character(len=16) :: rest

      symbol = '  '
      charge = 0
      valence = .false.
      tail = ''
      ok = .false.
      if (len_trim(text) < 1) return
      if (len_trim(text) >= 2) then
         two = normalise_symbol(text(1:2))
         if (table_index_of(two) > 0) then
            rest = lower_text(trim(adjustl(text(3:))))
            call parse_charge_suffix(rest, charge, valence, ok)
            if (ok) then
               symbol = two
               tail = rest
               ok = .true.
               return
            end if
         end if
      end if
      symbol = normalise_symbol(text(1:1))
      if (table_index_of(symbol) == 0) then
         symbol = '  '
         return
      end if
      rest = lower_text(trim(adjustl(text(2:))))
      call parse_charge_suffix(rest, charge, valence, ok)
      if (.not. ok) then
         symbol = '  '
         return
      end if
      tail = rest
      ok = .true.
   end subroutine split_species_name

   !> Parse the suffix of a species label: "", "val", "2-", "1+" or "+".
   pure subroutine parse_charge_suffix(text, charge, valence, ok)
      character(len=*), intent(in) :: text
      integer, intent(out) :: charge
      logical, intent(out) :: valence
      logical, intent(out) :: ok
      character(len=16) :: token
      integer :: i, n, digits

      charge = 0
      valence = .false.
      token = lower_text(text)
      n = len_trim(token)
      ok = .true.
      if (n == 0) return
      if (token(1:n) == 'val') then
         valence = .true.
         return
      end if
      digits = 0
      i = 1
      do while (i <= n)
         if (token(i:i) < '0' .or. token(i:i) > '9') exit
         digits = digits*10 + (iachar(token(i:i)) - iachar('0'))
         i = i + 1
      end do
      if (i /= n) then
         ok = .false.
         return
      end if
      if (token(i:i) /= '+' .and. token(i:i) /= '-') then
         ok = .false.
         return
      end if
      if (digits == 0 .and. i > 1) then
         ok = .false.
         return
      end if
      if (digits == 0) digits = 1
      charge = digits
      if (token(i:i) == '-') charge = -charge
   end subroutine parse_charge_suffix

   !> Canonical label of an element, charge and valence flag.
   pure function species_name(symbol, charge, valence) result(name)
      character(len=2), intent(in) :: symbol
      integer, intent(in) :: charge
      logical, intent(in) :: valence
      character(len=8) :: name
      character(len=8) :: digits

      if (valence) then
         name = trim(symbol)//'val'
      else if (charge == 0) then
         name = symbol
      else
         write (digits, '(i0)') abs(charge)
         name = trim(symbol)//trim(digits)//merge('+', '-', charge > 0)
      end if
   end function species_name

   !> Comma separated list of the tabulated states of one element.
   pure subroutine available_states(element, list)
      integer, intent(in) :: element
      character(len=*), intent(out) :: list
      integer :: idx

      list = ''
      do idx = 1, n_species
         if (species_element(idx) /= element) cycle
         if (len_trim(list) > 0) then
            list = trim(list)//', '//trim(species_label(idx))
         else
            list = species_label(idx)
         end if
      end do
   end subroutine available_states

end module sqc_elements
