! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Access to the periodic table and to X-ray / neutron scattering amplitudes.
module sqc_elements
   use sqc_kinds
   use sqc_element_data, only: n_elements, element_z, element_symbol, element_name, &
                               element_mass, element_neutron_b, element_has_neutron, &
                               element_xray_a, element_xray_b, element_xray_c, element_has_xray
   implicit none
   private

   public :: element_table_t, element_table, normalise_symbol

   !> Stateless lookup object for the bundled element data.
   type :: element_table_t
   contains
      procedure, nopass :: size => table_size
      procedure, nopass :: index_of => table_index_of
      procedure, nopass :: symbol_of => table_symbol_of
      procedure, nopass :: name_of => table_name_of
      procedure, nopass :: mass_of => table_mass_of
      procedure, nopass :: has_xray => table_has_xray
      procedure, nopass :: has_neutron => table_has_neutron
      procedure, nopass :: neutron_b => table_neutron_b
      procedure, nopass :: xray_f => table_xray_f
   end type element_table_t

   !> Shared instance used by the rest of the program.
   type(element_table_t), parameter :: element_table = element_table_t()

contains

   pure integer function table_size() result(n)
      n = n_elements
   end function table_size

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

   pure logical function table_has_xray(idx) result(found)
      integer, intent(in) :: idx
      found = idx >= 1 .and. idx <= n_elements .and. element_has_xray(idx)
   end function table_has_xray

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

   !> X-ray atomic form factor in electron units at momentum transfer q [1/A].
   !!
   !! f(q) = sum_i a_i exp(-b_i (q/4pi)^2) + c   (IT92 / Cromer-Mann)
   pure real(rk) function table_xray_f(idx, q) result(f)
      integer, intent(in) :: idx
      real(rk), intent(in) :: q
      real(rk) :: stol2
      integer :: i

      f = 0.0_rk
      if (idx < 1 .or. idx > n_elements) return
      stol2 = (q/(4.0_rk*acos(-1.0_rk)))**2
      f = element_xray_c(idx)
      do i = 1, 4
         f = f + element_xray_a(i, idx)*exp(-element_xray_b(i, idx)*stol2)
      end do
   end function table_xray_f

end module sqc_elements
