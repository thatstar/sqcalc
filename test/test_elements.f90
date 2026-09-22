! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Species labels: neutral atoms, ions, valence states and the aliases.
!!
!! The X-ray form factor rows are checked against the electron-count sum rule
!! f(0) = Z - q, which catches a mislabelled or mis-transcribed row without a
!! second copy of the table.  The lookup itself is exercised through the same
!! entry point the command line mapping uses.
program test_elements
   use sqc_kinds
   use sqc_elements, only: element_table
   use sqc_element_data, only: n_species, n_species_alias, species_element, &
                               species_charge, species_valence, species_label
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   implicit none

   character(len=256) :: message
   character(len=8) :: label
   integer :: idx, element, ierr, s, worst_label
   real(rk) :: f0, expected, worst

   ! --- neutral atoms -----------------------------------------------------
   call lookup('Si', idx, element, label, ierr, message)
   call require(ierr == 0 .and. trim(label) == 'Si', 'Si resolves to the neutral row', 0.0_rk)
   call require(element_table%species_charge_of(idx) == 0 .and. &
                .not. element_table%species_is_valence(idx), 'neutral Si carries no charge', 0.0_rk)
   call require(idx == element_table%neutral_species(element), &
                'neutral_species agrees with the label lookup', 0.0_rk)

   ! --- ions, case insensitive, with and without the explicit digit -------
   call lookup('si4+', idx, element, label, ierr, message)
   call require(ierr == 0 .and. trim(label) == 'Si4+', 'si4+ canonicalises to Si4+', 0.0_rk)
   call require(element_table%species_element_of(idx) == element_table%index_of('Si'), &
                'Si4+ belongs to the silicon element row', 0.0_rk)
   call lookup('O2-', idx, element, label, ierr, message)
   call require(ierr == 0 .and. element_table%species_charge_of(idx) == -2, &
                'O2- is the doubly charged anion', 0.0_rk)
   call require(abs(element_table%species_xray_f(idx, 0.0_rk) - 10.0_rk) < 1.0e-2_rk, &
                'O2- has ten electrons at q = 0', &
                abs(element_table%species_xray_f(idx, 0.0_rk) - 10.0_rk))

   ! --- valence states ----------------------------------------------------
   call lookup('Cval', idx, element, label, ierr, message)
   call require(ierr == 0 .and. element_table%species_is_valence(idx), &
                'Cval is a valence state', 0.0_rk)
   call require(abs(element_table%species_xray_f(idx, 0.0_rk) - 6.0_rk) < 1.0e-2_rk, &
                'a valence state keeps the neutral electron count', &
                abs(element_table%species_xray_f(idx, 0.0_rk) - 6.0_rk))

   ! --- aliases of the duplicate source rows ------------------------------
   call lookup('Siv', idx, element, label, ierr, message)
   call require(ierr == 0 .and. trim(label) == 'Si', 'Siv is an alias of Si', 0.0_rk)
   call lookup("H'", idx, element, label, ierr, message)
   call require(ierr == 0 .and. trim(label) == 'D', "H' is an alias of D", 0.0_rk)

   ! --- an element without an X-ray row is not an error of the label ------
   call lookup('Es', idx, element, label, ierr, message)
   call require(ierr == 0 .and. idx == 0, 'Es parses but has no X-ray species', 0.0_rk)
   call require(element == element_table%index_of('Es'), 'Es keeps its element index', 0.0_rk)

   ! --- the two error paths ----------------------------------------------
   call lookup('Xx', idx, element, label, ierr, message)
   call require(ierr == 1, 'an unknown element is reported as such', 0.0_rk)
   call lookup('O3-', idx, element, label, ierr, message)
   call require(ierr == 2, 'an unknown charge state is reported as such', 0.0_rk)
   call require(index(message, 'O2-') > 0, 'the error lists the available states', 0.0_rk)

   ! --- every row against the electron-count sum rule ---------------------
   worst = 0.0_rk
   worst_label = 0
   do s = 1, n_species
      f0 = element_table%species_xray_f(s, 0.0_rk)
      expected = real(element_table%z_of(species_element(s)), rk)
      if (.not. species_valence(s)) expected = expected - real(species_charge(s), rk)
      if (abs(f0 - expected) > worst) then
         worst = abs(f0 - expected)
         worst_label = s
      end if
   end do
   call require(worst < 0.1_rk, 'every species satisfies f(0) = Z - q', worst)
   write (output_unit, '(a,i0,a,f0.4)') '     largest deviation over ', n_species, &
      ' species: ', worst
   if (worst_label > 0) write (output_unit, '(a,a)') '     at ', trim(species_label(worst_label))

   write (output_unit, '(a,i0,a,i0,a)') 'ok elements: ', n_species, ' species, ', &
      n_species_alias, ' aliases'

contains

   !> Look a label up and keep the five outputs of the lookup together.
   subroutine lookup(text, sidx, elem, name, status, note)
      character(len=*), intent(in) :: text
      integer, intent(out) :: sidx, elem, status
      character(len=*), intent(out) :: name, note
      call element_table%find_species(text, sidx, elem, name, status, note)
   end subroutine lookup

   !> Fail loudly with the value that was out of tolerance.
   subroutine require(condition, label, value)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      real(rk), intent(in) :: value

      if (condition) return
      write (error_unit, '(a,a,a,es12.4)') 'FAIL: ', trim(label), ', deviation ', value
      stop 1
   end subroutine require

end program test_elements
