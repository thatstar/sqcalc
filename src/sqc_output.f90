! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Output container of the tables sqcalc writes.
!!
!! One vocabulary for the file formats, shared by the option parser, the
!! reciprocal and dynamic writers and the Debye pair histograms, so that
!! `--format` means the same thing everywhere and no module has to translate
!! another module's constant.
module sqc_output
   implicit none
   private

   public :: format_infer, format_text, format_hdf5, output_wants_hdf5

   !> format_infer (the default) reads the container off the file name.
   integer, parameter :: format_infer = -1
   integer, parameter :: format_text = 0
   integer, parameter :: format_hdf5 = 1

contains

   !> Does the table written to PATH go into an HDF5 container?
   !!
   !! FORMAT decides it when the user asked for one (`--format text|hdf5`);
   !! format_infer leaves it to the name, so a `.h5`/`.hdf5` suffix still
   !! selects HDF5 and every other name text.
   pure logical function output_wants_hdf5(format, path) result(hdf5)
      integer, intent(in) :: format
      character(len=*), intent(in) :: path

      select case (format)
      case (format_hdf5)
         hdf5 = .true.
      case (format_text)
         hdf5 = .false.
      case default
         hdf5 = index(path, '.h5') > 0 .or. index(path, '.hdf5') > 0
      end select
   end function output_wants_hdf5

end module sqc_output
