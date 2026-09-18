! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Dump the tables of an sqcalc HDF5 result file as text.
!!
!! Used by the test suite to compare the HDF5 output against the text output
!! without depending on python HDF5 bindings:
!!   h5read file.h5 grid     -> "# qx qy qz S(q)" table
!!   h5read file.h5 shell    -> "# q S(q)" table
!!   h5read file.h5 partials -> "# q Si-Si Si-O ..." table
!!   h5read file.h5 rdf      -> "# r g(r) Si-Si Si-O ..." table
program h5read
   use, intrinsic :: iso_fortran_env, only: int32, int64, real64, error_unit, output_unit
   use hdf5
   implicit none
   character(len=512) :: path, which
   integer(hid_t) :: file_id, group_id, subgroup_id
   integer :: hdferr, i, j, n, npairs
   real(real64), allocatable :: q(:), s(:), qx(:), qy(:), qz(:)
   integer(int32), allocatable :: h(:), k(:), l(:)
   real(real64), allocatable :: part(:, :), pv(:)
   character(len=18), allocatable :: dnames(:)

   call get_command_argument(1, path)
   call get_command_argument(2, which)
   if (len_trim(path) == 0 .or. len_trim(which) == 0) then
      write (error_unit, '(a)') 'usage: h5read FILE (grid|shell|partials|rdf)'
      stop 1
   end if
   call h5open_f(hdferr)
   call h5fopen_f(trim(path), H5F_ACC_RDONLY_F, file_id, hdferr)
   if (hdferr /= 0) then
      write (error_unit, '(a)') 'cannot open '//trim(path)
      stop 2
   end if

   select case (trim(which))
   case ('shell')
      call h5gopen_f(file_id, 'shell', group_id, hdferr)
      call read_real_1d(group_id, 'q', q, n)
      call read_real_1d(group_id, 'S', s, n)
      write (output_unit, '(a)') '# q S(q)'
      do i = 1, n
         write (output_unit, '(f14.6,2x,es20.12)') q(i), s(i)
      end do
      call h5gclose_f(group_id, hdferr)
   case ('grid')
      call h5gopen_f(file_id, 'grid', group_id, hdferr)
      call read_real_1d(group_id, 'qx', qx, n)
      call read_real_1d(group_id, 'qy', qy, n)
      call read_real_1d(group_id, 'qz', qz, n)
      call read_real_1d(group_id, 'S', s, n)
      write (output_unit, '(a)') '# qx qy qz S(q)'
      do i = 1, n
         write (output_unit, '(3(f14.8,2x),es20.12)') qx(i), qy(i), qz(i), s(i)
      end do
      call h5gclose_f(group_id, hdferr)
   case ('partials')
      ! partial structure factors of the shell table, one dataset per pair
      call h5gopen_f(file_id, 'shell', group_id, hdferr)
      call read_real_1d(group_id, 'q', q, n)
      call read_string_1d(group_id, 'pairs', dnames, npairs)
      allocate (part(n, npairs))
      do i = 1, npairs
         call h5gopen_f(group_id, 'S_partial', subgroup_id, hdferr)
         call read_real_1d(subgroup_id, trim(dnames(i)), pv, n)
         part(:, i) = pv
         call h5gclose_f(subgroup_id, hdferr)
      end do
      write (output_unit, '(a)', advance='no') '# q'
      do i = 1, npairs
         write (output_unit, '(a)', advance='no') ' '//trim(dnames(i))
      end do
      write (output_unit, '(a)') ''
      do i = 1, n
         write (output_unit, '(f14.6)', advance='no') q(i)
         do j = 1, npairs
            write (output_unit, '(2x,es20.12)', advance='no') part(i, j)
         end do
         write (output_unit, '(a)') ''
      end do
      call h5gclose_f(group_id, hdferr)
   case ('rdf')
      ! pair distribution functions of the Debye method: total g(r) plus one
      ! dataset per pair, under the same labels as the text file
      call h5gopen_f(file_id, 'rdf', group_id, hdferr)
      call read_real_1d(group_id, 'r', q, n)
      call read_real_1d(group_id, 'g', s, n)
      call read_string_1d(group_id, 'pairs', dnames, npairs)
      allocate (part(n, npairs))
      do i = 1, npairs
         call h5gopen_f(group_id, 'g_partial', subgroup_id, hdferr)
         call read_real_1d(subgroup_id, trim(dnames(i)), pv, n)
         part(:, i) = pv
         call h5gclose_f(subgroup_id, hdferr)
      end do
      write (output_unit, '(a)', advance='no') '# r g(r)'
      do i = 1, npairs
         write (output_unit, '(a)', advance='no') ' '//trim(dnames(i))
      end do
      write (output_unit, '(a)') ''
      do i = 1, n
         write (output_unit, '(f14.6,2x,es20.12)', advance='no') q(i), s(i)
         do j = 1, npairs
            write (output_unit, '(2x,es20.12)', advance='no') part(i, j)
         end do
         write (output_unit, '(a)') ''
      end do
      call h5gclose_f(group_id, hdferr)
   case default
      write (error_unit, '(a)') 'unknown table "'//trim(which)//'"'
      stop 3
   end select
   call h5fclose_f(file_id, hdferr)
   call h5close_f(hdferr)

contains

   subroutine read_real_1d(group, name, values, n)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      real(real64), allocatable, intent(out) :: values(:)
      integer, intent(out) :: n
      integer(hid_t) :: dset_id, space_id
      integer(hsize_t) :: dims(1), maxdims(1)
      integer :: hdferr

      call h5dopen_f(group, trim(name), dset_id, hdferr)
      call h5dget_space_f(dset_id, space_id, hdferr)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
      n = int(dims(1))
      allocate (values(n))
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, values, dims, hdferr)
      call h5sclose_f(space_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
   end subroutine read_real_1d

   !> Read a fixed length string dataset (1D).
   subroutine read_string_1d(group, name, values, n)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      character(len=18), allocatable, intent(out) :: values(:)
      integer, intent(out) :: n
      integer(hid_t) :: dset_id, space_id, type_id
      integer(hsize_t) :: dims(1), maxdims(1)
      integer :: hdferr

      call h5dopen_f(group, trim(name), dset_id, hdferr)
      call h5dget_space_f(dset_id, space_id, hdferr)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
      n = int(dims(1))
      allocate (values(n))
      call h5dget_type_f(dset_id, type_id, hdferr)
      call h5dread_f(dset_id, type_id, values, dims, hdferr)
      call h5tclose_f(type_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
   end subroutine read_string_1d

end program h5read
