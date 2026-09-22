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
!!   h5read file.h5 xrd      -> "# 2theta[deg] I" table
!!   h5read file.h5 partials -> "# q Si-Si Si-O ..." table
!!   h5read file.h5 rdf      -> "# r g(r) Si-Si Si-O ..." table
!!   h5read file.h5 s4       -> "# qx qy qz tau S4(q,t)" table
!!   h5read file.h5 chi4     -> "# tau Q(t) chi4(t)" table
!!   h5read file.h5 count G  -> the per-lag origin counts of the dynamic group G
!!   h5read file.h5 attr N   -> the value of the file attribute N
program h5read
   use, intrinsic :: iso_fortran_env, only: int32, int64, real64, error_unit, output_unit
   use hdf5
   implicit none
   character(len=512) :: path, which
   character(len=512) :: subject
   integer(hid_t) :: file_id, group_id, subgroup_id
   integer :: hdferr, i, j, n, npairs
   real(real64), allocatable :: q(:), s(:), qx(:), qy(:), qz(:)
   integer(int32), allocatable :: h(:), k(:), l(:)
   real(real64), allocatable :: part(:, :), pv(:)
   real(real64), allocatable :: part3(:, :, :)
   real(real64), allocatable :: q2(:, :), spec(:, :), spec2(:, :), axis(:)
   character(len=18), allocatable :: dnames(:)
   integer(int64), allocatable :: cnt(:, :)
   integer :: ncount1, ncount2
   integer :: nq, ncol, naxis, naxis_check, nq_check, npairs2, kk

   call get_command_argument(1, path)
   call get_command_argument(2, which)
   call get_command_argument(3, subject)
   if (len_trim(path) == 0 .or. len_trim(which) == 0) then
      write (error_unit, '(a)') 'usage: h5read FILE '// &
         '(grid|shell|partials|xrd|rdf|sqw|fqt|fqt_self|s4|chi4|pair_entropy|s2_accum|'// &
         'count GROUP|attr NAME)'
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
   case ('xrd')
      ! powder XRD pattern of the static methods
      call h5gopen_f(file_id, 'xrd', group_id, hdferr)
      call read_real_1d(group_id, 'two_theta', q, n)
      call read_real_1d(group_id, 'I', s, n)
      write (output_unit, '(a)') '# 2theta[deg] I'
      do i = 1, n
         write (output_unit, '(f12.4,2x,es20.12)') q(i), s(i)
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
   case ('sqw', 'fqt')
      ! dynamic structure factor (or F(q,t)): one row per (q, axis) sample,
      ! the same long table the text writer produces
      call h5gopen_f(file_id, trim(which), group_id, hdferr)
      call read_real_2d(group_id, 'q', q2, nq, ncol)
      call read_real_1d(group_id, axis_of(trim(which)), axis, naxis)
      if (trim(which) == 'sqw') then
         call read_real_2d(group_id, 'S', spec, nq, naxis_check)
         call read_string_1d(group_id, 'pairs', dnames, npairs)
         allocate (part3(nq, naxis, max(npairs, 1)))
         part3 = 0.0_real64
         do i = 1, npairs
            call h5gopen_f(group_id, 'S_partial', subgroup_id, hdferr)
            call read_real_2d(subgroup_id, trim(dnames(i)), spec2, nq_check, naxis_check)
            part3(:, :, i) = spec2
            call h5gclose_f(subgroup_id, hdferr)
            deallocate (spec2)
         end do
      else
         call read_real_2d(group_id, 'F', spec, nq, naxis_check)
         call read_string_1d(group_id, 'pairs', dnames, npairs)
         allocate (part3(nq, naxis, max(npairs, 1)))
         part3 = 0.0_real64
         do i = 1, npairs
            call h5gopen_f(group_id, 'F_partial', subgroup_id, hdferr)
            call read_real_2d(subgroup_id, trim(dnames(i)), spec2, nq_check, naxis_check)
            part3(:, :, i) = spec2
            call h5gclose_f(subgroup_id, hdferr)
            deallocate (spec2)
         end do
      end if
      write (output_unit, '(a)', advance='no') '# qx qy qz '//trim(axis_of(trim(which)))
      if (trim(which) == 'sqw') then
         write (output_unit, '(a)', advance='no') ' S(q,w)'
      else
         write (output_unit, '(a)', advance='no') ' F(q,t)'
      end if
      do i = 1, npairs
         write (output_unit, '(a)', advance='no') ' '//trim(dnames(i))
      end do
      write (output_unit, '(a)') ''
      do j = 1, nq
         do i = 1, naxis
            write (output_unit, '(3(f14.8,2x),f16.8,2x,es20.12)', advance='no') &
               q2(j, 1), q2(j, 2), q2(j, 3), axis(i), spec(j, i)
            do kk = 1, npairs
               write (output_unit, '(2x,es20.12)', advance='no') part3(j, i, kk)
            end do
            write (output_unit, '(a)') ''
         end do
      end do
      call h5gclose_f(group_id, hdferr)
   case ('s4')
      ! total four-point structure factor: one row per (q, lag) sample
      call h5gopen_f(file_id, 's4', group_id, hdferr)
      call read_real_2d(group_id, 'q', q2, nq, ncol)
      call read_real_1d(group_id, 'tau', axis, naxis)
      call read_real_2d(group_id, 'S4', spec, nq, naxis_check)
      write (output_unit, '(a)') '# qx qy qz tau S4(q,t)'
      do j = 1, nq
         do i = 1, naxis
            write (output_unit, '(3(f14.8,2x),f16.8,2x,es20.12)') &
               q2(j, 1), q2(j, 2), q2(j, 3), axis(i), spec(j, i)
         end do
      end do
      call h5gclose_f(group_id, hdferr)
   case ('chi4')
      ! average overlap Q(t) and dynamic susceptibility chi4(t)
      call h5gopen_f(file_id, 'chi4', group_id, hdferr)
      call read_real_1d(group_id, 'tau', axis, naxis)
      call read_real_1d(group_id, 'Q', q, naxis_check)
      call read_real_1d(group_id, 'chi4', s, naxis_check)
      write (output_unit, '(a)') '# tau Q(t) chi4(t)'
      do i = 1, naxis
         write (output_unit, '(f16.8,2x,es20.12,2x,es20.12)') axis(i), q(i), s(i)
      end do
      call h5gclose_f(group_id, hdferr)
   case ('fqt_self')
      ! self intermediate scattering function and its per-species parts
      call h5gopen_f(file_id, 'fqt_self', group_id, hdferr)
      call read_real_2d(group_id, 'q', q2, nq, ncol)
      call read_real_1d(group_id, 'tau', axis, naxis)
      call read_real_2d(group_id, 'F_s', spec, nq, naxis_check)
      call read_string_1d(group_id, 'pairs', dnames, npairs)
      allocate (part3(nq, naxis, max(npairs, 1)))
      part3 = 0.0_real64
      do i = 1, npairs
         call h5gopen_f(group_id, 'F_s_partial', subgroup_id, hdferr)
         call read_real_2d(subgroup_id, trim(dnames(i)), spec2, nq_check, naxis_check)
         part3(:, :, i) = spec2
         call h5gclose_f(subgroup_id, hdferr)
         deallocate (spec2)
      end do
      write (output_unit, '(a)', advance='no') '# qx qy qz tau F_s(q,t)'
      do i = 1, npairs
         write (output_unit, '(a)', advance='no') ' '//trim(dnames(i))
      end do
      write (output_unit, '(a)') ''
      do j = 1, nq
         do i = 1, naxis
            write (output_unit, '(3(f14.8,2x),f16.8,2x,es20.12)', advance='no') &
               q2(j, 1), q2(j, 2), q2(j, 3), axis(i), spec(j, i)
            do kk = 1, npairs
               write (output_unit, '(2x,es20.12)', advance='no') part3(j, i, kk)
            end do
            write (output_unit, '(a)') ''
         end do
      end do
      call h5gclose_f(group_id, hdferr)
   case ('pair_entropy')
      ! final pair entropy values: one numeric row, total then partials
      call h5gopen_f(file_id, 'pair_entropy', group_id, hdferr)
      call read_string_1d(group_id, 'pairs', dnames, npairs)
      call read_real_1d(group_id, 'S2', s, n)
      call read_real_1d(group_id, 'S2_total', q, nq_check)
      write (output_unit, '(a)', advance='no') '# S2_total'
      do i = 1, npairs
         write (output_unit, '(a)', advance='no') ' S2('//trim(dnames(i))//')'
      end do
      write (output_unit, '(a)') ''
      write (output_unit, '(es20.12)', advance='no') q(1)
      do i = 1, npairs
         write (output_unit, '(2x,es20.12)', advance='no') s(i)
      end do
      write (output_unit, '(a)') ''
      call h5gclose_f(group_id, hdferr)
   case ('s2_accum')
      ! accumulation curve: one row per r bin
      call h5gopen_f(file_id, 'pair_entropy', group_id, hdferr)
      call read_string_1d(group_id, 'pairs', dnames, npairs)
      call read_real_1d(group_id, 'r', q, n)
      call read_real_1d(group_id, 'S2_total_curve', s, naxis_check)
      call read_real_2d(group_id, 'S2_curve', spec, nq_check, npairs2)
      write (output_unit, '(a)', advance='no') '# r S2_total(r)'
      do i = 1, npairs
         write (output_unit, '(a)', advance='no') ' S2('//trim(dnames(i))//')(r)'
      end do
      write (output_unit, '(a)') ''
      do j = 1, n
         write (output_unit, '(f12.6,2x,es20.12)', advance='no') q(j), s(j)
         do i = 1, npairs
            write (output_unit, '(2x,es20.12)', advance='no') spec(j, i)
         end do
         write (output_unit, '(a)') ''
      end do
      call h5gclose_f(group_id, hdferr)
   case ('count')
      ! Origin count per (q, axis) sample of a dynamic group; the shell and
      ! chi4 groups need it to judge how much a long lag is worth.
      call h5gopen_f(file_id, trim(subject), group_id, hdferr)
      if (hdferr /= 0) then
         write (error_unit, '(a)') 'cannot open the group "'//trim(subject)//'"'
         stop 4
      end if
      call read_int_grid(group_id, 'count', cnt, ncount1, ncount2)
      write (output_unit, '(a,a,a)') '# ', trim(subject), '/count'
      do j = 1, ncount1
         write (output_unit, '(*(i0,1x))') (cnt(j, i), i = 1, ncount2)
      end do
      call h5gclose_f(group_id, hdferr)
   case ('attr')
      ! One file attribute, so that the suite can check the metadata a run
      ! writes next to its tables (the q sampling and, for a single lattice
      ! vector, the Miller indices, which |q| alone does not identify).
      call read_attr(file_id, trim(subject))
   case default
      write (error_unit, '(a)') 'unknown table "'//trim(which)//'"'
      stop 3
   end select
   call h5fclose_f(file_id, hdferr)
   call h5close_f(hdferr)

contains

   !> Name of the frequency (or time) axis dataset of a dynamic group.
   pure function axis_of(which) result(name)
      character(len=*), intent(in) :: which
      character(len=8) :: name
      if (trim(which) == 'sqw') then
         name = 'omega'
      else
         name = 'tau'
      end if
   end function axis_of

   !> Print one file attribute: an integer, a real or a string.
   !!
   !! The suite uses it to check the metadata a run writes next to its tables,
   !! in particular the q sampling and, for a single lattice vector, the Miller
   !! indices, which |q| alone does not identify.
   subroutine read_attr(loc, name)
      integer(hid_t), intent(in) :: loc
      character(len=*), intent(in) :: name
      integer(hid_t) :: attr_id, type_id
      integer :: hdferr, iclass
      integer(size_t) :: asize
      integer(int64) :: ibuf(1)
      real(real64) :: rbuf(1)
      character(len=:), allocatable :: sbuf(:)

      call h5aopen_f(loc, trim(name), attr_id, hdferr)
      if (hdferr /= 0) then
         write (error_unit, '(a)') 'no attribute "'//trim(name)//'"'
         stop 5
      end if
      call h5aget_type_f(attr_id, type_id, hdferr)
      call h5tget_class_f(type_id, iclass, hdferr)
      ! The HDF5 Fortran interface of this build exports the type classes as
      ! variables rather than parameters, so they cannot head a CASE.
      if (iclass == H5T_INTEGER_F) then
         call h5aread_f(attr_id, H5T_STD_I64LE, ibuf, [1_hsize_t], hdferr)
         write (output_unit, '(i0)') ibuf(1)
      else if (iclass == H5T_FLOAT_F) then
         call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, rbuf, [1_hsize_t], hdferr)
         write (output_unit, '(es20.12)') rbuf(1)
      else if (iclass == H5T_STRING_F) then
         call h5tget_size_f(type_id, asize, hdferr)
         allocate (character(len=max(int(asize), 1)) :: sbuf(1))
         call h5aread_f(attr_id, type_id, sbuf, [1_hsize_t], hdferr)
         write (output_unit, '(a)') trim(sbuf(1))
      else
         write (error_unit, '(a)') 'attribute "'//trim(name)//'" has an unsupported type'
         stop 6
      end if
      if (hdferr /= 0) then
         write (error_unit, '(a)') 'cannot read attribute "'//trim(name)//'"'
         stop 7
      end if
      call h5tclose_f(type_id, hdferr)
      call h5aclose_f(attr_id, hdferr)
   end subroutine read_attr

   !> Read a 2D real dataset.
   subroutine read_real_2d(group, name, values, n1, n2)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      real(real64), allocatable, intent(out) :: values(:, :)
      integer, intent(out) :: n1, n2
      integer(hid_t) :: dset_id, space_id
      integer(hsize_t) :: dims(2), maxdims(2)
      integer :: hdferr

      call h5dopen_f(group, trim(name), dset_id, hdferr)
      call h5dget_space_f(dset_id, space_id, hdferr)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
      n1 = int(dims(1))
      n2 = int(dims(2))
      allocate (values(n1, n2))
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, values, dims, hdferr)
      call h5sclose_f(space_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
   end subroutine read_real_2d

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

   !> Read the 64 bit origin counts of a group as rows.
   !!
   !! The coherent and S4 groups store one row per q mode, while chi4 is a
   !! scalar time series, so a 1D dataset is returned as a single row.
   subroutine read_int_grid(group, name, values, n1, n2)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      integer(int64), allocatable, intent(out) :: values(:, :)
      integer, intent(out) :: n1, n2
      integer(hid_t) :: dset_id, space_id
      integer(hsize_t) :: dims(2), maxdims(2)
      integer(hsize_t) :: dims1(1), maxdims1(1)
      integer :: hdferr, rank

      call h5dopen_f(group, trim(name), dset_id, hdferr)
      call h5dget_space_f(dset_id, space_id, hdferr)
      call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
      if (rank == 1) then
         call h5sget_simple_extent_dims_f(space_id, dims1, maxdims1, hdferr)
         n1 = 1
         n2 = int(dims1(1))
         allocate (values(n1, n2))
         call h5dread_f(dset_id, H5T_STD_I64LE, values(1, :), dims1, hdferr)
      else
         call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
         n1 = int(dims(1))
         n2 = int(dims(2))
         allocate (values(n1, n2))
         call h5dread_f(dset_id, H5T_STD_I64LE, values, dims, hdferr)
      end if
      call h5sclose_f(space_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
   end subroutine read_int_grid

   !> Read a fixed length string dataset (1D).
   subroutine read_string_1d(group, name, values, n)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      character(len=18), allocatable, intent(out) :: values(:)
      integer, intent(out) :: n
      integer(hid_t) :: dset_id, space_id, type_id
      integer(hsize_t) :: dims(1), maxdims(1)
      integer :: hdferr
      logical :: exists

      call h5lexists_f(group, trim(name), exists, hdferr)
      if (.not. exists) then
         ! optional dataset (for example the pair labels of a run written
         ! without partials): report an empty list instead of failing
         n = 0
         allocate (values(0))
         return
      end if
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
