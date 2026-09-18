!> HDF5 output for the structure factor results.
!!
!! The file holds the shell averaged table and, when requested, the full
!! reciprocal grid.  Datasets are plain one dimensional arrays (no dimension
!! ordering ambiguity), so any HDF5 reader (h5py, MATLAB, Julia, h5dump, ...)
!! sees them directly:
!!
!!   /shell/q        nq     doubles   shell centers [1/A]
!!   /shell/S        nq     doubles   S(q) averaged over the shell
!!   /shell/count    nq     int64     reciprocal lattice points per shell
!!   /grid/qx,qy,qz  nmodes doubles   reciprocal lattice vector [1/A]
!!   /grid/h,k,l     nmodes int32     Miller indices
!!   /grid/S         nmodes doubles   S(q) at that reciprocal lattice point
!!
!! plus file attributes describing the run (cell, weights, mapping, frames, ...).
module sqc_hdf5
   use, intrinsic :: iso_fortran_env, only: int32, int64, real64
   use, intrinsic :: iso_c_binding, only: c_loc, c_ptr
   use hdf5
   use sqc_kinds
   use sqc_structure, only: structure_factor_t
   use sqc_weights, only: weight_scheme_t
   implicit none
   private

   public :: hdf5_write_results, hdf5_support, hdf5_write_rdf

   !> HDF5 was compiled into this binary.
   logical, parameter :: hdf5_support = .true.

contains

   !> Write the shell table, the reciprocal grid and the run metadata.
   subroutine hdf5_write_results(path, method, scheme, natoms, cell, eps, device, &
                                 precision, ierr, message)
      character(len=*), intent(in) :: path
      class(structure_factor_t), intent(in) :: method
      type(weight_scheme_t), intent(in) :: scheme
      integer(ik), intent(in) :: natoms
      real(rk), intent(in) :: cell(3, 3)
      real(rk), intent(in) :: eps
      character(len=*), intent(in) :: device, precision
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: file_id, group_id, space_id, dset_id, attr_id, type_id
      integer(hsize_t) :: dims(1)
      integer :: hdferr, s
      integer(lk), allocatable :: counts(:)
      real(real64), allocatable :: qd(:), sd(:), qx(:), qy(:), qz(:)
      integer(int32), allocatable :: hh(:), kk(:), ll(:)
      real(real64) :: qc, cell_dbl(3, 3)
      character(len=256) :: mapping
      integer(int64) :: nframes
      integer :: nmodes_txt

      ierr = 0
      message = ''

      call h5open_f(hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5 library initialisation failed'
         return
      end if
      call h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'cannot create HDF5 file "'//trim(path)//'"'
         call h5close_f(hdferr)
         return
      end if

      ! --- attributes describing the run -------------------------------------
      call write_int_attr(file_id, 'nframes', int(method%nframes, int64), ierr, message)
      if (ierr /= 0) return
      call write_int_attr(file_id, 'natoms', int(natoms, int64), ierr, message)
      if (ierr /= 0) return
      call write_int_attr(file_id, 'nmodes', int(method%nmodes, int64), ierr, message)
      if (ierr /= 0) return
      call write_int_attr(file_id, 'nq', int(method%nq, int64), ierr, message)
      if (ierr /= 0) return
      call write_real_attr(file_id, 'qmin', real(method%qmin, real64), ierr, message)
      if (ierr /= 0) return
      call write_real_attr(file_id, 'qmax', real(method%qmax, real64), ierr, message)
      if (ierr /= 0) return
      call write_real_attr(file_id, 'eps', real(eps, real64), ierr, message)
      if (ierr /= 0) return
      cell_dbl = real(cell, real64)
      call h5screate_simple_f(2, [3_hsize_t, 3_hsize_t], space_id, hdferr)
      call check(hdferr, 'create cell attribute space')
      if (ierr /= 0) return
      call h5acreate_f(file_id, 'cell', H5T_NATIVE_DOUBLE, space_id, attr_id, hdferr)
      call check(hdferr, 'create cell attribute')
      if (ierr /= 0) return
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, cell_dbl, [3_hsize_t, 3_hsize_t], hdferr)
      call check(hdferr, 'write cell attribute')
      if (ierr /= 0) return
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      call write_string_attr(file_id, 'weight', trim(scheme%label()), ierr, message)
      if (ierr /= 0) return
      call write_string_attr(file_id, 'norm', norm_label(method%norm), ierr, message)
      if (ierr /= 0) return
      call write_string_attr(file_id, 'device', device, ierr, message)
      if (ierr /= 0) return
      call write_string_attr(file_id, 'precision', precision, ierr, message)
      if (ierr /= 0) return
      mapping = ''
      do s = 1, scheme%mapped_types()
         if (len_trim(scheme%symbols(s)) == 0) cycle
         if (len_trim(mapping) > 0) mapping = trim(mapping)//','
         write (mapping, '(a,i0,a,a)') trim(mapping), s, ':', trim(scheme%symbols(s))
      end do
      call write_string_attr(file_id, 'mapping', trim(mapping), ierr, message)
      if (ierr /= 0) return

      ! --- shell averaged table ---------------------------------------------
      call h5gcreate_f(file_id, 'shell', group_id, hdferr)
      call check(hdferr, 'create /shell group')
      if (ierr /= 0) return
      allocate (qd(method%nq), sd(method%nq), counts(method%nq))
      do s = 1, method%nq
         qd(s) = real(method%qmin + (real(s, rk) - 0.5_rk)*method%shell_dq, real64)
         sd(s) = real(method%shell_value(s), real64)
         counts(s) = method%shell_count(s)
      end do
      dims = [int(method%nq, hsize_t)]
      call write_dataset_f(group_id, 'q', H5T_NATIVE_DOUBLE, dims, qd, ierr, message)
      if (ierr /= 0) return
      call write_dataset_f(group_id, 'S', H5T_NATIVE_DOUBLE, dims, sd, ierr, message)
      if (ierr /= 0) return
      call write_dataset_i8_f(group_id, 'count', dims, counts, ierr, message)
      if (ierr /= 0) return
      deallocate (qd, sd, counts)
      call h5gclose_f(group_id, hdferr)

      ! --- reciprocal grid ---------------------------------------------------
      if (method%want_grid) then
         call h5gcreate_f(file_id, 'grid', group_id, hdferr)
         call check(hdferr, 'create /grid group')
         if (ierr /= 0) return
         nmodes_txt = int(method%nmodes)
         allocate (qx(method%nmodes), qy(method%nmodes), qz(method%nmodes), &
                   sd(method%nmodes), hh(method%nmodes), kk(method%nmodes), &
                   ll(method%nmodes))
         do s = 1, nmodes_txt
            qx(s) = real(method%qvec(1, s), real64)
            qy(s) = real(method%qvec(2, s), real64)
            qz(s) = real(method%qvec(3, s), real64)
            sd(s) = real(method%grid_value(int(s, lk)), real64)
            hh(s) = int(method%hkl(1, s), int32)
            kk(s) = int(method%hkl(2, s), int32)
            ll(s) = int(method%hkl(3, s), int32)
         end do
         dims = [int(method%nmodes, hsize_t)]
         call write_dataset_f(group_id, 'qx', H5T_NATIVE_DOUBLE, dims, qx, ierr, message)
         if (ierr /= 0) return
         call write_dataset_f(group_id, 'qy', H5T_NATIVE_DOUBLE, dims, qy, ierr, message)
         if (ierr /= 0) return
         call write_dataset_f(group_id, 'qz', H5T_NATIVE_DOUBLE, dims, qz, ierr, message)
         if (ierr /= 0) return
         call write_dataset_f(group_id, 'S', H5T_NATIVE_DOUBLE, dims, sd, ierr, message)
         if (ierr /= 0) return
         call write_dataset_i4_f(group_id, 'h', dims, hh, ierr, message)
         if (ierr /= 0) return
         call write_dataset_i4_f(group_id, 'k', dims, kk, ierr, message)
         if (ierr /= 0) return
         call write_dataset_i4_f(group_id, 'l', dims, ll, ierr, message)
         if (ierr /= 0) return
         deallocate (qx, qy, qz, sd, hh, kk, ll)
         call h5gclose_f(group_id, hdferr)
      end if

      call h5fclose_f(file_id, hdferr)
      call check(hdferr, 'close the file')
      call h5close_f(hdferr)

   contains

      !> Turn a failing HDF5 call into an sqcalc error.
      subroutine check(hdferr, what)
         integer, intent(in) :: hdferr
         character(len=*), intent(in) :: what
         if (hdferr /= 0) then
            ierr = 1
            write (message, '(a,a,a,i0,a)') 'HDF5: cannot ', trim(what), &
               ' (error ', hdferr, ')'
         end if
      end subroutine check
   end subroutine hdf5_write_results

   !> Text label of a normalization code (kept in sync with sqc_structure).
   pure function norm_label(norm) result(label)
      integer, intent(in) :: norm
      character(len=:), allocatable :: label
      select case (norm)
      case (2)
         label = 'self'
      case (3)
         label = 'n'
      case default
         label = 'mean'
      end select
   end function norm_label

   subroutine write_dataset_f(group, name, htype, dims, values, ierr, message)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      integer(hid_t), intent(in) :: htype
      integer(hsize_t), intent(in) :: dims(:)
      real(real64), intent(in) :: values(:)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: space_id, dset_id
      integer :: hdferr

      ierr = 0
      message = ''
      call h5screate_simple_f(size(dims), dims, space_id, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot create the dataspace of '//trim(name)
         return
      end if
      call h5dcreate_f(group, trim(name), htype, space_id, dset_id, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot create dataset '//trim(name)
         return
      end if
      call h5dwrite_f(dset_id, htype, values, dims, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot write dataset '//trim(name)
         return
      end if
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
   end subroutine write_dataset_f

   subroutine write_dataset_i4_f(group, name, dims, values, ierr, message)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      integer(hsize_t), intent(in) :: dims(:)
      integer(int32), intent(in) :: values(:)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: space_id, dset_id
      integer :: hdferr

      ierr = 0
      message = ''
      call h5screate_simple_f(size(dims), dims, space_id, hdferr)
      call h5dcreate_f(group, trim(name), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, values, dims, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot write dataset '//trim(name)
         return
      end if
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
   end subroutine write_dataset_i4_f

   subroutine write_dataset_i8_f(group, name, dims, values, ierr, message)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      integer(hsize_t), intent(in) :: dims(:)
      integer(int64), intent(in) :: values(:)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: space_id, dset_id
      integer :: hdferr

      ierr = 0
      message = ''
      call h5screate_simple_f(size(dims), dims, space_id, hdferr)
      call h5dcreate_f(group, trim(name), H5T_STD_I64LE, space_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_STD_I64LE, values, dims, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot write dataset '//trim(name)
         return
      end if
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
   end subroutine write_dataset_i8_f

   subroutine write_real_attr(obj, name, value, ierr, message)
      integer(hid_t), intent(in) :: obj
      character(len=*), intent(in) :: name
      real(real64), intent(in) :: value
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: space_id, attr_id
      integer :: hdferr
      real(real64) :: buffer(1)

      ierr = 0
      message = ''
      buffer(1) = value
      call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
      call h5acreate_f(obj, trim(name), H5T_NATIVE_DOUBLE, space_id, attr_id, hdferr)
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, buffer, [1_hsize_t], hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot write attribute '//trim(name)
         return
      end if
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(space_id, hdferr)
   end subroutine write_real_attr

   subroutine write_int_attr(obj, name, value, ierr, message)
      integer(hid_t), intent(in) :: obj
      character(len=*), intent(in) :: name
      integer(int64), intent(in) :: value
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: space_id, attr_id
      integer :: hdferr
      integer(int64) :: buffer(1)

      ierr = 0
      message = ''
      buffer(1) = value
      call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
      call h5acreate_f(obj, trim(name), H5T_STD_I64LE, space_id, attr_id, hdferr)
      call h5awrite_f(attr_id, H5T_STD_I64LE, buffer, [1_hsize_t], hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot write attribute '//trim(name)
         return
      end if
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(space_id, hdferr)
   end subroutine write_int_attr

   !> Fixed length string attribute (portable, readable by any HDF5 tool).
   subroutine write_string_attr(obj, name, text, ierr, message)
      integer(hid_t), intent(in) :: obj
      character(len=*), intent(in) :: name, text
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: space_id, attr_id, type_id
      integer(size_t) :: text_len
      integer :: hdferr
      character(len=len(text)) :: buffer(1)

      ierr = 0
      message = ''
      buffer(1) = text
      ! HDF5 rejects a zero length string type, so an empty value (for example
      ! the mapping attribute when no -m was given) is stored as a single blank.
      if (len(text) == 0) buffer(1) = ' '
      text_len = len(text)
      if (text_len == 0) text_len = 1_size_t
      call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
      call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, hdferr)
      call h5tset_size_f(type_id, text_len, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot create the string type of attribute '//trim(name)
         return
      end if
      call h5acreate_f(obj, trim(name), type_id, space_id, attr_id, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot create attribute '//trim(name)
         return
      end if
      call h5awrite_f(attr_id, type_id, buffer, [1_hsize_t], hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot write attribute '//trim(name)
         return
      end if
      call h5aclose_f(attr_id, hdferr)
      call h5tclose_f(type_id, hdferr)
      call h5sclose_f(space_id, hdferr)
   end subroutine write_string_attr

   !> Pair distribution functions: `/rdf/r`, `/rdf/g` (total), one dataset per
   !! pair under `/rdf/g/<label>` and the label list in `/rdf/pairs`.
   subroutine hdf5_write_rdf(path, r, g_total, g_partial, labels, ierr, message)
      character(len=*), intent(in) :: path
      real(rk), intent(in) :: r(:), g_total(:), g_partial(:, :, :)
      character(len=*), intent(in) :: labels(:)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: file_id, group_id, subgroup_id
      integer(hsize_t) :: dims(1), dims2(2)
      integer :: hdferr, ia, ib, k, p, ntypes, nbins, npair, idum
      real(real64), allocatable :: one_d(:)
      character(len=18), allocatable :: label_buf(:)

      ierr = 0
      message = ''
      ntypes = size(g_partial, 1)
      nbins = size(g_partial, 3)
      npair = size(labels)

      call h5open_f(hdferr)
      call h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, file_id, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'cannot create HDF5 file "'//trim(path)//'"'
         call h5close_f(idum)
         return
      end if
      call h5gcreate_f(file_id, 'rdf', group_id, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot create the /rdf group'
         call h5fclose_f(file_id, idum)
         return
      end if

      dims = [int(nbins, hsize_t)]
      allocate (one_d(nbins))
      one_d = real(r, real64)
      call write_dataset_f(group_id, 'r', H5T_NATIVE_DOUBLE, dims, one_d, ierr, message)
      if (ierr /= 0) return
      one_d = real(g_total, real64)
      call write_dataset_f(group_id, 'g', H5T_NATIVE_DOUBLE, dims, one_d, ierr, message)
      if (ierr /= 0) return

      ! labels of the pair list
      allocate (label_buf(npair))
      do p = 1, npair
         label_buf(p) = labels(p)
      end do
      dims2 = [int(npair, hsize_t), 1_hsize_t]
      call write_string_dataset(group_id, 'pairs', label_buf, ierr, message)
      if (ierr /= 0) return

      ! one dataset per partial pair, under /rdf/g/<label>
      call h5gcreate_f(group_id, 'g', subgroup_id, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot create the /rdf/g group'
         return
      end if
      p = 0
      do ia = 1, ntypes
         do ib = ia, ntypes
            p = p + 1
            do k = 1, nbins
               one_d(k) = real(g_partial(ia, ib, k), real64)
            end do
            call write_dataset_f(subgroup_id, trim(label_buf(p)), H5T_NATIVE_DOUBLE, &
                                 dims, one_d, ierr, message)
            if (ierr /= 0) return
         end do
      end do

      call h5gclose_f(subgroup_id, hdferr)
      call h5gclose_f(group_id, hdferr)
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)
      deallocate (one_d, label_buf)
   end subroutine hdf5_write_rdf

   !> Fixed length string dataset (1D).
   subroutine write_string_dataset(group, name, values, ierr, message)
      integer(hid_t), intent(in) :: group
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: values(:)
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(hid_t) :: space_id, dset_id, type_id
      integer(hsize_t) :: dims(1)
      integer(size_t) :: text_len
      integer :: hdferr
      character(len=18), allocatable :: buffer(:)

      ierr = 0
      message = ''
      allocate (buffer(size(values)))
      buffer = values
      text_len = len(buffer(1))
      dims = [int(size(values), hsize_t)]
      call h5screate_simple_f(1, dims, space_id, hdferr)
      call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, hdferr)
      call h5tset_size_f(type_id, text_len, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot create the string type of '//trim(name)
         return
      end if
      call h5dcreate_f(group, trim(name), type_id, space_id, dset_id, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot create dataset '//trim(name)
         return
      end if
      call h5dwrite_f(dset_id, type_id, buffer, dims, hdferr)
      if (hdferr /= 0) then
         ierr = 1
         message = 'HDF5: cannot write dataset '//trim(name)
         return
      end if
      call h5dclose_f(dset_id, hdferr)
      call h5tclose_f(type_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      deallocate (buffer)
   end subroutine write_string_dataset

end module sqc_hdf5
