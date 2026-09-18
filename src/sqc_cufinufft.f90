!> iso_c_binding wrapper for the cuFFT based cufinufft library.
!!
!! cufinufft is built as a separate library next to the CPU FINUFFT
!! (FINUFFT_USE_CUDA=ON).  Its C API mirrors the CPU guru interface but works on
!! device pointers, so this module also wraps the few CUDA runtime calls that
!! sqcalc needs to move data between host and device.
!!
!! The cufinufft_opts mirror must stay in sync with include/cufinufft_opts.h;
!! cufinufft_opts_is_consistent() checks the layout at run time.
module sqc_cufinufft
   use, intrinsic :: iso_c_binding
   implicit none
   private

   public :: cufinufft_opts_t, cufinufft_default_opts, cufinufft_makeplan, &
             cufinufft_setpts, cufinufft_execute, cufinufft_destroy, &
             cufinufft_opts_is_consistent, cufinufft_error_message, &
             cuda_success, cuda_malloc, cuda_free, cuda_upload_real, &
             cuda_download_real, cuda_upload_complex, cuda_download_complex, &
             cuda_device_count, cuda_error_message, bytes_per_double, &
             bytes_per_complex

   !> CUDA runtime return code for success.
   integer(c_int), parameter :: cuda_success = 0
   !> cudaMemcpyKind values.
   integer(c_int), parameter :: cuda_memcpy_host_to_device = 1
   integer(c_int), parameter :: cuda_memcpy_device_to_host = 2
   !> Storage sizes of the interoperable types used here.
   integer(c_size_t), parameter :: bytes_per_double = 8_c_size_t
   integer(c_size_t), parameter :: bytes_per_complex = 16_c_size_t

   !> C-compatible mirror of cufinufft_opts (FINUFFT v2.5.x, double precision).
   type, bind(c) :: cufinufft_opts_t
      real(c_double) :: upsampfac
      integer(c_int) :: gpu_method
      integer(c_int) :: gpu_sort
      integer(c_int) :: gpu_binsizex
      integer(c_int) :: gpu_binsizey
      integer(c_int) :: gpu_binsizez
      integer(c_int) :: gpu_obinsizex
      integer(c_int) :: gpu_obinsizey
      integer(c_int) :: gpu_obinsizez
      integer(c_int) :: gpu_maxsubprobsize
      integer(c_int) :: gpu_kerevalmeth
      integer(c_int) :: gpu_spreadinterponly
      integer(c_int) :: gpu_maxbatchsize
      integer(c_int) :: gpu_device_id
      type(c_ptr) :: gpu_stream
      integer(c_int) :: modeord
      integer(c_int) :: gpu_np
      integer(c_int) :: debug
   end type cufinufft_opts_t

   interface
      subroutine cufinufft_default_opts(opts) bind(c, name='cufinufft_default_opts')
         import :: cufinufft_opts_t
         type(cufinufft_opts_t), intent(out) :: opts
      end subroutine cufinufft_default_opts

      function cufinufft_makeplan(itype, ndim, n_modes, iflag, ntrans, eps, plan, opts) &
         bind(c, name='cufinufft_makeplan') result(ier)
         import :: c_int, c_int64_t, c_double, c_ptr, cufinufft_opts_t
         integer(c_int), value :: itype, ndim, iflag, ntrans
         integer(c_int64_t), intent(in) :: n_modes(*)
         real(c_double), value :: eps
         type(c_ptr), intent(out) :: plan
         type(cufinufft_opts_t), intent(inout) :: opts
         integer(c_int) :: ier
      end function cufinufft_makeplan

      function cufinufft_setpts(plan, m, dx, dy, dz, n, ds, dt, du) &
         bind(c, name='cufinufft_setpts') result(ier)
         import :: c_int, c_int64_t, c_ptr
         type(c_ptr), value :: plan, dx, dy, dz, ds, dt, du
         integer(c_int64_t), value :: m
         integer(c_int), value :: n
         integer(c_int) :: ier
      end function cufinufft_setpts

      function cufinufft_execute(plan, strengths, modes) bind(c, name='cufinufft_execute') &
         result(ier)
         import :: c_int, c_ptr
         type(c_ptr), value :: plan, strengths, modes
         integer(c_int) :: ier
      end function cufinufft_execute

      function cufinufft_destroy(plan) bind(c, name='cufinufft_destroy') result(ier)
         import :: c_int, c_ptr
         type(c_ptr), value :: plan
         integer(c_int) :: ier
      end function cufinufft_destroy

      ! --- CUDA runtime -----------------------------------------------------
      function cuda_malloc(devptr, nbytes) bind(c, name='cudaMalloc') result(istat)
         import :: c_ptr, c_size_t, c_int
         type(c_ptr), intent(out) :: devptr
         integer(c_size_t), value :: nbytes
         integer(c_int) :: istat
      end function cuda_malloc

      function cuda_free(devptr) bind(c, name='cudaFree') result(istat)
         import :: c_ptr, c_int
         type(c_ptr), value :: devptr
         integer(c_int) :: istat
      end function cuda_free

      function cuda_memcpy(dst, src, count, kind) bind(c, name='cudaMemcpy') result(istat)
         import :: c_ptr, c_size_t, c_int
         type(c_ptr), value :: dst, src
         integer(c_size_t), value :: count
         integer(c_int), value :: kind
         integer(c_int) :: istat
      end function cuda_memcpy

      function cuda_device_count(ndevices) bind(c, name='cudaGetDeviceCount') result(istat)
         import :: c_int
         integer(c_int), intent(out) :: ndevices
         integer(c_int) :: istat
      end function cuda_device_count

      function cuda_error_string(code) bind(c, name='cudaGetErrorString') result(text)
         import :: c_int, c_ptr
         integer(c_int), value :: code
         type(c_ptr) :: text
      end function cuda_error_string
   end interface

contains

   !> Copy a host real array to the device.
   integer(c_int) function cuda_upload_real(devptr, host) result(istat)
      type(c_ptr), value :: devptr
      real(c_double), target, contiguous, intent(in) :: host(:)

      istat = cuda_memcpy(devptr, c_loc(host(1)), int(size(host), c_size_t)*bytes_per_double, &
                          cuda_memcpy_host_to_device)
   end function cuda_upload_real

   !> Copy device data into a host real array.
   integer(c_int) function cuda_download_real(devptr, host) result(istat)
      type(c_ptr), value :: devptr
      real(c_double), target, contiguous, intent(out) :: host(:)

      istat = cuda_memcpy(c_loc(host(1)), devptr, int(size(host), c_size_t)*bytes_per_double, &
                          cuda_memcpy_device_to_host)
   end function cuda_download_real

   !> Copy a host complex array to the device.
   integer(c_int) function cuda_upload_complex(devptr, host) result(istat)
      type(c_ptr), value :: devptr
      complex(c_double_complex), target, contiguous, intent(in) :: host(:)

      istat = cuda_memcpy(devptr, c_loc(host(1)), int(size(host), c_size_t)*bytes_per_complex, &
                          cuda_memcpy_host_to_device)
   end function cuda_upload_complex

   !> Copy device data into a host complex array.
   integer(c_int) function cuda_download_complex(devptr, host) result(istat)
      type(c_ptr), value :: devptr
      complex(c_double_complex), target, contiguous, intent(out) :: host(:)

      istat = cuda_memcpy(c_loc(host(1)), devptr, int(size(host), c_size_t)*bytes_per_complex, &
                          cuda_memcpy_device_to_host)
   end function cuda_download_complex

   !> Human readable CUDA error text.
   function cuda_error_message(code) result(message)
      integer(c_int), intent(in) :: code
      character(len=:), allocatable :: message
      type(c_ptr) :: text
      character(kind=c_char), pointer :: buffer(:)
      integer :: n

      text = cuda_error_string(code)
      if (.not. c_associated(text)) then
         message = 'unknown CUDA error'
         return
      end if
      call c_f_pointer(text, buffer, [256])
      n = 0
      do while (n < size(buffer) .and. buffer(n + 1) /= c_null_char)
         n = n + 1
      end do
      allocate (character(len=max(n, 1)) :: message)
      do n = 1, len(message)
         message(n:n) = buffer(n)
      end do
   end function cuda_error_message

   !> Human readable cufinufft error code (the library returns FINUFFT error codes).
   function cufinufft_error_message(code) result(message)
      integer(c_int), intent(in) :: code
      character(len=:), allocatable :: message
      character(len=64) :: text

      write (text, '(a,i0)') 'cufinufft error code ', code
      message = trim(text)
   end function cufinufft_error_message

   !> Sanity check that our cufinufft_opts mirror matches the linked library.
   logical function cufinufft_opts_is_consistent() result(ok)
      type(cufinufft_opts_t) :: opts

      call cufinufft_default_opts(opts)
      ok = abs(opts%upsampfac) <= 1.0e-12_c_double .and. opts%gpu_method == 0 &
           .and. opts%gpu_sort == 1 .and. opts%gpu_kerevalmeth == 1 &
           .and. opts%gpu_maxsubprobsize == 1024 .and. opts%gpu_obinsizex == 0 &
           .and. opts%gpu_obinsizey == 0 .and. opts%gpu_obinsizez == 0 &
           .and. opts%gpu_binsizex == 0 .and. opts%gpu_binsizey == 0 &
           .and. opts%gpu_binsizez == 0 .and. opts%gpu_maxbatchsize == 0 &
           .and. opts%gpu_np == 0 .and. opts%debug == 0 .and. opts%modeord == 0 &
           .and. opts%gpu_device_id == 0 .and. opts%gpu_spreadinterponly == 0 &
           .and. .not. c_associated(opts%gpu_stream)
   end function cufinufft_opts_is_consistent

end module sqc_cufinufft
