!> GPU implementation of the structure factor, using cufinufft (cuFFT backend).
!!
!! The reciprocal grid, the normalization and the accumulators are exactly the
!! same as in the CPU NUFFT method (both extend structure_factor_t); the
!! difference is that the type-1 transforms run on the GPU.  Per frame we upload
!! the wrapped coordinates and the per-species strengths, run one batched
!! transform (one entry per species), copy the grid amplitudes back and do the
!! cheap combine/binning on the host.
module sqc_gpu
   use sqc_kinds
   use sqc_cell, only: two_pi
   use sqc_dump, only: frame_t
   use sqc_weights, only: weight_scheme_t
   use sqc_structure, only: structure_factor_t, sf_prepare_species
   use sqc_finufft, only: mode_cmcl, type1
   use sqc_cufinufft
   use, intrinsic :: iso_c_binding
   implicit none
   private

   public :: cufinufft_structure_factor_t

   !> Host staging buffers used by the device transfers.  They live at module
   !! scope because Fortran does not allow the TARGET attribute on components
   !! (c_loc needs it), and they are sized once per run for the single method
   !! instance that sqcalc creates.
   real(rk), allocatable, target :: host_x(:), host_y(:), host_z(:)
   complex(c_double_complex), allocatable, target :: host_c(:), host_fk(:)
   !> float32 staging buffers, used when the transform runs in single precision.
   real(c_float), allocatable, target :: host_x32(:), host_y32(:), host_z32(:)
   complex(c_float_complex), allocatable, target :: host_c32(:), host_fk32(:)
   !> |rho(q)|^2 of the kept modes for the current frame (parallel scratch).
   real(rk), allocatable :: mode_intensity(:)

   !> Structure factor evaluated with cufinufft on a CUDA device.
   type, extends(structure_factor_t) :: cufinufft_structure_factor_t
      !> cufinufft plan handle.
      type(c_ptr) :: plan = c_null_ptr
      !> Options passed to the plan (and the device id).
      type(cufinufft_opts_t) :: opts
      !> CUDA device to use.
      integer :: device = 0
      !> Evaluate the transform in single precision (cufinufftf) instead of double.
      logical :: single_precision = .false.
      !> Number of LAMMPS types present in the trajectory.
      integer :: nspecies = 0
      !> Species index of each atom and LAMMPS type id of each species.
      integer(ik), allocatable :: species_of(:)
      integer(ik), allocatable :: species_type(:)
      !> True when the amplitudes depend on q (X-ray).
      logical :: q_dependent = .false.
      !> Per-species amplitudes (constant, or per kept mode).
      real(rk), allocatable :: amp_const(:)
      real(rk), allocatable :: amp_table(:, :)
      !> Device buffers.
      type(c_ptr) :: d_x = c_null_ptr
      type(c_ptr) :: d_y = c_null_ptr
      type(c_ptr) :: d_z = c_null_ptr
      type(c_ptr) :: d_c = c_null_ptr
      type(c_ptr) :: d_fk = c_null_ptr
   contains
      procedure :: method_setup => gpu_setup
      procedure :: accumulate_frame => gpu_accumulate
      final :: gpu_finalize
   end type cufinufft_structure_factor_t

contains

   !> Bytes per coordinate value in device memory (float32 or float64).
   pure integer(c_size_t) function device_real_bytes(self) result(nbytes)
      class(cufinufft_structure_factor_t), intent(in) :: self
      nbytes = merge(bytes_per_real32, bytes_per_double, self%single_precision)
   end function device_real_bytes

   !> Bytes per transform value in device memory (complex float32 or float64).
   pure integer(c_size_t) function device_complex_bytes(self) result(nbytes)
      class(cufinufft_structure_factor_t), intent(in) :: self
      nbytes = merge(bytes_per_complex32, bytes_per_complex, self%single_precision)
   end function device_complex_bytes

   !> Allocate host and device buffers and create the cufinufft plan.
   subroutine gpu_setup(self, frame, scheme, ierr, message)
      class(cufinufft_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      integer(c_int64_t) :: n_modes(3)
      integer(c_int) :: istat, ier, ndevices

      ierr = 0
      message = ''

      if (.not. cufinufft_opts_is_consistent()) then
         ierr = 1
         message = 'the linked cufinufft uses an incompatible cufinufft_opts layout'
         return
      end if
      istat = cuda_device_count(ndevices)
      if (istat /= cuda_success .or. ndevices < 1) then
         ierr = 1
         message = 'no usable CUDA device ('//cuda_error_message(istat)//')'
         return
      end if
      if (self%device < 0 .or. self%device >= ndevices) then
         ierr = 1
         write (message, '(a,i0,a,i0,a)') 'CUDA device ', self%device, ' requested but ', &
            ndevices, ' device(s) available'
         return
      end if

      call sf_prepare_species(self, frame, scheme, self%species_of, self%species_type, &
                              self%amp_const, self%amp_table, self%q_dependent, ierr, message)
      if (ierr /= 0) return
      self%nspecies = size(self%species_type)

      ! Host staging buffers.
      allocate (host_x(frame%natoms), host_y(frame%natoms), host_z(frame%natoms))
      allocate (host_c(int(self%nspecies, lk)*int(frame%natoms, lk)))
      allocate (host_fk(int(self%nspecies, lk)*self%gridpoints))
      allocate (mode_intensity(self%nmodes))
      host_x = 0.0_rk
      host_y = 0.0_rk
      host_z = 0.0_rk
      host_c = (0.0_rk, 0.0_rk)
      host_fk = (0.0_rk, 0.0_rk)
      if (self%single_precision) then
         allocate (host_x32(frame%natoms), host_y32(frame%natoms), host_z32(frame%natoms))
         allocate (host_c32(int(self%nspecies, lk)*int(frame%natoms, lk)))
         allocate (host_fk32(int(self%nspecies, lk)*self%gridpoints))
         host_x32 = 0.0_c_float
         host_y32 = 0.0_c_float
         host_z32 = 0.0_c_float
         host_c32 = (0.0_c_float, 0.0_c_float)
         host_fk32 = (0.0_c_float, 0.0_c_float)
      end if

      ! Device buffers.
      istat = cuda_malloc(self%d_x, int(frame%natoms, c_size_t)*device_real_bytes(self))
      if (istat == cuda_success) &
         istat = cuda_malloc(self%d_y, int(frame%natoms, c_size_t)*device_real_bytes(self))
      if (istat == cuda_success) &
         istat = cuda_malloc(self%d_z, int(frame%natoms, c_size_t)*device_real_bytes(self))
      if (istat == cuda_success) &
         istat = cuda_malloc(self%d_c, int(self%nspecies, c_size_t) &
                             *int(frame%natoms, c_size_t)*device_complex_bytes(self))
      if (istat == cuda_success) &
         istat = cuda_malloc(self%d_fk, int(self%nspecies, c_size_t)*self%gridpoints &
                             *device_complex_bytes(self))
      if (istat /= cuda_success) then
         ierr = 1
         message = 'cannot allocate device memory ('//cuda_error_message(istat)//')'
         return
      end if

      call cufinufft_default_opts(self%opts)
      self%opts%modeord = mode_cmcl
      self%opts%gpu_device_id = int(self%device, c_int)
      n_modes = int(self%modes, c_int64_t)
      if (self%single_precision) then
         ier = cufinufftf_makeplan(type1, 3_c_int, n_modes, 1_c_int, int(self%nspecies, c_int), &
                                   real(self%eps, c_float), self%plan, self%opts)
      else
         ier = cufinufft_makeplan(type1, 3_c_int, n_modes, 1_c_int, int(self%nspecies, c_int), &
                                  real(self%eps, c_double), self%plan, self%opts)
      end if
      if (ier /= 0) then
         ierr = 1
         message = 'cufinufft failed to create a plan: '//cufinufft_error_message(ier)
         return
      end if
   end subroutine gpu_setup

   !> Run one frame through the GPU transform and accumulate |rho(q)|^2.
   subroutine gpu_accumulate(self, frame, scheme, ierr, message)
      class(cufinufft_structure_factor_t), intent(inout) :: self
      type(frame_t), intent(in) :: frame
      type(weight_scheme_t), intent(in) :: scheme
      integer, intent(out) :: ierr
      character(len=*), intent(out) :: message
      real(rk) :: s(3), amp, value, pi
      integer(c_int) :: istat, ier
      integer :: i, im, isp, sh, base
      integer(lk) :: g
      complex(c_double_complex) :: total

      ierr = 0
      message = ''
      pi = acos(-1.0_rk)

      ! Wrapped fractional coordinates scaled into [-pi, pi).
      do i = 1, frame%natoms
         s = frame%cell%fractional(frame%pos(:, i))
         s = s - floor(s)
         host_x(i) = two_pi*s(1) - pi
         host_y(i) = two_pi*s(2) - pi
         host_z(i) = two_pi*s(3) - pi
      end do
      if (self%single_precision) then
         host_x32 = real(host_x, c_float)
         host_y32 = real(host_y, c_float)
         host_z32 = real(host_z, c_float)
         istat = cuda_upload_real32(self%d_x, host_x32)
         if (istat == cuda_success) istat = cuda_upload_real32(self%d_y, host_y32)
         if (istat == cuda_success) istat = cuda_upload_real32(self%d_z, host_z32)
      else
         istat = cuda_upload_real(self%d_x, host_x)
         if (istat == cuda_success) istat = cuda_upload_real(self%d_y, host_y)
         if (istat == cuda_success) istat = cuda_upload_real(self%d_z, host_z)
      end if
      if (istat /= cuda_success) then
         ierr = 1
         message = 'cannot upload coordinates ('//cuda_error_message(istat)//')'
         return
      end if

      if (self%single_precision) then
         ier = cufinufftf_setpts(self%plan, int(frame%natoms, c_int64_t), self%d_x, self%d_y, &
                                 self%d_z, 0_c_int, c_null_ptr, c_null_ptr, c_null_ptr)
      else
         ier = cufinufft_setpts(self%plan, int(frame%natoms, c_int64_t), self%d_x, self%d_y, &
                                self%d_z, 0_c_int, c_null_ptr, c_null_ptr, c_null_ptr)
      end if
      if (ier /= 0) then
         ierr = 1
         message = 'cufinufft setpts failed: '//cufinufft_error_message(ier)
         return
      end if

      ! Strengths: one contiguous block of natoms entries per species.
      host_c = (0.0_rk, 0.0_rk)
      do isp = 1, self%nspecies
         base = int(isp - 1, ik)*frame%natoms
         do i = 1, frame%natoms
            if (self%species_of(i) == isp) host_c(base + i) = (1.0_rk, 0.0_rk)
         end do
      end do
      if (self%single_precision) then
         host_c32 = cmplx(host_c, kind=c_float_complex)
         istat = cuda_upload_complex32(self%d_c, host_c32)
      else
         istat = cuda_upload_complex(self%d_c, host_c)
      end if
      if (istat /= cuda_success) then
         ierr = 1
         message = 'cannot upload strengths ('//cuda_error_message(istat)//')'
         return
      end if

      if (self%single_precision) then
         ier = cufinufftf_execute(self%plan, self%d_c, self%d_fk)
      else
         ier = cufinufft_execute(self%plan, self%d_c, self%d_fk)
      end if
      if (ier /= 0) then
         ierr = 1
         message = 'cufinufft execute failed: '//cufinufft_error_message(ier)
         return
      end if
      if (self%single_precision) then
         istat = cuda_download_complex32(self%d_fk, host_fk32)
         if (istat == cuda_success) then
            ! Widen for the (double precision) combine and accumulation.
            !$omp parallel do schedule(static)
            do i = 1, int(size(host_fk, kind=lk))
               host_fk(i) = cmplx(real(host_fk32(i)), aimag(host_fk32(i)), &
                                  kind=c_double_complex)
            end do
            !$omp end parallel do
         end if
      else
         istat = cuda_download_complex(self%d_fk, host_fk)
      end if
      if (istat /= cuda_success) then
         ierr = 1
         message = 'cannot download the transform ('//cuda_error_message(istat)//')'
         return
      end if

      ! Combine the species with their amplitudes (parallel over modes: each
      ! mode needs its own species sum, but nothing is shared) and accumulate.
      !$omp parallel do schedule(static) private(im, isp, g, amp, total)
      do im = 1, int(self%nmodes)
         g = self%gidx(im)
         total = (0.0_rk, 0.0_rk)
         do isp = 1, self%nspecies
            if (self%q_dependent) then
               amp = self%amp_table(im, isp)
            else
               amp = self%amp_const(isp)
            end if
            total = total + cmplx(amp, 0.0_rk, c_double_complex) &
                            *host_fk(int(isp - 1, lk)*self%gridpoints + g)
         end do
         mode_intensity(im) = real(total*conjg(total), rk)
      end do
      !$omp end parallel do
      do im = 1, int(self%nmodes)
         value = mode_intensity(im)
         g = self%gidx(im)
         sh = self%shell(im)
         self%num(sh) = self%num(sh) + value
         if (self%want_grid) self%gnum(g) = self%gnum(g) + value
      end do
      self%nframes = self%nframes + 1
   end subroutine gpu_accumulate

   !> Release device memory and the plan.
   subroutine gpu_finalize(self)
      type(cufinufft_structure_factor_t), intent(inout) :: self
      integer(c_int) :: ier

      if (c_associated(self%plan)) then
         ier = cufinufft_destroy(self%plan)
         self%plan = c_null_ptr
      end if
      if (c_associated(self%d_fk)) then
         ier = cuda_free(self%d_fk)
         self%d_fk = c_null_ptr
      end if
      if (c_associated(self%d_c)) then
         ier = cuda_free(self%d_c)
         self%d_c = c_null_ptr
      end if
      if (c_associated(self%d_z)) then
         ier = cuda_free(self%d_z)
         self%d_z = c_null_ptr
      end if
      if (c_associated(self%d_y)) then
         ier = cuda_free(self%d_y)
         self%d_y = c_null_ptr
      end if
      if (c_associated(self%d_x)) then
         ier = cuda_free(self%d_x)
         self%d_x = c_null_ptr
      end if
      if (allocated(host_x)) deallocate (host_x)
      if (allocated(host_y)) deallocate (host_y)
      if (allocated(host_z)) deallocate (host_z)
      if (allocated(host_c)) deallocate (host_c)
      if (allocated(host_fk)) deallocate (host_fk)
      if (allocated(host_x32)) deallocate (host_x32)
      if (allocated(host_y32)) deallocate (host_y32)
      if (allocated(host_z32)) deallocate (host_z32)
      if (allocated(host_c32)) deallocate (host_c32)
      if (allocated(host_fk32)) deallocate (host_fk32)
      if (allocated(mode_intensity)) deallocate (mode_intensity)
   end subroutine gpu_finalize

end module sqc_gpu
