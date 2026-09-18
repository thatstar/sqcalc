! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Thin iso_c_binding wrapper around the double-precision FINUFFT C API.
!!
!! We bind the C API (finufft.h) rather than the F77 translation layer so that
!! sqcalc works with any FINUFFT build.  The finufft_opts layout below must stay
!! in sync with include/finufft_opts.h; finufft_opts_is_consistent() verifies the
!! layout at run time against the values written by finufft_default_opts().
module sqc_finufft
   use, intrinsic :: iso_c_binding
   use sqc_kinds, only: rk
   implicit none
   private

   public :: finufft_opts_t, finufft_default_opts, finufft_makeplan, finufft_setpts, &
             finufft_execute, finufft_destroy, finufft_opts_is_consistent, &
             mode_cmcl, mode_fft, type1

   !> Mode ordering options of finufft_opts.
   integer, parameter :: mode_cmcl = 0, mode_fft = 1
   !> Transform types supported by the guru interface.
   integer, parameter :: type1 = 1, type2 = 2, type3 = 3

   !> C-compatible mirror of finufft_opts (double precision).
   !!
   !! Field list taken from the vendored include/finufft_opts.h.  Note that the
   !! development branch of FINUFFT has one more field
   !! (allow_eps_too_small) which the v2.5.1 release does not have, so this
   !! mirror must be checked against the vendored header whenever FINUFFT is
   !! updated - finufft_opts_is_consistent() does that at run time.
   type, bind(c) :: finufft_opts_t
      integer(c_int) :: modeord
      integer(c_int) :: spreadinterponly
      integer(c_int) :: debug
      integer(c_int) :: spread_debug
      integer(c_int) :: showwarn
      integer(c_int) :: nthreads
      integer(c_int) :: fftw
      integer(c_int) :: spread_sort
      integer(c_int) :: spread_kerevalmeth
      integer(c_int) :: spread_kerpad
      real(c_double) :: upsampfac
      integer(c_int) :: spread_thread
      integer(c_int) :: maxbatchsize
      integer(c_int) :: spread_nthr_atomic
      integer(c_int) :: spread_max_sp_size
      integer(c_int) :: spread_kerformula
      type(c_funptr) :: fftw_lock_fun
      type(c_funptr) :: fftw_unlock_fun
      type(c_ptr) :: fftw_lock_data
   end type finufft_opts_t

   interface
      subroutine finufft_default_opts(opts) bind(c, name='finufft_default_opts')
         import :: finufft_opts_t
         type(finufft_opts_t), intent(out) :: opts
      end subroutine finufft_default_opts

      function finufft_makeplan(itype, ndim, n_modes, iflag, ntrans, tol, plan, opts) &
         bind(c, name='finufft_makeplan') result(ier)
         import :: c_int, c_int64_t, c_double, c_ptr, finufft_opts_t
         integer(c_int), value :: itype, ndim, iflag, ntrans
         integer(c_int64_t), intent(in) :: n_modes(*)
         real(c_double), value :: tol
         type(c_ptr), intent(out) :: plan
         type(finufft_opts_t), intent(in) :: opts
         integer(c_int) :: ier
      end function finufft_makeplan

      function finufft_setpts(plan, nj, xj, yj, zj, nk, s, t, u) &
         bind(c, name='finufft_setpts') result(ier)
         import :: c_int, c_int64_t, c_double, c_ptr
         type(c_ptr), value :: plan
         integer(c_int64_t), value :: nj, nk
         real(c_double), intent(in) :: xj(*), yj(*), zj(*), s(*), t(*), u(*)
         integer(c_int) :: ier
      end function finufft_setpts

      function finufft_execute(plan, weights, result) bind(c, name='finufft_execute') &
         result(ier)
         import :: c_int, c_double_complex, c_ptr
         type(c_ptr), value :: plan
         complex(c_double_complex), intent(in) :: weights(*)
         complex(c_double_complex), intent(inout) :: result(*)
         integer(c_int) :: ier
      end function finufft_execute

      function finufft_destroy(plan) bind(c, name='finufft_destroy') result(ier)
         import :: c_int, c_ptr
         type(c_ptr), value :: plan
         integer(c_int) :: ier
      end function finufft_destroy
   end interface

contains

   !> Sanity check that our finufft_opts mirror matches the linked library.
   !!
   !! Catches a stale static type definition (for example after updating the
   !! vendored FINUFFT) before it silently corrupts an actual transform.
   logical function finufft_opts_is_consistent() result(ok)
      type(finufft_opts_t) :: opts

      call finufft_default_opts(opts)
      ok = opts%modeord == 0 .and. opts%spreadinterponly == 0 &
           .and. opts%debug == 0 .and. opts%spread_debug == 0 .and. opts%showwarn == 1 &
           .and. opts%nthreads == 0 .and. (opts%fftw == 64 .or. opts%fftw == -1) &
           .and. opts%spread_sort == 2 .and. abs(opts%upsampfac) <= 1.0e-12_c_double &
           .and. opts%maxbatchsize == 0 .and. opts%spread_nthr_atomic == -1 &
           .and. opts%spread_max_sp_size == 0 .and. opts%spread_kerformula == 0 &
           .and. .not. c_associated(opts%fftw_lock_fun) &
           .and. .not. c_associated(opts%fftw_unlock_fun) &
           .and. .not. c_associated(opts%fftw_lock_data)
   end function finufft_opts_is_consistent

end module sqc_finufft
