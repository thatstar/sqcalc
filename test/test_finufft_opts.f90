!> Verify that the Fortran mirror of finufft_opts matches the linked library.
program test_finufft_opts
   use sqc_finufft, only: finufft_opts_is_consistent
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   implicit none

   if (.not. finufft_opts_is_consistent()) then
      write (error_unit, '(a)') 'FAIL: finufft_opts layout does not match the linked FINUFFT'
      error stop 1
   end if
   write (output_unit, '(a)') 'PASS: finufft_opts layout matches the linked FINUFFT'
end program test_finufft_opts
