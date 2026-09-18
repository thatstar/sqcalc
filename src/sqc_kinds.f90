!> Kind parameters used throughout sqcalc.
module sqc_kinds
   use, intrinsic :: iso_fortran_env, only: real64, int32, int64
   implicit none
   private

   !> Working real precision (double).
   integer, parameter, public :: rk = real64
   !> Default integer kind.
   integer, parameter, public :: ik = int32
   !> Kind used for array sizes / counts that may exceed 2^31.
   integer, parameter, public :: lk = int64

end module sqc_kinds
