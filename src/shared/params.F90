MODULE PARAMS
#ifndef NDEBUG
  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE, INTRINSIC :: IEEE_FEATURES
#endif
  USE, INTRINSIC :: ISO_C_BINDING
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
  USE VN_TYPES_F
  IMPLICIT NONE

  INTEGER, PARAMETER :: UINP = INPUT_UNIT
  INTEGER, PARAMETER :: UOUT = OUTPUT_UNIT
  INTEGER, PARAMETER :: ULOG = ERROR_UNIT

  ! Max file name length.
  INTEGER, PARAMETER :: FNL = 252

#ifdef MAX_CORES_PER_RUN
  INTEGER, PARAMETER :: MAXCPR = MAX_CORES_PER_RUN
#else
  INTEGER, PARAMETER :: MAXCPR = 96 ! max. number of 2nd-level threads
#endif

#ifdef MAX_THREADS_PER_CORE
  INTEGER, PARAMETER :: MAXTPC = MAX_THREADS_PER_CORE
#else
  INTEGER, PARAMETER :: MAXTPC = 64 ! max. number of 1st-level threads
#endif

#include "vn_params.F90"
  REAL(KIND=DWP), PARAMETER :: D_CS_PI_4 = SQRT(D_TWO) / D_TWO
  REAL(KIND=DWP), PARAMETER :: D_ALPHA = (D_ONE + SQRT(17.0_DWP)) / 8.0_DWP

CONTAINS
END MODULE PARAMS
