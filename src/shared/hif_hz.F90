MODULE HIF_HZ
  USE BLAS_UTILS
  USE DTYPES
  IMPLICIT NONE

CONTAINS

#ifdef MKL_DIRECT_CALL
#include "mkl_direct_call.fi"
#endif

#include "hif.F90"
#include "hz.F90"
END MODULE HIF_HZ
