MODULE BLAS_UTILS
#ifdef USE_MKL
#ifdef USE_INTEL
#ifdef USE_X200
  USE IFCORE
#endif
  USE MKL_SERVICE
#endif
#endif
  USE OMP_LIB
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION BLAS_PREPARE()
    IMPLICIT NONE
#ifdef USE_MKL
#ifndef USE_INTEL
#ifndef MKL_NEST_SEQ
    EXTERNAL :: MKL_SET_DYNAMIC
#endif
#else
#ifdef USE_X200
    IF (FOR_GET_HBW_AVAILABILITY() .NE. FOR_K_HBW_AVAILABLE) THEN
       BLAS_PREPARE = FOR_SET_FASTMEM_POLICY(FOR_K_FASTMEM_RETRY_WARN)
       BLAS_PREPARE = MKL_SET_MEMORY_LIMIT(MKL_MEM_MCDRAM, 0)
    END IF
#endif
#endif
#endif
    CALL OMP_SET_NESTED(.TRUE._c_int)
    CALL OMP_SET_DYNAMIC(.FALSE._c_int)
#ifdef USE_MKL
#ifndef MKL_NEST_SEQ
    CALL MKL_SET_DYNAMIC(1_c_int)
#endif
#endif
    BLAS_PREPARE = INT(OMP_GET_MAX_ACTIVE_LEVELS())
  END FUNCTION BLAS_PREPARE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION BLAS_SET_NUM_THREADS(NT)
#ifdef USE_MKL
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NT
#ifdef MKL_NEST_SEQ
    BLAS_SET_NUM_THREADS = MIN(MAX(0,NT),1)
#else
#ifdef USE_GNU
    INTEGER(c_int), EXTERNAL :: MKL_SET_NUM_THREADS_LOCAL
#endif

    BLAS_SET_NUM_THREADS = INT(MKL_SET_NUM_THREADS_LOCAL(INT(NT,c_int)))
#endif
#else
    USE OMP_LIB
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT

    BLAS_SET_NUM_THREADS = INT(OMP_GET_NUM_THREADS())
    CALL OMP_SET_NUM_THREADS(INT(NT,c_int))
#endif
  END FUNCTION BLAS_SET_NUM_THREADS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE BLAS_ZROTM_RR(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    DOUBLE PRECISION, INTENT(IN) :: C1, C2
    DOUBLE PRECISION, INTENT(IN) :: S1
    DOUBLE PRECISION, INTENT(IN) :: S2
    DOUBLE COMPLEX, INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_RR

  PURE SUBROUTINE BLAS_ZROTM_RC(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    DOUBLE PRECISION, INTENT(IN) :: C1, C2
    DOUBLE PRECISION, INTENT(IN) :: S1
    DOUBLE COMPLEX, INTENT(IN) :: S2
    DOUBLE COMPLEX, INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_RC

  PURE SUBROUTINE BLAS_ZROTM_CR(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    DOUBLE PRECISION, INTENT(IN) :: C1, C2
    DOUBLE COMPLEX, INTENT(IN) :: S1
    DOUBLE PRECISION, INTENT(IN) :: S2
    DOUBLE COMPLEX, INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_CR

  PURE SUBROUTINE BLAS_ZROTM_CC(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    DOUBLE PRECISION, INTENT(IN) :: C1, C2
    DOUBLE COMPLEX, INTENT(IN) :: S1
    DOUBLE COMPLEX, INTENT(IN) :: S2
    DOUBLE COMPLEX, INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_CC

  PURE SUBROUTINE BLAS_ZROTM(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    DOUBLE PRECISION, INTENT(IN) :: C1, C2
    DOUBLE COMPLEX, INTENT(IN) :: S1, S2
    DOUBLE COMPLEX, INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    IF (AIMAG(S1) .EQ. D_ZERO) THEN
       ! Rx
       IF (AIMAG(S2) .EQ. D_ZERO) THEN
          ! RR
          CALL BLAS_ZROTM_RR(SIDE, N, ZX, INCX, ZY, INCY, C1, DBLE(S1), DBLE(S2), C2, INFO)
       ELSE
          ! RC
          CALL BLAS_ZROTM_RC(SIDE, N, ZX, INCX, ZY, INCY, C1, DBLE(S1), S2, C2, INFO)
       END IF
    ELSE
       ! Cx
       IF (AIMAG(S2) .EQ. D_ZERO) THEN
          ! CR
          CALL BLAS_ZROTM_CR(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, DBLE(S2), C2, INFO)
       ELSE
          ! CC
          CALL BLAS_ZROTM_CC(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
       END IF
    END IF
  END SUBROUTINE BLAS_ZROTM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE BLAS_ZVROTM(N, ZX, ZY, C1, S1, S2, C2)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: C1, C2
    DOUBLE COMPLEX, INTENT(IN) :: S1, S2
    DOUBLE COMPLEX, INTENT(INOUT) :: ZX(*), ZY(*)

    DOUBLE COMPLEX :: W(DSIMDL), Z(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: W,Z

    DOUBLE PRECISION ReW(DSIMDL),ImW(DSIMDL), ReZ(DSIMDL),ImZ(DSIMDL), Re1(DSIMDL),Im1(DSIMDL), Re2(DSIMDL),Im2(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: ReW,ImW, ReZ,ImZ, Re1,Im1, Re2,Im2
    DOUBLE PRECISION :: aC1,aC2, ReS1,ImS1, ReS2,ImS2
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: aC1,aC2, ReS1,ImS1, ReS2,ImS2
    INTEGER :: I, J, K

    !DIR$ ASSUME_ALIGNED ZX:ALIGNB
    !DIR$ ASSUME_ALIGNED ZY:ALIGNB

    IF (N .LE. 0) RETURN

    aC1 = C1
    aC2 = C2
    ReS1 = DBLE(S1)
    ImS1 = AIMAG(S1)
    ReS2 = DBLE(S2)
    ImS2 = AIMAG(S2)

    DO I = 1, N, DSIMDL
       K = MIN(DSIMDL, N-(I-1))
       !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
       DO J = 1, K
          W(J) = ZX(I+(J-1))
          Z(J) = ZY(I+(J-1))
          ReW(J) = DBLE(W(J))
          ImW(J) = AIMAG(W(J))
          ReZ(J) = DBLE(Z(J))
          ImZ(J) = AIMAG(Z(J))
          Re1(J) = ReW(J) * aC1 + ReZ(J) * ReS1 - ImZ(J) * ImS1
          Im1(J) = ImW(J) * aC1 + ImZ(J) * ReS1 + ReZ(J) * ImS1
          Re2(J) = ReW(J) * ReS2 - ImW(J) * ImS2 + ReZ(J) * aC2
          Im2(J) = ReW(J) * ImS2 + ImW(J) * ReS2 + ImZ(J) * aC2
          ZX(I+(J-1)) = DCMPLX(Re1(J), Im1(J))
          ZY(I+(J-1)) = DCMPLX(Re2(J), Im2(J))
       END DO
    END DO
  END SUBROUTINE BLAS_ZVROTM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE BLAS_ZSWAPC(N, ZX, INCX, ZY, INCY)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: INCX, INCY, N
    DOUBLE COMPLEX, INTENT(INOUT) :: ZX(*), ZY(*)

    DOUBLE COMPLEX :: ZTEMP
    INTEGER :: I, IX, IY

    IF (N .LE. 0) RETURN
    IF ((INCX .EQ. 1) .AND. (INCY .EQ. 1)) THEN
       DO I = 1, N
          ZTEMP = DCONJG(ZX(I))
          ZX(I) = DCONJG(ZY(I))
          ZY(I) = ZTEMP
       END DO
    ELSE
       I = 1 - N
       IF (INCX .LT. 0) THEN
          IX = I*INCX + 1
       ELSE
          IX = 1
       END IF
       IF (INCY .LT. 0) THEN
          IY = I*INCY + 1
       ELSE
          IY = 1
       END IF
       DO I = 1, N
          ZTEMP = DCONJG(ZX(IX))
          ZX(IX) = DCONJG(ZY(IY))
          ZY(IY) = ZTEMP
          IX = IX + INCX
          IY = IY + INCY
       END DO
    END IF
  END SUBROUTINE BLAS_ZSWAPC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE BLAS_UTILS
