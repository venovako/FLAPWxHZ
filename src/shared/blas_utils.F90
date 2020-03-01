MODULE BLAS_UTILS
#ifdef USE_MKL
#ifdef USE_INTEL
#ifdef USE_X200
  USE IFCORE
#endif
  USE MKL_SERVICE
#endif
#endif
  USE PARAMS
  USE VN_BLAS_F
#ifdef HAVE_IMAGINARY
  USE VN_IMAGINARY_F
#endif
  USE OMP_LIB
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BLAS_ALLOC(PTR, AL, SZ, INFO)
    IMPLICIT NONE

    INTERFACE
       FUNCTION C_ALIGNED_MALLOC(PTR, AL, SZ) BIND(C,NAME='posix_memalign')
         USE, INTRINSIC :: ISO_C_BINDING
         TYPE(c_ptr), INTENT(OUT), TARGET :: PTR
         INTEGER(c_size_t), INTENT(IN), VALUE :: AL, SZ
         INTEGER(c_int) :: C_ALIGNED_MALLOC
       END FUNCTION C_ALIGNED_MALLOC
    END INTERFACE
#ifndef NDEBUG
    INTERFACE
       FUNCTION C_MEMSET(S, C, N) BIND(C,NAME='memset')
         USE, INTRINSIC :: ISO_C_BINDING
         TYPE(c_ptr), INTENT(IN), VALUE :: S
         INTEGER(c_int), INTENT(IN), VALUE :: C
         INTEGER(c_size_t), INTENT(IN), VALUE :: N
         INTEGER(c_intptr_t) :: C_MEMSET
       END FUNCTION C_MEMSET
    END INTERFACE
#endif
    TYPE(c_ptr), INTENT(OUT) :: PTR
    INTEGER, INTENT(IN) :: AL, SZ
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: MAL

    IF (AL .LE. 0) THEN
       INFO = -2
    ELSE IF (SZ .LE. 0) THEN
       INFO = -3
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    ! align to at least the page size...
    MAL = MAX(AL,PGSIZB)
    INFO = INT(C_ALIGNED_MALLOC(PTR, INT(MAL,c_size_t), INT(SZ,c_size_t)))
#ifndef NDEBUG
    ! fill the data with -1/qNaNs
    IF (INFO .EQ. 0) THEN
       MAL = INT(C_MEMSET(PTR, INT(-1,c_int), INT(SZ,c_size_t)))
       IF (MAL .EQ. 0) INFO = -1
    END IF
#endif
  END SUBROUTINE BLAS_ALLOC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BLAS_FREE(PTR)
    IMPLICIT NONE

    INTERFACE
       SUBROUTINE C_ALIGNED_FREE(PTR) BIND(C,NAME='free')
         USE, INTRINSIC :: ISO_C_BINDING
         TYPE(c_ptr), VALUE :: PTR
       END SUBROUTINE C_ALIGNED_FREE
    END INTERFACE

    TYPE(c_ptr), INTENT(IN) :: PTR

    CALL C_ALIGNED_FREE(PTR)
  END SUBROUTINE BLAS_FREE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION BLAS_PREPARE()
    IMPLICIT NONE
#ifdef USE_MKL
#ifdef USE_INTEL
#ifdef USE_X200
    IF (FOR_GET_HBW_AVAILABILITY() .NE. FOR_K_HBW_AVAILABLE) THEN
       BLAS_PREPARE = FOR_SET_FASTMEM_POLICY(FOR_K_FASTMEM_RETRY_WARN)
       BLAS_PREPARE = MKL_SET_MEMORY_LIMIT(MKL_MEM_MCDRAM, 0)
    END IF
#endif
#endif
#endif
    BLAS_PREPARE = VN_BLAS_PREPARE()
  END FUNCTION BLAS_PREPARE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION BLAS_SET_NUM_THREADS(NT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT
    BLAS_SET_NUM_THREADS = VN_BLAS_SET_NUM_THREADS(NT)
  END FUNCTION BLAS_SET_NUM_THREADS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE BLAS_ZROTM_RR(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    REAL(KIND=DWP), INTENT(IN) :: S1
    REAL(KIND=DWP), INTENT(IN) :: S2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_RR

  PURE SUBROUTINE BLAS_ZROTM_RC(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    REAL(KIND=DWP), INTENT(IN) :: S1
    COMPLEX(KIND=DWP), INTENT(IN) :: S2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_RC

  PURE SUBROUTINE BLAS_ZROTM_CR(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    COMPLEX(KIND=DWP), INTENT(IN) :: S1
    REAL(KIND=DWP), INTENT(IN) :: S2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_CR

  PURE SUBROUTINE BLAS_ZROTM_CC(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    COMPLEX(KIND=DWP), INTENT(IN) :: S1
    COMPLEX(KIND=DWP), INTENT(IN) :: S2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_CC

#ifdef HAVE_IMAGINARY

  PURE SUBROUTINE BLAS_ZROTM_RI(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, DS2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    REAL(KIND=DWP), INTENT(IN) :: S1
    REAL(KIND=DWP), INTENT(IN) :: DS2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    TYPE(ZIMAGINARY) :: S2
    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
    S2%J = DS2
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_RI

  PURE SUBROUTINE BLAS_ZROTM_IR(SIDE, N, ZX, INCX, ZY, INCY, C1, DS1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    REAL(KIND=DWP), INTENT(IN) :: DS1
    REAL(KIND=DWP), INTENT(IN) :: S2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    TYPE(ZIMAGINARY) :: S1
    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
    S1%J = DS1
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_IR

  PURE SUBROUTINE BLAS_ZROTM_CI(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, DS2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    COMPLEX(KIND=DWP), INTENT(IN) :: S1
    REAL(KIND=DWP), INTENT(IN) :: DS2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    TYPE(ZIMAGINARY) :: S2
    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
    S2%J = DS2
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_CI

  PURE SUBROUTINE BLAS_ZROTM_IC(SIDE, N, ZX, INCX, ZY, INCY, C1, DS1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    REAL(KIND=DWP), INTENT(IN) :: DS1
    COMPLEX(KIND=DWP), INTENT(IN) :: S2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    TYPE(ZIMAGINARY) :: S1
    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
    S1%J = DS1
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_IC

  PURE SUBROUTINE BLAS_ZROTM_II(SIDE, N, ZX, INCX, ZY, INCY, C1, DS1, DS2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    REAL(KIND=DWP), INTENT(IN) :: DS1
    REAL(KIND=DWP), INTENT(IN) :: DS2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    TYPE(ZIMAGINARY) :: S1
    TYPE(ZIMAGINARY) :: S2
    LOGICAL :: RIGHT
    COMPLEX(KIND=DWP) :: W, Z
    INTEGER :: I, KX, KY
    S1%J = DS1
    S2%J = DS2
#include "blas_zrotm_xx.F90"
  END SUBROUTINE BLAS_ZROTM_II

#endif

  PURE SUBROUTINE BLAS_ZROTM(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, S2, C2, INFO)
    IMPLICIT NONE

    CHARACTER, INTENT(IN) :: SIDE
    INTEGER, INTENT(IN) :: INCX, INCY, N
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    COMPLEX(KIND=DWP), INTENT(IN) :: S1, S2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)
    INTEGER, INTENT(OUT) :: INFO

    IF (AIMAG(S1) .EQ. D_ZERO) THEN
       ! Rx
       IF (AIMAG(S2) .EQ. D_ZERO) THEN
          ! RR
          CALL BLAS_ZROTM_RR(SIDE, N, ZX, INCX, ZY, INCY, C1, REAL(S1), REAL(S2), C2, INFO)
#ifdef HAVE_IMAGINARY
       ELSE IF (REAL(S2) .EQ. D_ZERO) THEN
          ! RI
          CALL BLAS_ZROTM_RI(SIDE, N, ZX, INCX, ZY, INCY, C1, REAL(S1), AIMAG(S2), C2, INFO)
#endif
       ELSE
          ! RC
          CALL BLAS_ZROTM_RC(SIDE, N, ZX, INCX, ZY, INCY, C1, REAL(S1), S2, C2, INFO)
       END IF
#ifdef HAVE_IMAGINARY
    ELSE IF (REAL(S1) .EQ. D_ZERO) THEN
       ! Ix
       IF (AIMAG(S2) .EQ. D_ZERO) THEN
          ! IR
          CALL BLAS_ZROTM_IR(SIDE, N, ZX, INCX, ZY, INCY, C1, AIMAG(S1), REAL(S2), C2, INFO)
       ELSE IF (REAL(S2) .EQ. D_ZERO) THEN
          ! II
          CALL BLAS_ZROTM_II(SIDE, N, ZX, INCX, ZY, INCY, C1, AIMAG(S1), AIMAG(S2), C2, INFO)
       ELSE
          ! IC
          CALL BLAS_ZROTM_IC(SIDE, N, ZX, INCX, ZY, INCY, C1, AIMAG(S1), S2, C2, INFO)
       END IF
#endif
    ELSE
       ! Cx
       IF (AIMAG(S2) .EQ. D_ZERO) THEN
          ! CR
          CALL BLAS_ZROTM_CR(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, REAL(S2), C2, INFO)
#ifdef HAVE_IMAGINARY
       ELSE IF (REAL(S2) .EQ. D_ZERO) THEN
          ! CI
          CALL BLAS_ZROTM_CI(SIDE, N, ZX, INCX, ZY, INCY, C1, S1, AIMAG(S2), C2, INFO)
#endif
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
    REAL(KIND=DWP), INTENT(IN) :: C1, C2
    COMPLEX(KIND=DWP), INTENT(IN) :: S1, S2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)

    COMPLEX(KIND=DWP) :: W(DSIMDL), Z(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: W,Z

    REAL(KIND=DWP) :: ReW(DSIMDL),ImW(DSIMDL), ReZ(DSIMDL),ImZ(DSIMDL), Re1(DSIMDL),Im1(DSIMDL), Re2(DSIMDL),Im2(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: ReW,ImW, ReZ,ImZ, Re1,Im1, Re2,Im2
    REAL(KIND=DWP) :: aC1,aC2, ReS1,ImS1, ReS2,ImS2
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: aC1,aC2, ReS1,ImS1, ReS2,ImS2
    INTEGER :: I, J, K

    !DIR$ ASSUME_ALIGNED ZX:ALIGNB
    !DIR$ ASSUME_ALIGNED ZY:ALIGNB

    IF (N .LE. 0) RETURN

    aC1 = C1
    aC2 = C2
    ReS1 = REAL(S1)
    ImS1 = AIMAG(S1)
    ReS2 = REAL(S2)
    ImS2 = AIMAG(S2)

    DO I = 1, N, DSIMDL
       K = MIN(DSIMDL, N-(I-1))
       !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
       DO J = 1, K
          W(J) = ZX(I+(J-1))
          Z(J) = ZY(I+(J-1))
          ReW(J) = REAL(W(J))
          ImW(J) = AIMAG(W(J))
          ReZ(J) = REAL(Z(J))
          ImZ(J) = AIMAG(Z(J))
          Re1(J) = ReW(J) * aC1 + ReZ(J) * ReS1 - ImZ(J) * ImS1
          Im1(J) = ImW(J) * aC1 + ImZ(J) * ReS1 + ReZ(J) * ImS1
          Re2(J) = ReW(J) * ReS2 - ImW(J) * ImS2 + ReZ(J) * aC2
          Im2(J) = ReW(J) * ImS2 + ImW(J) * ReS2 + ImZ(J) * aC2
          ZX(I+(J-1)) = CMPLX(Re1(J), Im1(J), DWP)
          ZY(I+(J-1)) = CMPLX(Re2(J), Im2(J), DWP)
       END DO
    END DO
  END SUBROUTINE BLAS_ZVROTM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE BLAS_ZSWAPC(N, ZX, INCX, ZY, INCY)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: INCX, INCY, N
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZX(*), ZY(*)

    COMPLEX(KIND=DWP) :: ZTEMP
    INTEGER :: I, IX, IY

    IF (N .LE. 0) RETURN
    IF ((INCX .EQ. 1) .AND. (INCY .EQ. 1)) THEN
       DO I = 1, N
          ZTEMP = CONJG(ZX(I))
          ZX(I) = CONJG(ZY(I))
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
          ZTEMP = CONJG(ZX(IX))
          ZX(IX) = CONJG(ZY(IY))
          ZY(IY) = ZTEMP
          IX = IX + INCX
          IY = IY + INCY
       END DO
    END IF
  END SUBROUTINE BLAS_ZSWAPC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE INTEGER FUNCTION BLAS_IZAMAX(N, ZX, INCX)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: INCX, N
    COMPLEX(KIND=DWP), INTENT(IN) :: ZX(*)

    REAL(KIND=DWP) :: DMAX, DTMP
    INTEGER :: I, IX

    BLAS_IZAMAX = 0
    IF ((N .LT. 1) .OR. (INCX .LE. 0)) RETURN
    BLAS_IZAMAX = 1
    IF (N .EQ. 1) RETURN

    IF (INCX .EQ. 1) THEN
       DMAX = ABS(ZX(1))
       DO I = 2, N
          DTMP = ABS(ZX(I))
          IF (DTMP .GT. DMAX) THEN
             BLAS_IZAMAX = I
             DMAX = DTMP
          END IF
       END DO
    ELSE
       IX = 1
       DMAX = ABS(ZX(1))
       DO I = 2, N
          IX = IX + INCX
          DTMP = ABS(ZX(IX))
          IF (DTMP .GT. DMAX) THEN
             BLAS_IZAMAX = I
             DMAX = DTMP
          END IF
       END DO
    END IF
  END FUNCTION BLAS_IZAMAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE BLAS_ZVSCAL(M, DX, ZY)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M
    REAL(KIND=DWP), INTENT(IN) :: DX(M)
    COMPLEX(KIND=DWP), INTENT(INOUT) :: ZY(M)

    INTEGER :: I
    !DIR$ ASSUME_ALIGNED DX:ALIGNB

    !DIR$ VECTOR ALWAYS ASSERT
    DO I = 1, M
       ZY(I) = DX(I) * ZY(I)
    END DO
  END SUBROUTINE BLAS_ZVSCAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE BLAS_UTILS
