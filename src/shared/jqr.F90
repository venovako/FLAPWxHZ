MODULE JQR
  USE BLAS_UTILS
  IMPLICIT NONE

CONTAINS

  ! [ A B ]
  ! [ C D ]

  PURE FUNCTION ZDET2(A, B, C, D)
    IMPLICIT NONE

    COMPLEX(KIND=DWP), INTENT(IN) :: A, B, C, D
    COMPLEX(KIND=DWP) :: ZDET2

    ZDET2 = A * D - B * C
  END FUNCTION ZDET2

  PURE FUNCTION DDET2(A, B, D)
    IMPLICIT NONE

    REAL(KIND=DWP), INTENT(IN) :: A, D
    COMPLEX(KIND=DWP), INTENT(IN) :: B
    REAL(KIND=DWP) :: DDET2

    DDET2 = A * D - (REAL(B) * REAL(B) + AIMAG(B) * AIMAG(B))
  END FUNCTION DDET2

  PURE SUBROUTINE ZHINV2(A, B, D, AA, BB, DD)
    IMPLICIT NONE

    REAL(KIND=DWP), INTENT(IN) :: A, D
    COMPLEX(KIND=DWP), INTENT(IN) :: B
    REAL(KIND=DWP), INTENT(OUT) :: AA, DD
    COMPLEX(KIND=DWP), INTENT(OUT) :: BB

    REAL(KIND=DWP) :: RDET2

    RDET2 = D_ONE / DDET2(A, B, D)
    AA = RDET2 * D
    DD = RDET2 * A
    BB = (-RDET2) * B
  END SUBROUTINE ZHINV2

  PURE SUBROUTINE ZKSQRT2(A, B, D, J, AA, BB, CC, DD, INFO)
    IMPLICIT NONE

    REAL(KIND=DWP), INTENT(IN) :: A, D
    COMPLEX(KIND=DWP), INTENT(IN) :: B
    INTEGER, INTENT(IN) :: J
    COMPLEX(KIND=DWP), INTENT(OUT) :: AA, BB, CC, DD
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: K2SQRTMDET, K2TRACE, DENOM

    !        [  A B ]   [ J  0 ]
    ! K**2 = [ ~B D ] * [ 0 -J ]
    AA = CMPLX(A * J, D_ZERO, DWP)
    BB = B * (-J)
    CC = CONJG(B) * J
    DD = CMPLX(D * (-J), D_ZERO, DWP)

    K2SQRTMDET = SQRT(-DDET2(A, B, D))
    K2TRACE = REAL(AA) + REAL(DD)
    DENOM = K2TRACE + D_TWO * K2SQRTMDET

    IF (DENOM .EQ. D_ZERO) THEN
       INFO = 0
    ELSE IF (DENOM .LT. D_ZERO) THEN
       INFO = -1
       DENOM = D_MONE / SQRT(-DENOM)
       AA = CMPLX(D_ZERO, (REAL(AA) + K2SQRTMDET) * DENOM, DWP)
       BB = CMPLX(-AIMAG(BB) * DENOM, REAL(BB) * DENOM, DWP)
       CC = CMPLX(-AIMAG(CC) * DENOM, REAL(CC) * DENOM, DWP)
       DD = CMPLX(D_ZERO, (REAL(DD) + K2SQRTMDET) * DENOM, DWP)
    ELSE ! DENOM .GT. D_ZERO
       INFO = 1
       DENOM = D_ONE / SQRT(DENOM)
       AA = CMPLX((REAL(AA) + K2SQRTMDET) * DENOM, D_ZERO, DWP)
       BB = BB * DENOM
       CC = CC * DENOM
       DD = CMPLX((REAL(DD) + K2SQRTMDET) * DENOM, D_ZERO, DWP)
    END IF
  END SUBROUTINE ZKSQRT2

  ! ZX^H JJ ZY
  PURE FUNCTION ZJDOT(M, ZX, ZY, JJ)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, JJ(M)
    COMPLEX(KIND=DWP), INTENT(IN) :: ZX(M), ZY(M)
    COMPLEX(KIND=DWP) :: ZJDOT
    
    COMPLEX(KIND=DWP) :: X(DSIMDL), Y(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: X, Y
    REAL(KIND=DWP) :: DR(DSIMDL), DI(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: DR, DI
    INTEGER :: I, J, K

    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
    X = Z_ZERO
    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
    Y = Z_ZERO
    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
    DR = D_ZERO
    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
    DI = D_ZERO

    DO I = 1, M, DSIMDL
       K = MIN(DSIMDL, M-(I-1))
       !DIR$ VECTOR ALWAYS ASSERT
       DO J = 1, K
          X(J) = ZX(I+(J-1))
          Y(J) = ZY(I+(J-1))
          DR(J) = DR(J) + JJ(I+(J-1)) * (REAL(X(J))*REAL(Y(J)) + AIMAG(X(J))*AIMAG(Y(J)))
          DI(J) = DI(J) + JJ(I+(J-1)) * (REAL(X(J))*AIMAG(Y(J))- AIMAG(X(J))*REAL(Y(J)))
       END DO
    END DO

    ZJDOT = CMPLX(SUM(DR), SUM(DI), DWP)
  END FUNCTION ZJDOT

  ! ZZ^H JJ ZZ
  PURE FUNCTION DJNRM2(M, ZZ, JJ)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, JJ(M)
    COMPLEX(KIND=DWP), INTENT(IN) :: ZZ(M)
    REAL(KIND=DWP) :: DJNRM2

    COMPLEX(KIND=DWP) :: Z(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: Z
    REAL(KIND=DWP) :: D(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: D
    INTEGER :: I, J, K

    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
    Z = Z_ZERO
    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
    D = D_ZERO

    DO I = 1, M, DSIMDL
       K = MIN(DSIMDL, M-(I-1))
       !DIR$ VECTOR ALWAYS ASSERT
       DO J = 1, K
          Z(J) = ZZ(I+(J-1))
          D(J) = D(J) + JJ(I+(J-1)) * (REAL(Z(J))*REAL(Z(J)) + AIMAG(Z(J))*AIMAG(Z(J)))
       END DO
    END DO

    DJNRM2 = SUM(D)
  END FUNCTION DJNRM2

  PURE SUBROUTINE ZMAXJ(M, ZZ, JJ, IDXS, VALS)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, JJ(M)
    COMPLEX(KIND=DWP), INTENT(IN) :: ZZ(M)
    INTEGER, INTENT(OUT) :: IDXS(3)
    REAL(KIND=DWP), INTENT(OUT) :: VALS(3)

    INTEGER :: I, J
    REAL(KIND=DWP) :: A

    IDXS = 0
    VALS = D_MONE

    DO I = 1, M
       A = ABS(ZZ(I))
       J = JJ(I) + 2
       IF (A .GT. VALS(J)) THEN
          IDXS(J) = I
          VALS(J) = A
       END IF
    END DO
  END SUBROUTINE ZMAXJ

  SUBROUTINE ZJH(M, N, A, LDA, JJ, CNRMJ, TPC, T, LDT, INFO)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N, LDA, TPC, LDT
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    INTEGER, INTENT(INOUT) :: JJ(M)
    REAL(KIND=DWP), INTENT(INOUT) :: CNRMJ
    COMPLEX(KIND=DWP), INTENT(OUT) :: T(LDT,N)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: IDXS(3), I, J, K
    REAL(KIND=DWP) :: VALS(3), R, AG1, FCT
    COMPLEX(KIND=DWP) :: EIA, SCL, F1

    EXTERNAL :: ZAXPY, ZSWAP

    IF (M .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (LDA .LT. M) THEN
       INFO = -4
    ELSE IF (TPC .LT. 1) THEN
       INFO = -7
    ELSE IF (LDT .LT. M) THEN
       INFO = -9
    ELSE
       INFO = 0
    END IF

    IF (INFO .NE. 0) RETURN
    IF (M .EQ. 0) RETURN
    IF (N .EQ. 0) RETURN

    CALL ZMAXJ(M, A(1,1), JJ, IDXS, VALS)

    IF (CNRMJ .EQ. D_ZERO) THEN
       K = 2
    ELSE IF (CNRMJ .LT. D_ZERO) THEN
       K = 1
    ELSE ! CNRMJ .GT. D_ZERO
       K = 3
    END IF

    INFO = IDXS(K)
    IF (IDXS(K) .GT. 1) THEN
       CALL ZSWAP(N, A(1,1), LDA, A(IDXS(K),1), LDA)
       IF (JJ(IDXS(K)) .NE. JJ(1)) THEN
          I = JJ(1)
          JJ(1) = JJ(IDXS(K))
          JJ(IDXS(K)) = I
       END IF
    END IF
    IF (K .EQ. 2) THEN
       !DIR$ VECTOR ALWAYS ASSERT
       DO I = 1, M
          T(I,1) = A(I,1)
          A(I,1) = Z_ZERO
       END DO
       IF (IDXS(K) .EQ. 0) INFO = 1
       RETURN
    END IF

    R = SQRT(ABS(CNRMJ))
    AG1 = ABS(A(1,1))
    IF (AG1 .EQ. D_ZERO) THEN
       EIA = Z_ONE
    ELSE
       EIA = A(1,1) / AG1
    END IF
    F1 = R * EIA

    IF (JJ(1) .EQ. 1) THEN
       !DIR$ FMA
       FCT = CNRMJ + AG1 * R
    ELSE IF (JJ(1) .EQ. -1) THEN
       !DIR$ FMA
       FCT = CNRMJ - AG1 * R
    ELSE ! |JJ(1)| .NE. 1
       INFO = -5
       RETURN
    END IF
    FCT = D_MONE / FCT
    IF (ABS(FCT) .GT. HUGE(FCT)) THEN
       INFO = -6
       RETURN
    END IF

    T(1,1) = F1 + A(1,1)
    A(1,1) = -F1
    !DIR$ VECTOR ALWAYS ASSERT
    DO I = 2, M
       T(I,1) = A(I,1)
       A(I,1) = Z_ZERO
    END DO

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(J,K,SCL) SHARED(A,T,JJ,FCT,M,N,TPC)
    K = BLAS_SET_NUM_THREADS(TPC)
    !$OMP DO
    DO J = 2, N
       SCL = FCT * ZJDOT(M, T(1,1), A(1,J), JJ)
       IF (SCL .NE. Z_ZERO) CALL ZAXPY(M, SCL, T(1,1), 1, A(1,J), 1)
    END DO
    !$OMP END DO
    K = BLAS_SET_NUM_THREADS(K)
    !$OMP END PARALLEL

    CNRMJ = FCT
  END SUBROUTINE ZJH

  ! A => R
  ! T: J-Householder reflector generators
  ! INFO < 0: error; else, # of pairs of the pivot columns
  ! FCT: FCTs from ZJH
  ! ROW(I) = INFO(ZJH) + ROW(I)
  SUBROUTINE ZJR(M, N, A, LDA, JJ, T, LDT, TPC, P, FCT, ROW, WORK, INFO)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N, LDA, LDT, TPC
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    INTEGER, INTENT(INOUT) :: JJ(M)
    COMPLEX(KIND=DWP), INTENT(OUT) :: T(LDT,N)
    REAL(KIND=DWP), INTENT(OUT) :: FCT(N), WORK(N)
    INTEGER, INTENT(OUT) :: P(N), ROW(N), INFO

    INTEGER :: I, J, K, S
    REAL(KIND=DWP) :: V, LAM, SIG, CS1
    COMPLEX(KIND=DWP) :: Z, SN1

    INTEGER, EXTERNAL :: IDAMAX
    EXTERNAL :: ZROT, ZSWAP, ZLAEV2

    IF (M .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (M .LT. N) THEN
       INFO = -3
    ELSE IF (LDA .LT. M) THEN
       INFO = -4
    ELSE IF (LDT .LT. M) THEN
       INFO = -7
    ELSE IF (TPC .LT. 1) THEN
       INFO = -8
    ELSE
       INFO = 0
    END IF

    IF (INFO .NE. 0) RETURN
    IF (M .EQ. 0) RETURN
    IF (N .EQ. 0) RETURN

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) SHARED(N,P,ROW,FCT,WORK,T)
    DO J = 1, N
       P(J) = J
       ROW(J) = J - 1
       FCT(J) = D_ZERO
       WORK(J) = D_ZERO
       IF (ROW(J) .GE. 1) THEN
          !DIR$ VECTOR ALWAYS
          DO I = 1, ROW(J)
             T(I,J) = Z_ZERO
          END DO
       END IF
    END DO
    !$OMP END PARALLEL DO

    K = 1
    DO WHILE (K .LE. N)
       S = 1

       ! ...PIVOTING...

       ! See Algorithm C in Sect. 3.3 (p. 169) in: J. R. Bunch and L. Kaufman,
       ! Some Stable Methods for Calculating Inertia and Solving Symmetric Linear Systems.
       ! Mathematics of Computation, Vol. 31, No. 137. (Jan., 1977), pp. 163-179.

       ! diagonal pivoting

       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J) SHARED(A,JJ,FCT,M,N,K)
       DO J = K, N
          FCT(J) = DJNRM2(M-(K-1), A(K,J), JJ(K))
       END DO
       !$OMP END PARALLEL DO

       I = IDAMAX(N-(K-1), FCT(K), 1) + (K-1)
       IF (I .GT. K) THEN
          CALL ZSWAP(M, A(1,K), 1, A(1,I), 1)
          CS1 = FCT(K)
          FCT(K) = FCT(I)
          FCT(I) = CS1
          J = P(K)
          P(K) = P(I)
          P(I) = J
       END IF

       V = ABS(FCT(K))
       IF (K .GE. N) GOTO 1

       ! Bunch-Kaufman-Parlett pivoting

       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,Z) SHARED(A,T,JJ,FCT,ROW,WORK,M,N,K)
       DO J = K+1, N
          Z = ZJDOT(M-(K-1), A(K,K), A(K,J), JJ(K))
          WORK(J) = ABS(Z)
          T(K,J) = Z
       END DO
       !$OMP END PARALLEL DO

       I = IDAMAX(N-K, WORK(K+1), 1) + K
       LAM = WORK(I)
       CS1 = D_ALPHA * LAM
       IF (V .GE. CS1) GOTO 1

       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J) SHARED(A,JJ,WORK,M,N,K,I)
       DO J = K, N
          IF (J .EQ. I) THEN
             WORK(J) = D_MZERO
          ELSE
             WORK(J) = ABS(ZJDOT(M-(K-1), A(K,I), A(K,J), JJ(K)))
          END IF
       END DO
       !$OMP END PARALLEL DO

       J = IDAMAX(N-(K-1), WORK(K), 1) + (K-1)
       SIG = WORK(J)
       CS1 = CS1 * LAM
       IF ((V * SIG) .GE. CS1) GOTO 1

       IF (ABS(FCT(I)) .GE. (D_ALPHA * SIG)) THEN
          CALL ZSWAP(M, A(1,K), 1, A(1,I), 1)
          CS1 = FCT(K)
          FCT(K) = FCT(I)
          FCT(I) = CS1
          J = P(K)
          P(K) = P(I)
          P(I) = J
          GOTO 1
       ELSE ! 2x2 pivot
          S = 2
          INFO = INFO + 1
          Z = T(K,I)
          IF (I .NE. (K+1)) THEN
             CALL ZSWAP(M, A(1,K+1), 1, A(1,I), 1)
             CS1 = FCT(K+1)
             FCT(K+1) = FCT(I)
             FCT(I) = CS1
             J = P(K+1)
             P(K+1) = P(I)
             P(I) = J
             I = K + 1
          END IF
       END IF

       ! ...J-HOUSEHOLDERs...

1      IF (S .EQ. 1) THEN
          CALL ZJH(M-(K-1), N-(K-1), A(K,K), LDA, JJ(K), FCT(K), TPC, T(K,K), LDT, J)
          IF (J .LE. 0) THEN
             WRITE (ULOG,'(A,I2)') 'ZJR: ZJH=', J
             INFO = -6
             RETURN
          END IF
          ROW(K) = ROW(K) + J
       ELSE ! S .EQ. 2
          WORK(1) = FCT(K)
          WORK(2) = FCT(I)
          ! [ CS1  CONJG(SN1) ] [    A     B ] [ CS1 -CONJG(SN1) ] = [ RT1  0  ]
          ! [-SN1     CS1     ] [ CONJG(B) C ] [ SN1     CS1     ]   [  0  RT2 ]
          CALL ZLAEV2(WORK(1), Z, WORK(2), LAM, SIG, CS1, SN1)
          IF (SIG .EQ. D_ZERO) THEN
             INFO = -5
             RETURN
          END IF

          CALL ZROT(M-(K-1), A(K,K), 1, A(K,I), 1, CS1, SN1)

          CALL ZJH(M-(K-1), N-(K-1), A(K,K), LDA, JJ(K), LAM, TPC, T(K,K), LDT, J)
          IF (J .LE. 0) THEN
             WRITE (ULOG,'(A,I2)') 'ZJR: ZJH=', J
             INFO = -6
             RETURN
          END IF
          FCT(K) = LAM
          ROW(K) = ROW(K) + J

          CALL ZJH(M-(I-1), N-(I-1), A(I,I), LDA, JJ(I), SIG, TPC, T(I,I), LDT, J)
          IF (J .LE. 0) THEN
             WRITE (ULOG,'(A,I2)') 'ZJR: ZJH=', J
             INFO = -6
             RETURN
          END IF
          FCT(I) = SIG
          ROW(I) = ROW(I) + J

          CALL ZROT(2, A(K,K), 1, A(K,I), 1, CS1, -SN1)
       END IF

       K = K + S
    END DO
  END SUBROUTINE ZJR

  ! B(I,J) = A(I,P(J))
  SUBROUTINE ZCPIVCP(M, N, A, LDA, B, LDB, P, INFO)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N, LDA, LDB, P(N)
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    COMPLEX(KIND=DWP), INTENT(OUT) :: B(LDB,N)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: I, J, K

    IF (M .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (LDA .LT. M) THEN
       INFO = -4
    ELSE IF (LDB .LT. M) THEN
       INFO = -6
    ELSE
       INFO = 0
    END IF

    IF (INFO .NE. 0) RETURN
    IF (M .EQ. 0) RETURN
    IF (N .EQ. 0) RETURN

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) SHARED(M,N,A,B,P)
    DO J = 1, N
       K = P(J)
       !DIR$ VECTOR ALWAYS ASSERT
       DO I = 1, M
          B(I,J) = A(I,K)
       END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE ZCPIVCP

  SUBROUTINE ZQRF(M, N, A, LDA, INFO)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N, LDA
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    INTEGER, INTENT(OUT) :: INFO

    COMPLEX(KIND=DWP), ALLOCATABLE :: TAU(:), WORK(:)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: TAU, WORK

    INTEGER :: K, LWORK
    COMPLEX(KIND=DWP) :: WORK1(1)

    EXTERNAL :: ZGEQRF, ZLASET

    IF (M .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (N .GT. M) THEN
       INFO = -2
    ELSE IF (LDA .LT. M) THEN
       INFO = -4
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    K = MIN(M,N)
    ALLOCATE(TAU(K))

    ! A = Q*R
    WORK1 = Z_ZERO
    LWORK = -1
    CALL ZGEQRF(M, N, A, LDA, TAU, WORK1, LWORK, INFO)
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I2)') 'ZGEQRF(workspace query): ', INFO
       RETURN
    END IF
    LWORK = CEILING(REAL(WORK1(1), DWP))
    ALLOCATE(WORK(LWORK))
    CALL ZGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
    DEALLOCATE(WORK)
    DEALLOCATE(TAU)
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I2)') 'ZGEQR: ', INFO
       RETURN
    END IF

    ! set everything below the diagonal of A to 0
    CALL ZLASET('L', M-1, N, Z_ZERO, Z_ZERO, A(2,1), LDA)
  END SUBROUTINE ZQRF

  SUBROUTINE ZTSR(M, N, A, LDA, INFO)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N, LDA
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    INTEGER, INTENT(OUT) :: INFO

    COMPLEX(KIND=DWP), ALLOCATABLE :: T(:), WORK(:)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: T, WORK

    INTEGER :: K, TSIZE, LWORK
    COMPLEX(KIND=DWP) :: T5(5), WORK1(1)

    EXTERNAL :: ZGEQR, ZLASET

    IF (M .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (N .GT. M) THEN
       INFO = -2
    ELSE IF (LDA .LT. M) THEN
       INFO = -4
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    K = MIN(M,N)

    ! A = Q*R
    T5 = Z_ZERO
    TSIZE = -1
    WORK1 = Z_ZERO
    LWORK = -1
    CALL ZGEQR(M, N, A, LDA, T5, TSIZE, WORK1, LWORK, INFO)
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I2)') 'ZGEQR(workspace query): ', INFO
       RETURN
    END IF
    TSIZE = CEILING(REAL(T5(1), DWP))
    ALLOCATE(T(TSIZE))
    LWORK = CEILING(REAL(WORK1(1), DWP))
    ALLOCATE(WORK(LWORK))
    CALL ZGEQR(M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO)
    DEALLOCATE(WORK)
    DEALLOCATE(T)
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I2)') 'ZGEQR: ', INFO
       RETURN
    END IF

    ! set everything below the diagonal of A to 0
    CALL ZLASET('L', M-1, N, Z_ZERO, Z_ZERO, A(2,1), LDA)
  END SUBROUTINE ZTSR

  SUBROUTINE ZTSQR(M, N, A, LDA, Q, LDQ, INFO)
    IMPLICIT NONE

    CHARACTER, PARAMETER :: SIDE = 'L'
    CHARACTER, PARAMETER :: TRANS = 'N'

    INTEGER, INTENT(IN) :: M, N, LDA, LDQ
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    COMPLEX(KIND=DWP), INTENT(OUT) :: Q(LDQ,N)
    INTEGER, INTENT(OUT) :: INFO

    COMPLEX(KIND=DWP), ALLOCATABLE :: T(:), WORK(:)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: T, WORK

    INTEGER :: K, TSIZE, LWORK, ITMP
    COMPLEX(KIND=DWP) :: T5(5), WORK1(1)

    EXTERNAL :: ZGEQR, ZGEMQR, ZLASET

    IF (M .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (N .GT. M) THEN
       INFO = -2
    ELSE IF (LDA .LT. M) THEN
       INFO = -4
    ELSE IF (LDQ .LT. M) THEN
       INFO = -6
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    K = MIN(M,N)

    ! A = Q*R
    T5 = Z_ZERO
    TSIZE = -1
    WORK1 = Z_ZERO
    LWORK = -1
    CALL ZGEQR(M, N, A, LDA, T5, TSIZE, WORK1, LWORK, INFO)
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I2)') 'ZGEQR(workspace query): ', INFO
       RETURN
    END IF
    TSIZE = CEILING(REAL(T5(1), DWP))
    ALLOCATE(T(TSIZE))
    LWORK = CEILING(REAL(WORK1(1), DWP))
    ALLOCATE(WORK(LWORK))
    CALL ZGEQR(M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO)
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I2)') 'ZGEQR: ', INFO
       DEALLOCATE(WORK)
       DEALLOCATE(T)
       RETURN
    END IF

    ! compute Q
    ITMP = LWORK
    WORK1 = Z_ZERO
    LWORK = -1
    CALL ZGEMQR(SIDE, TRANS, M, N, K, A, LDA, T, TSIZE, Q, LDQ, WORK1, LWORK, INFO)
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I3)') 'ZGEMQR(workspace query): ', INFO
       DEALLOCATE(WORK)
       DEALLOCATE(T)
       RETURN
    END IF
    LWORK = CEILING(REAL(WORK1(1), DWP))
    IF (LWORK .GT. ITMP) THEN
       DEALLOCATE(WORK)
       ALLOCATE(WORK(LWORK))
    ELSE
       LWORK = ITMP
    END IF
    CALL ZLASET('A', M, N, Z_ZERO, Z_ONE, Q, LDQ) ! Q = I
    CALL ZGEMQR(SIDE, TRANS, M, N, K, A, LDA, T, TSIZE, Q, LDQ, WORK, LWORK, INFO)
    DEALLOCATE(WORK)
    DEALLOCATE(T)
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I3)') 'ZGEMQR: ', INFO
       RETURN
    END IF

    ! set everything below the diagonal of A to 0
    ITMP = M - 1
    CALL ZLASET('L', ITMP, N, Z_ZERO, Z_ZERO, A(2,1), LDA)
  END SUBROUTINE ZTSQR

END MODULE JQR
