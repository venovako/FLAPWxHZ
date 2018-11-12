MODULE JQR
  USE BLAS_UTILS
  IMPLICIT NONE

CONTAINS

  ! [ A B ]
  ! [ C D ]

  PURE FUNCTION ZDET2(A, B, C, D)
    IMPLICIT NONE

    DOUBLE COMPLEX, INTENT(IN) :: A, B, C, D
    DOUBLE COMPLEX :: ZDET2

    ZDET2 = A * D - B * C
  END FUNCTION ZDET2

  PURE FUNCTION DDET2(A, B, D)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: A, D
    DOUBLE COMPLEX, INTENT(IN) :: B
    DOUBLE PRECISION :: DDET2

    DDET2 = A * D - (DBLE(B) * DBLE(B) + AIMAG(B) * AIMAG(B))
  END FUNCTION DDET2

  PURE SUBROUTINE ZHINV2(A, B, D, AA, BB, DD)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: A, D
    DOUBLE COMPLEX, INTENT(IN) :: B
    DOUBLE PRECISION, INTENT(OUT) :: AA, DD
    DOUBLE COMPLEX, INTENT(OUT) :: BB

    DOUBLE PRECISION :: RDET2

    RDET2 = D_ONE / DDET2(A, B, D)
    AA = RDET2 * D
    DD = RDET2 * A
    BB = (-RDET2) * B
  END SUBROUTINE ZHINV2

  PURE SUBROUTINE ZKSQRT2(A, B, D, J, AA, BB, CC, DD, INFO)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: A, D
    DOUBLE COMPLEX, INTENT(IN) :: B
    INTEGER, INTENT(IN) :: J
    DOUBLE COMPLEX, INTENT(OUT) :: AA, BB, CC, DD
    INTEGER, INTENT(OUT) :: INFO

    DOUBLE PRECISION :: K2SQRTMDET, K2TRACE, DENOM

    !        [  A B ]   [ J  0 ]
    ! K**2 = [ ~B D ] * [ 0 -J ]
    AA = DCMPLX(A * J, D_ZERO)
    BB = B * (-J)
    CC = DCONJG(B) * J
    DD = DCMPLX(D * (-J), D_ZERO)

    K2SQRTMDET = SQRT(-DDET2(A, B, D))
    K2TRACE = DBLE(AA) + DBLE(DD)
    DENOM = K2TRACE + D_TWO * K2SQRTMDET

    IF (DENOM .EQ. D_ZERO) THEN
       INFO = 0
    ELSE IF (DENOM .LT. D_ZERO) THEN
       INFO = -1
       DENOM = D_MONE / SQRT(-DENOM)
       AA = DCMPLX(D_ZERO, (DBLE(AA) + K2SQRTMDET) * DENOM)
       BB = DCMPLX(-AIMAG(BB) * DENOM, DBLE(BB) * DENOM)
       CC = DCMPLX(-AIMAG(CC) * DENOM, DBLE(CC) * DENOM)
       DD = DCMPLX(D_ZERO, (DBLE(DD) + K2SQRTMDET) * DENOM)
    ELSE ! DENOM .GT. D_ZERO
       INFO = 1
       DENOM = D_ONE / SQRT(DENOM)
       AA = DCMPLX((DBLE(AA) + K2SQRTMDET) * DENOM, D_ZERO)
       BB = BB * DENOM
       CC = CC * DENOM
       DD = DCMPLX((DBLE(DD) + K2SQRTMDET) * DENOM, D_ZERO)
    END IF
  END SUBROUTINE ZKSQRT2

  ! ZX^H JJ ZY
  PURE FUNCTION ZJDOT(M, ZX, ZY, JJ)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, JJ(M)
    DOUBLE COMPLEX, INTENT(IN) :: ZX(M), ZY(M)
    DOUBLE COMPLEX :: ZJDOT
    
    DOUBLE COMPLEX :: X(DSIMDL), Y(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: X, Y
    DOUBLE PRECISION :: DR(DSIMDL), DI(DSIMDL)
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
          DR(J) = DR(J) + JJ(I+(J-1)) * (DBLE(X(J))*DBLE(Y(J)) + AIMAG(X(J))*AIMAG(Y(J)))
          DI(J) = DI(J) + JJ(I+(J-1)) * (DBLE(X(J))*AIMAG(Y(J))- AIMAG(X(J))*DBLE(Y(J)))
       END DO
    END DO

    ZJDOT = DCMPLX(SUM(DR), SUM(DI))
  END FUNCTION ZJDOT

  ! ZZ^H JJ ZZ
  PURE FUNCTION DJNRM2(M, ZZ, JJ)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, JJ(M)
    DOUBLE COMPLEX, INTENT(IN) :: ZZ(M)
    DOUBLE PRECISION :: DJNRM2

    DOUBLE COMPLEX :: Z(DSIMDL)
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: Z
    DOUBLE PRECISION :: D(DSIMDL)
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
          D(J) = D(J) + JJ(I+(J-1)) * (DBLE(Z(J))*DBLE(Z(J)) + AIMAG(Z(J))*AIMAG(Z(J)))
       END DO
    END DO

    DJNRM2 = SUM(D)
  END FUNCTION DJNRM2

  SUBROUTINE CNRMJS(M, N, A, LDA, JJ, FCT, INFO)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N, LDA, JJ(M)
    DOUBLE COMPLEX, INTENT(IN) :: A(LDA,N)
    DOUBLE PRECISION, INTENT(OUT) :: FCT(N)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: J

    IF (M .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (LDA .LT. M) THEN
       INFO = -4
    ELSE
       INFO = 0
    END IF

    IF (INFO .NE. 0) RETURN
    IF (M .EQ. 0) RETURN
    IF (N .EQ. 0) RETURN

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J) SHARED(A,JJ,FCT,M,N)
    DO J = 1, N
       FCT(J) = DJNRM2(M, A(1,J), JJ)
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE CNRMJS

  PURE SUBROUTINE ZMAXJ(M, ZZ, JJ, IDXS, VALS)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, JJ(M)
    DOUBLE COMPLEX, INTENT(IN) :: ZZ(M)
    INTEGER, INTENT(OUT) :: IDXS(3)
    DOUBLE PRECISION, INTENT(OUT) :: VALS(3)

    INTEGER :: I, J
    DOUBLE PRECISION :: A

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

  SUBROUTINE ZJH(M, N, A, LDA, JJ, CNRMJ, T, LDT, INFO)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N, LDA, LDT
    DOUBLE COMPLEX, INTENT(INOUT) :: A(LDA,N)
    INTEGER, INTENT(INOUT) :: JJ(M)
    DOUBLE PRECISION, INTENT(INOUT) :: CNRMJ
    DOUBLE COMPLEX, INTENT(OUT) :: T(LDT,N)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: IDXS(3), I, J, K
    DOUBLE PRECISION :: VALS(3), F, FCT
    DOUBLE COMPLEX :: SCL

    EXTERNAL :: ZAXPY, ZSWAP

    IF (M .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (LDA .LT. M) THEN
       INFO = -4
    ELSE IF (LDT .LT. M) THEN
       INFO = -8
    ELSE
       INFO = 0
    END IF

    IF (INFO .NE. 0) RETURN
    IF (M .EQ. 0) RETURN
    IF (N .EQ. 0) RETURN

    CALL ZMAXJ(M, A(1,1), JJ, IDXS, VALS)

    IF (CNRMJ .EQ. D_ZERO) THEN
       K = 2
       !DIR$ VECTOR ALWAYS ASSERT
       DO I = 1, M
          T(I,1) = Z_ZERO
       END DO
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
       IF (IDXS(K) .EQ. 0) INFO = 1
       RETURN
    END IF

    ! f**H JJ f == g**H JJ g
    F = -SIGN(SQRT(ABS(CNRMJ)), DBLE(A(1,1)))
    IF (JJ(1) .EQ. 1) THEN
       !DIR$ FMA
       FCT = CNRMJ - F * DBLE(A(1,1))
    ELSE IF (JJ(1) .EQ. -1) THEN
       !DIR$ FMA
       FCT = CNRMJ + F * DBLE(A(1,1))
    ELSE ! ABS(JJ(1)) .NE. 1
       STOP 'ZJH: JJ(1)'
    END IF
    FCT = D_MONE / FCT
    IF (ABS(FCT) .GT. HUGE(FCT)) THEN
       INFO = 0
       RETURN
    END IF

    T(1,1) = DCMPLX(F - DBLE(A(1,1)), -AIMAG(A(1,1)))
    A(1,1) = DCMPLX(F, D_ZERO)
    !DIR$ VECTOR ALWAYS ASSERT
    DO I = 2, M
       T(I,1) = -A(I,1)
       A(I,1) = Z_ZERO
    END DO

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(J,K,SCL) SHARED(A,T,JJ,M,N,FCT)
    K = BLAS_SET_NUM_THREADS(1)
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
  SUBROUTINE ZJQR(M, N, A, LDA, JJ, T, LDT, P, FCT, ROW, WORK, INFO)
    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: ALPHA = SCALE(D_ONE + SQRT(17.0D0), -3)
    
    INTEGER, INTENT(IN) :: M, N, LDA, LDT
    DOUBLE COMPLEX, INTENT(INOUT) :: A(LDA,N)
    INTEGER, INTENT(INOUT) :: JJ(M)
    DOUBLE COMPLEX, INTENT(OUT) :: T(LDT,N)
    DOUBLE PRECISION, INTENT(OUT) :: FCT(N), WORK(N)
    INTEGER, INTENT(OUT) :: P(N), ROW(N), INFO

    INTEGER :: I, J, K, S
    DOUBLE PRECISION :: V, LAM, SIG

    EXTERNAL :: ZSWAP, ZLASET

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
    ELSE
       INFO = 0
    END IF

    IF (INFO .NE. 0) RETURN
    IF (M .EQ. 0) RETURN
    IF (N .EQ. 0) RETURN

    IF (N .GT. 1) CALL ZLASET('U', N-1, N-1, Z_ZERO, Z_ZERO, T(1,2), LDT)

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I) SHARED(P,ROW,N)
    DO I = 1, N
       P(I) = I
       ROW(I) = I - 1
    END DO
    !$OMP END PARALLEL DO

    K = 1
    DO WHILE (K .LE. N)
       S = 1

       ! ...COLUMN J-NORMs...

       CALL CNRMJS(M-(K-1), N-(K-1), A(K,K), LDA, JJ(K), FCT(K), I)
       IF (I .NE. 0) THEN
          WRITE (ULOG,'(A,I2)') 'ZJQR: CNRMSJ=', I
          INFO = -5
          RETURN
       END IF

       ! ...PIVOTING...

       V = ABS(FCT(K))

       ! diagonal pivoting

       ! I = K
       ! DO J = K+1, N
       !    WORK(1) = ABS(FCT(J))
       !    IF (WORK(1) .GT. V) THEN
       !       I = J
       !       V = WORK(1)
       !    END IF
       ! END DO

       ! IF (I .GT. K) THEN
       !    CALL ZSWAP(M, A(1,K), 1, A(1,I), 1)
       !    WORK(1) = FCT(K)
       !    FCT(K) = FCT(I)
       !    FCT(I) = WORK(1)
       !    J = P(K)
       !    P(K) = P(I)
       !    P(I) = J
       ! END IF

       ! Bunch-Kaufman-Parlett pivoting

       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J) SHARED(A,JJ,WORK,M,N,K)
       DO J = K+1, N
          WORK(J) = ABS(ZJDOT(M-(K-1), A(K,K), A(K,J), JJ(K)))
       END DO
       !$OMP END PARALLEL DO

       I = K
       LAM = D_MONE
       DO J = K+1, N
          IF (WORK(J) .GT. LAM) THEN
             I = J
             LAM = WORK(J)
          END IF
       END DO
       WORK(I) = ALPHA * LAM
       IF (V .GE. WORK(I)) GOTO 1

       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J) SHARED(A,JJ,WORK,M,N,K,I)
       DO J = K, N
          IF (J .NE. I) WORK(J) = ABS(ZJDOT(M-(K-1), A(K,I), A(K,J), JJ(K)))
       END DO
       !$OMP END PARALLEL DO

       SIG = D_MONE
       DO J = K, N
          IF ((J .NE. I) .AND. (WORK(J) .GT. SIG)) SIG = WORK(J)
       END DO
       WORK(I) = WORK(I) * LAM
       IF ((V * SIG) .GE. WORK(I)) GOTO 1

       IF (ABS(FCT(I)) .GT. (ALPHA * SIG)) THEN
          CALL ZSWAP(M, A(1,K), 1, A(1,I), 1)
          WORK(1) = FCT(K)
          FCT(K) = FCT(I)
          FCT(I) = WORK(1)
          J = P(K)
          P(K) = P(I)
          P(I) = J
          GOTO 1
       ELSE IF (I .NE. (K+1)) THEN
          CALL ZSWAP(M, A(1,K+1), 1, A(1,I), 1)
          WORK(1) = FCT(K+1)
          FCT(K+1) = FCT(I)
          FCT(I) = WORK(1)
          J = P(K+1)
          P(K+1) = P(I)
          P(I) = J
       END IF

       S = 2
       INFO = INFO + 1

       ! ...J-HOUSEHOLDERs...

1      IF (S .EQ. 1) THEN
          IF (FCT(K) .EQ. D_ZERO) THEN
             I = -6
          ELSE
             CALL ZJH(M-(K-1), N-(K-1), A(K,K), LDA, JJ(K), FCT(K), T(K,K), LDT, I)
          END IF
          IF (I .LE. 0) THEN
             WRITE (ULOG,'(A,I2)') 'ZJQR: ZJH=', I
             INFO = -6
             RETURN
          END IF
          ROW(K) = ROW(K) + I
       ELSE ! S .EQ. 2
          ! ...TODO...
          CONTINUE
       END IF

       K = K + S
    END DO
  END SUBROUTINE ZJQR

  ! B(I,J) = A(I,P(J))
  SUBROUTINE ZCPIVCP(M, N, A, LDA, B, LDB, P, INFO)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: M, N, LDA, LDB, P(N)
    DOUBLE COMPLEX, INTENT(IN) :: A(LDA,N)
    DOUBLE COMPLEX, INTENT(OUT) :: B(LDB,N)
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

  SUBROUTINE ZTSR(M, N, A, LDA, INFO)
    IMPLICIT NONE

    CHARACTER, PARAMETER :: SIDE = 'L'
    CHARACTER, PARAMETER :: TRANS = 'N'

    INTEGER, INTENT(IN) :: M, N, LDA
    DOUBLE COMPLEX, INTENT(INOUT) :: A(LDA,N)
    INTEGER, INTENT(OUT) :: INFO

    DOUBLE COMPLEX, ALLOCATABLE :: T(:), WORK(:)
#ifdef USE_X200
    !DIR$ ATTRIBUTES MEMKIND:HBW, ALIGN:ALIGNB :: T, WORK
#else
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: T, WORK
#endif

    INTEGER :: K, TSIZE, LWORK
    DOUBLE COMPLEX :: T5(5), WORK1(1)

    EXTERNAL :: ZGEQR, ZGEMQR, ZLASET

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
    TSIZE = CEILING(DBLE(T5(1)))
    ALLOCATE(T(TSIZE))
    LWORK = CEILING(DBLE(WORK1(1)))
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
    DOUBLE COMPLEX, INTENT(INOUT) :: A(LDA,N)
    DOUBLE COMPLEX, INTENT(OUT) :: Q(LDQ,N)
    INTEGER, INTENT(OUT) :: INFO

    DOUBLE COMPLEX, ALLOCATABLE :: T(:), WORK(:)
#ifdef USE_X200
    !DIR$ ATTRIBUTES MEMKIND:HBW, ALIGN:ALIGNB :: T, WORK
#else
    !DIR$ ATTRIBUTES ALIGN:ALIGNB :: T, WORK
#endif

    INTEGER :: K, TSIZE, LWORK, ITMP
    DOUBLE COMPLEX :: T5(5), WORK1(1)

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
    TSIZE = CEILING(DBLE(T5(1)))
    ALLOCATE(T(TSIZE))
    LWORK = CEILING(DBLE(WORK1(1)))
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
    LWORK = CEILING(DBLE(WORK1(1)))
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
