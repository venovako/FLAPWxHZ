SUBROUTINE BIO_READ_ALL(FN, M, N, YY, LDY, WW, LDW, JJ, INFO)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: FN
  INTEGER, INTENT(IN) :: M, N, LDY, LDW
  DOUBLE COMPLEX, INTENT(OUT), TARGET :: YY(LDY,N), WW(LDW,N)
  INTEGER, INTENT(OUT), TARGET :: JJ(M)
  INTEGER, INTENT(OUT) :: INFO

  INTEGER :: SZ(3), FD(3), I, J, TN, NT

  SZ = -1
  FD = -1

  IF (LEN_TRIM(FN) .LE. 0) THEN
     INFO = -1
  ELSE IF (M .LT. 0) THEN
     INFO = -2
  ELSE IF (N .LT. 0) THEN
     INFO = -3
  ELSE IF (N .GT. M) THEN
     INFO = -3
  ELSE IF (LDY .LT. M) THEN
     INFO = -5
  ELSE IF (LDW .LT. M) THEN
     INFO = -7
  ELSE
     INFO = 0
  END IF
  IF (INFO .NE. 0) RETURN

  CALL BOPEN_RO((TRIM(FN)//c_char_'.YY'), SZ(1), FD(1))
  IF (FD(1) .LT. 0) THEN
     INFO = 1
     GOTO 1
  END IF

  CALL BOPEN_RO((TRIM(FN)//c_char_'.WW'), SZ(2), FD(2))
  IF (FD(2) .LT. 0) THEN
     INFO = 2
     GOTO 1
  END IF

  CALL BOPEN_RO((TRIM(FN)//c_char_'.JJ'), SZ(3), FD(3))
  IF (FD(3) .LT. 0) THEN
     INFO = 3
     GOTO 1
  END IF

  SZ(1) = M * C_SIZEOF(Z_ZERO)
  SZ(2) = M * C_SIZEOF(Z_ZERO)
  SZ(3) =     C_SIZEOF(0)

  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(I,J,TN,NT) SHARED(YY,WW,JJ,M,N,FD,SZ) REDUCTION(MAX:INFO)
  INFO = 0
  !$OMP DO
  DO J = 1, N
     I = BREAD(FD(1), C_LOC(YY(1,J)), SZ(1), (J-1) * SZ(1))
     IF (I .NE. SZ(1)) INFO = MAX(INFO,1)
     I = BREAD(FD(2), C_LOC(WW(1,J)), SZ(2), (J-1) * SZ(2))
     IF (I .NE. SZ(2)) INFO = MAX(INFO,2)
  END DO
  !$OMP END DO
  TN = OMP_GET_THREAD_NUM()
  NT = OMP_GET_NUM_THREADS()
  IF (TN .EQ. 0) THEN
     I = (M / NT) + MOD(M,NT)
     J = 0
  ELSE ! TN .GT. 0
     I = M / NT
     J = TN * I + MOD(M,NT)
  END IF
  IF (I .GT. 0) THEN
     I = I * SZ(3)
     J = BREAD(FD(3), C_LOC(JJ(J+1)), I, J * SZ(3))
     IF (J .NE. I) INFO = MAX(INFO,3)
  END IF
  !$OMP END PARALLEL

  SZ(1) = SZ(1) * N
  SZ(2) = SZ(2) * N
  SZ(3) = SZ(3) * M

1 CALL BCLOSEN(FD, 3)
END SUBROUTINE BIO_READ_ALL

SUBROUTINE BIO_WRITE_ALL(FN, N, Y, LDY, W, LDW, J8, P, ROW, INFO)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: FN
  INTEGER, INTENT(IN) :: N, LDY, LDW
  INTEGER, INTENT(IN), TARGET :: J8(N), P(N), ROW(N)
  DOUBLE COMPLEX, INTENT(IN), TARGET :: Y(LDY,N), W(LDW,N)
  INTEGER, INTENT(OUT) :: INFO

  INTEGER :: SZ(5), FD(5), I, J, TN, NT

  SZ = -1
  FD = -1

  IF (LEN_TRIM(FN) .LE. 0) THEN
     INFO = -1
  ELSE IF (N .LT. 0) THEN
     INFO = -2
  ELSE IF (LDY .LT. N) THEN
     INFO = -4
  ELSE IF (LDW .LT. N) THEN
     INFO = -6
  ELSE
     INFO = 0
  END IF
  IF (INFO .NE. 0) RETURN

  SZ(1) = N * N * C_SIZEOF(Z_ZERO)
  SZ(2) = N * N * C_SIZEOF(Z_ZERO)
  SZ(3) =     N * C_SIZEOF(0)
  SZ(4) =     N * C_SIZEOF(0)
  SZ(5) =     N * C_SIZEOF(0)

  CALL BOPEN_RW((TRIM(FN)//c_char_'.Y'), SZ(1), FD(1))
  IF (FD(1) .LT. 0) THEN
     INFO = 1
     GOTO 2
  END IF

  CALL BOPEN_RW((TRIM(FN)//c_char_'.W'), SZ(2), FD(2))
  IF (FD(2) .LT. 0) THEN
     INFO = 2
     GOTO 2
  END IF

  CALL BOPEN_RW((TRIM(FN)//c_char_'.J'), SZ(3), FD(3))
  IF (FD(3) .LT. 0) THEN
     INFO = 3
     GOTO 2
  END IF

  CALL BOPEN_RW((TRIM(FN)//c_char_'.P'), SZ(4), FD(4))
  IF (FD(4) .LT. 0) THEN
     INFO = 4
     GOTO 2
  END IF

  CALL BOPEN_RW((TRIM(FN)//c_char_'.O'), SZ(5), FD(5))
  IF (FD(5) .LT. 0) THEN
     INFO = 5
     GOTO 2
  END IF

  SZ(1) = N * C_SIZEOF(Z_ZERO) ! Y
  SZ(2) = N * C_SIZEOF(Z_ZERO) ! W
  SZ(3) =     C_SIZEOF(0)      ! J
  SZ(4) =     C_SIZEOF(0)      ! P
  SZ(5) =     C_SIZEOF(0)      ! O

  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(I,J,TN,NT) SHARED(Y,W,J8,P,ROW,N,FD,SZ) REDUCTION(MAX:INFO)
  INFO = 0
  !$OMP DO
  DO J = 1, N
     I = BWRITE(FD(1), C_LOC(Y(1,J)), SZ(1), (J-1) * SZ(1))
     IF (I .NE. SZ(1)) INFO = MAX(INFO,1)
     I = BWRITE(FD(2), C_LOC(W(1,J)), SZ(2), (J-1) * SZ(2))
     IF (I .NE. SZ(2)) INFO = MAX(INFO,2)
  END DO
  !$OMP END DO
  TN = OMP_GET_THREAD_NUM()
  NT = OMP_GET_NUM_THREADS()
  IF (TN .EQ. 0) THEN
     I = (N / NT) + MOD(N,NT)
     J = 0
  ELSE ! TN .GT. 0
     I = N / NT
     J = TN * I + MOD(N,NT)
  END IF
  IF (I .GT. 0) THEN       
     IF (BWRITE(FD(3), C_LOC(J8(J+1)), I * SZ(3), J * SZ(3)) .NE. (I * SZ(3))) INFO = MAX(INFO,3)
     IF (BWRITE(FD(4), C_LOC(P(J+1)), I * SZ(4), J * SZ(4)) .NE. (I * SZ(4))) INFO = MAX(INFO,4)
     IF (BWRITE(FD(5), C_LOC(ROW(J+1)), I * SZ(5), J * SZ(5)) .NE. (I * SZ(5))) INFO = MAX(INFO,5)
  END IF
  !$OMP END PARALLEL

  SZ(1) = SZ(1) * N
  SZ(2) = SZ(2) * N
  SZ(3) = SZ(3) * N
  SZ(4) = SZ(4) * N
  SZ(5) = SZ(5) * N

2 CALL BCLOSEN(FD, 5)
END SUBROUTINE BIO_WRITE_ALL
