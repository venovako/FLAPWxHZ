  SUBROUTINE BOPEN_YWJ_RO(FN, M, N, SZ, FD, INFO)
    IMPLICIT NONE
    CHARACTER(LEN=*,KIND=c_char), INTENT(IN) :: FN
    INTEGER, INTENT(IN) :: M, N
    INTEGER, INTENT(OUT) :: SZ(3), FD(3), INFO

    INTEGER :: EXPTSZ(3), DIFFSZ(3)

    SZ = -1
    FD = -1

    INFO = 0
    IF (N .LT. 0) INFO = -3
    IF (M .LT. 0) INFO = -2
    IF (INFO .NE. 0) RETURN

    EXPTSZ(1) = M * N * C_SIZEOF(Z_ZERO)
    EXPTSZ(2) = M * N * C_SIZEOF(Z_ZERO)
    EXPTSZ(3) = M     * C_SIZEOF(0)

    DIFFSZ = 0

    CALL BOPEN_RO((TRIM(FN)//c_char_'.Y'), SZ(1), FD(1))
    IF (FD(1) .LT. 0) THEN
       INFO = 1
       GOTO 1
    END IF
    DIFFSZ(1) = SZ(1) - EXPTSZ(1)
    IF (DIFFSZ(1) .NE. 0) THEN
       INFO = 1
       GOTO 1
    END IF
 
    CALL BOPEN_RO((TRIM(FN)//c_char_'.W'), SZ(2), FD(2))
    IF (FD(2) .LT. 0) THEN
       INFO = 2
       GOTO 1
    END IF
    DIFFSZ(2) = SZ(2) - EXPTSZ(2)
    IF (DIFFSZ(2) .NE. 0) THEN
       INFO = 2
       GOTO 1
    END IF

    CALL BOPEN_RO((TRIM(FN)//c_char_'.J'), SZ(3), FD(3))
    IF (FD(3) .LT. 0) THEN
       INFO = 3
       GOTO 1
    END IF
    DIFFSZ(3) = SZ(3) - EXPTSZ(3)
    IF (DIFFSZ(3) .NE. 0) THEN
       INFO = 3
       GOTO 1
    END IF

    RETURN

1   SZ = DIFFSZ
  END SUBROUTINE BOPEN_YWJ_RO

  SUBROUTINE BREAD_YW(FD, Y, W, M, IFCOL, NCOLS, SZ, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FD(2), M, IFCOL, NCOLS
    COMPLEX(KIND=DWP), INTENT(OUT), TARGET :: Y(M,NCOLS)
    COMPLEX(KIND=DWP), INTENT(OUT), TARGET :: W(M,NCOLS)
    INTEGER, INTENT(OUT) :: SZ(2), INFO

    INFO = 0

    SZ(1) = BREAD(FD(1), C_LOC(Y), C_SIZEOF(Y), (IFCOL-1) * M * C_SIZEOF(Z_ZERO))
    IF (SZ(1) .NE. C_SIZEOF(Y)) INFO = 1

    SZ(2) = BREAD(FD(2), C_LOC(W), C_SIZEOF(W), (IFCOL-1) * M * C_SIZEOF(Z_ZERO))
    IF (SZ(2) .NE. C_SIZEOF(W)) INFO = 2
  END SUBROUTINE BREAD_YW

  SUBROUTINE BREAD_J(FD, J, M, SZ, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FD, M
    INTEGER, INTENT(OUT), TARGET :: J(M)
    INTEGER, INTENT(OUT) :: SZ, INFO

    INFO = 0

    SZ = BREAD(FD, C_LOC(J), C_SIZEOF(J), 0)
    IF (SZ .NE. C_SIZEOF(J)) INFO = 1
  END SUBROUTINE BREAD_J

  SUBROUTINE BOPEN_EZS_RW(FN, M, N, SZ, FD, INFO)
    IMPLICIT NONE
    CHARACTER(LEN=*,KIND=c_char), INTENT(IN) :: FN
    INTEGER, INTENT(IN) :: M, N
    INTEGER, INTENT(OUT) :: SZ(9), FD(9), INFO

    SZ = -1
    FD = -1

    INFO = 0
    IF (N .LT. 0) INFO = -3
    IF (M .LT. 0) INFO = -2
    IF (INFO .NE. 0) RETURN

    SZ(1) =     N * C_SIZEOF(D_ZERO) !  E
    SZ(2) = N * N * C_SIZEOF(Z_ZERO) !  Z

    SZ(3) =     N * C_SIZEOF(D_ZERO) ! SS
    SZ(4) = M * N * C_SIZEOF(Z_ZERO) ! YU
    SZ(5) =     N * C_SIZEOF(D_ZERO) ! EY
    SZ(6) =     N * C_SIZEOF(D_ZERO) ! SY
    SZ(7) = M * N * C_SIZEOF(Z_ZERO) ! WV
    SZ(8) =     N * C_SIZEOF(D_ZERO) ! EW
    SZ(9) =     N * C_SIZEOF(D_ZERO) ! SW

    CALL BOPEN_RW((TRIM(FN)//c_char_'.E'), SZ(1), FD(1))
    IF (FD(1) .LT. 0) THEN
       INFO = 1
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.Z'), SZ(2), FD(2))
    IF (FD(2) .LT. 0) THEN
       INFO = 2
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.SS'), SZ(3), FD(3))
    IF (FD(3) .LT. 0) THEN
       INFO = 3
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.YU'), SZ(4), FD(4))
    IF (FD(4) .LT. 0) THEN
       INFO = 4
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.EY'), SZ(5), FD(5))
    IF (FD(5) .LT. 0) THEN
       INFO = 5
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.SY'), SZ(6), FD(6))
    IF (FD(6) .LT. 0) THEN
       INFO = 6
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.WV'), SZ(7), FD(7))
    IF (FD(7) .LT. 0) THEN
       INFO = 7
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.EW'), SZ(8), FD(8))
    IF (FD(8) .LT. 0) THEN
       INFO = 8
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.SW'), SZ(9), FD(9))
    IF (FD(9) .LT. 0) THEN
       INFO = 9
       RETURN
    END IF
  END SUBROUTINE BOPEN_EZS_RW

  SUBROUTINE BWRITE_EZS(FD, E, Z, SS, YU, EY, SY, WV, EW, SW, M, N, IFCOL, NCOLS, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FD(9), M, N, IFCOL, NCOLS
    REAL(KIND=DWP), INTENT(IN), TARGET :: E(NCOLS)
    COMPLEX(KIND=DWP), INTENT(IN), TARGET :: Z(N,NCOLS)
    REAL(KIND=DWP), INTENT(IN), TARGET :: SS(NCOLS)
    REAL(KIND=DWP), INTENT(IN), TARGET :: EY(NCOLS)
    REAL(KIND=DWP), INTENT(IN), TARGET :: SY(NCOLS)
    REAL(KIND=DWP), INTENT(IN), TARGET :: EW(NCOLS)
    REAL(KIND=DWP), INTENT(IN), TARGET :: SW(NCOLS)
    COMPLEX(KIND=DWP), INTENT(IN), TARGET :: YU(M,NCOLS)
    COMPLEX(KIND=DWP), INTENT(IN), TARGET :: WV(M,NCOLS)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER(c_size_t) :: SZ(9), S

    SZ(1) =     NCOLS * C_SIZEOF(D_ZERO) !  E
    SZ(2) = N * NCOLS * C_SIZEOF(Z_ZERO) !  Z

    SZ(3) =     NCOLS * C_SIZEOF(D_ZERO) ! SS
    SZ(4) = M * NCOLS * C_SIZEOF(Z_ZERO) ! YU
    SZ(5) =     NCOLS * C_SIZEOF(D_ZERO) ! EY
    SZ(6) =     NCOLS * C_SIZEOF(D_ZERO) ! SY
    SZ(7) = M * NCOLS * C_SIZEOF(Z_ZERO) ! WV
    SZ(8) =     NCOLS * C_SIZEOF(D_ZERO) ! EW
    SZ(9) =     NCOLS * C_SIZEOF(D_ZERO) ! SW

    INFO = 0

    S = BWRITE(FD(1), C_LOC(E), C_SIZEOF(E), (IFCOL-1) * C_SIZEOF(D_ZERO))
    IF (S .NE. SZ(1)) INFO = 1
    S = BWRITE(FD(2), C_LOC(Z), C_SIZEOF(Z), (IFCOL-1) * N * C_SIZEOF(Z_ZERO))
    IF (S .NE. SZ(2)) INFO = 2

    S = BWRITE(FD(3), C_LOC(SS), C_SIZEOF(SS), (IFCOL-1) * C_SIZEOF(D_ZERO))
    IF (S .NE. SZ(3)) INFO = 3
    S = BWRITE(FD(4), C_LOC(YU), C_SIZEOF(YU), (IFCOL-1) * M * C_SIZEOF(Z_ZERO))
    IF (S .NE. SZ(4)) INFO = 4
    S = BWRITE(FD(5), C_LOC(EY), C_SIZEOF(EY), (IFCOL-1) * C_SIZEOF(D_ZERO))
    IF (S .NE. SZ(5)) INFO = 5
    S = BWRITE(FD(6), C_LOC(SY), C_SIZEOF(SY), (IFCOL-1) * C_SIZEOF(D_ZERO))
    IF (S .NE. SZ(6)) INFO = 6
    S = BWRITE(FD(7), C_LOC(WV), C_SIZEOF(WV), (IFCOL-1) * M * C_SIZEOF(Z_ZERO))
    IF (S .NE. SZ(7)) INFO = 7
    S = BWRITE(FD(8), C_LOC(EW), C_SIZEOF(EW), (IFCOL-1) * C_SIZEOF(D_ZERO))
    IF (S .NE. SZ(8)) INFO = 8
    S = BWRITE(FD(9), C_LOC(SW), C_SIZEOF(SW), (IFCOL-1) * C_SIZEOF(D_ZERO))
    IF (S .NE. SZ(9)) INFO = 9
  END SUBROUTINE BWRITE_EZS
