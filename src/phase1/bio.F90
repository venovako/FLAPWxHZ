  SUBROUTINE BOPEN_XTU_RO(FN, L, a, G, SZ, FD, INFO)
    IMPLICIT NONE
    CHARACTER(LEN=*,KIND=c_char), INTENT(IN) :: FN
    INTEGER, INTENT(IN) :: L, a, G
    INTEGER, INTENT(OUT) :: SZ(3), FD(3), INFO

    INTEGER :: EXPTSZ(3), DIFFSZ(3)

    SZ = -1
    FD = -1

    INFO = 0
    IF (G .LT. 0) INFO = -4
    IF (a .LT. 0) INFO = -3
    IF (L .LT. 0) INFO = -2
    IF (INFO .NE. 0) RETURN

    EXPTSZ(1) = (2 * L) * (G)     * a * C_SIZEOF(Z_ZERO)
    EXPTSZ(2) = (2 * L) * (2 * L) * a * C_SIZEOF(Z_ZERO)
    EXPTSZ(3) = (L)               * a * C_SIZEOF(D_ZERO)

    DIFFSZ = 0

    CALL BOPEN_RO((TRIM(FN)//c_char_'.X'), SZ(1), FD(1))
    IF (FD(1) .LT. 0) THEN
       INFO = 1
       GOTO 1
    END IF
    DIFFSZ(1) = SZ(1) - EXPTSZ(1)
    IF (DIFFSZ(1) .NE. 0) THEN
       INFO = 1
       GOTO 1
    END IF
 
    CALL BOPEN_RO((TRIM(FN)//c_char_'.T'), SZ(2), FD(2))
    IF (FD(2) .LT. 0) THEN
       INFO = 2
       GOTO 1
    END IF
    DIFFSZ(2) = SZ(2) - EXPTSZ(2)
    IF (DIFFSZ(2) .NE. 0) THEN
       INFO = 2
       GOTO 1
    END IF

    CALL BOPEN_RO((TRIM(FN)//c_char_'.U'), SZ(3), FD(3))
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
  END SUBROUTINE BOPEN_XTU_RO

  SUBROUTINE BREAD_XTU(FD, X, T, U, L, G, ATOM, SZ, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FD(3), L, G, ATOM
    DOUBLE COMPLEX, INTENT(OUT), TARGET :: X(2*L,G)
    DOUBLE COMPLEX, INTENT(OUT), TARGET :: T(2*L,2*L)
    DOUBLE PRECISION, INTENT(OUT), TARGET :: U(L)
    INTEGER, INTENT(OUT) :: SZ(3), INFO

    INFO = 0

    SZ(1) = BREAD(FD(1), C_LOC(X), C_SIZEOF(X), ATOM * C_SIZEOF(X))
    IF (SZ(1) .NE. C_SIZEOF(X)) INFO = 1

    SZ(2) = BREAD(FD(2), C_LOC(T), C_SIZEOF(T), ATOM * C_SIZEOF(T))
    IF (SZ(2) .NE. C_SIZEOF(T)) INFO = 2

    SZ(3) = BREAD(FD(3), C_LOC(U), C_SIZEOF(U), ATOM * C_SIZEOF(U))
    IF (SZ(3) .NE. C_SIZEOF(U)) INFO = 3
  END SUBROUTINE BREAD_XTU

  SUBROUTINE BOPEN_YWJ_RW(FN, L, a, G, SZ, FD, INFO)
    IMPLICIT NONE
    CHARACTER(LEN=*,KIND=c_char), INTENT(IN) :: FN
    INTEGER, INTENT(IN) :: L, a, G
    INTEGER, INTENT(OUT) :: SZ(3), FD(3), INFO

    SZ = -1
    FD = -1

    INFO = 0
    IF (G .LT. 0) INFO = -4
    IF (a .LT. 0) INFO = -3
    IF (L .LT. 0) INFO = -2
    IF (INFO .NE. 0) RETURN

    SZ(1) = (2 * L) * (G) * a * C_SIZEOF(Z_ZERO)
    SZ(2) = (2 * L) * (G) * a * C_SIZEOF(Z_ZERO)
    SZ(3) = (2 * L)       * a * C_SIZEOF(0)

    CALL BOPEN_RW((TRIM(FN)//c_char_'.YY'), SZ(1), FD(1))
    IF (FD(1) .LT. 0) THEN
       INFO = 1
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.WW'), SZ(2), FD(2))
    IF (FD(2) .LT. 0) THEN
       INFO = 2
       RETURN
    END IF

    CALL BOPEN_RW((TRIM(FN)//c_char_'.JJ'), SZ(3), FD(3))
    IF (FD(3) .LT. 0) THEN
       INFO = 3
       RETURN
    END IF
  END SUBROUTINE BOPEN_YWJ_RW

  SUBROUTINE BWRITE_YWJ(FD, Y, W, J, L, G, ATOM, SZ, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FD(3), L, G, ATOM
    DOUBLE COMPLEX, INTENT(IN), TARGET :: Y(2*L,G)
    DOUBLE COMPLEX, INTENT(IN), TARGET :: W(2*L,G)
    INTEGER, INTENT(IN), TARGET :: J(2*L)
    INTEGER, INTENT(OUT) :: SZ(3), INFO

    INFO = 0

    SZ(1) = BWRITE(FD(1), C_LOC(Y), C_SIZEOF(Y), ATOM * C_SIZEOF(Y))
    IF (SZ(1) .NE. C_SIZEOF(Y)) INFO = 1

    SZ(2) = BWRITE(FD(2), C_LOC(W), C_SIZEOF(W), ATOM * C_SIZEOF(W))
    IF (SZ(2) .NE. C_SIZEOF(W)) INFO = 2

    SZ(3) = BWRITE(FD(3), C_LOC(J), C_SIZEOF(J), ATOM * C_SIZEOF(J))
    IF (SZ(3) .NE. C_SIZEOF(J)) INFO = 3
  END SUBROUTINE BWRITE_YWJ
