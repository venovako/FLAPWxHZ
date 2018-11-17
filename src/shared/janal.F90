SUBROUTINE JANAL(M, J, NPLUS, JNBLKS, JNSTIX,JNLENS)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: M, J(M)
  INTEGER, INTENT(OUT) :: NPLUS, JNBLKS, JNSTIX(*),JNLENS(*)

  INTEGER :: I

#ifndef NDEBUG
  IF (M .LT. 0) STOP 'JANAL: M < 0'
#endif

  ! Analyse the non-unity blocks of J.
  NPLUS = 0
  JNBLKS = 0

  I = 1
  DO WHILE (I .LE. M)
     IF (J(I) .EQ. 1) THEN
        NPLUS = NPLUS + 1
        I = I + 1
     ELSE IF (J(I) .EQ. -1) THEN
        JNBLKS = JNBLKS + 1
        JNSTIX(JNBLKS) = I
        JNLENS(JNBLKS) = 1
        I = I + 1
        DO WHILE (I .LE. M)
           IF (J(I) .EQ. -1) THEN
              JNLENS(JNBLKS) = JNLENS(JNBLKS) + 1
              I = I + 1
           ELSE
              EXIT
           END IF
        END DO
     ELSE ! |J(I)| .NE. 1
        STOP 'JANAL: J not supported'
     END IF
  END DO
END SUBROUTINE JANAL
