PROGRAM JSTRAT_TEST
  USE JSTRAT_F
  IMPLICIT NONE

  INTEGER :: JS(8)
  INTEGER :: ID, N, I, J, INFO
  INTEGER :: ARR(2,2,128), COMM(2,128), PAIR(2,128)

  WRITE (*,*) 'JS='
  READ (*,*) ID
  WRITE (*,*) 'N='
  READ (*,*) N

  CALL JSTRAT_INIT(JS, ID, N, INFO)
  WRITE (*,*) 'JSTRAT_INIT = ', INFO

  DO I = 1, N
     CALL JSTRAT_NEXT_WC(JS, ARR, INFO)
     WRITE (*,*) 'JSTRAT_NEXT_WC(', I, ')=', INFO
     CALL JSTRAT_UNPACK_WC(JS, ARR, PAIR, COMM, INFO)
     WRITE (*,*) 'JSTRAT_UNPACK_WC(', I, ')=', INFO
     DO J = 1, N/2
        IF ((PAIR(1,J) .LT. 1) .OR. (PAIR(1,J) .GT. N)) THEN
           WRITE (*,*) J
           STOP 'A'
        END IF
        IF ((PAIR(2,J) .LT. 1) .OR. (PAIR(2,J) .GT. N)) THEN
           WRITE (*,*) J, PAIR(2,J)
           STOP 'B'
        END IF
        IF ((ABS(COMM(1,J)) .LT. 1) .OR. (ABS(COMM(1,J)) .GT. N/2)) THEN
           WRITE (*,*) J, COMM(1,J)
           STOP 'C'
        END IF
        IF ((ABS(COMM(2,J)) .LT. 1) .OR. (ABS(COMM(2,J)) .GT. N/2)) THEN
           WRITE (*,*) J, COMM(2,J)
           STOP 'D'
        END IF
     END DO
  END DO
  CALL JSTRAT_FREE(JS)
END PROGRAM JSTRAT_TEST
