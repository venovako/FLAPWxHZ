MODULE TIMER
  USE VN_TIMER_F
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: DUS2S = 1D-6
  DOUBLE PRECISION, PARAMETER :: DNS2S = 1D-9

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION GET_THREAD_NS()
    IMPLICIT NONE
    GET_THREAD_NS = INT(VN_GET_THREAD_NS())
  END FUNCTION GET_THREAD_NS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION GET_SYS_US()
    IMPLICIT NONE
    GET_SYS_US = INT(VN_GET_SYS_US())
  END FUNCTION GET_SYS_US

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE DOUBLE PRECISION FUNCTION TIMER2DBLE(CLK)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: CLK(3)

    TIMER2DBLE = (CLK(2) - CLK(1)) / DBLE(CLK(3))
  END FUNCTION TIMER2DBLE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TIMER_START(CLK)
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: CLK(3)

    CALL SYSTEM_CLOCK(CLK(1), CLK(3))
  END SUBROUTINE TIMER_START

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TIMER_STOP(CLK)
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: CLK(3)

    CALL SYSTEM_CLOCK(CLK(2))
  END SUBROUTINE TIMER_STOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TIMER_PRINT(CLK)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: CLK(3)

    INTEGER :: C, Q, R

    C = CLK(2) - CLK(1)
    Q = C / CLK(3)
    R = MOD(C, CLK(3))

    SELECT CASE (CLK(3))
    CASE (1000)       ! ms
       WRITE (*,1) Q, R
    CASE (1000000)    ! us
       WRITE (*,2) Q, R
    CASE (1000000000) ! ns
       WRITE (*,3) Q, R
    CASE DEFAULT      ! other scale
       WRITE (*,4) Q, R, CLK(3)
    END SELECT

1   FORMAT(I5,'.',I3.3,' s')
2   FORMAT(I5,'.',I6.6,' s')
3   FORMAT(I5,'.',I9.9,' s')
4   FORMAT(I5,'+',I12,'/',I12,' s')
  END SUBROUTINE TIMER_PRINT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE TIMER
