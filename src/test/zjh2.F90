PROGRAM ZJHT
  USE JQR
  IMPLICIT NONE

  INTEGER, PARAMETER :: M = 4, N = 3
  INTEGER :: JJ(M), OJ(M), P(N), ROW(N), INFO, I, J
  REAL(KIND=DWP) :: FCT(N), WORK(N), MM(N,N)
  COMPLEX(KIND=DWP) :: A(M,N), AA(M,N), AH(N,M), T(M,N), MM1(N,N), MM2(N,N)

  REAL(KIND=DWP), EXTERNAL :: DZNRM2

  JJ(1) =  1
  JJ(2) =  1
  JJ(3) = -1
  JJ(4) = -1

  A(1,1) = CMPLX( 1.0E+5_DWP, D_ZERO, DWP)
  A(2,1) = CMPLX( 1.0E-4_DWP, D_ZERO, DWP)
  A(3,1) = CMPLX( 1.0E-2_DWP, D_ZERO, DWP)
  A(4,1) = CMPLX( 1.0E-4_DWP, D_ZERO, DWP)

  A(1,2) = CMPLX( 1.0E-4_DWP, D_ZERO, DWP)
  A(2,2) = CMPLX(-1.0E-5_DWP, D_ZERO, DWP)
  A(3,2) = CMPLX( 1.0E-4_DWP, D_ZERO, DWP)
  A(4,2) = CMPLX(-1.0E-5_DWP, D_ZERO, DWP)

  A(1,3) = CMPLX( 1.0E-4_DWP, D_ZERO, DWP)
  A(2,3) = CMPLX( 1.0E-4_DWP, D_ZERO, DWP)
  A(3,3) = CMPLX( 1.0E-4_DWP, D_ZERO, DWP)
  A(4,3) = CMPLX( D_ZERO,     D_ZERO, DWP)
  AA = A
  OJ = JJ
  T = Z_ZERO
  AH = Z_ZERO
  CALL ZJR(M, N, A, M, JJ, T, M, 1, P, FCT, ROW, WORK, INFO)
  IF (INFO .LT. 0) THEN
     WRITE (UOUT,'(A,I2)') 'INFO: ', INFO
  ELSE
     DO I = 1, N
        WRITE (UOUT,'(A,I1,A,I2)') 'J(', I, ')=', JJ(I)
     END DO
     DO I = N+1, M
        WRITE (UOUT,'(A,I1,A,I2)') 'JJ(', I, ')=', JJ(I)
     END DO
     DO J = 1, N
        WRITE (UOUT,'(2(A,I1))') 'P(', J, ')=', P(J)
     END DO
     DO I = 1, N
        WRITE (UOUT,'(2(A,I1))') 'O(', I, ')=', ROW(I)
     END DO
     CALL WRITE_MTX_3x3(M, A, 'R')
  END IF
CONTAINS

  SUBROUTINE WRITE_MTX_3x3(LDA, A, L)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LDA
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,3)
    CHARACTER(LEN=*), INTENT(IN) :: L

    WRITE (UOUT,'(A)') TRIM(L)
    WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
         '(',REAL(A(1,1)),',',AIMAG(A(1,1)),') ', &
         '(',REAL(A(1,2)),',',AIMAG(A(1,2)),') ', &
         '(',REAL(A(1,3)),',',AIMAG(A(1,3)),')'
    WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
         '(',REAL(A(2,1)),',',AIMAG(A(2,1)),') ', &
         '(',REAL(A(2,2)),',',AIMAG(A(2,2)),') ', &
         '(',REAL(A(2,3)),',',AIMAG(A(2,3)),')'
    WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
         '(',REAL(A(3,1)),',',AIMAG(A(3,1)),') ', &
         '(',REAL(A(3,2)),',',AIMAG(A(3,2)),') ', &
         '(',REAL(A(3,3)),',',AIMAG(A(3,3)),')'
  END SUBROUTINE WRITE_MTX_3x3
END PROGRAM ZJHT
