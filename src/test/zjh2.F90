PROGRAM ZJHT
  USE JQR
  IMPLICIT NONE

  INTEGER, PARAMETER :: M = 4, N = 3
  INTEGER :: JJ(M), OJ(M), P(N), ROW(N), INFO, I, J
  REAL(KIND=DWP) :: FCT(N), WORK(N), MM(N,N)
  COMPLEX(KIND=DWP) :: A(M,N), AA(M,N), AH(N,M), T(M,N), MM1(N,N), MM2(N,N)

  REAL(KIND=DWP), EXTERNAL :: DNRM2

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
     WRITE (UOUT,'(4(A,I2))') 'OJ: ', OJ(1), ',', OJ(2), ',', OJ(3), ',', OJ(4)
     WRITE (UOUT,'(4(A,I2))') 'JJ: ', JJ(1), ',', JJ(2), ',', JJ(3), ';', JJ(4)
     WRITE (UOUT,'(3(A,I2))') ' P: ', P(1), ',', P(2), ',', P(3)
     WRITE (UOUT,'(3(A,I2))') ' O: ', ROW(1), ',', ROW(2), ',', ROW(3)
     CALL WRITE_MTX_3x3(M, A, 'R')

     DO J = 1, N
        DO I = 1, M
           T(I,J) = OJ(I) * AA(I,P(J))
           AH(J,I) = CONJG(AA(I,P(J)))
        END DO
     END DO
     MM1 = MATMUL(AH, T)

     DO J = 1, N
        DO I = 1, M
           T(I,J) = JJ(I) * A(I,J)
           AH(J,I) = CONJG(A(I,J))
        END DO
     END DO
     MM2 = MATMUL(AH, T)

     DO J = 1, N
        DO I = 1, N
           MM(I,J) = ABS(MM1(I,J) - MM2(I,J))
        END DO
     END DO
     WRITE (UOUT,'(A,I2,A,ES25.17E3)') 'INFO=', INFO, '; ||A^H J A - R^H J R||_F =', DNRM2(N*N, MM, 1)
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
