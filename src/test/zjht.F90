PROGRAM ZJHT
  USE JQR
  IMPLICIT NONE

  INTEGER, PARAMETER :: M = 4, N = 3
  INTEGER :: JJ(M), OJ(M), P(N), ROW(N), INFO, I, J, K, J1, J2, J3, J4
  REAL(KIND=DWP) :: FCT(N), WORK(N), MM(N,N)
  COMPLEX(KIND=DWP) :: A(M,N), AA(M,N), AH(N,M), T(M,N), MM1(N,N), MM2(N,N)

  K = 1

  DO J1 = -1, 1, 2
     DO J2 = -1, 1, 2
        DO J3 = -1, 1, 2
           DO J4 = -1, 1, 2
              JJ(1) = J1
              JJ(2) = J2
              JJ(3) = J3
              JJ(4) = J4

              A(1,1) = CMPLX(D_ONE, D_ONE, DWP)
              A(2,1) = CMPLX(D_ONE, D_MONE, DWP)
              A(3,1) = CMPLX(D_TWO, D_ONE, DWP)
              A(4,1) = CMPLX(D_MTWO, D_ONE, DWP)

              A(1,2) = CMPLX(D_MONE, D_ONE, DWP)
              A(2,2) = CMPLX(D_MONE, D_MONE, DWP)
              A(3,2) = CMPLX(D_MTWO, D_ONE, DWP)
              A(4,2) = CMPLX(D_TWO, D_MONE, DWP)

              A(1,3) = CMPLX(D_ZERO, D_ONE, DWP)
              A(2,3) = CMPLX(D_ONE, D_ZERO, DWP)
              A(3,3) = CMPLX(D_ZERO, D_MTWO, DWP)
              A(4,3) = CMPLX(D_MTWO, D_MONE, DWP)
              AA = A
              OJ = JJ
              ! CALL WRITE_MTX_3x3(MM1, 'A^H J A')
              T = Z_ZERO
              AH = Z_ZERO
              CALL ZJQR(M, N, A, M, JJ, T, M, 1, P, FCT, ROW, WORK, INFO)
              IF (INFO .LT. 0) THEN
                 WRITE (UOUT,'(I2,A,I2)') K, ': ', INFO
              ELSE
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
                 MM2 = MATMUL(AH,T)
                 DO J = 1, N
                    DO I = 1, N
                       MM(I,J) = ABS(MM1(I,J) - MM2(I,J))
                    END DO
                 END DO
                 WRITE (UOUT,'(2(I2,A),ES25.17E3)') K, ',', INFO, ': ||A^H J A - R^H J R||_F =', NORM2(MM)
              END IF
              K = K + 1
           END DO
        END DO
     END DO
  END DO

! CONTAINS

!   SUBROUTINE WRITE_MTX_3x3(A, L)
!     IMPLICIT NONE
!     COMPLEX(KIND=DWP), INTENT(IN) :: A(3,3)
!     CHARACTER(LEN=*), INTENT(IN) :: L

!     WRITE (UOUT,'(A)') TRIM(L)
!     WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
!          '(',REAL(A(1,1)),',',AIMAG(A(1,1)),') ', &
!          '(',REAL(A(1,2)),',',AIMAG(A(1,2)),') ', &
!          '(',REAL(A(1,3)),',',AIMAG(A(1,3)),')'
!     WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
!          '(',REAL(A(2,1)),',',AIMAG(A(2,1)),') ', &
!          '(',REAL(A(2,2)),',',AIMAG(A(2,2)),') ', &
!          '(',REAL(A(2,3)),',',AIMAG(A(2,3)),')'
!     WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
!          '(',REAL(A(3,1)),',',AIMAG(A(3,1)),') ', &
!          '(',REAL(A(3,2)),',',AIMAG(A(3,2)),') ', &
!          '(',REAL(A(3,3)),',',AIMAG(A(3,3)),')'
!   END SUBROUTINE WRITE_MTX_3x3

END PROGRAM ZJHT
