PROGRAM ZJHT
  USE JQR
  IMPLICIT NONE

  INTEGER, PARAMETER :: M = 4, N = 3
  INTEGER :: JJ(M), OJ(M), P(N), ROW(N), INFO, I, J, K, J1, J2, J3, J4
  DOUBLE PRECISION :: FCT(N), WORK(N), MM(N,N)
  DOUBLE COMPLEX :: A(M,N), AA(M,N), AH(N,M), T(M,N), MM1(N,N), MM2(N,N)

  K = 1

  DO J1 = -1, 1, 2
     DO J2 = -1, 1, 2
        DO J3 = -1, 1, 2
           DO J4 = -1, 1, 2
              JJ(1) = J1
              JJ(2) = J2
              JJ(3) = J3
              JJ(4) = J4

              A(1,1) = DCMPLX(D_ONE, D_ONE)
              A(2,1) = DCMPLX(D_ONE, D_MONE)
              A(3,1) = DCMPLX(D_TWO, D_ONE)
              A(4,1) = DCMPLX(D_MTWO, D_ONE)

              A(1,2) = DCMPLX(D_MONE, D_ONE)
              A(2,2) = DCMPLX(D_MONE, D_MONE)
              A(3,2) = DCMPLX(D_MTWO, D_ONE)
              A(4,2) = DCMPLX(D_TWO, D_MONE)

              A(1,3) = DCMPLX(D_ZERO, D_ONE)
              A(2,3) = DCMPLX(D_ONE, D_ZERO)
              A(3,3) = DCMPLX(D_ZERO, D_MTWO)
              A(4,3) = DCMPLX(D_MTWO, D_MONE)
              AA = A
              OJ = JJ
              ! CALL WRITE_MTX_3x3(MM1, 'A^H J A')
              T = Z_ZERO
              AH = Z_ZERO
              CALL ZJQR(M, N, A, M, JJ, T, M, P, FCT, ROW, WORK, INFO)
              IF (INFO .LT. 0) THEN
                 WRITE (UOUT,'(I2,A,I2)') K, ': ', INFO
              ELSE
                 DO J = 1, N
                    DO I = 1, M
                       T(I,J) = OJ(I) * AA(I,P(J))
                       AH(J,I) = DCONJG(AA(I,P(J)))
                    END DO
                 END DO
                 MM1 = MATMUL(AH, T)
                 DO J = 1, N
                    DO I = 1, M
                       T(I,J) = JJ(I) * A(I,J)
                       AH(J,I) = DCONJG(A(I,J))
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

CONTAINS

  SUBROUTINE WRITE_MTX_3x3(A, L)
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(IN) :: A(3,3)
    CHARACTER(LEN=*), INTENT(IN) :: L

    WRITE (UOUT,'(A)') TRIM(L)
    WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
         '(',DBLE(A(1,1)),',',AIMAG(A(1,1)),') ', &
         '(',DBLE(A(1,2)),',',AIMAG(A(1,2)),') ', &
         '(',DBLE(A(1,3)),',',AIMAG(A(1,3)),')'
    WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
         '(',DBLE(A(2,1)),',',AIMAG(A(2,1)),') ', &
         '(',DBLE(A(2,2)),',',AIMAG(A(2,2)),') ', &
         '(',DBLE(A(2,3)),',',AIMAG(A(2,3)),')'
    WRITE (UOUT,'(3(A,ES25.17E3,A,ES25.17E3,A))') &
         '(',DBLE(A(3,1)),',',AIMAG(A(3,1)),') ', &
         '(',DBLE(A(3,2)),',',AIMAG(A(3,2)),') ', &
         '(',DBLE(A(3,3)),',',AIMAG(A(3,3)),')'
  END SUBROUTINE WRITE_MTX_3x3

END PROGRAM ZJHT
