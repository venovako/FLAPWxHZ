PROGRAM ROT_TEST
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: D_ZERO = 0.0D0, D_ONE = 1.0D0, D_TWO = 2.0D0
  DOUBLE COMPLEX, PARAMETER :: Z_ZERO = (D_ZERO,D_ZERO), Z_ONE = (D_ONE,D_ZERO)

  DOUBLE PRECISION :: A_11, A_22
  DOUBLE COMPLEX :: A_12, B_12

  DOUBLE COMPLEX :: A(2,2), B(2,2), F(2,2), G(2,2), H(2,2), W(2,2)

  DOUBLE PRECISION :: BB, U, V, T, E, TG, SG, CG, S, T2T,C2T,S2T
  DOUBLE PRECISION :: CPHI, CPSI
  DOUBLE COMPLEX :: UV, EIASPHI, EMIBSPSI

  INTEGER :: I, J, K

  WRITE (*,'(A)',ADVANCE='NO') 'A_11= '
  READ (*,*) A_11
  A(1,1) = A_11

  WRITE (*,'(A)',ADVANCE='NO') 'A_22= '
  READ (*,*) A_22
  A(2,2) = A_22

  WRITE (*,'(A)',ADVANCE='NO') 'Re(A_12)= '
  READ (*,*) U

  WRITE (*,'(A)',ADVANCE='NO') 'Im(A_12)= '
  READ (*,*) V

  A_12 = DCMPLX(U,V)
  A(1,2) = A_12
  A(2,1) = DCONJG(A_12)

  CALL WRITE_MTX_2x2(A, 'A =')

  B(1,1) = Z_ONE
  B(2,2) = Z_ONE

  WRITE (*,'(A)',ADVANCE='NO') 'Re(B_12)= '
  READ (*,*) U

  WRITE (*,'(A)',ADVANCE='NO') 'Im(B_12)= '
  READ (*,*) V

  B_12 = DCMPLX(U,V)
  B(1,2) = B_12
  B(2,1) = DCONJG(B_12)

  CALL WRITE_MTX_2x2(B, 'B =')

  ! compute the rotation F

  E = A_22 - A_11
  
  BB = ABS(B_12)
  ! T = \sqrt(1 - x^2), hope for FMA
  T = SQRT(D_ONE - BB * BB)

  ! UV = (B^*_{12} / |B_{12}|) * A_{12}
  UV = (DCONJG(B_12) / BB) * A_12
  U = DBLE(UV)
  V = AIMAG(UV)

  IF ((E .EQ. D_ZERO) .AND. (V .EQ. D_ZERO)) THEN
     CG = D_ONE / SQRT(D_TWO)
     SG = -CG
     F(1,1) = CG
     F(2,1) = (DCONJG(B_12) / BB) * -SG
     F(1,2) = (B_12 / BB) * SG
     F(2,2) = CG
  ELSE
     S = SIGN(D_ONE, E)
     TG = 2 * V / E
     CG = D_ONE / SQRT(D_ONE + TG * TG)
     IF (CG .NE. D_ZERO) THEN
        SG = TG * CG
     ELSE
        SG = SIGN(D_ONE, TG)
     END IF

     T2T = (S * (2*U - (A_11 + A_22)*BB)) / (T * SQRT(E*E + 4*V*V))
     C2T = D_ONE / SQRT(D_ONE + T2T * T2T)
     IF (C2T .NE. D_ZERO) THEN
        S2T = T2T * C2T
     ELSE
        S2T = SIGN(D_ONE, T2T)
     END IF

     CPHI = SQRT((D_ONE + BB*S2T + T*C2T * CG) / 2)
     CPSI = SQRT((D_ONE - BB*S2T + T*C2T * CG) / 2)

     EIASPHI = ((B_12 / BB) * CPSI) * &
          (DCMPLX(S2T - BB, T * SG * C2T) / (D_ONE - BB * S2T + T * CG * C2T))
     EMIBSPSI = ((DCONJG(B_12) / BB) * CPHI) * &
          (DCMPLX(S2T + BB, -T * SG * C2T) / (D_ONE + BB * S2T + T * CG * C2T))

     F(1,1) =  CPHI / T
     F(2,1) = -EMIBSPSI / T
     F(1,2) =  EIASPHI / T
     F(2,2) =  CPSI / T
  END IF

  CALL WRITE_MTX_2x2(F, 'F =')

  ! G = F^H
  G(1,1) = F(1,1)
  G(2,1) = DCONJG(F(1,2))
  G(1,2) = DCONJG(F(2,1))
  G(2,2) = F(2,2)

  H = Z_ZERO
  W = Z_ZERO

  DO I = 1, 2
     DO J = 1, 2
        DO K = 1, 2
           W(I,J) = W(I,J) + A(I,K) * F(K,J)
        END DO
     END DO
  END DO
  DO I = 1, 2
     DO J = 1, 2
        DO K = 1, 2
           H(I,J) = H(I,J) + G(I,K) * W(K,J)
        END DO
     END DO
  END DO
  CALL WRITE_MTX_2x2(H, 'F^H A F =')

  H = Z_ZERO
  W = Z_ZERO

  DO I = 1, 2
     DO J = 1, 2
        DO K = 1, 2
           W(I,J) = W(I,J) + B(I,K) * F(K,J)
        END DO
     END DO
  END DO
  DO I = 1, 2
     DO J = 1, 2
        DO K = 1, 2
           H(I,J) = H(I,J) + G(I,K) * W(K,J)
        END DO
     END DO
  END DO
  CALL WRITE_MTX_2x2(H, 'F^H B F =')

CONTAINS

  SUBROUTINE WRITE_MTX_2x2(A, L)
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(IN) :: A(2,2)
    CHARACTER(LEN=*), INTENT(IN) :: L

    WRITE (*,'(A)') TRIM(L)
    WRITE (*,'(2(A,ES25.17E3,A,ES25.17E3,A))') '(',DBLE(A(1,1)),',',AIMAG(A(1,1)),') ', '(',DBLE(A(1,2)),',',AIMAG(A(1,2)),')'
    WRITE (*,'(2(A,ES25.17E3,A,ES25.17E3,A))') '(',DBLE(A(2,1)),',',AIMAG(A(2,1)),') ', '(',DBLE(A(2,2)),',',AIMAG(A(2,2)),')'
  END SUBROUTINE WRITE_MTX_2x2

END PROGRAM ROT_TEST
