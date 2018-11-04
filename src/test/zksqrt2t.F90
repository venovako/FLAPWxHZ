PROGRAM ZKSQRT2T
  USE JQR
  IMPLICIT NONE

  DOUBLE PRECISION :: A, D, ReB, ImB
  INTEGER :: I, J, K, INFO
  DOUBLE COMPLEX :: B, AA, BB, CC, DD, K2(2,2), SQ(2,2)

  WRITE (*,'(A)',ADVANCE='NO') 'A  = '
  READ (*,*) A
  WRITE (*,'(A)',ADVANCE='NO') 'ReB= '
  READ (*,*) ReB
  WRITE (*,'(A)',ADVANCE='NO') 'ImB= '
  READ (*,*) ImB
  WRITE (*,'(A)',ADVANCE='NO') 'D  = '
  READ (*,*) D
  WRITE (*,'(A)',ADVANCE='NO') 'J  = '
  READ (*,*) J

  B = DCMPLX(ReB, ImB)
  K2(1,1) = DCMPLX(A * J, D_ZERO)
  K2(1,2) = B * (-J)
  K2(2,1) = DCONJG(B) * J
  K2(2,2) = DCMPLX(D * (-J), D_ZERO)
  CALL WRITE_MTX_2x2(K2, 'K**2')

  CALL ZKSQRT2(A, B, D, J, AA, BB, CC, DD, INFO)
  SQ(1,1) = AA
  SQ(1,2) = BB
  SQ(2,1) = CC
  SQ(2,2) = DD
  WRITE (*,'(I2)',ADVANCE='NO') INFO
  CALL WRITE_MTX_2x2(SQ, ' K')

  SQ = K2 - MATMUL(SQ, SQ)
  CALL WRITE_MTX_2x2(SQ, 'K2 - SQ*SQ')

CONTAINS

  SUBROUTINE WRITE_MTX_2x2(A, L)
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(IN) :: A(2,2)
    CHARACTER(LEN=*), INTENT(IN) :: L

    WRITE (*,'(A)') TRIM(L)
    WRITE (*,'(2(A,ES25.17E3,A,ES25.17E3,A))') '(',DBLE(A(1,1)),',',AIMAG(A(1,1)),') ', '(',DBLE(A(1,2)),',',AIMAG(A(1,2)),')'
    WRITE (*,'(2(A,ES25.17E3,A,ES25.17E3,A))') '(',DBLE(A(2,1)),',',AIMAG(A(2,1)),') ', '(',DBLE(A(2,2)),',',AIMAG(A(2,2)),')'
  END SUBROUTINE WRITE_MTX_2x2

END PROGRAM ZKSQRT2T
