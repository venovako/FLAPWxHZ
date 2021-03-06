! Prepare and factorize.
SUBROUTINE ZPREP_BLKS(M,NCB,K, Y,YU,LDY, W,LDW, JNSTIX,JNLENS,JNBLKS,NPLUS,&
     IPIV,JVEC,IPL,INVP, BNPLUS, BH,BS,BZ,LDB, INFO2)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: M,NCB,K, LDY,LDW,LDB
  INTEGER, INTENT(IN) :: JNBLKS,JNSTIX(JNBLKS),JNLENS(JNBLKS),NPLUS
  INTEGER, INTENT(OUT) :: IPIV(2*K),JVEC(2*K),IPL(K),INVP(K), BNPLUS, INFO2(2)
  COMPLEX(KIND=DWP), INTENT(IN) :: Y(LDY,NCB),W(LDW,NCB)
  COMPLEX(KIND=DWP), INTENT(OUT) :: YU(LDY,NCB), BH(LDB,K),BS(LDB,K),BZ(LDB,K)

  INTEGER :: B, I, INFO
  INTEGER :: BNRANK, BNPOS, BN2PIV

  EXTERNAL :: ZLASET

  INFO2 = 0

  ! bordering possibilities
  ! NCP == NCQ, NCP+NCQ=NCB == K
  ! NCP == NCQ+1, NCP+NCQ=NCB == K-1
  ! NCP == NCQ, NCP+NCQ=NCB == K-2
  B = K - NCB
  IF ((B .LT. 0) .OR. (B .GT. 2)) STOP 'ZHZL2: |K - NCP+NCQ| .GT. 2'

  ! Prepare BH = F^H J F.
  ! Compute into BZ...
  IF (NPLUS .LT. M) THEN
     CALL ZJMUL(M,NCB, Y,YU,LDY, JNSTIX,JNLENS,JNBLKS, BZ,LDB)
  ELSE ! J = I
     CALL ZIMUL(M,NCB, Y,LDY, BZ,LDB)
  END IF
  ! ...factorize...
  CALL HIF_ZHEBPC(NCB, BZ,LDB, BNRANK, BNPOS, BN2PIV, IPIV, JVEC, INFO)
  IF (INFO .EQ. 0) THEN
     CALL HIF_JPART(NCB, BNRANK, BZ,LDB, JVEC, BNPOS, IPL, INVP)
     BNPLUS = BNPOS
  ELSE
     INFO2(1) = INFO
     RETURN
  END IF
  ! ...and BH = BZ^H
  CALL ZCOPYH(NCB, BZ,LDB, BH,LDB)
  IF (B .EQ. 1) THEN
     DO I = 1, NCB
        BH(I,NCB+1) = Z_ZERO
     END DO
     BH(NCB+1,NCB+1) = Z_ONE
     DO I = NCB, 1, -1
        BH(NCB+1,I) = Z_ZERO
     END DO
  END IF

  ! Prepare BS = G^H G.
  ! Compute into BZ...
  CALL ZIMUL(M,NCB, W,LDW, BZ,LDB)
  ! Call the Hermitian indefinite factorization on a positive definite block
  ! to get the diagonal pivoting, and to verify that the block is pos. def.
  ! and of full column rank; if it is not, fail immediately.
  CALL HIF_ZHEBPC(NCB, BZ,LDB, BNRANK, BNPOS, BN2PIV, IPIV(K), JVEC(K), INFO)
  IF (INFO .EQ. 0) THEN
     IF ((BNPOS .NE. BNRANK) .OR. (BN2PIV .NE. 0)) THEN
        WRITE (ULOG,'(4(A,I11))') 'NC(BS)=',NCB,', NRANK(BS)=',BNRANK,', NPLUS(BS)=',BNPOS,', N2PIV(BS)=',BN2PIV
        STOP 'ZPREP_BLKS: indefinite BS'
     END IF
  ELSE
     INFO2(2) = INFO
     RETURN
  END IF
  ! ...and BS = BZ^H
  CALL ZCOPYH(NCB, BZ,LDB, BS,LDB)
  IF (B .EQ. 1) THEN
     DO I = 1, NCB
        BS(I,NCB+1) = Z_ZERO
     END DO
     BS(NCB+1,NCB+1) = Z_ONE
     DO I = NCB, 1, -1
        BS(NCB+1,I) = Z_ZERO
     END DO
  END IF

  ! Prepare BZ = I.
  IF (B .EQ. 2) THEN
     I = NCB
  ELSE
     I = K
  END IF
  CALL ZLASET('A', I,I, Z_ZERO, Z_ONE, BZ,LDB)
END SUBROUTINE ZPREP_BLKS
