PROGRAM PHASE1
  USE BINIO
  USE DTYPES
  USE HIF_HZ
  USE TIMER
  USE OMP_LIB
  IMPLICIT NONE

  CHARACTER(LEN=FNL,KIND=c_char) :: FN
  INTEGER :: L, a, G
  ! Cores Per Run and Threads Per Core
  INTEGER :: CPR, TPC
  INTEGER :: INFO, ATOM
  INTEGER :: SZ(6), FD(6)
  ! timing
  INTEGER :: CLK(3)
  DOUBLE PRECISION :: GTIMES(3)

  TYPE(Z1MEM), ALLOCATABLE :: BUF(:)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: BUF

  CPR = MAX(INT(OMP_GET_MAX_THREADS()),1)
  CALL READCL(FN, L, a, G, TPC, INFO)
  IF (INFO .NE. 0) THEN
     IF (INFO .LT. 0) THEN
        WRITE (ULOG,'(A,I2)') 'Cannot read argument', -INFO
     ELSE
        WRITE (ULOG,'(A,I2)') 'Illegal value of argument', INFO
     END IF
     STOP 'phase1.exe FN L a G TPC'
#ifndef NDEBUG
  ELSE ! INFO .EQ. 0
     WRITE (ULOG,'(2(A),3(I11),I3,I2)') 'phase1.exe ', TRIM(FN), L, a, G, CPR, TPC
#endif
  END IF

  CALL BOPEN_XTU_RO(FN, L, a, G, SZ, FD, INFO)
#ifndef NDEBUG
  WRITE (ULOG,'(A,I20,A)') 'X file has ', SZ(1), ' B.'
#endif
  IF (INFO .EQ. 1) STOP 'X file cannot be opened'
#ifndef NDEBUG
  WRITE (ULOG,'(A,I20,A)') 'T file has ', SZ(2), ' B.'
#endif
  IF (INFO .EQ. 2) STOP 'T file cannot be opened'
#ifndef NDEBUG
  WRITE (ULOG,'(A,I20,A)') 'U file has ', SZ(3), ' B.'
#endif
  IF (INFO .EQ. 3) STOP 'U file cannot be opened'
  IF (INFO .NE. 0) STOP 'BOPEN_XTU_RO: error'

  CALL BOPEN_YWJ_RW(FN, L, a, G, SZ(4), FD(4), INFO)
  IF (INFO .EQ. 1) STOP 'Y file cannot be opened for writing'
  IF (INFO .EQ. 2) STOP 'W file cannot be opened for writing'
  IF (INFO .EQ. 3) STOP 'J file cannot be opened for writing'
  IF (INFO .NE. 0) STOP 'BOPEN_YWJ_RW: error'

  ATOM = 0
  INFO = 0
  CLK = 0
  GTIMES = D_MZERO

  ALLOCATE(BUF(a),STAT=INFO)
  IF (INFO .GT. 0) STOP 'Error allocating BUF'
  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ATOM) SHARED(BUF,a)
  DO ATOM = 1, a
     BUF(ATOM)%X => NULL()
     BUF(ATOM)%T => NULL()
     BUF(ATOM)%Y => NULL()
     BUF(ATOM)%U => NULL()
     BUF(ATOM)%IPIV => NULL()
     BUF(ATOM)%JVEC => NULL()
     BUF(ATOM)%IPL => NULL()
     BUF(ATOM)%INVP => NULL()
  END DO
  !$OMP END PARALLEL DO

  INFO = BLAS_PREPARE()

  WRITE (UOUT,'(A)') '"ATOM","NRANK","NPLUS","N2PIV","INFO","TIME1","TIME2","TIME3"'
  CALL TIMER_START(CLK)
  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ATOM) SHARED(FD,BUF,L,a,G,TPC) REDUCTION(MAX:INFO)
  DO ATOM = 1, a
     CALL PAR_WORK(FD, BUF(ATOM), L, G, (ATOM-1), TPC, INFO)
  END DO
  !$OMP END PARALLEL DO
  CALL TIMER_STOP(CLK)

  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ATOM) SHARED(BUF,a)
  DO ATOM = 1, a
     IF (ASSOCIATED(BUF(ATOM)%INVP)) THEN
        DEALLOCATE(BUF(ATOM)%INVP)
        BUF(ATOM)%INVP => NULL()
     END IF
     IF (ASSOCIATED(BUF(ATOM)%IPL)) THEN
        DEALLOCATE(BUF(ATOM)%IPL)
        BUF(ATOM)%IPL => NULL()
     END IF
     IF (ASSOCIATED(BUF(ATOM)%JVEC)) THEN
        DEALLOCATE(BUF(ATOM)%JVEC)
        BUF(ATOM)%JVEC => NULL()
     END IF
     IF (ASSOCIATED(BUF(ATOM)%IPIV)) THEN
        DEALLOCATE(BUF(ATOM)%IPIV)
        BUF(ATOM)%IPIV => NULL()
     END IF
     IF (ASSOCIATED(BUF(ATOM)%U)) THEN
        DEALLOCATE(BUF(ATOM)%U)
        BUF(ATOM)%U => NULL()
     END IF
     IF (ASSOCIATED(BUF(ATOM)%Y)) THEN
        DEALLOCATE(BUF(ATOM)%Y)
        BUF(ATOM)%Y => NULL()
     END IF
     IF (ASSOCIATED(BUF(ATOM)%T)) THEN
        DEALLOCATE(BUF(ATOM)%T)
        BUF(ATOM)%T => NULL()
     END IF
     IF (ASSOCIATED(BUF(ATOM)%X)) THEN
        DEALLOCATE(BUF(ATOM)%X)
        BUF(ATOM)%X => NULL()
     END IF
  END DO
  !$OMP END PARALLEL DO
  DEALLOCATE(BUF)

  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11,A,I1)') 'PAR_WORK: atom ', (INFO/4), ' error ', MOD(INFO,4)
  ELSE ! INFO .EQ. 0
     GTIMES(1) = TIMER2DBLE(CLK)
     WRITE (UOUT,'(5(I11,A),2(F11.6,A),F11.6)') &
          -a,',', L,',', G,',', CPR,',', TPC,',', GTIMES(1),',', GTIMES(2),',', GTIMES(3)
  END IF

  CALL BCLOSEN(FD(4), 3)
  IF (FD(6) .NE. 0) STOP 'J file cannot be closed after writing'
  IF (FD(5) .NE. 0) STOP 'W file cannot be closed after writing'
  IF (FD(4) .NE. 0) STOP 'Y file cannot be closed after writing'

  CALL BCLOSEN(FD(1), 3)
  IF (FD(3) .NE. 0) STOP 'U file cannot be closed'
  IF (FD(2) .NE. 0) STOP 'T file cannot be closed'
  IF (FD(1) .NE. 0) STOP 'X file cannot be closed'

CONTAINS
#include "readcl.F90"
#include "bio.F90"
#include "par_work.F90"
END PROGRAM PHASE1
