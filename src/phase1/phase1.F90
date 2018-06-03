PROGRAM PHASE1
  USE BINIO
  USE HIF_HZ
  USE TIMER
  IMPLICIT NONE

  INTEGER, PARAMETER :: FNL = 20

  CHARACTER(LEN=FNL,KIND=c_char) :: FN
  INTEGER :: L, a, G
  ! Cores Per Run and Threads Per Core
  INTEGER :: CPR, TPC
  INTEGER :: INFO, ATOM
  INTEGER :: SZ(6), FD(6)
  ! timing
  INTEGER :: CLK(3)
  DOUBLE PRECISION :: GTIMES(3)

  CALL READCL(FN, L, a, G, CPR, TPC, INFO)
  IF (INFO .EQ. 0) THEN
     WRITE (ULOG,'(2(A),3(I11),I3,I2)') 'phase1.exe ', TRIM(FN), L, a, G, CPR, TPC
  ELSE
     IF (INFO .LT. 0) THEN
        WRITE (ULOG,'(A,I2)') 'Cannot read argument', -INFO
     ELSE
        WRITE (ULOG,'(A,I2)') 'Illegal value of argument', INFO
     END IF
     STOP 'phase1.exe FN L a G CPR TPC'
  END IF

  CALL BOPEN_XTU_RO(FN, L, a, G, SZ, FD, INFO)
  WRITE (ULOG,'(A,I20,A)') 'X file has ', SZ(1), ' B.'
  IF (INFO .EQ. 1) STOP 'X file cannot be opened'
  WRITE (ULOG,'(A,I20,A)') 'T file has ', SZ(2), ' B.'
  IF (INFO .EQ. 2) STOP 'T file cannot be opened'
  WRITE (ULOG,'(A,I20,A)') 'U file has ', SZ(3), ' B.'
  IF (INFO .EQ. 3) STOP 'U file cannot be opened'
  IF (INFO .NE. 0) STOP 'BOPEN_XTU_RO: error'

  CALL BOPEN_YWJ_RW(FN, L, a, G, SZ(4), FD(4), INFO)
  IF (INFO .EQ. 1) STOP 'Y file cannot be opened for writing'
  IF (INFO .EQ. 2) STOP 'W file cannot be opened for writing'
  IF (INFO .EQ. 3) STOP 'J file cannot be opened for writing'
  IF (INFO .NE. 0) STOP 'BOPEN_YWJ_RW: error'

  INFO = BLAS_PREPARE()
  ! TODO: print some warning...
#ifndef MKL_NEST_SEQ
  TPC = MAX(1,TPC)
  IF (INFO .LE. 1) TPC = 1
#else
  TPC = MAX(0,TPC)
  IF (TPC .NE. 0) TPC = 0
#endif
  ATOM = 0
  INFO = 0
  CLK = 0
  GTIMES = -D_ZERO

  WRITE (UOUT,'(A)') '"ATOM","NRANK","NPLUS","N2PIV","INFO","TIME1","TIME2","TIME3"'
  CALL TIMER_START(CLK)
  !$OMP PARALLEL DO SHARED(FD,L,a,G,TPC) PRIVATE(ATOM,INFO) NUM_THREADS(CPR) PROC_BIND(SPREAD)
  DO ATOM = 1, a
     CALL PAR_WORK(FD, L, G, (ATOM-1), TPC, INFO)
  END DO
  !$OMP END PARALLEL DO
  CALL TIMER_STOP(CLK)
  GTIMES(1) = TIMER2DBLE(CLK)
  WRITE (UOUT,'(5(I11,A),2(F11.6,A),F11.6)') &
       -a,',', L,',', G,',', CPR,',', TPC,',', GTIMES(1),',', GTIMES(2),',', GTIMES(3)

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
#include "par_work.F90"
END PROGRAM PHASE1
