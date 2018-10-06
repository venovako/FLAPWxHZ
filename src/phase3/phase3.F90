PROGRAM PHASE3
  USE BINIO
  USE HIF_HZ
  USE TIMER
  IMPLICIT NONE

  INTEGER, PARAMETER :: FNL = 20

  CHARACTER(LEN=FNL,KIND=c_char) :: FN
  INTEGER :: M, N
  ! 2nd-level and 1st-level threads
  INTEGER :: CPR, TPC
  INTEGER :: JSTRAT(2), NSWP(2)
  INTEGER :: SZ(12), FD(12)
  TYPE(ZSHENT) :: PSHBUF(MAXCPR)
  INTEGER :: PSTATS(8)
  DOUBLE PRECISION :: MXTIME
  INTEGER :: INFO

  CALL READCL(FN, M, N, CPR, TPC, JSTRAT, NSWP, INFO)
  IF (INFO .NE. 0) THEN
     IF (INFO .LT. 0) THEN
        WRITE (ULOG,'(A,I2)') 'Cannot read argument', -INFO
     ELSE
        WRITE (ULOG,'(A,I2)') 'Illegal value of argument', INFO
     END IF
  END IF

  CALL BOPEN_YWJ_RO(FN, M, N, SZ, FD, INFO)
  IF (INFO .EQ. 1) STOP 'Y file cannot be opened'
  IF (INFO .EQ. 2) STOP 'W file cannot be opened'
  IF (INFO .EQ. 3) STOP 'J file cannot be opened'
  IF (INFO .NE. 0) STOP 'BOPEN_YWJ_RO: error'

  CALL BOPEN_EZS_RW(FN, M, N, SZ(4), FD(4), INFO)
  IF (INFO .EQ. 1) STOP ' E file cannot be opened for writing'
  IF (INFO .EQ. 2) STOP ' Z file cannot be opened for writing'
  IF (INFO .EQ. 3) STOP 'SS file cannot be opened for writing'
  IF (INFO .EQ. 4) STOP 'YU file cannot be opened for writing'
  IF (INFO .EQ. 5) STOP 'EY file cannot be opened for writing'
  IF (INFO .EQ. 6) STOP 'SY file cannot be opened for writing'
  IF (INFO .EQ. 7) STOP 'WV file cannot be opened for writing'
  IF (INFO .EQ. 8) STOP 'EW file cannot be opened for writing'
  IF (INFO .EQ. 9) STOP 'SW file cannot be opened for writing'
  IF (INFO .NE. 0) STOP 'BOPEN_ZES_RW: error'

  INFO = BLAS_PREPARE()
#ifndef MKL_NEST_SEQ
  TPC = MAX(1,TPC)
  IF (INFO .LE. 1) TPC = 1
#else
  TPC = MAX(0,TPC)
  IF (TPC .NE. 0) TPC = 0
#endif
  PSTATS = 0
  MXTIME = D_ZERO
  INFO = 0
#ifndef MKL_NEST_SEQ
  !$OMP PARALLEL SHARED(FD,M,N,CPR,TPC,JSTRAT,NSWP,PSHBUF,PSTATS) NUM_THREADS(CPR) PROC_BIND(SPREAD) REDUCTION(MAX:MXTIME,INFO)
#else
  !$OMP PARALLEL SHARED(FD,M,N,CPR,JSTRAT,NSWP,PSHBUF,PSTATS) NUM_THREADS(CPR) PROC_BIND(SPREAD) REDUCTION(MAX:MXTIME,INFO)
#endif
#ifndef MKL_NEST_SEQ
  CALL PAR_WORK(FD, M,N, CPR,TPC, JSTRAT,NSWP, PSHBUF,PSTATS, MXTIME,INFO)
#else
  CALL PAR_WORK(FD, M,N, CPR, JSTRAT,NSWP, PSHBUF,PSTATS, MXTIME,INFO)
#endif
  !$OMP END PARALLEL
  WRITE (UOUT,'(A)') &
       '"MKL","FN","M","N","CPR","TPC","JS1","NS1","JS2","NS2","INFO","TIME","SWP","T1","T2","T3","ALLROT","BIGROT"'
  ! SEQ: use sequential MKL and OMP_PLACES=CORES (at most 64 cores on our Xeon Phi)
  ! PAR: use multi-threaded MKL and OMP_PLACES=THREADS (at most 64*4=256 threads on our Xeon Phi)
#ifdef MKL_NEST_SEQ
#ifdef NDEBUG
  WRITE (UOUT,'(A)',ADVANCE='NO') 'SEQ,'
#else
  WRITE (UOUT,'(A)',ADVANCE='NO') 'seq,'
#endif
#else
#ifdef NDEBUG
  WRITE (UOUT,'(A)',ADVANCE='NO') 'PAR,'
#else
  WRITE (UOUT,'(A)',ADVANCE='NO') 'par,'
#endif
#endif
  WRITE (UOUT,'(2(A),2(I11,A),2(I3,A),2(I2,A,I3,A))',ADVANCE='NO') TRIM(FN),',', M,',',N,',', CPR,',',TPC,',', &
       JSTRAT(1),',',NSWP(1),',', JSTRAT(2),',',NSWP(2),','
  WRITE (UOUT,'(I3,A,F12.6,A,I3)',ADVANCE='NO') INFO,',',MXTIME,',',PSTATS(8)
  WRITE (UOUT,'(3(A,F11.6))',ADVANCE='NO') ',',(PSTATS(1)*DNS2S),',',(PSTATS(2)*DNS2S),',',(PSTATS(3)*DNS2S)
  WRITE (UOUT,'(2(A,I20))') ',',PSTATS(6),',',PSTATS(7)

  CALL BCLOSEN(FD(4), 9)
  IF (FD(12) .NE. 0) STOP 'SW file cannot be closed after writing'
  IF (FD(11) .NE. 0) STOP 'EW file cannot be closed after writing'
  IF (FD(10) .NE. 0) STOP 'WV file cannot be closed after writing'
  IF (FD( 9) .NE. 0) STOP 'SY file cannot be closed after writing'
  IF (FD( 8) .NE. 0) STOP 'EY file cannot be closed after writing'
  IF (FD( 7) .NE. 0) STOP 'YU file cannot be closed after writing'
  IF (FD( 6) .NE. 0) STOP 'SS file cannot be closed after writing'
  IF (FD( 5) .NE. 0) STOP ' Z file cannot be closed after writing'
  IF (FD( 4) .NE. 0) STOP ' E file cannot be closed after writing'

  CALL BCLOSEN(FD(1), 3)
  IF (FD(3) .NE. 0) STOP 'U file cannot be closed'
  IF (FD(2) .NE. 0) STOP 'T file cannot be closed'
  IF (FD(1) .NE. 0) STOP 'X file cannot be closed'

CONTAINS
#include "readcl.F90"
#include "bio.F90"
#include "par_work.F90"
END PROGRAM PHASE3
