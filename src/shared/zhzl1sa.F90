! L1 double complex HZ (parallel, vectorized).
SUBROUTINE ZHZL1SA(M,N, H,LDH, JVEC, S,LDS, Z,LDZ, JS,JSPAIR, NSWP,CPR,TPC,&
     EE,EY,EW, SY,SW,SS, NROT,INFO)
#ifndef NDEBUG
  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE, INTRINSIC :: IEEE_FEATURES
  USE, INTRINSIC :: ISO_C_BINDING
#endif
  USE OMP_LIB
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: M,N, LDH,LDS,LDZ, JVEC(M), JS(JSMLEX),JSPAIR(2,JS(JSMLEX),JS(JSMLEX-1)), NSWP,CPR,TPC
  DOUBLE COMPLEX, INTENT(INOUT) :: H(LDH,N),S(LDS,N),Z(LDZ,N)
  DOUBLE PRECISION, INTENT(OUT) :: EE(N),EY(N),EW(N), SY(N),SW(N),SS(N)
  INTEGER, INTENT(OUT) :: NROT(2),INFO

  ! vector variables

  INTEGER :: HZ(ISIMDL,CPR)
  DOUBLE PRECISION :: DHZ(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: HZ, DHZ

  DOUBLE PRECISION :: RE_H_PP(DSIMDL,CPR)
  DOUBLE PRECISION :: RE_H_QQ(DSIMDL,CPR)
  DOUBLE PRECISION :: RE_H_PQ(DSIMDL,CPR)
  DOUBLE PRECISION :: IM_H_PQ(DSIMDL,CPR)
  DOUBLE PRECISION :: AV_H_PQ(DSIMDL,CPR)
  DOUBLE PRECISION :: CA_H_PQ(DSIMDL,CPR)
  DOUBLE PRECISION :: SA_H_PQ(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: RE_H_PP,RE_H_QQ,RE_H_PQ,IM_H_PQ, AV_H_PQ,CA_H_PQ,SA_H_PQ

  DOUBLE PRECISION :: RE_S_PP(DSIMDL,CPR)
  DOUBLE PRECISION :: RE_S_QQ(DSIMDL,CPR)
  DOUBLE PRECISION :: RE_S_PQ(DSIMDL,CPR)
  DOUBLE PRECISION :: IM_S_PQ(DSIMDL,CPR)
  DOUBLE PRECISION :: AV_S_PQ(DSIMDL,CPR)
  DOUBLE PRECISION :: CA_S_PQ(DSIMDL,CPR)
  DOUBLE PRECISION :: SA_S_PQ(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: RE_S_PP,RE_S_QQ,RE_S_PQ,IM_S_PQ, AV_S_PQ,CA_S_PQ,SA_S_PQ

  DOUBLE PRECISION :: T(DSIMDL,CPR)
  DOUBLE PRECISION :: U(DSIMDL,CPR)
  DOUBLE PRECISION :: V(DSIMDL,CPR)
  DOUBLE PRECISION :: E(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: T,U,V,E

  DOUBLE PRECISION :: TG(DSIMDL,CPR)
  DOUBLE PRECISION :: CG(DSIMDL,CPR)
  DOUBLE PRECISION :: SG(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: TG,CG,SG

  DOUBLE PRECISION :: T2T(DSIMDL,CPR)
  DOUBLE PRECISION :: C2T(DSIMDL,CPR)
  DOUBLE PRECISION :: S2T(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: T2T,C2T,S2T

  DOUBLE PRECISION :: CPHI(DSIMDL,CPR)
  DOUBLE PRECISION :: CPSI(DSIMDL,CPR)
  DOUBLE PRECISION :: RE_ASPHI(DSIMDL,CPR)
  DOUBLE PRECISION :: IM_ASPHI(DSIMDL,CPR)
  DOUBLE PRECISION :: RE_MBSPSI(DSIMDL,CPR)
  DOUBLE PRECISION :: IM_MBSPSI(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: CPHI,CPSI, RE_ASPHI,IM_ASPHI, RE_MBSPSI,IM_MBSPSI

  DOUBLE COMPLEX :: ZTMP1(DSIMDL,CPR)
  DOUBLE COMPLEX :: ZTMP2(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: ZTMP1,ZTMP2

  DOUBLE PRECISION :: DTMP1(DSIMDL,CPR)
  DOUBLE PRECISION :: DTMP2(DSIMDL,CPR)
  DOUBLE PRECISION :: DTMP3(DSIMDL,CPR)
  DOUBLE PRECISION :: DTMP4(DSIMDL,CPR)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: DTMP1,DTMP2,DTMP3,DTMP4

  INTEGER :: SNROT(2), LNROT(2)
  DOUBLE PRECISION :: DTOL
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: SNROT,LNROT, DTOL

  INTEGER :: NSTEPS, NPAIRS
  INTEGER :: PPV, VPS
  INTEGER :: SWEEP, STEP, VEC
  INTEGER :: PIX,PAIR, P,Q,R, I,J,L
#ifndef NDEBUG
  LOGICAL(c_int) :: LFHALT(5)
#endif

  EXTERNAL :: ZDSCAL, ZLASET

  !DIR$ ASSUME_ALIGNED JS:ALIGNB
  !DIR$ ASSUME_ALIGNED JSPAIR:ALIGNB
  !DIR$ ASSUME_ALIGNED H:ALIGNB
  !DIR$ ASSUME_ALIGNED S:ALIGNB
  !DIR$ ASSUME_ALIGNED Z:ALIGNB
  !DIR$ ASSUME_ALIGNED EE:ALIGNB
  !DIR$ ASSUME_ALIGNED EY:ALIGNB
  !DIR$ ASSUME_ALIGNED EW:ALIGNB
  !DIR$ ASSUME_ALIGNED SY:ALIGNB
  !DIR$ ASSUME_ALIGNED SY:ALIGNB
  !DIR$ ASSUME_ALIGNED SS:ALIGNB

  IF (M .LT. 0) THEN
     INFO = -1
  ELSE IF (N .LT. 0) THEN
     INFO = -2
  ELSE IF (N .GT. M) THEN
     INFO = -2
  ELSE IF (MOD(N,2) .NE. 0) THEN
     INFO = -2
  ELSE IF (LDH .LT. M) THEN
     INFO = -4
  ELSE IF (MOD(LDH,ZALIGN) .NE. 0) THEN
     INFO = -4
  ELSE IF (LDS .LT. M) THEN
     INFO = -7
  ELSE IF (MOD(LDS,ZALIGN) .NE. 0) THEN
     INFO = -7
  ELSE IF (LDZ .LT. N) THEN
     INFO = -9
  ELSE IF (MOD(LDZ,ZALIGN) .NE. 0) THEN
     INFO = -9
  ELSE IF (NSWP .LT. 0) THEN
     INFO = -12
  ELSE IF (CPR .LT. 1) THEN
     INFO = -13
  ELSE IF (CPR .GT. MAXCPR) THEN
     INFO = -13
  ELSE IF (TPC .LT. 1) THEN
     INFO = -14
  ELSE IF (TPC .GT. MAXTPC) THEN
     INFO = -14
  ELSE
     INFO = 0
  END IF
  IF (INFO .NE. 0) RETURN

  !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
  NROT = 0
  IF (N .EQ. 0) RETURN

#ifndef NDEBUG
  DO L = 1, 5
     CALL IEEE_GET_HALTING_MODE(IEEE_ALL(L), LFHALT(L))
  END DO
  CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .FALSE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .FALSE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_UNDERFLOW, .FALSE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_INEXACT, .FALSE._c_int)
#endif

  DTOL = SCALE(SQRT(DBLE(M)), -53)

  NSTEPS = JS(JSMLEX-1)
  NPAIRS = JS(JSMLEX)
  ! pairs per vector
  PPV = MIN(NPAIRS, DSIMDL)
  ! vectors per step
  VPS = (NPAIRS + (PPV - 1)) / PPV

  CALL ZLASET('A', N, N, Z_ZERO, Z_ONE, Z, LDZ)

  DO SWEEP = 1, NSWP
     !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
     SNROT = 0
     DO STEP = 1, NSTEPS
        !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
        LNROT = 0
        !$OMP  PARALLEL DEFAULT(NONE) REDUCTION(+:LNROT) PRIVATE(VEC,I,J,L,P,Q,R,PIX,PAIR) &
        !$OMP& SHARED(M,N,H,LDH,S,LDS,Z,LDZ,JVEC,DTOL,NPAIRS,JS,JSPAIR,STEP,PPV,VPS, HZ,DHZ, &
        !$OMP& RE_H_PP,RE_H_QQ,RE_H_PQ,IM_H_PQ,AV_H_PQ,CA_H_PQ,SA_H_PQ,RE_S_PP,RE_S_QQ,RE_S_PQ,IM_S_PQ,AV_S_PQ,CA_S_PQ,SA_S_PQ,&
        !$OMP& T,U,V,E,TG,CG,SG,T2T,C2T,S2T,CPHI,CPSI,RE_ASPHI,IM_ASPHI,RE_MBSPSI,IM_MBSPSI,ZTMP1,ZTMP2,DTMP1,DTMP2,DTMP3,DTMP4)
        !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
        LNROT = 0
        !$OMP DO
        DO VEC = 1, VPS
           R = INT(OMP_GET_THREAD_NUM()) + 1

           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, ISIMDL
              HZ(I,R) = 0
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              DHZ(I,R) = D_ZERO
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              RE_S_PP(I,R) = D_ONE
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              RE_S_QQ(I,R) = D_ONE
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              RE_S_PQ(I,R) = D_ZERO
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              IM_S_PQ(I,R) = D_ZERO
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              AV_S_PQ(I,R) = D_ZERO
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              CA_S_PQ(I,R) = D_ONE
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              SA_S_PQ(I,R) = D_ONE
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              RE_H_PP(I,R) = D_ONE
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              RE_H_QQ(I,R) = D_ONE
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              RE_H_PQ(I,R) = D_ZERO
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              IM_H_PQ(I,R) = D_ZERO
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              AV_H_PQ(I,R) = D_ZERO
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              CA_H_PQ(I,R) = D_ONE
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              SA_H_PQ(I,R) = D_ONE
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              TG(I,R) = D_ZERO
           END DO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO I = 1, DSIMDL
              SG(I,R) = D_ONE
           END DO

           ! compute the dot products

           DO PIX = 1, PPV
              ! ``global'' pair index
              PAIR = (VEC - 1) * PPV + PIX
              IF (PAIR .LE. NPAIRS) THEN
                 P = JSPAIR(1,PAIR,STEP)
                 Q = JSPAIR(2,PAIR,STEP)
                 ! ...dot products...

                 ! S

                 RE_S_PP(PIX,R) = D_ZERO
                 RE_S_QQ(PIX,R) = D_ZERO
                 RE_S_PQ(PIX,R) = D_ZERO
                 IM_S_PQ(PIX,R) = D_ZERO

                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    DTMP1(I,R) = D_ZERO ! RE_S_PP
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    DTMP2(I,R) = D_ZERO ! RE_S_QQ
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    DTMP3(I,R) = D_ZERO ! RE_S_PQ
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    DTMP4(I,R) = D_ZERO ! IM_S_PQ
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    ZTMP1(I,R) = Z_ZERO ! ZP
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    ZTMP2(I,R) = Z_ZERO ! ZQ
                 END DO

                 DO I = 1, M, DSIMDL
                    L = MIN(DSIMDL, M-(I-1))
                    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                    DO J = 1, L
                       ZTMP1(J,R) = S(I+(J-1),P)
                       ZTMP2(J,R) = S(I+(J-1),Q)
                       !DBLE(DCONJG(ZTMP1(J))*ZTMP1(J))
                       DTMP1(J,R) = DTMP1(J,R) + (DBLE(ZTMP1(J,R))*DBLE(ZTMP1(J,R)) + AIMAG(ZTMP1(J,R))*AIMAG(ZTMP1(J,R)))
                       !DBLE(DCONJG(ZTMP2(J))*ZTMP2(J))
                       DTMP2(J,R) = DTMP2(J,R) + (DBLE(ZTMP2(J,R))*DBLE(ZTMP2(J,R)) + AIMAG(ZTMP2(J,R))*AIMAG(ZTMP2(J,R)))
                       ! += DCONJG(ZTMP1(J)) * ZTMP2(J)
                       DTMP3(J,R) = DTMP3(J,R) + (DBLE(ZTMP1(J,R))*DBLE(ZTMP2(J,R)) + AIMAG(ZTMP1(J,R))*AIMAG(ZTMP2(J,R)))
                       DTMP4(J,R) = DTMP4(J,R) + (DBLE(ZTMP1(J,R))*AIMAG(ZTMP2(J,R))- AIMAG(ZTMP1(J,R))*DBLE(ZTMP2(J,R)))
                    END DO
                 END DO

                 RE_S_PP(PIX,R) = SUM(DTMP1(:,R))
                 RE_S_QQ(PIX,R) = SUM(DTMP2(:,R))
                 RE_S_PQ(PIX,R) = SUM(DTMP3(:,R))
                 IM_S_PQ(PIX,R) = SUM(DTMP4(:,R))

                 IF (RE_S_PP(PIX,R) .NE. RE_S_PP(PIX,R)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(S_pp)'
                 ELSE IF (RE_S_PP(PIX,R) .LE. D_ZERO) THEN
                    ! should never happen
                    STOP 'ZHZL1: S_pp .LE. 0'
                 ELSE IF (RE_S_PP(PIX,R) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    ! A joint prescaling of H and S needed...
                    STOP 'ZHZL1: Infinity(S_pp)'
                 ELSE IF (RE_S_PP(PIX,R) .NE. D_ONE) THEN
                    RE_S_PP(PIX,R) = D_ONE / SQRT(RE_S_PP(PIX,R))
                 END IF

                 IF (RE_S_QQ(PIX,R) .NE. RE_S_QQ(PIX,R)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(S_qq)'
                 ELSE IF (RE_S_QQ(PIX,R) .LE. D_ZERO) THEN
                    ! should never happen
                    STOP 'ZHZL1: S_qq .LE. 0'
                 ELSE IF (RE_S_QQ(PIX,R) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    ! A joint prescaling of H and S needed...
                    STOP 'ZHZL1: Infinity(S_qq)'
                 ELSE IF (RE_S_QQ(PIX,R) .NE. D_ONE) THEN
                    RE_S_QQ(PIX,R) = D_ONE / SQRT(RE_S_QQ(PIX,R))
                 END IF

                 ! H

                 RE_H_PP(PIX,R) = D_ZERO
                 RE_H_QQ(PIX,R) = D_ZERO
                 RE_H_PQ(PIX,R) = D_ZERO
                 IM_H_PQ(PIX,R) = D_ZERO

                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    DTMP1(I,R) = D_ZERO ! RE_H_PP
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    DTMP2(I,R) = D_ZERO ! RE_H_QQ
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    DTMP3(I,R) = D_ZERO ! RE_H_PQ
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    DTMP4(I,R) = D_ZERO ! IM_H_PQ
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    ZTMP1(I,R) = Z_ZERO ! ZP
                 END DO
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO I = 1, DSIMDL
                    ZTMP2(I,R) = Z_ZERO ! ZQ
                 END DO

                 DO I = 1, M, DSIMDL
                    L = MIN(DSIMDL, M-(I-1))
                    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                    DO J = 1, L
                       ZTMP1(J,R) = H(I+(J-1),P)
                       ZTMP2(J,R) = H(I+(J-1),Q)
                       DTMP1(J,R) = DTMP1(J,R) + JVEC(I+(J-1)) * (DBLE(ZTMP1(J,R))*DBLE(ZTMP1(J,R)) + AIMAG(ZTMP1(J,R))*AIMAG(ZTMP1(J,R)))
                       DTMP2(J,R) = DTMP2(J,R) + JVEC(I+(J-1)) * (DBLE(ZTMP2(J,R))*DBLE(ZTMP2(J,R)) + AIMAG(ZTMP2(J,R))*AIMAG(ZTMP2(J,R)))
                       DTMP3(J,R) = DTMP3(J,R) + JVEC(I+(J-1)) * (DBLE(ZTMP1(J,R))*DBLE(ZTMP2(J,R)) + AIMAG(ZTMP1(J,R))*AIMAG(ZTMP2(J,R)))
                       DTMP4(J,R) = DTMP4(J,R) + JVEC(I+(J-1)) * (DBLE(ZTMP1(J,R))*AIMAG(ZTMP2(J,R))- AIMAG(ZTMP1(J,R))*DBLE(ZTMP2(J,R)))
                    END DO
                 END DO

                 RE_H_PP(PIX,R) = SUM(DTMP1(:,R))
                 RE_H_QQ(PIX,R) = SUM(DTMP2(:,R))
                 RE_H_PQ(PIX,R) = SUM(DTMP3(:,R))
                 IM_H_PQ(PIX,R) = SUM(DTMP4(:,R))

                 IF (RE_H_PP(PIX,R) .NE. RE_H_PP(PIX,R)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(H_pp)'
                 ELSE IF (RE_H_PP(PIX,R) .EQ. D_ZERO) THEN
                    ! should never happen
                    STOP 'ZHZL1: H_pp .EQ. 0'
                 ELSE IF (RE_H_PP(PIX,R) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    ! A joint prescaling of H and S needed...
                    STOP 'ZHZL1: Infinity(H_pp)'
                 END IF

                 IF (RE_H_QQ(PIX,R) .NE. RE_H_QQ(PIX,R)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(H_qq)'
                 ELSE IF (RE_H_QQ(PIX,R) .EQ. D_ZERO) THEN
                    ! should never happen
                    STOP 'ZHZL1: H_qq .EQ. 0'
                 ELSE IF (RE_H_QQ(PIX,R) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    ! A joint prescaling of H and S needed...
                    STOP 'ZHZL1: Infinity(H_qq)'
                 END IF

                 IF (RE_H_PQ(PIX,R) .NE. RE_H_PQ(PIX,R)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(Re(H_pq))'
                 ELSE IF (ABS(RE_H_PQ(PIX,R)) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    STOP 'ZHZL1: Infinity(|Re(H_pq)|)'
                 END IF

                 IF (IM_H_PQ(PIX,R) .NE. IM_H_PQ(PIX,R)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(Im(H_pq))'
                 ELSE IF (ABS(IM_H_PQ(PIX,R)) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    STOP 'ZHZL1: Infinity(|Im(H_pq)|)'
                 END IF
              END IF
           END DO

           ! compute the transformation for a pair corresponding to the vector lane
  
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO PIX = 1, DSIMDL ! PPV
              ! compute the scales
              DTMP1(PIX,R) = RE_S_PP(PIX,R) * RE_S_PP(PIX,R)
              DTMP2(PIX,R) = RE_S_QQ(PIX,R) * RE_S_QQ(PIX,R)
              DTMP3(PIX,R) = RE_S_PP(PIX,R) * RE_S_QQ(PIX,R)
              ! scale H
              RE_H_PP(PIX,R) = RE_H_PP(PIX,R) * DTMP1(PIX,R)
              RE_H_QQ(PIX,R) = RE_H_QQ(PIX,R) * DTMP2(PIX,R)
              RE_H_PQ(PIX,R) = RE_H_PQ(PIX,R) * DTMP3(PIX,R)
              IM_H_PQ(PIX,R) = IM_H_PQ(PIX,R) * DTMP3(PIX,R)
              ! scale S
              RE_S_PQ(PIX,R) = RE_S_PQ(PIX,R) * DTMP3(PIX,R)
              IM_S_PQ(PIX,R) = IM_S_PQ(PIX,R) * DTMP3(PIX,R)
              ! compute ABS
              AV_H_PQ(PIX,R) = HYPOT(RE_H_PQ(PIX,R), IM_H_PQ(PIX,R))
              AV_S_PQ(PIX,R) = HYPOT(RE_S_PQ(PIX,R), IM_S_PQ(PIX,R))
              ! rotate or not
              ! 1 if H has to be rotated, else 0
              DTMP3(PIX,R) = SCALE(SIGN(D_ONE, AV_H_PQ(PIX,R) - SQRT(RE_H_PP(PIX,R))*SQRT(RE_H_QQ(PIX,R))*DTOL) + D_ONE, -1)
              ! 1 if S has to be rotated, else 0
              DTMP4(PIX,R) = SCALE(SIGN(D_ONE, AV_S_PQ(PIX,R) - DTOL) + D_ONE, -1)
              ! 1 if either H or S have to be rotated, else 0
              DHZ(PIX,R) = MAX(DTMP3(PIX,R), DTMP4(PIX,R))
              HZ(PIX,R) = INT(DHZ(PIX,R))
           END DO

           IF (MAXVAL(AV_H_PQ(:,R)) .GT. HUGE(D_ZERO)) STOP 'ZHZL1: |H_pq| overflow.'

           J = SUM(HZ(:,R))
           IF (J .EQ. 0) CYCLE
#ifndef NDEBUG
           IF (J .LT. 0) STOP 'ZHZL1: LNROT < 0'
           IF (J .GT. PPV) STOP 'ZHZL1: LNROT > PPV'
#endif
           LNROT(1) = LNROT(1) + J

           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO PIX = 1, DSIMDL ! PPV
              ! get the polar form
              DTMP1(PIX,R) = D_ONE / AV_H_PQ(PIX,R)
              DTMP2(PIX,R) = D_ONE / AV_S_PQ(PIX,R)
              ! expect 0/0 := NaN, 0*Inf := NaN, MIN(x,NaN) == x
              CA_H_PQ(PIX,R) = MIN(RE_H_PQ(PIX,R) * DTMP1(PIX,R), CA_H_PQ(PIX,R))
              SA_H_PQ(PIX,R) = MIN(IM_H_PQ(PIX,R) * DTMP1(PIX,R), SA_H_PQ(PIX,R)) * DTMP3(PIX,R)
              CA_S_PQ(PIX,R) = MIN(RE_S_PQ(PIX,R) * DTMP2(PIX,R), CA_S_PQ(PIX,R))
              SA_S_PQ(PIX,R) = MIN(IM_S_PQ(PIX,R) * DTMP2(PIX,R), SA_S_PQ(PIX,R)) * DTMP4(PIX,R)
              ! compute the temps
              T(PIX,R) = SQRT(D_ONE - AV_S_PQ(PIX,R) * AV_S_PQ(PIX,R))
              U(PIX,R) = CA_S_PQ(PIX,R) * RE_H_PQ(PIX,R) + SA_S_PQ(PIX,R) * IM_H_PQ(PIX,R)
              V(PIX,R) = CA_S_PQ(PIX,R) * IM_H_PQ(PIX,R) - SA_S_PQ(PIX,R) * RE_H_PQ(PIX,R)
              E(PIX,R) = RE_H_QQ(PIX,R) - RE_H_PP(PIX,R)
              DHZ(PIX,R) = SIGN(D_ONE, E(PIX,R))
              ! compute fns of \gamma
              DTMP1(PIX,R) = SCALE(V(PIX,R) / E(PIX,R), 1)
              ! V==0 & E==0 ==> avoid NaN(TG), set TG=0 (expect MAX(x,NaN) == x)
              TG(PIX,R) = SIGN(MAX(ABS(DTMP1(PIX,R)), TG(PIX,R)), DTMP1(PIX,R))
              CG(PIX,R) = D_ONE / SQRT(D_ONE + TG(PIX,R) * TG(PIX,R))
              ! beware of Inf(TG)
              SG(PIX,R) = SIGN(MIN(TG(PIX,R) * CG(PIX,R), SG(PIX,R)), TG(PIX,R))
              ! compute fns of 2\vartheta
              T2T(PIX,R) = (DHZ(PIX,R) * (SCALE(U(PIX,R), 1) - (RE_H_PP(PIX,R) + RE_H_QQ(PIX,R)) * AV_S_PQ(PIX,R))) / &
                   (T(PIX,R) * SQRT(E(PIX,R) * E(PIX,R) + SCALE(V(PIX,R) * V(PIX,R), 2)))
              C2T(PIX,R) = D_ONE / SQRT(D_ONE + T2T(PIX,R) * T2T(PIX,R))
              S2T(PIX,R) = T2T(PIX,R) * C2T(PIX,R)
              DTMP1(PIX,R) = D_ONE + T(PIX,R) * C2T(PIX,R) * CG(PIX,R)
              DHZ(PIX,R) = AV_S_PQ(PIX,R) * S2T(PIX,R)
              DTMP2(PIX,R) = DTMP1(PIX,R) - DHZ(PIX,R)
              DTMP1(PIX,R) = DTMP1(PIX,R) + DHZ(PIX,R)
              ! compute the transformation
              CPHI(PIX,R) = SQRT(SCALE(DTMP1(PIX,R), -1))
              CPSI(PIX,R) = SQRT(SCALE(DTMP2(PIX,R), -1))
              ! for big/small rot
              DHZ(PIX,R) = (D_ONE - CPHI(PIX,R)) + (D_ONE - CPSI(PIX,R))
              DTMP3(PIX,R) = S2T(PIX,R) - AV_S_PQ(PIX,R)
              DTMP4(PIX,R) = T(PIX,R) * SG(PIX,R) * C2T(PIX,R)
              TG(PIX,R) = -SA_S_PQ(PIX,R) * DTMP4(PIX,R)
              T2T(PIX,R) = CA_S_PQ(PIX,R) * DTMP4(PIX,R)
              T(PIX,R) = D_ONE / T(PIX,R)
              RE_ASPHI(PIX,R) = CA_S_PQ(PIX,R) * DTMP3(PIX,R) + TG(PIX,R)
              IM_ASPHI(PIX,R) = T2T(PIX,R) + SA_S_PQ(PIX,R) * DTMP3(PIX,R)
              DTMP4(PIX,R) = (CPSI(PIX,R) * RE_S_PP(PIX,R) * T(PIX,R)) / DTMP2(PIX,R)
              RE_ASPHI(PIX,R) = RE_ASPHI(PIX,R) * DTMP4(PIX,R)
              IM_ASPHI(PIX,R) = IM_ASPHI(PIX,R) * DTMP4(PIX,R)
              DTMP3(PIX,R) = S2T(PIX,R) + AV_S_PQ(PIX,R)
              RE_MBSPSI(PIX,R) = CA_S_PQ(PIX,R) * DTMP3(PIX,R) + TG(PIX,R)
              IM_MBSPSI(PIX,R) = T2T(PIX,R) + SA_S_PQ(PIX,R) * DTMP3(PIX,R)
              DTMP4(PIX,R) = (CPHI(PIX,R) * RE_S_QQ(PIX,R) * T(PIX,R)) / DTMP1(PIX,R)
              RE_MBSPSI(PIX,R) = -RE_MBSPSI(PIX,R) * DTMP4(PIX,R)
              IM_MBSPSI(PIX,R) = IM_MBSPSI(PIX,R) * DTMP4(PIX,R)
              CPHI(PIX,R) = CPHI(PIX,R) * RE_S_PP(PIX,R) * T(PIX,R)
              CPSI(PIX,R) = CPSI(PIX,R) * RE_S_QQ(PIX,R) * T(PIX,R)
           END DO

           ! apply the transformations

           DO PIX = 1, PPV
              IF (HZ(PIX,R) .EQ. 0) CYCLE
              IF (DHZ(PIX,R) .GT. D_ZERO) LNROT(2) = LNROT(2) + 1
              ! ``global'' pair index
              PAIR = (VEC - 1) * PPV + PIX
              IF (PAIR .LE. NPAIRS) THEN
                 P = JSPAIR(1,PAIR,STEP)
                 Q = JSPAIR(2,PAIR,STEP)
                 ! ...transform...
                 IF (.NOT. (CPHI(PIX,R) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: F_11 overflow or NaN.'
                 DTMP1(PIX,R) = CPHI(PIX,R)
                 IF (.NOT. (ABS(RE_MBSPSI(PIX,R)) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: |Re(F_21)| overflow or NaN.'
                 IF (.NOT. (ABS(IM_MBSPSI(PIX,R)) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: |Im(F_21)| overflow or NaN.'
                 ZTMP1(PIX,R) = DCMPLX(RE_MBSPSI(PIX,R), IM_MBSPSI(PIX,R))
                 IF (.NOT. (ABS(RE_ASPHI(PIX,R)) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: |Re(F_12)| overflow or NaN.'
                 IF (.NOT. (ABS(IM_ASPHI(PIX,R)) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: |Im(F_12)| overflow or NaN.'
                 ZTMP2(PIX,R) = DCMPLX(RE_ASPHI(PIX,R), IM_ASPHI(PIX,R))
                 IF (.NOT. (CPSI(PIX,R) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: F_22 overflow or NaN.'
                 DTMP2(PIX,R) = CPSI(PIX,R)
                 CALL BLAS_ZVROTM(M, H(1,P), H(1,Q), DTMP1(PIX,R), ZTMP1(PIX,R), ZTMP2(PIX,R), DTMP2(PIX,R))
                 CALL BLAS_ZVROTM(M, S(1,P), S(1,Q), DTMP1(PIX,R), ZTMP1(PIX,R), ZTMP2(PIX,R), DTMP2(PIX,R))
                 CALL BLAS_ZVROTM(N, Z(1,P), Z(1,Q), DTMP1(PIX,R), ZTMP1(PIX,R), ZTMP2(PIX,R), DTMP2(PIX,R))
              END IF
           END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
        !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
        SNROT = SNROT + LNROT
     END DO
     WRITE (ULOG,'(I3,A,I20,A,I20)') SWEEP,',',SNROT(1),',',SNROT(2)
     !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
     NROT = NROT + SNROT
     IF (SNROT(2) .EQ. 0) EXIT
  END DO
  INFO = SWEEP

  ! Scaling of H,S,Z.

  IF (NROT(1) .GT. 0) THEN
     !$OMP  PARALLEL DEFAULT(NONE) PRIVATE(I,J,L,P,Q,R,DTOL) &
     !$OMP& SHARED(M,N,H,LDH,S,LDS,Z,LDZ,JVEC,TPC,EE,EY,EW,SY,SW,SS, DTMP1,DTMP2,ZTMP1,ZTMP2)
     Q = BLAS_SET_NUM_THREADS(TPC)
     !$OMP DO
     DO J = 1, N
        R = INT(OMP_GET_THREAD_NUM()) + 1

        !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
        DO I = 1, DSIMDL
           DTMP1(I,R) = D_ZERO
        END DO

        DO I = 1, M, DSIMDL
           P = MIN(DSIMDL, M-(I-1))
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO L = 1, P
              ZTMP1(L,R) = H(I+(L-1),J)
              DTMP1(L,R) = DTMP1(L,R) + JVEC(I+(L-1))*(DBLE(ZTMP1(L,R))*DBLE(ZTMP1(L,R)) + AIMAG(ZTMP1(L,R))*AIMAG(ZTMP1(L,R)))
           END DO
        END DO
        EY(J) = SUM(DTMP1(:,R))
        IF (EY(J) .NE. D_ZERO) THEN
           DTOL = ABS(EY(J))
           IF (DTOL .NE. D_ONE) THEN
              SY(J) = SQRT(DTOL)
              DTOL = D_ONE / SY(J)
              CALL ZDSCAL(M, DTOL, H(1,J), 1)
           ELSE
              SY(J) = D_ONE
           END IF
        ELSE
           SY(J) = D_ZERO
        END IF

        !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
        DO I = 1, DSIMDL
           DTMP2(I,R) = D_ZERO
        END DO

        DO I = 1, M, DSIMDL
           P = MIN(DSIMDL, M-(I-1))
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO L = 1, P
              ZTMP2(L,R) = S(I+(L-1),J)
              DTMP2(L,R) = DTMP2(L,R) + (DBLE(ZTMP2(L,R))*DBLE(ZTMP2(L,R)) + AIMAG(ZTMP2(L,R))*AIMAG(ZTMP2(L,R)))
           END DO
        END DO
        EW(J) = SUM(DTMP2(:,R))
        IF (EW(J) .NE. D_ZERO) THEN
           DTOL = EW(J)
           IF (DTOL .NE. D_ONE) THEN
              SW(J) = SQRT(DTOL)
              DTOL = D_ONE / SW(J)
              CALL ZDSCAL(M, DTOL, S(1,J), 1)
           ELSE
              SW(J) = D_ONE
           END IF
        ELSE
           SW(J) = D_ZERO
        END IF

        EE(J) = EY(J) / EW(J)
        SS(J) = HYPOT(SY(J), SW(J))
        IF (SS(J) .NE. D_ONE) THEN
           ! underflow
           IF (SS(J) .LT. TINY(D_ZERO)) STOP 'ZHZL1: Scale of Z underflows.'
           ! overflow
           IF (SS(J) .GT. HUGE(D_ZERO)) STOP 'ZHZL1: Scale of Z overflows.'
           DTOL = D_ONE / SS(J)
           IF (DTOL .LT. TINY(D_ZERO)) THEN
              ! underflow
              DTOL = SS(J)
              SY(J) = SY(J) / DTOL
              SW(J) = SW(J) / DTOL
              DO I = 1, N, DSIMDL
                 P = MIN(DSIMDL, N-(I-1))
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DO L = 1, P
                    Z(I+(L-1),J) = Z(I+(L-1),J) / DTOL
                 END DO
              END DO
           ELSE
              SY(J) = SY(J) * DTOL
              SW(J) = SW(J) * DTOL
              CALL ZDSCAL(N, DTOL, Z(1,J), 1)
           END IF
        END IF
     END DO
     !$OMP END DO
     Q = BLAS_SET_NUM_THREADS(Q)
     !$OMP END PARALLEL
  END IF

#ifndef NDEBUG
  DO L = 1, 5
     CALL IEEE_SET_HALTING_MODE(IEEE_ALL(L), LFHALT(L))
  END DO
#endif
END SUBROUTINE ZHZL1SA
