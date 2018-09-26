! L1 double complex HZ (parallel, vectorized).
SUBROUTINE ZHZL1P(M,N, H,LDH, JVEC, S,LDS, Z,LDZ, JS,JSPAIR, NSWP,CPR,&
     EE,EY,EW, SY,SW,SS, NROT,INFO)
#ifndef NDEBUG
  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE, INTRINSIC :: IEEE_FEATURES
#endif
  USE JACSTR
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: M,N, LDH,LDS,LDZ, JVEC(M), JS(JSMLEX),JSPAIR(2,JS(JSMLEX),JS(JSMLEX-1)), NSWP,CPR
  DOUBLE COMPLEX, INTENT(INOUT) :: H(LDH,N),S(LDS,N),Z(LDZ,N)
  DOUBLE PRECISION, INTENT(OUT) :: EE(N),EY(N),EW(N), SY(N),SW(N),SS(N)
  INTEGER, INTENT(OUT) :: NROT(2),INFO

  INTEGER :: NSTEPS, NPAIRS
  INTEGER :: PPV, VPS
  INTEGER :: SWEEP, STEP, VEC, PIX, PAIR
  INTEGER :: P, Q, I, J, L

  INTEGER :: SNROT(2)
  DOUBLE PRECISION :: DTOL
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: SNROT, DTOL

  ! vector variables

  INTEGER :: HZ(ISIMDL)
  DOUBLE PRECISION :: DHZ(DSIMDL)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: HZ, DHZ

  DOUBLE PRECISION :: RE_H_PP(DSIMDL)
  DOUBLE PRECISION :: RE_H_QQ(DSIMDL)
  DOUBLE PRECISION :: RE_H_PQ(DSIMDL)
  DOUBLE PRECISION :: IM_H_PQ(DSIMDL)
  DOUBLE PRECISION :: AV_H_PQ(DSIMDL)
  DOUBLE PRECISION :: CA_H_PQ(DSIMDL)
  DOUBLE PRECISION :: SA_H_PQ(DSIMDL)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: RE_H_PP,RE_H_QQ,RE_H_PQ,IM_H_PQ, AV_H_PQ,CA_H_PQ,SA_H_PQ

  DOUBLE PRECISION :: RE_S_PP(DSIMDL)
  DOUBLE PRECISION :: RE_S_QQ(DSIMDL)
  DOUBLE PRECISION :: RE_S_PQ(DSIMDL)
  DOUBLE PRECISION :: IM_S_PQ(DSIMDL)
  DOUBLE PRECISION :: AV_S_PQ(DSIMDL)
  DOUBLE PRECISION :: CA_S_PQ(DSIMDL)
  DOUBLE PRECISION :: SA_S_PQ(DSIMDL)  
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: RE_S_PP,RE_S_QQ,RE_S_PQ,IM_S_PQ, AV_S_PQ,CA_S_PQ,SA_S_PQ

  DOUBLE PRECISION :: T(DSIMDL)
  DOUBLE PRECISION :: U(DSIMDL)
  DOUBLE PRECISION :: V(DSIMDL)
  DOUBLE PRECISION :: E(DSIMDL)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: T,U,V,E

  DOUBLE PRECISION :: TG(DSIMDL)
  DOUBLE PRECISION :: CG(DSIMDL)
  DOUBLE PRECISION :: SG(DSIMDL)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: TG,CG,SG

  DOUBLE PRECISION :: T2T(DSIMDL)
  DOUBLE PRECISION :: C2T(DSIMDL)
  DOUBLE PRECISION :: S2T(DSIMDL)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: T2T,C2T,S2T

  DOUBLE PRECISION :: CPHI(DSIMDL)
  DOUBLE PRECISION :: CPSI(DSIMDL)
  DOUBLE PRECISION :: RE_ASPHI(DSIMDL)
  DOUBLE PRECISION :: IM_ASPHI(DSIMDL)
  DOUBLE PRECISION :: RE_MBSPSI(DSIMDL)
  DOUBLE PRECISION :: IM_MBSPSI(DSIMDL)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: CPHI,CPSI, RE_ASPHI,IM_ASPHI, RE_MBSPSI,IM_MBSPSI

  DOUBLE COMPLEX :: ZTMP1(DSIMDL)
  DOUBLE COMPLEX :: ZTMP2(DSIMDL)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: ZTMP1,ZTMP2

  DOUBLE PRECISION :: DTMP1(DSIMDL)
  DOUBLE PRECISION :: DTMP2(DSIMDL)
  DOUBLE PRECISION :: DTMP3(DSIMDL)
  DOUBLE PRECISION :: DTMP4(DSIMDL)
  !DIR$ ATTRIBUTES ALIGN:ALIGNB :: DTMP1,DTMP2,DTMP3,DTMP4

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
#ifdef NDEBUG
  !DIR$ ASSUME (M .GE. 0)
  !DIR$ ASSUME (N .GE. 0)
  !DIR$ ASSUME (M .GE. N)
  !DIR$ ASSUME (MOD(N,2) .EQ. 0)
  !DIR$ ASSUME (MOD(LDH,ZALIGN) .EQ. 0)
  !DIR$ ASSUME (MOD(LDS,ZALIGN) .EQ. 0)
  !DIR$ ASSUME (MOD(LDZ,ZALIGN) .EQ. 0)
  !DIR$ ASSUME (NSWP .GE. 0)
  !DIR$ ASSUME (CPR .GE. 1)
#else
  IF (M .LT. 0) STOP 'ZHZL1: M .LT. 0'
  IF (N .LT. 0) STOP 'ZHZL1: N .LT. 0'
  IF (M .LT. N) STOP 'ZHZL1: M .LT. N'
  IF (MOD(N,2) .NE. 0) STOP 'ZHZL1: MOD(N,2) .NE. 0'
  IF (MOD(LDH,ZALIGN) .NE. 0) STOP 'ZHZL1: MOD(LDH,ZALIGN) .NE. 0'
  IF (MOD(LDS,ZALIGN) .NE. 0) STOP 'ZHZL1: MOD(LDS,ZALIGN) .NE. 0'
  IF (MOD(LDZ,ZALIGN) .NE. 0) STOP 'ZHZL1: MOD(LDZ,ZALIGN) .NE. 0'
  IF (NSWP .LT. 0) STOP 'ZHZL1: NSWP .LT. 0'
  IF (CPR .LT. 1) STOP 'ZHZL1: CPR .LT. 1'
#endif

  INFO = 0
  !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
  NROT = 0
  IF (N .LE. 0) RETURN

#ifndef NDEBUG
  CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .FALSE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .FALSE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE._c_int)
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
        !$OMP  PARALLEL DO DEFAULT(NONE) NUM_THREADS(CPR) PROC_BIND(SPREAD)              &
        !$OMP& SHARED(M,N,JSPAIR,NPAIRS,STEP,PPV,VPS,H,S,Z,JVEC,DTOL)                    &
        !$OMP& PRIVATE(VEC,PIX,PAIR, P,Q, I,J,L,                                         &
        !$OMP& RE_H_PP,RE_H_QQ,RE_H_PQ,IM_H_PQ, RE_S_PP,RE_S_QQ,RE_S_PQ,IM_S_PQ,         &
        !$OMP& HZ,DHZ, AV_H_PQ,CA_H_PQ,SA_H_PQ, AV_S_PQ,CA_S_PQ,SA_S_PQ, T,U,V,E,        &
        !$OMP& TG,CG,SG, T2T,C2T,S2T, CPHI,CPSI, RE_ASPHI,IM_ASPHI, RE_MBSPSI,IM_MBSPSI, &
        !$OMP& ZTMP1,ZTMP2, DTMP1,DTMP2,DTMP3,DTMP4) REDUCTION(+:SNROT)
        DO VEC = 1, VPS
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           HZ = 0
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DHZ = D_ZERO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           RE_S_PP = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           RE_S_QQ = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           RE_S_PQ = D_ZERO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED  
           IM_S_PQ = D_ZERO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           AV_S_PQ = D_ZERO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           CA_S_PQ = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           SA_S_PQ = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           RE_H_PP = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           RE_H_QQ = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           RE_H_PQ = D_ZERO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           IM_H_PQ = D_ZERO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           AV_H_PQ = D_ZERO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           CA_H_PQ = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           SA_H_PQ = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           TG = D_ZERO
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           SG = D_ONE

           ! compute the dot products

           DO PIX = 1, PPV
              ! ``global'' pair index
              PAIR = (VEC - 1) * PPV + PIX
              IF (PAIR .LE. NPAIRS) THEN
                 P = JSPAIR(1,PAIR,STEP)
                 Q = JSPAIR(2,PAIR,STEP)
                 ! ...dot products...

                 ! S

                 RE_S_PP(PIX) = D_ZERO
                 RE_S_QQ(PIX) = D_ZERO
                 RE_S_PQ(PIX) = D_ZERO
                 IM_S_PQ(PIX) = D_ZERO

                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DTMP1 = D_ZERO ! RE_S_PP
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DTMP2 = D_ZERO ! RE_S_QQ
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DTMP3 = D_ZERO ! RE_S_PQ
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DTMP4 = D_ZERO ! IM_S_PQ
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 ZTMP1 = Z_ZERO ! ZP
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 ZTMP2 = Z_ZERO ! ZQ

                 DO I = 1, M, DSIMDL
                    L = MIN(DSIMDL, M-(I-1))
                    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                    DO J = 1, L
                       ZTMP1(J) = S(I+(J-1),P)
                       ZTMP2(J) = S(I+(J-1),Q)
                       !DBLE(DCONJG(ZTMP1(J))*ZTMP1(J))
                       DTMP1(J) = DTMP1(J) + (DBLE(ZTMP1(J))*DBLE(ZTMP1(J)) + AIMAG(ZTMP1(J))*AIMAG(ZTMP1(J)))
                       !DBLE(DCONJG(ZTMP2(J))*ZTMP2(J))
                       DTMP2(J) = DTMP2(J) + (DBLE(ZTMP2(J))*DBLE(ZTMP2(J)) + AIMAG(ZTMP2(J))*AIMAG(ZTMP2(J)))
                       ! += DCONJG(ZTMP1(J)) * ZTMP2(J)
                       DTMP3(J) = DTMP3(J) + (DBLE(ZTMP1(J))*DBLE(ZTMP2(J)) + AIMAG(ZTMP1(J))*AIMAG(ZTMP2(J)))
                       DTMP4(J) = DTMP4(J) + (DBLE(ZTMP1(J))*AIMAG(ZTMP2(J))- AIMAG(ZTMP1(J))*DBLE(ZTMP2(J)))
                    END DO
                 END DO

                 RE_S_PP(PIX) = SUM(DTMP1)
                 RE_S_QQ(PIX) = SUM(DTMP2)
                 RE_S_PQ(PIX) = SUM(DTMP3)
                 IM_S_PQ(PIX) = SUM(DTMP4)

                 IF (RE_S_PP(PIX) .NE. RE_S_PP(PIX)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(S_pp)'
                 ELSE IF (RE_S_PP(PIX) .LE. D_ZERO) THEN
                    ! should never happen
                    STOP 'ZHZL1: S_pp .LE. 0'
                 ELSE IF (RE_S_PP(PIX) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    ! A joint prescaling of H and S needed...
                    STOP 'ZHZL1: Infinity(S_pp)'
                 ELSE IF (RE_S_PP(PIX) .NE. D_ONE) THEN
                    RE_S_PP(PIX) = D_ONE / SQRT(RE_S_PP(PIX))
                 END IF

                 IF (RE_S_QQ(PIX) .NE. RE_S_QQ(PIX)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(S_qq)'
                 ELSE IF (RE_S_QQ(PIX) .LE. D_ZERO) THEN
                    ! should never happen
                    STOP 'ZHZL1: S_qq .LE. 0'
                 ELSE IF (RE_S_QQ(PIX) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    ! A joint prescaling of H and S needed...
                    STOP 'ZHZL1: Infinity(S_qq)'
                 ELSE IF (RE_S_QQ(PIX) .NE. D_ONE) THEN
                    RE_S_QQ(PIX) = D_ONE / SQRT(RE_S_QQ(PIX))
                 END IF
                 
                 ! H

                 RE_H_PP(PIX) = D_ZERO
                 RE_H_QQ(PIX) = D_ZERO
                 RE_H_PQ(PIX) = D_ZERO
                 IM_H_PQ(PIX) = D_ZERO

                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DTMP1 = D_ZERO ! RE_H_PP
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DTMP2 = D_ZERO ! RE_H_QQ
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DTMP3 = D_ZERO ! RE_H_PQ
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 DTMP4 = D_ZERO ! IM_H_PQ
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 ZTMP1 = Z_ZERO ! ZP
                 !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                 ZTMP2 = Z_ZERO ! ZQ

                 DO I = 1, M, DSIMDL
                    L = MIN(DSIMDL, M-(I-1))
                    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                    DO J = 1, L
                       ZTMP1(J) = H(I+(J-1),P)
                       ZTMP2(J) = H(I+(J-1),Q)
                       DTMP1(J) = DTMP1(J) + JVEC(M) * (DBLE(ZTMP1(J))*DBLE(ZTMP1(J)) + AIMAG(ZTMP1(J))*AIMAG(ZTMP1(J)))
                       DTMP2(J) = DTMP2(J) + JVEC(M) * (DBLE(ZTMP2(J))*DBLE(ZTMP2(J)) + AIMAG(ZTMP2(J))*AIMAG(ZTMP2(J)))
                       DTMP3(J) = DTMP3(J) + JVEC(M) * (DBLE(ZTMP1(J))*DBLE(ZTMP2(J)) + AIMAG(ZTMP1(J))*AIMAG(ZTMP2(J)))
                       DTMP4(J) = DTMP4(J) + JVEC(M) * (DBLE(ZTMP1(J))*AIMAG(ZTMP2(J))- AIMAG(ZTMP1(J))*DBLE(ZTMP2(J)))
                    END DO
                 END DO

                 RE_H_PP(PIX) = SUM(DTMP1)
                 RE_H_QQ(PIX) = SUM(DTMP2)
                 RE_H_PQ(PIX) = SUM(DTMP3)
                 IM_H_PQ(PIX) = SUM(DTMP4)

                 IF (RE_H_PP(PIX) .NE. RE_H_PP(PIX)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(H_pp)'
                 ELSE IF (RE_H_PP(PIX) .EQ. D_ZERO) THEN
                    ! should never happen
                    STOP 'ZHZL1: H_pp .EQ. 0'
                 ELSE IF (RE_H_PP(PIX) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    ! A joint prescaling of H and S needed...
                    STOP 'ZHZL1: Infinity(H_pp)'
                 END IF

                 IF (RE_H_QQ(PIX) .NE. RE_H_QQ(PIX)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(H_qq)'
                 ELSE IF (RE_H_QQ(PIX) .EQ. D_ZERO) THEN
                    ! should never happen
                    STOP 'ZHZL1: H_qq .EQ. 0'
                 ELSE IF (RE_H_QQ(PIX) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    ! A joint prescaling of H and S needed...
                    STOP 'ZHZL1: Infinity(H_qq)'
                 END IF

                 IF (RE_H_PQ(PIX) .NE. RE_H_PQ(PIX)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(Re(H_pq))'
                 ELSE IF (ABS(RE_H_PQ(PIX)) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    STOP 'ZHZL1: Infinity(|Re(H_pq)|)'
                 END IF

                 IF (IM_H_PQ(PIX) .NE. IM_H_PQ(PIX)) THEN
                    ! NaN
                    STOP 'ZHZL1: NaN(Im(H_pq))'
                 ELSE IF (ABS(IM_H_PQ(PIX)) .GT. HUGE(D_ZERO)) THEN
                    ! overflow
                    STOP 'ZHZL1: Infinity(|Im(H_pq)|)'
                 END IF
              END IF
           END DO

           ! compute the transformation for a pair corresponding to the vector lane

           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO PIX = 1, DSIMDL ! PPV
              ! compute the scales
              DTMP1(PIX) = RE_S_PP(PIX) * RE_S_PP(PIX)
              DTMP2(PIX) = RE_S_QQ(PIX) * RE_S_QQ(PIX)
              DTMP3(PIX) = RE_S_PP(PIX) * RE_S_QQ(PIX)
              ! scale H
              RE_H_PP(PIX) = RE_H_PP(PIX) * DTMP1(PIX)
              RE_H_QQ(PIX) = RE_H_QQ(PIX) * DTMP2(PIX)
              RE_H_PQ(PIX) = RE_H_PQ(PIX) * DTMP3(PIX)
              IM_H_PQ(PIX) = IM_H_PQ(PIX) * DTMP3(PIX)
              ! scale S
              RE_S_PQ(PIX) = RE_S_PQ(PIX) * DTMP3(PIX)
              IM_S_PQ(PIX) = IM_S_PQ(PIX) * DTMP3(PIX)
              ! compute ABS
              AV_H_PQ(PIX) = HYPOT(RE_H_PQ(PIX), IM_H_PQ(PIX))
              AV_S_PQ(PIX) = HYPOT(RE_S_PQ(PIX), IM_S_PQ(PIX))
              ! rotate or not
              ! 1 if H has to be rotated, else 0
              DTMP3(PIX) = SCALE(SIGN(D_ONE, AV_H_PQ(PIX) - SQRT(RE_H_PP(PIX)) * SQRT(RE_H_QQ(PIX)) * DTOL) + D_ONE, -1)
              ! 1 if S has to be rotated, else 0
              DTMP4(PIX) = SCALE(SIGN(D_ONE, AV_S_PQ(PIX) - DTOL) + D_ONE, -1)
              ! 1 if either H or S have to be rotated, else 0
              DHZ(PIX) = MAX(DTMP3(PIX), DTMP4(PIX))
              HZ(PIX) = INT(DHZ(PIX))
           END DO

           IF (MAXVAL(AV_H_PQ) .GT. HUGE(D_ZERO)) STOP 'ZHZL1: |H_pq| overflow.'

           J = SUM(HZ)
           IF (J .EQ. 0) CYCLE
#ifndef NDEBUG
           IF (J .LT. 0) STOP 'ZHZL1: SNROT < 0'
           IF (J .GT. PPV) STOP 'ZHZL1: SNROT > PPV'
#endif
           SNROT(1) = SNROT(1) + J

           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO PIX = 1, DSIMDL ! PPV
              ! get the polar form
              DTMP1(PIX) = D_ONE / AV_H_PQ(PIX)
              DTMP2(PIX) = D_ONE / AV_S_PQ(PIX)
              ! expect 0/0 := NaN, 0*Inf := NaN, MIN(x,NaN) == x
              CA_H_PQ(PIX) = MIN(RE_H_PQ(PIX) * DTMP1(PIX), CA_H_PQ(PIX))
              SA_H_PQ(PIX) = MIN(IM_H_PQ(PIX) * DTMP1(PIX), SA_H_PQ(PIX)) * DTMP3(PIX)
              CA_S_PQ(PIX) = MIN(RE_S_PQ(PIX) * DTMP2(PIX), CA_S_PQ(PIX))
              SA_S_PQ(PIX) = MIN(IM_S_PQ(PIX) * DTMP2(PIX), SA_S_PQ(PIX)) * DTMP4(PIX)
              ! compute the temps
              T(PIX) = SQRT(D_ONE - AV_S_PQ(PIX) * AV_S_PQ(PIX))
              U(PIX) = CA_S_PQ(PIX) * RE_H_PQ(PIX) + SA_S_PQ(PIX) * IM_H_PQ(PIX)
              V(PIX) = CA_S_PQ(PIX) * IM_H_PQ(PIX) - SA_S_PQ(PIX) * RE_H_PQ(PIX)
              E(PIX) = RE_H_QQ(PIX) - RE_H_PP(PIX)
              DHZ(PIX) = SIGN(D_ONE, E(PIX))
              ! compute fns of \gamma
              DTMP1(PIX) = SCALE(V(PIX) / E(PIX), 1)
              ! V==0 & E==0 ==> avoid NaN(TG), set TG=0 (expect MAX(x,NaN) == x)
              TG(PIX) = SIGN(MAX(ABS(DTMP1(PIX)), TG(PIX)), DTMP1(PIX))
              CG(PIX) = D_ONE / SQRT(D_ONE + TG(PIX) * TG(PIX))
              ! beware of Inf(TG)
              SG(PIX) = SIGN(MIN(TG(PIX) * CG(PIX), SG(PIX)), TG(PIX))
              ! compute fns of 2\vartheta
              T2T(PIX) = (DHZ(PIX) * (SCALE(U(PIX), 1) - (RE_H_PP(PIX) + RE_H_QQ(PIX)) * AV_S_PQ(PIX))) / &
                   (T(PIX) * SQRT(E(PIX) * E(PIX) + SCALE(V(PIX) * V(PIX), 2)))
              C2T(PIX) = D_ONE / SQRT(D_ONE + T2T(PIX) * T2T(PIX))
              S2T(PIX) = T2T(PIX) * C2T(PIX)
              DTMP1(PIX) = D_ONE + T(PIX) * C2T(PIX) * CG(PIX)
              DHZ(PIX) = AV_S_PQ(PIX) * S2T(PIX)
              DTMP2(PIX) = DTMP1(PIX) - DHZ(PIX)
              DTMP1(PIX) = DTMP1(PIX) + DHZ(PIX)
              ! compute the transformation
              CPHI(PIX) = SQRT(SCALE(DTMP1(PIX), -1))
              CPSI(PIX) = SQRT(SCALE(DTMP2(PIX), -1))
              ! for big/small rot
              DHZ(PIX) = (D_ONE - CPHI(PIX)) + (D_ONE - CPSI(PIX))
              DTMP3(PIX) = S2T(PIX) - AV_S_PQ(PIX)
              DTMP4(PIX) = T(PIX) * SG(PIX) * C2T(PIX)
              TG(PIX) = -SA_S_PQ(PIX) * DTMP4(PIX)
              T2T(PIX) = CA_S_PQ(PIX) * DTMP4(PIX)
              T(PIX) = D_ONE / T(PIX)
              RE_ASPHI(PIX) = CA_S_PQ(PIX) * DTMP3(PIX) + TG(PIX)
              IM_ASPHI(PIX) = T2T(PIX) + SA_S_PQ(PIX) * DTMP3(PIX)
              DTMP4(PIX) = (CPSI(PIX) * RE_S_PP(PIX) * T(PIX)) / DTMP2(PIX)
              RE_ASPHI(PIX) = RE_ASPHI(PIX) * DTMP4(PIX)
              IM_ASPHI(PIX) = IM_ASPHI(PIX) * DTMP4(PIX)
              DTMP3(PIX) = S2T(PIX) + AV_S_PQ(PIX)
              RE_MBSPSI(PIX) = CA_S_PQ(PIX) * DTMP3(PIX) + TG(PIX)
              IM_MBSPSI(PIX) = T2T(PIX) + SA_S_PQ(PIX) * DTMP3(PIX)
              DTMP4(PIX) = (CPHI(PIX) * RE_S_QQ(PIX) * T(PIX)) / DTMP1(PIX)
              RE_MBSPSI(PIX) = -RE_MBSPSI(PIX) * DTMP4(PIX)
              IM_MBSPSI(PIX) = IM_MBSPSI(PIX) * DTMP4(PIX)
              CPHI(PIX) = CPHI(PIX) * RE_S_PP(PIX) * T(PIX)
              CPSI(PIX) = CPSI(PIX) * RE_S_QQ(PIX) * T(PIX)
           END DO
           
           ! apply the transformations

           DO PIX = 1, PPV
              IF (HZ(PIX) .EQ. 0) CYCLE
              IF (DHZ(PIX) .GT. D_ZERO) SNROT(2) = SNROT(2) + 1
              ! ``global'' pair index
              PAIR = (VEC - 1) * PPV + PIX
              IF (PAIR .LE. NPAIRS) THEN
                 P = JSPAIR(1,PAIR,STEP)
                 Q = JSPAIR(2,PAIR,STEP)
                 ! ...transform...
                 IF (.NOT. (CPHI(PIX) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: F_11 overflow or NaN.'
                 DTMP1(PIX) = CPHI(PIX)
                 IF (.NOT. (ABS(RE_MBSPSI(PIX)) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: |Re(F_21)| overflow or NaN.'
                 IF (.NOT. (ABS(IM_MBSPSI(PIX)) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: |Im(F_21)| overflow or NaN.'
                 ZTMP1(PIX) = DCMPLX(RE_MBSPSI(PIX), IM_MBSPSI(PIX))
                 IF (.NOT. (ABS(RE_ASPHI(PIX)) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: |Re(F_12)| overflow or NaN.'
                 IF (.NOT. (ABS(IM_ASPHI(PIX)) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: |Im(F_12)| overflow or NaN.'
                 ZTMP2(PIX) = DCMPLX(RE_ASPHI(PIX), IM_ASPHI(PIX))
                 IF (.NOT. (CPSI(PIX) .LE. HUGE(D_ZERO))) STOP 'ZHZL1: F_22 overflow or NaN.'
                 DTMP2(PIX) = CPSI(PIX)
                 CALL BLAS_ZVROTM(M, H(1,P), H(1,Q), DTMP1(PIX), ZTMP1(PIX), ZTMP2(PIX), DTMP2(PIX))
                 CALL BLAS_ZVROTM(M, S(1,P), S(1,Q), DTMP1(PIX), ZTMP1(PIX), ZTMP2(PIX), DTMP2(PIX))
                 CALL BLAS_ZVROTM(N, Z(1,P), Z(1,Q), DTMP1(PIX), ZTMP1(PIX), ZTMP2(PIX), DTMP2(PIX))
              END IF
           END DO
        END DO
        !$OMP END PARALLEL DO
     END DO
     WRITE (ULOG,'(I3,A,I20,A,I20)') SWEEP,',',SNROT(1),',',SNROT(2)
     IF (SNROT(2) .EQ. 0) EXIT
     !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
     DO I = 1, 2
        NROT(I) = NROT(I) + SNROT(I)
     END DO
  END DO

  INFO = SWEEP

  ! Scaling of H,S,Z.

  IF (NROT(1) .GT. 0) THEN
     !$OMP  PARALLEL DO DEFAULT(NONE) NUM_THREADS(CPR) PROC_BIND(SPREAD) &
     !$OMP& SHARED(M,N,H,S,Z,JVEC) PRIVATE(I,J,L,P,DTOL,EE,EY,EW,SY,SW,SS,DTMP1,DTMP2,ZTMP1,ZTMP2)
     DO J = 1, N
        !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
        DTMP1 = D_ZERO

        DO I = 1, M, DSIMDL
           P = MIN(DSIMDL, M-(I-1))
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO L = 1, P
              ZTMP1(L) = H(I+(L-1),J)
              DTMP1(L) = DTMP1(L) + JVEC(M) * (DBLE(ZTMP1(L))*DBLE(ZTMP1(L)) + AIMAG(ZTMP1(L))*AIMAG(ZTMP1(L)))
           END DO
        END DO
        EY(J) = SUM(DTMP1)
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
        DTMP2 = D_ZERO

        DO I = 1, M, DSIMDL
           P = MIN(DSIMDL, M-(I-1))
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           DO L = 1, P
              ZTMP2(L) = S(I+(L-1),J)
              DTMP2(L) = DTMP2(L) + (DBLE(ZTMP2(L))*DBLE(ZTMP2(L)) + AIMAG(ZTMP2(L))*AIMAG(ZTMP2(L)))
           END DO
        END DO
        EW(J) = SUM(DTMP2)
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
              !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
              DO I = 1, N
                 Z(I,J) = Z(I,J) / DTOL
              END DO
           ELSE
              CALL ZDSCAL(N, DTOL, Z(1,J), 1)
           END IF
        END IF
     END DO
     !$OMP END PARALLEL DO
  END IF

#ifndef NDEBUG
  CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .TRUE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .TRUE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .TRUE._c_int)
#endif
END SUBROUTINE ZHZL1P
