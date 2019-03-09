! L1 double complex HZ.
SUBROUTINE ZHZL1(K, BH,NPLUS, BS,BZ, LDB, JS,JSPAIR, NSWP, NROT,INFO)
#ifndef NDEBUG
  USE, INTRINSIC :: IEEE_ARITHMETIC
  USE, INTRINSIC :: IEEE_FEATURES
#endif
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: K,NPLUS,LDB, JS(JSMLEX),JSPAIR(2,JS(JSMLEX),JS(JSMLEX-1)), NSWP
  DOUBLE COMPLEX, INTENT(INOUT) :: BH(LDB,K),BS(LDB,K),BZ(LDB,K)
  INTEGER, INTENT(OUT) :: NROT(2),INFO

  INTEGER :: NSTEPS, NPAIRS
  INTEGER :: PPV, VPS
  INTEGER :: SWEEP, STEP, VEC, PIX, PAIR
  INTEGER :: P, Q, I, J, L

  INTEGER :: SNROT(2)
  DOUBLE PRECISION :: DTOL, DSCL(3)
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

  EXTERNAL :: ZDSCAL

  !DIR$ ASSUME_ALIGNED JS:ALIGNB
  !DIR$ ASSUME_ALIGNED JSPAIR:ALIGNB
  !DIR$ ASSUME_ALIGNED BH:ALIGNB
  !DIR$ ASSUME_ALIGNED BS:ALIGNB
  !DIR$ ASSUME_ALIGNED BZ:ALIGNB
#ifdef NDEBUG
  !DIR$ ASSUME (K .GE. 0)
  !DIR$ ASSUME (MOD(K,2) .EQ. 0)
  !DIR$ ASSUME (NPLUS .GE. 0)
  !DIR$ ASSUME (NPLUS .LE. K)
  !DIR$ ASSUME (MOD(LDB,ZALIGN) .EQ. 0)
  !DIR$ ASSUME (NSWP .GE. 0)
#else
  IF (K .LT. 0) STOP 'ZHZL1: K .LT. 0'
  IF (MOD(K,2) .NE. 0) STOP 'ZHZL1: MOD(K,2) .NE. 0'
  IF (NPLUS .LT. 0) STOP 'ZHZL1: NPLUS .LT. 0'
  IF (NPLUS .GT. K) STOP 'ZHZL1: NPLUS .GT. K'
  IF (MOD(LDB,ZALIGN) .NE. 0) STOP 'ZHZL1: MOD(LDB,ZALIGN) .NE. 0'
  IF (NSWP .LT. 0) STOP 'ZHZL1: NSWP .LT. 0'
#endif

  INFO = 0
  !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
  NROT = 0
  IF (K .LE. 0) RETURN

#ifndef NDEBUG
  CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .FALSE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .FALSE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE._c_int)
#endif

  DTOL = SCALE(SQRT(DBLE(K)), -53)

  NSTEPS = JS(JSMLEX-1)
  NPAIRS = JS(JSMLEX)
  ! pairs per vector
  PPV = MIN(NPAIRS, DSIMDL)
  ! vectors per step
  VPS = (NPAIRS + (PPV - 1)) / PPV

  DO SWEEP = 1, NSWP
     !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
     SNROT = 0
     DO STEP = 1, NSTEPS
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
           SG = D_ONE
           !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
           S2T = D_ONE

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

                 DO I = 1, K, DSIMDL
                    L = MIN(DSIMDL, K-(I-1))
                    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                    DO J = 1, L
                       ZTMP1(J) = BS(I+(J-1),P)
                       ZTMP2(J) = BS(I+(J-1),Q)
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
                    ! A joint prescaling of BH and BS needed...
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
                    ! A joint prescaling of BH and BS needed...
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

                 DO I = 1, NPLUS, DSIMDL
                    L = MIN(DSIMDL, NPLUS-(I-1))
                    !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
                    DO J = 1, L
                       ZTMP1(J) = BH(I+(J-1),P)
                       ZTMP2(J) = BH(I+(J-1),Q)
                       !DBLE(DCONJG(ZTMP1(J))*ZTMP1(J))
                       DTMP1(J) = DTMP1(J) + (DBLE(ZTMP1(J))*DBLE(ZTMP1(J)) + AIMAG(ZTMP1(J))*AIMAG(ZTMP1(J)))
                       !DBLE(DCONJG(ZTMP2(J))*ZTMP2(J))
                       DTMP2(J) = DTMP2(J) + (DBLE(ZTMP2(J))*DBLE(ZTMP2(J)) + AIMAG(ZTMP2(J))*AIMAG(ZTMP2(J)))
                       ! += DCONJG(ZTMP1(J)) * ZTMP2(J)
                       DTMP3(J) = DTMP3(J) + (DBLE(ZTMP1(J))*DBLE(ZTMP2(J)) + AIMAG(ZTMP1(J))*AIMAG(ZTMP2(J)))
                       DTMP4(J) = DTMP4(J) + (DBLE(ZTMP1(J))*AIMAG(ZTMP2(J))- AIMAG(ZTMP1(J))*DBLE(ZTMP2(J)))
                    END DO
                 END DO

                 DO I = NPLUS+1, K, DSIMDL
                    L = MIN(DSIMDL, K-(I-1))
                    !DIR$ VECTOR ALWAYS ASSERT
                    DO J = 1, L
                       ZTMP1(J) = BH(I+(J-1),P)
                       ZTMP2(J) = BH(I+(J-1),Q)
                       !DBLE(DCONJG(ZTMP1(J))*ZTMP1(J))
                       DTMP1(J) = DTMP1(J) - (DBLE(ZTMP1(J))*DBLE(ZTMP1(J)) + AIMAG(ZTMP1(J))*AIMAG(ZTMP1(J)))
                       !DBLE(DCONJG(ZTMP2(J))*ZTMP2(J))
                       DTMP2(J) = DTMP2(J) - (DBLE(ZTMP2(J))*DBLE(ZTMP2(J)) + AIMAG(ZTMP2(J))*AIMAG(ZTMP2(J)))
                       ! -= DCONJG(ZTMP1(J)) * ZTMP2(J)
                       DTMP3(J) = DTMP3(J) - (DBLE(ZTMP1(J))*DBLE(ZTMP2(J)) + AIMAG(ZTMP1(J))*AIMAG(ZTMP2(J)))
                       DTMP4(J) = DTMP4(J) - (DBLE(ZTMP1(J))*AIMAG(ZTMP2(J))- AIMAG(ZTMP1(J))*DBLE(ZTMP2(J)))
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
                    ! A joint prescaling of BH and BS needed...
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
                    ! A joint prescaling of BH and BS needed...
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

           J = 0
           DO PIX = 1, DSIMDL
              IF (PIX .GT. PPV) THEN
                 HZ(PIX) = 0
              ELSE
                 J = J + HZ(PIX)
              END IF
           END DO
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
              ! V==0 & E==0 ==> NaN(TG)
              DTMP1(PIX) = D_ONE - MAX(V(PIX) / V(PIX), D_ZERO)
              DTMP2(PIX) = D_ONE - MAX(E(PIX) / E(PIX), D_ZERO)
              HZ(PIX) = HZ(PIX) + INT(SCALE(DTMP1(PIX) * DTMP2(PIX), 1))
              ! compute fns of \gamma
              TG(PIX) = SCALE(V(PIX) / E(PIX), 1)
              CG(PIX) = D_ONE / SQRT(D_ONE + TG(PIX) * TG(PIX))
              ! beware of Inf(TG), expect MIN(x,NaN) == x
              SG(PIX) = SIGN(MIN(TG(PIX) * CG(PIX), SG(PIX)), TG(PIX))
              ! compute fns of 2\vartheta
              DHZ(PIX) = SIGN(D_ONE, E(PIX))
              T2T(PIX) = (DHZ(PIX) * (SCALE(U(PIX), 1) - (RE_H_PP(PIX) + RE_H_QQ(PIX)) * AV_S_PQ(PIX))) / &
                   (T(PIX) * SQRT(E(PIX) * E(PIX) + SCALE(V(PIX) * V(PIX), 2)))
              C2T(PIX) = D_ONE / SQRT(D_ONE + T2T(PIX) * T2T(PIX))
              S2T(PIX) = SIGN(MIN(T2T(PIX) * C2T(PIX), S2T(PIX)), T2T(PIX))
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
              IF (MOD(HZ(PIX),2) .EQ. 0) CYCLE
              ! ``global'' pair index
              PAIR = (VEC - 1) * PPV + PIX
              IF (PAIR .LE. NPAIRS) THEN
                 P = JSPAIR(1,PAIR,STEP)
                 Q = JSPAIR(2,PAIR,STEP)
                 ! ...transform...
                 IF (HZ(PIX) .EQ. 3) THEN
                    CPHI(PIX) = D_CS_PI_4 * RE_S_PP(PIX)
                    RE_MBSPSI(PIX) = CA_S_PQ(PIX) * D_CS_PI_4
                    IM_MBSPSI(PIX) = -SA_S_PQ(PIX) * D_CS_PI_4
                    RE_ASPHI(PIX) = -RE_MBSPSI(PIX) * RE_S_PP(PIX)
                    IM_ASPHI(PIX) = IM_MBSPSI(PIX) * RE_S_PP(PIX)
                    RE_MBSPSI(PIX) = RE_MBSPSI(PIX) * RE_S_QQ(PIX)
                    IM_MBSPSI(PIX) = IM_MBSPSI(PIX) * RE_S_QQ(PIX)
                    CPSI(PIX) = D_CS_PI_4 * RE_S_QQ(PIX)

                    DTMP1(PIX) = D_ONE / SQRT(D_ONE + AV_S_PQ(PIX))
                    DTMP2(PIX) = D_ONE / SQRT(D_ONE - AV_S_PQ(PIX))
                    CPHI(PIX) = CPHI(PIX) * DTMP1(PIX)
                    RE_MBSPSI(PIX) = RE_MBSPSI(PIX) * DTMP1(PIX)
                    IM_MBSPSI(PIX) = IM_MBSPSI(PIX) * DTMP1(PIX)
                    RE_ASPHI(PIX) = RE_ASPHI(PIX) * DTMP2(PIX)
                    IM_ASPHI(PIX) = IM_ASPHI(PIX) * DTMP2(PIX)
                    CPSI(PIX) = CPSI(PIX) * DTMP2(PIX)

                    DHZ(PIX) = D_TWO - SQRT(D_TWO)
                 END IF
                 ! form the transformation
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
                 IF (DHZ(PIX) .GT. D_ZERO) SNROT(2) = SNROT(2) + 1
                 ! apply the transformation
                 IF ((DTMP1(PIX) .NE. D_ONE) .OR. (DTMP2(PIX) .NE. D_ONE) .OR. &
                      (ZTMP1(PIX) .NE. Z_ZERO) .OR. (ZTMP2(PIX) .NE. Z_ZERO)) THEN
                    CALL BLAS_ZVROTM(K, BH(1,P), BH(1,Q), DTMP1(PIX), ZTMP1(PIX), ZTMP2(PIX), DTMP2(PIX))
                    CALL BLAS_ZVROTM(K, BS(1,P), BS(1,Q), DTMP1(PIX), ZTMP1(PIX), ZTMP2(PIX), DTMP2(PIX))
                    CALL BLAS_ZVROTM(K, BZ(1,P), BZ(1,Q), DTMP1(PIX), ZTMP1(PIX), ZTMP2(PIX), DTMP2(PIX))
                 ELSE ! identity
                    SNROT(1) = SNROT(1) - 1
                 END IF
              END IF
           END DO
        END DO
     END DO
     IF (SNROT(1) .EQ. 0) EXIT
     !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
     DO I = 1, 2
        NROT(I) = NROT(I) + SNROT(I)
     END DO
  END DO

  INFO = SWEEP

  ! Scaling of Z.

  DO J = 1, K
     ! compute the J-norm
     !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
     DTMP1 = D_ZERO

     DO I = 1, NPLUS, DSIMDL
        P = MIN(DSIMDL, NPLUS-(I-1))
        !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
        DO L = 1, P
           ZTMP1(L) = BH(I+(L-1),J)
           !DBLE(DCONJG(ZTMP1(L))*ZTMP1(L))
           DTMP1(L) = DTMP1(L) + (DBLE(ZTMP1(L))*DBLE(ZTMP1(L)) + AIMAG(ZTMP1(L))*AIMAG(ZTMP1(L)))
        END DO
     END DO
     DO I = NPLUS+1, K, DSIMDL
        P = MIN(DSIMDL, K-(I-1))
        !DIR$ VECTOR ALWAYS ASSERT
        DO L = 1, P
           ZTMP1(L) = BH(I+(L-1),J)
           !DBLE(DCONJG(ZTMP1(L))*ZTMP1(L))
           DTMP1(L) = DTMP1(L) - (DBLE(ZTMP1(L))*DBLE(ZTMP1(L)) + AIMAG(ZTMP1(L))*AIMAG(ZTMP1(L)))
        END DO
     END DO
     DTOL = SUM(DTMP1)
     IF (DTOL .NE. D_ZERO) THEN
        DTOL = ABS(DTOL)
        IF (DTOL .NE. D_ONE) THEN
           DSCL(1) = SQRT(DTOL)
        ELSE
           DSCL(1) = D_ONE
        END IF
     ELSE
        DSCL(1) = D_ZERO
     END IF

     ! compute the norm
     !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
     DTMP2 = D_ZERO

     DO I = 1, K, DSIMDL
        P = MIN(DSIMDL, K-(I-1))
        !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
        DO L = 1, P
           ZTMP2(L) = BS(I+(L-1),J)
           !DBLE(DCONJG(ZTMP2(L))*ZTMP2(L))
           DTMP2(L) = DTMP2(L) + (DBLE(ZTMP2(L))*DBLE(ZTMP2(L)) + AIMAG(ZTMP2(L))*AIMAG(ZTMP2(L)))
        END DO
     END DO
     DTOL = SUM(DTMP2)
     IF (DTOL .NE. D_ZERO) THEN
        IF (DTOL .NE. D_ONE) THEN
           DSCL(2) = SQRT(DTOL)
        ELSE
           DSCL(2) = D_ONE
        END IF
     ELSE
        DSCL(2) = D_ZERO
     END IF

     DSCL(3) = HYPOT(DSCL(1), DSCL(2))
     IF (DSCL(3) .NE. D_ONE) THEN
        ! underflow
        IF (DSCL(3) .LT. TINY(D_ZERO)) STOP 'ZHZL1: Scale of Z underflows.'
        ! overflow
        IF (DSCL(3) .GT. HUGE(D_ZERO)) STOP 'ZHZL1: Scale of Z overflows.'
        DTOL = D_ONE / DSCL(3)
        IF (DTOL .LT. TINY(D_ZERO)) THEN
           ! underflow
           DTOL = DSCL(3)
           DO I = 1, K, DSIMDL
              P = MIN(DSIMDL, K-(I-1))
              !DIR$ VECTOR ALWAYS ASSERT,ALIGNED
              DO L = 1, P
                 BZ(I+(L-1),J) = BZ(I+(L-1),J) / DTOL
              END DO
           END DO
        ELSE
           CALL ZDSCAL(K, DTOL, BZ(1,J), 1)
        END IF
     END IF
  END DO

#ifndef NDEBUG
  CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .TRUE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .TRUE._c_int)
  CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .TRUE._c_int)
#endif
END SUBROUTINE ZHZL1
