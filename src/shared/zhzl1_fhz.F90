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
