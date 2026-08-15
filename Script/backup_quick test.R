# back up Quick test

# Create function to calculate the final parameters
Calc_final_params <- function(parm.chem, parm.phys){
  with(as.list(c(parm.chem, parm.phys)), {
    # Physiological ----
    
    VIL <- VILc * BW                # L, Volume of intestinal lumen
    VI <- VIc * BW                  # L, Volume of intestine
    SA_SI <- SA_SIc * BW
    
    VL <- VLc * BW                  # L, Volume of liver
    VL_ic <- VL_icc*VL              # L, Volume liver intracellular                 
    VL_ec <- VL_ecc*VL              # L, Volume liver extracellular
    
    VK <- VKc*BW                    # L, Volume of kidney
    VPT <- VPTc*VK                  # L, Volume of proximal tubule
    VPTT <- VPT*VPTTc               # L, Volume of proximal tubule tissue
    VPTL <- VPT*VPTLc               # L, Volume of proximal tubule lumen
    VRKL <- (VK - VPT)*VRKLc        # L, Volume of rest of kidney lumen
    VRKT <- VK - VPTT - VPTL - VRKL # L, Volume of rest of kidney tissue, used in the model
    
    VA <- VAc * BW/0.9 + MAc/0.9    # L, Volume of adipose
    
    VP <- VPc * (1-Hct) * BW        # L, Volume of plasma
    
    VTotc <- 0.96
    VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
    
    VR <- VTot - (VI + VL + VK + VA + VP)             # L, Volume of the lumped rest compartment
    
    QTotc <- 0.988
    QI <- QIc/QTotc * QC            # L/d, Intestinal
    QL <- QLc/QTotc * QC            # L/d, Liver
    QK <- QKc/QTotc * QC            # L/d, Kidney
    QA <- QAc/QTotc * QC            # L/d, Adipose 
    QR <- QC - (QA + QI + QK + QL)  # L/d, Rest
    
    QUr <- QUrc * BW                # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
    GFR <- GFR
    QT <- QTc*GFR                   # L/d, Proximal tubule fluid flow
    
    
    # Toxicokinetic ----
    
    # affinities towards tissue constituents
    
    k_ML <- logD_MembraneLipidWater      # neutral phospholipids:water partition coefficient
    k_SP <- logD_StructuralProteinWater  # structural protein:water partition coefficient
    k_ALB <- logD_AlbuminWater           # albumin:water partition coefficient
    k_SL <- logD_StorageLipidWater       # structural lipids:water partition coefficient
    k_GLOB <- logD_GlobulinWater         # globulin:water partition coefficient
    
    
    # albumin concentrations in different matrices; Akihiro Tojo and Satoshi Kinugasa 2012 doi:10.1155/2012/481520
    Calb_P <- 43.74214         # g/L plasma, calculated based on lifestage equations
    Calb_PTL <- 14.4e-3        # g/L proximal tubule
    Calb_IL <- 0.007e-3        # g/L intestinal lumen
    Calb_L_ec <- 0.0012e-3     # g/L liver extracellular 
    
    # calculate albumin ratio; based on the Poulin and Theil 2009, below Table 6
    R_PTL <- Calb_P/Calb_PTL   # plasma to proximal tubule lumen albumin ratio
    R_L_ec <- Calb_P/Calb_L_ec # plasma to liver albumin ratio
    R_IL <- Calb_P/Calb_IL     # plasma to intestinal lumen albumin ratio
    
    
    ## fraction unbound ----
    # calculate fraction unbound; based on the equation from https://doi.org/10.1021/acs.est.3c07415
    
    DAlb <- 1.36                 # kg/L density of albumin
    D_ALB <- k_ALB*DAlb          # L_w/L_hsa; distribution coefficient
    FV_Alb <- Calb_P*10^-3         # kg/kg of albumin in serum
    
    DGlob <- 1.20                # kg/L density of globulin
    D_Glob <- k_GLOB*DGlob       # L_w/L_hsa; distribution coefficient
    FV_Glob <- 0.031             # kg/kg of globulin in serum
    
    FV_W <- 1 - FV_Alb - FV_Glob # serum mass balance
    
    # fraction unbound in plasma
    fup <- 1/(1 + (D_ALB*(FV_Alb/FV_W)) + (D_Glob*(FV_Glob/FV_W)))
    
    # calculate fraction unbound in different matrices; based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012; adapted to not account for fraction unionised as OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
    fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup))   # proximal tubule lumen
    fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # liver intracellular space
    fu_In  <- fu_Lic                            # intestine intracellular space
    
    
    ## partition coefficients ---- 
    #  code is a slightly modified version from Utsey et al. 2020 https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_Schmitt.R
    
    Physio.Tissues <- Physio.data %>% filter(!Tissue %in% c("Blood", "Plasma")) 
    
    Kp <- (Physio.Tissues$f_W +
             (k_ML*Physio.Tissues$f_ML) +
             (k_SL*Physio.Tissues$f_SL) +
             (k_SP*Physio.Tissues$f_SP) +
             (k_ALB*Physio.Tissues$f_ALB) 
    ) * fup 
    
    Names <- Physio.Tissues$Tissue %>% substr(1,2)
    Names <- paste("Kp",Names ,sep="")
    Kp.df <- as.data.frame(list(Tissue = Names)) %>%
      mutate(Value = Kp) %>% 
      pivot_wider(names_from = Tissue, values_from = Value)
    
    
    # calculate Plasma/Rest of the body partition coefficient
    # KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
    V_skinFraction <- 0.038
    V_lungFraction <- 0.007
    V_stomachFraction <- 0.002
    V_brainFraction <- 0.022
    V_heartFraction <- 0.004
    V_muscleFraction <- 0.292
    V_spleenFraction <- 0.002
    V_reproFraction <- 0.002
    V_boneFraction <- 0.029
    
    KpRe <- (Kp.df$KpSk * V_skinFraction +
               Kp.df$KpLu * V_lungFraction +
               Kp.df$KpSt * V_stomachFraction +
               Kp.df$KpBr * V_brainFraction +
               Kp.df$KpHe * V_heartFraction +
               Kp.df$KpMu * V_muscleFraction +
               Kp.df$KpSp * V_spleenFraction +
               Kp.df$KpGo * V_reproFraction +
               Kp.df$KpBo * V_boneFraction) / (
                 V_skinFraction +
                   V_lungFraction +
                   V_stomachFraction +
                   V_brainFraction +
                   V_heartFraction +
                   V_muscleFraction +
                   V_spleenFraction +
                   V_reproFraction +
                   V_boneFraction)  
    Kp.df$KpRe <- KpRe
    
    # final partition coefficients
    PI <- Kp.df$KpIn   # Intestine
    PL <- Kp.df$KpLi   # Liver
    PK <- Kp.df$KpKi   # Kidney
    PA <- Kp.df$KpAd   # Adipose
    PR <- Kp.df$KpRe   # Rest
    
    ## clearances ----
    
    CL_IL <- (Papp_SI*SA_SI*1e-3)    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/d -> L/d /1000)
    
    CL_OATP2B1_I <- CL_OATP2B1c*SF_OATP2B1_I*VI # L/d
    CF_I <- QI/CL_OATP2B1_I
    CF_I <- ifelse(CF_I < 1, CF_I, 1)
    CL_OATP2B1_I <- CL_OATP2B1_I * CF_I 
    
    CL_OATP1B1 <- CL_OATP1B1c*SF_OATP1B1*VL_ec   # L/d
    CL_OATP1B3 <- CL_OATP1B3c*SF_OATP1B3*VL_ec   # L/d
    CL_OATP2B1_L <- CL_OATP2B1c*SF_OATP2B1*VL_ec # L/d
    CL_NTCP <- CL_NTCPc*SF_NTCP*VL_ec            # L/d    
    CF_L_in <- (QI + QL)/(CL_OATP1B1 + CL_OATP1B3 + CL_OATP2B1_L + CL_NTCP)
    CF_L_in <- ifelse(CF_L_in < 1, CF_L_in, 1)
    
    CL_OATP1B1 <- CL_OATP1B1 * CF_L_in
    CL_OATP1B3 <- CL_OATP1B3 * CF_L_in
    CL_OATP2B1_L <- CL_OATP2B1_L * CF_L_in
    CL_NTCP <- CL_NTCP * CF_L_in
    
    CL_BSEP <- CL_BSEPc*SF_BSEP*VL_ic  # L/d
    CF_L_out <- (QI + QL)/(CL_BSEP)
    CF_L_out <- ifelse(CF_L_out < 1, CF_L_out, 1)
    CL_BSEP <- CL_BSEP * CF_L_out
    
    CL_OAT3 <- CL_OAT3c*SF_OAT3*VPT  # L/d 
    CL_OAT1 <- CL_OAT1c*SF_OAT1*VPT  # L/d 
    CF_K_in <- QK/(CL_OAT1 + CL_OAT3)
    CF_K_in <- ifelse(CF_K_in < 1, CF_K_in, 1)
    CL_OAT1 <- CL_OAT1 * CF_K_in
    CL_OAT3 <- CL_OAT3 * CF_K_in
    
    CL_OAT4 <- CL_OAT4c*SF_OAT4*VPT  # L/d 
    CF_K_out <- GFR/CL_OAT4
    CF_K_out <- ifelse(CF_K_out < 1, CF_K_out, 1)
    CL_OAT4 <- CL_OAT4 * CF_K_out
    
    # Clearance via menstruation
    CL_menses <- CL_menses  # L/d
    
    
    return(list(
      BW = BW,
      VIL = VIL,
      VI = VI,
      SA_SI = SA_SI, 
      VL = VL,
      VL_ic = VL_ic,
      VL_ec = VL_ec,
      VK = VK,
      VPT = VPT,
      VPTT = VPTT,
      VPTL = VPTL, 
      VRKL = VRKL, 
      VRKT = VRKT,
      VA = VA,
      VP = VP,
      VR = VR, 
      QC = QC,
      QI = QI,
      QL = QL,
      QK = QK,
      QA = QA,
      QR = QR,
      QUr = QUr,
      GFR = GFR,
      QT = QT,
      tco = tco,
      fup = fup,
      fu_PTL = fu_PTL,
      fu_Lic = fu_Lic,
      fu_In = fu_In,
      PI = PI,
      PL = PL,
      PK = PK,
      PA = PA,
      PR = PR,
      CF_I = CF_I,
      CF_L_in = CF_L_in, 
      CF_L_out = CF_L_out,
      CF_K_in = CF_K_in, 
      CF_K_out = CF_K_out,
      CL_IL = CL_IL,
      CL_OATP2B1_I = CL_OATP2B1_I,
      CL_OATP1B1 = CL_OATP1B1,
      CL_OATP1B3 = CL_OATP1B3,
      CL_OATP2B1_L = CL_OATP2B1_L,
      CL_NTCP = CL_NTCP,
      CL_BSEP = CL_BSEP,
      CL_OAT3 = CL_OAT3,
      CL_OAT1 = CL_OAT1,
      CL_OAT4 = CL_OAT4,
      CL_menses = CL_menses
    ))
  })
}

# Create function to run the PBK model
PBK.model <- function(t, state, pars){
  with(as.list(c(state, pars)), {
    
    ## Concentrations -------------------------
    
    CIL <- AIL/VIL               # ug/L, Intestinal lumen
    CI <- AI/VI                  # ug/L, Intestine 
    CVI <- CI/PI                 # ug/L, Intestine venous
    
    CL_ec <- AL_ec/VL_ec         # ug/L, Liver extracellular
    CL_ic <- AL_ic/VL_ic         # ug/L, Liver intracellular
    CVL_ec <- CL_ec/PL           # ug/L, Liver extracellular venous
    
    CPTT <- APTT/VPTT            # ug/L, Kidney, proximal tubule tissue
    CPTL <- APTL/VPTL            # ug/L, Kidney, proximal tubule lumen
    CRKT <- ARKT/VRKT            # ug/L, Rest of kidney tissue
    CRKL <- ARKL/VRKL            # ug/L, Rest of kidney lumen
    CVRKT <- CRKT/PK             # ug/L, Rest of kidney venous
    
    CA <- AA/VA                  # ug/L, Adipose
    CVA <- CA/PA                 # ug/L, Adipose venous
    
    CR <- AR/VR                  # ug/L, Rest
    CVR <- CR/PR                 # ug/L, Rest venous
    
    CP <- AP/VP                  # ug/L, Plasma
    
    
    ## Differential equations -------------------------
    
    dAIL <- - tco*AIL - CL_IL*CIL +
      - CL_OATP2B1_I*CL_ic*fu_Lic +
      + CL_BSEP*CL_ic*fu_Lic          # ug/d, Intestine lumen
    
    dAI <- QI*(CP - CVI) + CL_IL*CIL + 
      + CL_OATP2B1_I*CL_ic*fu_Lic #CL_OATP2B1_I*CI*fu_In                     # ug/d, Intestinal
    
    dAFe <-  tco*AIL                                       # ug/d, Feces
    
    
    dAL_ec <- + QI*CVI + QL*CP - (QI+QL)*CVL_ec + 
      - CL_OATP2B1_L*CL_ec*fup +
      - CL_OATP1B1*CL_ec*fup +
      - CL_OATP1B3*CL_ec*fup +          # ug/d, Liver extracellular space (vascular + interstitial space)
      - CL_NTCP*CL_ec*fup
    
    dAL_ic <- CL_OATP2B1_L*CL_ec*fup +
      + CL_OATP1B1*CL_ec*fup +
      + CL_OATP1B3*CL_ec*fup +
      + CL_NTCP*CL_ec*fup +
      - CL_BSEP*CL_ic*fu_Lic             # ug/d, Liver intracellular space
    
    
    dAPTT <- QK*(CP - CPTT) + 
      + CL_OAT1*fup*CP +
      + CL_OAT3*fup*CP +
      + CL_OAT4*CPTL*fu_PTL # ug/d, Proximal tubule tissue 
    
    dAPTL <- + fup*GFR*CP - QT*CPTL +
      - CL_OAT4*CPTL*fu_PTL   # ug/d, Proximal tubule lumen    
    
    dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
    
    dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
    
    dAUr <- QUr*CRKL                                         # ug/d, Urine
    
    
    dAA <- QA*(CP-CVA)                                      # ug/d, Adipose
    
    
    dAR <- QR*(CP-CVR)                                      # ug/d, Rest
    
    dAP <- - (QI + QL + QA + QR + QK)*CP - fup*GFR*CP +     # ug/d, Arterial Plasma
      + (QL+QI)*CVL_ec + QK*CVRKT + QA*CVA + QR*CVR +
      - CL_OAT1*fup*CP - CL_OAT3*fup*CP +
      - CL_menses*CP # ug/d, Venous Plasma
    
    dAMp <- + CL_menses*CP
    
    # Mass Balance
    Atot <- AIL +
      AI + AFe + 
      AL_ec + AL_ic +
      APTT + APTL + ARKT + ARKL + AUr +
      AA + 
      AR +
      AP +
      AMp
    
    MB <- 0.1 - Atot + 1    
    
    # End
    
    list(c(dAIL,
           dAI, 
           dAFe,
           dAL_ec,
           dAL_ic,
           dAPTT,
           dAPTL,
           dARKT,
           dARKL,
           dAUr,
           dAA,
           dAR, 
           dAP,
           dAMp
    ), 
    c(CIL = CIL,
      CI = CI, CVI = CVI,
      CL_ec = CL_ec,
      CL_ic = CL_ic,
      CVL_ec = CVL_ec,
      CPTT = CPTT,
      CPTL = CPTL,
      CRKT = CRKT,
      CRKL = CRKL,
      CA = CA, CVA = CVA,
      CR = CR, CVR = CVR,
      CP = CP
    ))
  })
}