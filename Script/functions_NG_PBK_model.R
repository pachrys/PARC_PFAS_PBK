# File containing model functions

# function to calculate the final parameters for Morris test (combines physchem and physiological parameters)
Morris_Calc_final_params <- function(parms){
  with(as.list(c(parms)), {
    # Physiological ----
    
    VI_L <- VILc * BW                # L, Volume of intestinal lumen
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
    QK <- QKc/QTotc * QC             # L/d, Kidney
    QA <- QAc/QTotc * QC             # L/d, Adipose 
    QR <- QC - (QA + QI + QK + QL)   # L/d, Rest
    
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
    Calb_P <- 40         # (35-50) g/L plasma
    Calb_PTL <- 0.84     #  (0.01-1.12) g/L proximal tubule lumen
    Calb_L_ec <- 18.60      #18.60   # (12-40) g/L liver extracellular  
    Calb_PT_ec <- 18.60     #18.60   # (12-40) g/L liver extracellular  
    
    
    # calculate albumin ratio; based on the Poulin and Theil 2009, below Table 6
    PMR_PTL <- Calb_P/Calb_PTL   # plasma to proximal tubule lumen albumin ratio
    PMR_L_ec <- Calb_P/Calb_L_ec # plasma to liver albumin ratio
    PMR_PT_ec <- Calb_P/Calb_PT_ec # plasma to liver albumin ratio
    
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
    
    # calculate fraction unbound in different matrices; 
    # based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012; adapted to not account for fraction unionised as OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
    fu_PTL <- PMR_PTL*fup/(1 + ((PMR_PTL-1)*fup))   # proximal tubule lumen
    fu_L_ec <- PMR_L_ec*fup/(1 + ((PMR_L_ec-1)*fup)) # liver extracellular space
    fu_PT_ec <- PMR_PT_ec*fup/(1 + ((PMR_PT_ec-1)*fup)) # liver extracellular space
    
    
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
    
    # intestine
    Clear_apparent <- (Papp_SI*SA_SI*1e-3)     # L/d, passive absorption (1e-3 to convert (cm^3)ml to L)
    
    Clear_OATP2B1_I <- Clear_OATP2B1c*SF_OATP2B1_I*VI_L  # L/d, active uptake, it is assume0d that the fraction unbound is 1 
    
    CF_I <- QI/(Clear_apparent+Clear_OATP2B1_I)
    CF_I_corrected <- ifelse(CF_I < 1, CF_I, 1) # correction factor in case the active clearance is higher than the blood flow
    
    Clear_uptake_I <- (Clear_apparent+Clear_OATP2B1_I) #* CF_I_corrected
    
    # liver
    Clear_OATP1B1 <- Clear_OATP1B1c*SF_OATP1B1*VL_ic*fu_L_ec   # L/d
    Clear_OATP1B3 <- Clear_OATP1B3c*SF_OATP1B3*VL_ic*fu_L_ec   # L/d
    Clear_OATP2B1_L <- Clear_OATP2B1c*SF_OATP2B1*VL_ic*fu_L_ec  # L/d
    Clear_NTCP <- Clear_NTCPc*SF_NTCP*VL_ic*fu_L_ec            # L/d    
    Clear_BSEP <- Clear_BSEPc*SF_BSEP*VL_ic  # L/d, it is assumed that fraction unbound is 1 as binding proteins (albumin in this case) are in the interstitium so not in the cells
    
    CF_L_uptake <- (QI + QL)/(Clear_OATP1B1 + 
                                Clear_OATP1B3 + 
                                Clear_OATP2B1_L + 
                                Clear_NTCP)
    CF_L_uptake_corrected <- ifelse(CF_L_uptake < 1, CF_L_uptake, 1) # correction factor in case the active clearance is higher than the blood flow
    
    Clear_uptake_L <- (Clear_OATP1B1 +
                         Clear_OATP1B3 + 
                         Clear_OATP2B1_L +
                         Clear_NTCP) #* CF_L_uptake_corrected
    
    CF_L_excretion <- (QI + QL)/(Clear_BSEP)
    CF_L_excretion_corrected <- ifelse(CF_L_excretion < 1, CF_L_excretion, 1) # correction factor in case the active clearance is higher than the blood flow
    
    Clear_excretion_L <- Clear_BSEP #* CF_L_excretion_corrected
    
    # kidney
    Clear_OAT3 <- Clear_OAT3c*SF_OAT3*VPT*fu_PT_ec  # L/d 
    Clear_OAT1 <- Clear_OAT1c*SF_OAT1*VPT*fu_PT_ec  # L/d 
    Clear_OAT4 <- Clear_OAT4c*SF_OAT4*VPT*fu_PTL  # L/d 
    
    # CF_K_excretion <- QK/(Clear_OAT1 + Clear_OAT3)
    # CF_K_excretion_corrected <- ifelse(CF_K_excretion < 1, CF_K_excretion, 1)
    # Clear_K_excretion <- CF_K_excretion_corrected * (Clear_OAT1 + Clear_OAT3)
    
    # calculate net reabsorption (difference between active elimination and reabsorption)
    Clear_net_reabsorption_unscaled <- Clear_OAT4 - (Clear_OAT1 + Clear_OAT3)
    
    CF_K_reabsorption <- QK/Clear_net_reabsorption_unscaled
    CF_K_reabsorption_corrected <- ifelse(CF_K_reabsorption < 1, CF_K_reabsorption, 1)
    
    Clear_net_reabsorption <- Clear_net_reabsorption_unscaled #* CF_K_reabsorption_corrected # this already corrects for the fraction unbound in the proximal tubule lumen 
    
    
    
    # Clearance via menstruation
    Qmenses <- Qmenses  # L/d
    
    
    return(list(
      BW = BW, 
      VI_L = VI_L,
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
      fu_L_ec = fu_L_ec,
      fu_PT_ec = fu_PT_ec,
      
      PI = PI,
      PL = PL,
      PK = PK,
      PA = PA,
      PR = PR,
      
      SF_OATP2B1_I = SF_OATP2B1_I,
      SF_OATP2B1 = SF_OATP2B1,
      SF_OATP1B3 = SF_OATP1B3,
      SF_OATP1B1 = SF_OATP1B1,
      SF_NTCP = SF_NTCP,
      SF_BSEP = SF_BSEP,
      SF_OAT1 = SF_OAT1,
      SF_OAT3 = SF_OAT3,
      SF_OAT4 = SF_OAT4,
      
      Clear_apparent = Clear_apparent,
      Clear_OATP2B1_I = Clear_OATP2B1_I,
      CF_I = CF_I_corrected,
      Clear_uptake_I = Clear_uptake_I,
      
      Clear_OATP1B1 = Clear_OATP1B1,
      Clear_OATP1B3 = Clear_OATP1B3, 
      Clear_OATP2B1_L = Clear_OATP2B1_L,
      Clear_NTCP = Clear_NTCP,
      Clear_BSEP = Clear_BSEP,
      CF_L_uptake = CF_L_uptake_corrected,
      Clear_uptake_L = Clear_uptake_L,
      CF_L_excretion = CF_L_excretion_corrected,
      Clear_excretion_L = Clear_excretion_L,
      
      Clear_OAT3 = Clear_OAT3,
      Clear_OAT1 = Clear_OAT1,
      Clear_OAT4 = Clear_OAT4,
      Clear_net_reabsorption_unscaled = Clear_net_reabsorption_unscaled,
      CF_K_reabsorption = CF_K_reabsorption_corrected,
      Clear_net_reabsorption = Clear_net_reabsorption,
      
      Qmenses = Qmenses
    ))
  })
}

# function to calculate the final parameters
Calc_final_params <- function(parm.chem, parm.phys){
  with(as.list(c(parm.chem, parm.phys)), {
    # Physiological ----
    
    VI_L <- VILc * BW                # L, Volume of intestinal lumen
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
    QK <- QKc/QTotc * QC             # L/d, Kidney
    QA <- QAc/QTotc * QC             # L/d, Adipose 
    QR <- QC - (QA + QI + QK + QL)   # L/d, Rest
    
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
    Calb_P <- 40         # (35-50) g/L plasma
    Calb_PTL <- 0.84     #  (0.01-1.12) g/L proximal tubule lumen
    Calb_L_ec <- 18.60      #18.60   # (12-40) g/L liver extracellular  
    Calb_PT_ec <- 18.60     #18.60   # (12-40) g/L liver extracellular  
    
    
    # calculate albumin ratio; based on the Poulin and Theil 2009, below Table 6
    PMR_PTL <- Calb_P/Calb_PTL   # plasma to proximal tubule lumen albumin ratio
    PMR_L_ec <- Calb_P/Calb_L_ec # plasma to liver albumin ratio
    PMR_PT_ec <- Calb_P/Calb_PT_ec # plasma to liver albumin ratio

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
    
    # calculate fraction unbound in different matrices; 
    # based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012; adapted to not account for fraction unionised as OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
    fu_PTL <- PMR_PTL*fup/(1 + ((PMR_PTL-1)*fup))   # proximal tubule lumen
    fu_L_ec <- PMR_L_ec*fup/(1 + ((PMR_L_ec-1)*fup)) # liver extracellular space
    fu_PT_ec <- PMR_PT_ec*fup/(1 + ((PMR_PT_ec-1)*fup)) # liver extracellular space
    
    
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
    
    # intestine
    Clear_apparent <- (Papp_SI*SA_SI*1e-3)     # L/d, passive absorption (1e-3 to convert (cm^3)ml to L)
    
    Clear_OATP2B1_I <- Clear_OATP2B1c*SF_OATP2B1_I*VI_L  # L/d, active uptake, it is assume0d that the fraction unbound is 1 
    
    CF_I <- QI/(Clear_apparent+Clear_OATP2B1_I)
    CF_I_corrected <- ifelse(CF_I < 1, CF_I, 1) # correction factor in case the active clearance is higher than the blood flow
    
    Clear_uptake_I <- (Clear_apparent+Clear_OATP2B1_I) #* CF_I_corrected
    
    # liver
    Clear_OATP1B1 <- Clear_OATP1B1c*SF_OATP1B1*VL_ic*fu_L_ec   # L/d
    Clear_OATP1B3 <- Clear_OATP1B3c*SF_OATP1B3*VL_ic*fu_L_ec   # L/d
    Clear_OATP2B1_L <- Clear_OATP2B1c*SF_OATP2B1*VL_ic*fu_L_ec  # L/d
    Clear_NTCP <- Clear_NTCPc*SF_NTCP*VL_ic*fu_L_ec            # L/d    
    Clear_BSEP <- Clear_BSEPc*SF_BSEP*VL_ic  # L/d, it is assumed that fraction unbound is 1 as binding proteins (albumin in this case) are in the interstitium so not in the cells
    
    CF_L_uptake <- (QI + QL)/(Clear_OATP1B1 + 
                             Clear_OATP1B3 + 
                             Clear_OATP2B1_L + 
                             Clear_NTCP)
    CF_L_uptake_corrected <- ifelse(CF_L_uptake < 1, CF_L_uptake, 1) # correction factor in case the active clearance is higher than the blood flow
    
    Clear_uptake_L <- (Clear_OATP1B1 +
                         Clear_OATP1B3 + 
                         Clear_OATP2B1_L +
                         Clear_NTCP) #* CF_L_uptake_corrected
    
    CF_L_excretion <- (QI + QL)/(Clear_BSEP)
    CF_L_excretion_corrected <- ifelse(CF_L_excretion < 1, CF_L_excretion, 1) # correction factor in case the active clearance is higher than the blood flow
    
    Clear_excretion_L <- Clear_BSEP #* CF_L_excretion_corrected
    
    # kidney
    Clear_OAT3 <- Clear_OAT3c*SF_OAT3*VPT*fu_PT_ec  # L/d 
    Clear_OAT1 <- Clear_OAT1c*SF_OAT1*VPT*fu_PT_ec  # L/d 
    Clear_OAT4 <- Clear_OAT4c*SF_OAT4*VPT*fu_PTL  # L/d 
    
    # CF_K_excretion <- QK/(Clear_OAT1 + Clear_OAT3)
    # CF_K_excretion_corrected <- ifelse(CF_K_excretion < 1, CF_K_excretion, 1)
    # Clear_K_excretion <- CF_K_excretion_corrected * (Clear_OAT1 + Clear_OAT3)
    
    # calculate net reabsorption (difference between active elimination and reabsorption)
    Clear_net_reabsorption_unscaled <- Clear_OAT4 - (Clear_OAT1 + Clear_OAT3)
    
    CF_K_reabsorption <- QK/Clear_net_reabsorption_unscaled
    CF_K_reabsorption_corrected <- ifelse(CF_K_reabsorption < 1, CF_K_reabsorption, 1)
    
    Clear_net_reabsorption <- Clear_net_reabsorption_unscaled #* CF_K_reabsorption_corrected # this already corrects for the fraction unbound in the proximal tubule lumen 
    
    
    
    # Clearance via menstruation
    Qmenses <- Qmenses  # L/d
    
    
    return(list(
      BW = BW, 
      VI_L = VI_L,
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
      fu_L_ec = fu_L_ec,
      fu_PT_ec = fu_PT_ec,
      
      PI = PI,
      PL = PL,
      PK = PK,
      PA = PA,
      PR = PR,
      
      SF_OATP2B1_I = SF_OATP2B1_I,
      SF_OATP2B1 = SF_OATP2B1,
      SF_OATP1B3 = SF_OATP1B3,
      SF_OATP1B1 = SF_OATP1B1,
      SF_NTCP = SF_NTCP,
      SF_BSEP = SF_BSEP,
      SF_OAT1 = SF_OAT1,
      SF_OAT3 = SF_OAT3,
      SF_OAT4 = SF_OAT4,
      
      Clear_apparent = Clear_apparent,
      Clear_OATP2B1_I = Clear_OATP2B1_I,
      CF_I = CF_I_corrected,
      Clear_uptake_I = Clear_uptake_I,
      
      Clear_OATP1B1 = Clear_OATP1B1,
      Clear_OATP1B3 = Clear_OATP1B3, 
      Clear_OATP2B1_L = Clear_OATP2B1_L,
      Clear_NTCP = Clear_NTCP,
      Clear_BSEP = Clear_BSEP,
      CF_L_uptake = CF_L_uptake_corrected,
      Clear_uptake_L = Clear_uptake_L,
      CF_L_excretion = CF_L_excretion_corrected,
      Clear_excretion_L = Clear_excretion_L,
      
      Clear_OAT3 = Clear_OAT3,
      Clear_OAT1 = Clear_OAT1,
      Clear_OAT4 = Clear_OAT4,
      Clear_net_reabsorption_unscaled = Clear_net_reabsorption_unscaled,
      CF_K_reabsorption = CF_K_reabsorption_corrected,
      Clear_net_reabsorption = Clear_net_reabsorption,
      
      Qmenses = Qmenses
    ))
  })
}

# function to run the PBK model
PBK.model <- function(t, state, pars){
  with(as.list(c(state, pars)), {
    
    # # to help with solver stability... force positive values and remove nan and inf
    state <- pmax(state, 0)
    state[is.nan(state) | is.infinite(state)] <- 0
    
    ## Concentrations -------------------------
    
    CIL <- AIL/VI_L               # ug/L, Intestinal lumen
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
    
    dAIL <- - tco*AIL +
      - Clear_uptake_I*CIL +
      + Clear_excretion_L*CL_ic          # ug/d, Intestine lumen
    
    dAI <- QI*(CP - CVI) +
      + Clear_uptake_I*CIL                      # ug/d, Intestinal
    
    dAFe <- + tco*AIL                                       # ug/d, Feces
    
    dAL_ec <- + QI*CVI + QL*CP - (QI+QL)*CVL_ec + 
      - Clear_uptake_L*CL_ec          # ug/d, Liver extracellular space (vascular + interstitial space)

    dAL_ic <- + Clear_uptake_L*CL_ec +
      - Clear_excretion_L*CL_ic             # ug/d, Liver intracellular space
    
    dAPTT <- + (QK - (fup*GFR))*CP - QK*CPTT + 
      + Clear_net_reabsorption*CPTL # ug/d, Proximal tubule tissue 
    
    dAPTL <- +(fup*GFR*CP) - QT*CPTL +
      - Clear_net_reabsorption*CPTL   # ug/d, Proximal tubule lumen    
    
    dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
    
    dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
    
    dAUr <- QUr*CRKL                                         # ug/d, Urine
    
    dAA <- QA*(CP-CVA)                                      # ug/d, Adipose
    
    dAR <- QR*(CP-CVR)                                      # ug/d, Rest
    
    dAP <- - (QI + QL + QA + QR + (QK- (fup*GFR)) + (fup*GFR))*CP +     # ug/d, Arterial Plasma
      + (QL+QI)*CVL_ec + QK*CVRKT + QA*CVA + QR*CVR + # ug/d, Venous Plasma
      - Qmenses*CP # plasma clearance during menstruation
    
    dAMp <- + Qmenses*CP
    
    # Mass Balance
    Atot <- AIL +
      AI + AFe + 
      AL_ec + AL_ic +
      APTT + APTL + ARKT + ARKL + AUr +
      AA + 
      AR +
      AP +
      AMp
    
    dInMB <- 0 # used for mass balance calculations
    
    AIn <- InMB
    MB <- AIn - Atot + 1    
    
    # End
    
    # if (AIL < 1e-08) dAIL <- max(dAIL, 0)
    # if (AI < 1e-08) dAI <- max(dAI, 0)
    # if (AFe < 1e-08) dAFe <- max(dAFe, 0)
    # if (AL_ec < 1e-08) dAL_ec <- max(dAL_ec, 0)
    # if (AL_ic < 1e-08) dAL_ic <- max(dAL_ic, 0)
    # if (APTT < 1e-08) dAPTT <- max(dAPTT, 0)
    # if (APTL < 1e-08) dAPTL <- max(dAPTL, 0)
    # if (ARKT < 1e-08) dARKT <- max(dARKT, 0)
    # if (ARKL < 1e-08) dARKL <- max(dARKL, 0)
    # if (AUr < 1e-08) dAUr <- max(dAUr, 0)
    # if (AA < 1e-08) dAA <- max(dAA, 0)
    # if (AR < 1e-08) dAR <- max(dAR, 0)
    # if (AP < 1e-08) dAP <- max(dAP, 0)
    # if (AMp < 1e-08) dAMp <- max(dAMp, 0)
    # if (InMB < 1e-08) dInMB <- max(dInMB, 0)
    
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
           dAMp,
           dInMB
    ), 
    c(AIn = AIn,
      MB = MB,
      CIL = CIL, CI = CI, CVI = CVI,
      CL_ec = CL_ec, CVL_ec = CVL_ec, CL_ic = CL_ic,
      CPTT = CPTT, CPTL = CPTL, CRKT = CRKT, CVRKT = CVRKT, CRKL = CRKL,
      CA = CA, CVA = CVA,
      CR = CR, CVR = CVR,
      CP = CP
    ))
  })
}
