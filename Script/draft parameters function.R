# Calc parameters function

Calc_final_params <- function(parm.MC, parm.phys){
  
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
  
  # neutral phospholipids:water partition coefficient
  k_ML <- PhysChemTK_params %>% filter(Parameter == "logD_MembraneLipidWater") %>% pull(Value) %>% as.numeric() |> exp()
  
  # structural protein:water partition coefficient
  k_SP <- PhysChemTK_params %>% filter(Parameter == "logD_StructuralProteinWater") %>% pull(Value) %>% as.numeric() |> exp()
  
  # albumin:water partition coefficient
  k_ALB <- PhysChemTK_params %>% filter(Parameter == "logD_AlbuminWater") %>% pull(Value) %>% as.numeric() |> exp()
  
  # structural lipids:water partition coefficient
  k_SL <- PhysChemTK_params %>% filter(Parameter == "logD_StorageLipidWater") %>% pull(Value) %>% as.numeric() |> exp()
  
  # globulin:water partition coefficient
  k_GLOB <- PhysChemTK_params %>% filter(Parameter == "logD_GlobulinWater") %>% pull(Value) %>% as.numeric() |> exp()
  
  
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
  
  Physio.data <- Tissue.c %>% 
    filter(Species == "Human") %>%
    filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
    select(- c(Comment, Species)) %>% 
    mutate(across(-Tissue, as.numeric))
  
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
  PIc <- Kp.df$KpIn   # Intestine
  PLc <- Kp.df$KpLi   # Liver
  PKc <- Kp.df$KpKi   # Kidney
  PAc <- Kp.df$KpAd   # Adipose
  PRc <- Kp.df$KpRe   # Rest
  
  ## clearances ----
  
  # un scaled clearances
  Papp_SI <- PhysChemTK_params %>% filter(Parameter == "Papp_SI") %>% pull(Value) %>% as.numeric() |> exp() # cm/s
  CL_OATP2B1c <- PhysChemTK_params %>% filter(Parameter == "CL_OATP2B1c") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  CL_OATP1B1c <- PhysChemTK_params %>% filter(Parameter == "CL_OATP1B1c") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  CL_OATP1B3c <- PhysChemTK_params %>% filter(Parameter == "CL_OATP1B3c") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  CL_OATP1B3c <- PhysChemTK_params %>% filter(Parameter == "CL_OATP1B3c") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  CL_NTCPc <- PhysChemTK_params %>% filter(Parameter == "CL_NTCPc") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  CL_BSEPc <- PhysChemTK_params %>% filter(Parameter == "CL_BSEPc") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  CL_OAT1c <- PhysChemTK_params %>% filter(Parameter == "CL_OAT1c") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  CL_OAT3c <- PhysChemTK_params %>% filter(Parameter == "CL_OAT3c") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  CL_OAT4c <- PhysChemTK_params %>% filter(Parameter == "CL_OAT4c") %>% pull(Value) %>% as.numeric() # L/d/mg transporter
  
  # scaling factors
  SF_OATP2B1_I <- PhysChemTK_params %>% filter(Parameter == "SF_OATP2B1_I") %>% pull(Value) %>% as.numeric() # mg transporter/kg tissue
  SF_OATP1B1 <- PhysChemTK_params %>% filter(Parameter == "SF_OATP1B1") %>% pull(Value) %>% as.numeric()
  SF_OATP1B3 <- PhysChemTK_params %>% filter(Parameter == "SF_OATP1B3") %>% pull(Value) %>% as.numeric()
  SF_OATP2B1 <- PhysChemTK_params %>% filter(Parameter == "SF_OATP2B1") %>% pull(Value) %>% as.numeric() # mg transporter/kg tissue
  SF_NTCP <- PhysChemTK_params %>% filter(Parameter == "SF_NTCP") %>% pull(Value) %>% as.numeric() # mg transporter/kg tissue
  SF_BSEP <- PhysChemTK_params %>% filter(Parameter == "SF_BSEP") %>% pull(Value) %>% as.numeric() # mg transporter/kg tissue
  SF_OAT1 <- PhysChemTK_params %>% filter(Parameter == "SF_OAT1") %>% pull(Value) %>% as.numeric() # mg transporter/kg tissue
  SF_OAT3 <- PhysChemTK_params %>% filter(Parameter == "SF_OAT1") %>% pull(Value) %>% as.numeric() # mg transporter/kg tissue
  
  
  print("Intestine - passive absorption")
  CL_IL <- (Papp_SI*SA_SI*1e-3)    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/d -> L/d /1000)
  
  CL_OATP2B1_I <- CL_OATP2B1c*SF_OATP2B1_I*VI # L/d
  print("Intestine - active absorption")
  print(QI)
  print(CL_OATP2B1_I)
  
  print("Liver uptake")
  CL_OATP1B1 <- CL_OATP1B1c*SF_OATP1B1*VL_ec   # L/d
  CL_OATP1B3 <- CL_OATP1B3c*SF_OATP1B3*VL_ec   # L/d
  CL_OATP2B1_L <- CL_OATP2B1c*SF_OATP2B1*VL_ec # L/d
  CL_NTCP <- CL_NTCPc*SF_NTCP*VL_ec            # L/d    
  print(QI+QL)
  print(CL_OATP1B1 + CL_OATP1B3 + CL_OATP2B1_L)
  
  print("Liver excretion")
  CL_BSEP <- CL_BSEPc*SF_BSEP*VL_ic  # L/d
  print(QI + QL)
  print(CL_BSEP)
  
  print("Kidney uptake")
  CL_OAT3 <- CL_OAT3c*SF_OAT3*VPT  # L/d 
  CL_OAT1 <- CL_OAT1c*SF_OAT1*VPT  # L/d 
  print(QK)
  print(CL_OAT1 + CL_OAT3)
  
  print("Kidney out")
  CL_OAT4 <- CL_OAT4c*SF_OAT4*VPT  # L/d 
  print(GFR)
  print(CL_OAT4)
  
  # Clearance via menstruation
  CL_menses <- CL_menses  # L/d
  
  
  return(list(
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
    QI = QI,
    QL = QL,
    QK = QK,
    QA = QA,
    QR = QR,
    QUr = QUr,
    GFR = GFR,
    QT = QT,
    fup = fup,
    fu_PTL = fu_PTL,
    fu_Lic = fu_Lic,
    fu_In = fu_In,
    PI = PI,
    PL = PL,
    PK = PK,
    PA = PA,
    PR = PR,
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
  
}
