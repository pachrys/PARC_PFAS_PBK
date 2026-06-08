# --------------------------------------------------------------------------- #
# SCRIPT FOR PBK MODEL FOR ORAL EXPOSURE
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

# ORAL ####
ORAL_PBK_RUN <- function(y, parms, times){ # Input for ode
  
  PBK.model <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      
      ### Physiological ----
      
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
      
      VP <- VPc * (1-Hct) * BW                  # L, Volume of plasma
      
      VTotc <- 0.96
      VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
      
      VR <- VTot - (VI + VL + VK + VA + VP)             # L, Volume of the lumped rest compartment
      
      QTotc <- 0.988
      QI <- QIc/QTotc * QC                  # L/d, Intestinal
      QL <- QLc/QTotc * QC                  # L/d, Liver
      QK <- QKc/QTotc * QC                  # L/d, Kidney
      QA <- QAc/QTotc * QC                  # L/d, Adipose 
      QR <- QC - (QA + QI + QK + QL)  # L/d, Rest
      
      QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
      GFR <- GFR
      QT <- QTc*GFR                      # L/d, Proximal tubule fluid flow
      
      tco <- tco                    # /d, Bowel residence times in the colon
      
      # Fraction unbound in plasma
      fup <- fup #0.061/100 #
      
      PI <- PIc * fup     # Intestinal
      PL <- PLc * fup     # Liver
      PK <- PKc * fup     # Kidney
      PA <- PAc * fup     # Adipose
      PR <- PRc * fup     # Rest
      
      # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
      # Equation was adapted to not account for fraction unionised
      # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
      fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
      fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
      fu_In  <- fu_Lic # Intestine intracellular space
      
      
      ### Kinetic ----
      
      # Clearances
      CL_IL <- (Papp_SI*SA_SI*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24)
      message(glue("Abs {CL_IL}"))
      
      CL_OATP2B1_I <- CL_OATP2B1c*60*24*SF_OATP2B1_I*VI #L/d
      message(glue("OATP2B1_I {CL_OATP2B1_I}"))
      
      CL_OATP1B1 <- CL_OATP1B1c*60*24*SF_OATP1B1*VL_ec #L/d
      message(glue("OATP1B1 {CL_OATP1B1}"))
      
      CL_OATP1B3 <- CL_OATP1B3c*60*24*SF_OATP1B3*VL_ec #L/d
      message(glue("OATP1B3 {CL_OATP1B3}"))
      
      CL_OATP2B1_L <- CL_OATP2B1c*60*24*SF_OATP2B1_L*VL_ec #L/d
      message(glue("OATP2B1_L {CL_OATP2B1_L}"))
      
      CL_NTCP <- CL_NTCPc*60*24*SF_NTCP*VL_ec # L/d    
      message(glue("NTCP {CL_NTCP}"))
      
      CL_BSEP <- CL_BSEPc*60*24*SF_BSEP*VL_ic         # L/d
      message(glue("BSEP {CL_BSEP}"))
      
      CL_OAT4 <- CL_OAT4c*60*24*SF_OAT4*VPT           # L/d (min -> d; /mg protein 
      message(glue("OAT4 {CL_OAT4}"))
      
      CL_OAT3 <- CL_OAT3c*60*24*SF_OAT3*VPT           # L/d (min -> d)
      message(glue("OAT3 {CL_OAT3}"))
      
      CL_OAT1 <- CL_OAT1c*60*24*SF_OAT1*VPT           # L/d (min -> d)
      message(glue("OAT1 {CL_OAT1}"))
      
      # Clearance via menstruation
      CL_menses <- CL_menses     # L/d, plasma menstrual loss
      
      ## Dose -------------------------
      
      if(t<expSTOP){DoseOn=1} else{DoseOn=0}
      
      ## Oral exposure ##
      DOral <- expOral*BW*DoseOn         # ug, PFOA oral dose
      OralD <- DOral
      
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
      
      dOD <-  OralD - OD             # ug/d, Oral dose input OralD
      
      dAIL <- + OD - tco*AIL - CL_IL*CIL +
        - CL_OATP2B1_I*CIL + #*fu_In +
        + CL_BSEP*CL_ic*fu_Lic          # ug/d, Intestine lumen
      
      dAI <- QI*(CP - CVI) + CL_IL*CIL + 
        + CL_OATP2B1_I*CIL #*fu_In                     # ug/d, Intestinal
      
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
      Atot <- OD +
        AIL + AI + AFe + 
        AL_ec + AL_ic +
        APTT + APTL + ARKT + ARKL + AUr +
        AA + 
        AR +
        AP +
        AMp
      
      dAin <- OralD #1.702140 # to be used if repeated exposure
      MB <- Ain - Atot + 1    # to be used if repeated exposure
      # MB <- DOral - Atot + 1
      
      # End
      
      list(c(dOD,
             dAIL,
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
             dAin
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
        CP = CP,
        Atot = Atot, 
        MB = MB
      )
      )
    })
  }
  
  output_PBK <-lsoda(y = y, 
                      times = times, 
                      func = PBK.model,
                      parms = parms)
  
  return(as.data.frame(output_PBK))
}
