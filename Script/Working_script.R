# --------------------------------------------------------------------------- #
# Script for running the Morris Test
# CP, 30-06-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Load packages
library(here)
library(ggplot2)
library(deSolve)
library(tidyverse)
library(sensitivity)
library(ggrepel)
library(scales)
library(patchwork)
library(pracma)
library(PKNCA)
library(glue)


# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()



CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'), #element_blank(), #
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), 
    plot.title = element_text(hjust = 0, size = 35, face = "bold"),
    plot.title.position = "plot", 
    
    panel.background = element_blank(), #rect(fill = 'white', color = NA)
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_line(color = "grey90",linewidth = 0.1),  
    
    axis.title = element_text(size = 35), #, face = "bold"
    axis.text = element_text(size = 35),
    # axis.text.y = element_blank(),
    # axis.text.x = element_text(size = 64, face = "bold"),
    axis.line = element_blank(), #element_line(linewidth = 0.5),
    axis.ticks = element_line(color = "grey90", linewidth = 0.1),
    axis.ticks.length = unit(0.05, "cm"),
    
    legend.title = element_text(size = 35), 
    legend.text = element_text(size = 30),
    # legend.margin = margin(0, 0, 0, 0, "cm"),
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "cm"), 
    legend.spacing.x = unit(0.2, "cm"),
    legend.position = "right", # legend.position = c(0.9, 0.9),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
    legend.margin = margin(1, 1, 1, 1),
    
    strip.text.x = element_text(size = 35)
  )


# Set output storage directory
OUTPUT <- here("Output", "Testing clearance values and model stability") 
dir.create(OUTPUT, recursive = TRUE)

# Input ----
load(here("Input", "LSA_initial_params.RData"))

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
    
    # Intestine
    # print("Intestine")
    CL_OATP2B1_I <- CL_OATP2B1c*60*24*SF_OATP2B1_I*VI #L/d
    
    CF_I <- QI/CL_OATP2B1_I
    CF_I <- ifelse(CF_I < 1, CF_I, 1)
    # print("Correction factor")
    # print(CF_I)
    CL_OATP2B1_I <- CL_OATP2B1_I * CF_I 
    # print(QI)
    # print("OATP2B1_I")
    # print(CL_OATP2B1_I)
    # 
    
    # Liver in
    # print("Liver in")
    CL_OATP1B1 <- CL_OATP1B1c*60*24*SF_OATP1B1*VL_ec #L/d
    CL_OATP1B3 <- CL_OATP1B3c*60*24*SF_OATP1B3*VL_ec #L/d
    CL_OATP2B1_L <- CL_OATP2B1c*60*24*SF_OATP2B1*VL_ec #L/d
    CL_NTCP <- CL_NTCPc*60*24*SF_NTCP*VL_ec # L/d   
    
    CF_L_in <- (QI + QL)/(CL_OATP1B1 + CL_OATP1B3 + CL_OATP2B1_L + CL_NTCP)
    CF_L_in <- ifelse(CF_L_in < 1, CF_L_in, 1)
    
    # print("Correction factor")
    # print(CF_L_in)
    CL_OATP1B1 <- CL_OATP1B1 * CF_L_in
    CL_OATP1B3 <- CL_OATP1B3 * CF_L_in
    CL_OATP2B1_L <- CL_OATP2B1_L * CF_L_in
    CL_NTCP <- CL_NTCP * CF_L_in
    
    # print(QI + QL)
    # print(CL_OATP1B1 + CL_OATP1B3 + CL_OATP2B1_L + CL_NTCP)
    # print("1B1")
    # print(CL_OATP1B1)
    # print("1B3")
    # print(CL_OATP1B3)
    # print("2B1")
    # print(CL_OATP2B1_L)
    # print("NTCP")
    # print(CL_NTCP)
    
    # print("Liver out")
    CL_BSEP <- CL_BSEPc*60*24*SF_BSEP*VL_ic         # L/d
    CF_L_out <- (QI + QL)/(CL_BSEP)
    CF_L_out <- ifelse(CF_L_out < 1, CF_L_out, 1)
    
    # print("Correction factor")
    # print(CF_L_out)
    
    CL_BSEP <- CL_BSEP * CF_L_out
    # print(QI + QL)
    # print("BSEP")
    # print(CL_BSEP)
    
    # print("Kidney in")
    CL_OAT3 <- CL_OAT3c*60*24*SF_OAT3*VPT           # L/d (min -> d)
    CL_OAT1 <- CL_OAT1c*60*24*SF_OAT1*VPT           # L/d (min -> d)
    
    CF_K_in <- QK/(CL_OAT1 + CL_OAT3)
    CF_K_in <- ifelse(CF_K_in < 1, CF_K_in, 1)
    
    # print("Correction factor")
    # print(CF_K_in)
    
    CL_OAT1 <- CL_OAT1 * CF_K_in
    CL_OAT3 <- CL_OAT3 * CF_K_in
    
    # print(QK)
    # print("CL_OAT1")
    # print(CL_OAT1)
    # print("CL_OAT3")
    # print(CL_OAT3)
    
    # print("Kidney out")
    CL_OAT4 <- CL_OAT4c*60*24*SF_OAT4*VPT           # L/d (min -> d)
    
    CF_K_out <- GFR/CL_OAT4
    CF_K_out <- ifelse(CF_K_out < 1, CF_K_out, 1)
    # print("Correction factor")
    # print(CF_K_out)
    CL_OAT4 <- CL_OAT4 * CF_K_out
    # print(GFR)
    # print("CL_OAT4")
    # print(CL_OAT4)
    
    
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
    ),
    c(QK = QK,
      CL_OAT1 = CL_OAT1, 
      CL_OAT3 = CL_OAT3,
      CF_K_in = CF_K_in,
      GFR = GFR,
      CL_OAT4 = CL_OAT4,
      CF_K_out = CF_K_out,
      QH = QL + QI,
      CL_OATP1B1 = CL_OATP1B1,
      CL_OATP1B3 = CL_OATP1B3,
      CL_OATP2B1_L = CL_OATP2B1_L, 
      CL_NTCP = CL_NTCP,
      CF_L_in = CF_L_in,
      CL_BSEP = CL_BSEP,
      CF_L_out = CF_L_out,
      QI = QI,
      CL_OATP2B1_I = CL_OATP2B1_I,
      CF_I = CF_I,
      CL_IL = CL_IL)
    )
  })
}

PFAS_list <- c("ADONA", 
               "HFPO-DA",
               "PFBA",
               "PFBS",
               "PFDA",
               "PFDoDA",
               "PFHpA",
               "PFHpS",
               "PFHxA",
               "PFHxS",
               "PFNA",
               "PFOA",
               "PFOS",
               "PFPeA",
               "PFPeS",
               "PFTrDA",
               "PFUnDA")

nPFAS <- length(PFAS_list)

# Run ----


Predictions_collect <- list()
HL_collect <- list()
Cmax_collect <- list()

for (i in 1:nPFAS) {
  
  # i = 1
  PFAS <- PFAS_list[i]
  message(glue("{PFAS}----------------------------"))
  
  pars <- ModelParams[[PFAS]]
  
  A_init <- c(OD = 0,
              AIL = 0.1,
              AI = 0, 
              AFe = 0,
              AL_ec = 0,
              AL_ic = 0,
              APTT = 0,
              APTL = 0,
              ARKT = 0,
              ARKL = 0,
              AUr = 0,
              AA = 0,
              AR = 0,
              AP = 0,
              AMp = 0,
              Ain = 0)
  
  times = seq(0, 5*356, by = 1)
  
  PBK.out <- lsoda(y = A_init,
                   times = times,
                   func = PBK.model, 
                   parms = pars)
  
  PBK_OUTPUT <- as.data.frame(PBK.out) %>%
    mutate(CP = ifelse(CP < 0, 0, CP)) %>% 
    mutate(PFAS = PFAS) %>% 
    na.omit()

  Cmax <- round(max(PBK_OUTPUT$CP), 4) # ng/ml
  message(glue("-Cmax {Cmax} ng/ml--"))

  Tmax <- filter(PBK_OUTPUT, CP == max(PBK_OUTPUT$CP)) %>% pull(time)
  Tmax <- round(Tmax / 365, 4)
  message(glue("-Tmax {Tmax} years--"))

  half_life <- pk.calc.half.life(
    PBK_OUTPUT$CP,
    PBK_OUTPUT$time / 365
  )
  
  half_life <- round(half_life$half.life, 4)
  message(glue("-HL {half_life} years--"))
  
  print(ggplot(PBK_OUTPUT) +
    geom_path(aes(time/365, CP), size = 0.5) +
      labs(title = PFAS, 
           x = "Time (years)", 
           y = "Plasma concentration (ng/ml)") +
      theme_minimal() +
      theme(panel.border = element_rect(color = "black", 
                                        fill = NA, 
                                        linewidth = 0.5))
    )
  
  Predictions_collect[[i]] <- PBK_OUTPUT
  HL_collect[i] <- half_life
  Cmax_collect[i] <- Cmax
  
}

Predictions <- do.call(rbind, Predictions_collect)
TK_parameters <- Predictions %>% select(QK:PFAS) %>% distinct()
# HL_collect <- unlist(HL_collect)
# Cmax_collect <- unlist(Cmax_collect)


names(Predictions_collect) = names(HL_collect) = names(Cmax_collect) = PFAS_list

Results <- NULL
Results$Predictions = Predictions_collect
Results$HL = HL_collect
Results$Cmax = Cmax_collect

save(Results, 
     file = here(OUTPUT, "Results_mean_values.RData"))
