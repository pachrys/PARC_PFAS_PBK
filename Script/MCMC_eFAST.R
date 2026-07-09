# --------------------------------------------------------------------------- #
# Script for running the MCMC using eFAST
# Chrysanthi Pachoulide, June 2026
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment


# Load packages
library(here)
library(ggplot2)
library(deSolve)
library(tidyverse)
library(sensitivity)
library(glue)
library(readxl)
library(PKNCA)
library(parallel)
library(pracma)
library(viridis)

# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'), #element_blank(), #
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), 
    plot.title = element_text(hjust = 0, size = 35, face = "bold"),
    plot.title.position = "plot", 
    
    panel.background = element_blank(), #rect(fill = 'white', color = NA)
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_line(color = "grey90",linewidth = 0.1),  
    
    axis.title = element_text(size = 35), #, face = "bold"
    axis.text = element_text(size = 35),
    # axis.text.y = element_blank(),
    axis.text.x = element_text(size = 35, face = "bold", angle = 30),
    axis.line = element_blank(), #element_line(linewidth = 0.5),
    axis.ticks = element_line(color = "grey90", linewidth = 0.1),
    axis.ticks.length = unit(0.05, "cm"),
    
    legend.title = element_blank(), #element_text(size = 35), 
    legend.text = element_text(size = 30),
    # legend.margin = margin(0, 0, 0, 0, "cm"),
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "cm"), 
    legend.spacing.x = unit(0.2, "cm"),
    legend.position = "right", # legend.position = c(0.9, 0.9),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
    legend.margin = margin(1, 1, 1, 1),
    
    strip.text.x = element_blank()
    # strip.text.x = element_text(size = 35)
  )

theme_quick <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'), 
    plot.title = element_text(hjust = 0, face = "bold"),
    plot.title.position = "plot", 
    
    panel.background = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_line(color = "grey90"),  
    
    axis.line = element_blank(), 
    
    legend.position = "right"
  )

# Set output storage directory
OUTPUT <- here("Output", "MonteCarlo_1000") #format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
# dir.create(OUTPUT, recursive = TRUE)

# Load data

## PFAS

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

## Tissue constituents for partition coefficient calculation
Physio.data <- read_csv(here("Input", "TissueComposition.csv")) %>% 
  filter(Species == "Human") %>%
  filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
  select(- c(Comment, Species)) %>% 
  mutate(across(-Tissue, as.numeric))

## All PBK parameters
load(here("Input", "LSA_initial_params.RData"))

## MC parameters
MC_TK_params <- read_xlsx(here("Input", "MC_TK_Parameters.xlsx")) |>
  mutate(Predicted = ifelse(Specific == "calculated", "yes", Predicted))|>
  filter(Predicted == "yes") |>
  select(PFAS, Parameter, Value, SD, Distribution)
  # mutate(SD = ifelse(Distribution == "lognormal", SD-0.2, SD)) #SD = ifelse((Distribution == "uniform") & (Value - SD <=0), Value - (Value/2), SD),


# Create function to calculate the final parameters
Calc_final_params <- function(parm.MC, parm.phys){
  with(as.list(c(parm.MC, parm.phys)), {
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

# Create Monte Carlo function
run_Monte_Carlo <- function(parm.MC, parm.phys){
  
  A_init <- c(AIL = 0.1,
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
              AMp = 0)
  
  pars <- Calc_final_params(parm.MC = parm.MC, 
                            parm.phys = parm.phys)
  # save parameters to check when analysing results
  run_pars <- as.data.frame(pars)
  PBK.out <- lsoda(y = A_init,
                   times = seq(0, 1825, by = 10),
                   func = PBK.model, 
                   parms = pars)
  
  PBK_OUTPUT <- as.data.frame(PBK.out) %>% 
    mutate(CP = ifelse(CP < 0-1000, 0, CP)) 
  
  PBK_OUTPUT <- cbind(PBK_OUTPUT, run_pars)
  
 return(PBK_OUTPUT)
  
}

Monte_Carlo_PFAS_design <- list()
Monte_Carlo_PFAS_RESULTS <- list()
Monte_Carlo_PFAS_y <- list()
Monte_Carlo_PFAS_savey <- list()

for (i in 1:nPFAS) {


  
  PFAS_i <- PFAS_list[i]#"PFOA" 
  message(glue("{PFAS_i}------------------------------------"))
  
  # physiological parameters; they are kept fix
  parm.phys <- ModelParams[[PFAS_i]]
  
  
  PhysChemTK_params <- MC_TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(!PFAS) 
  
  parm.MC <- PhysChemTK_params |> 
    mutate(meanlog = ifelse(Distribution == "lognormal", Value, NA),
           sdlog = ifelse(Distribution == "lognormal", SD, NA),
           min = ifelse(Distribution == "uniform", Value - SD, NA),
           max = ifelse(Distribution == "uniform", Value + SD, NA),
           Distribution = ifelse(Distribution == "lognormal", "qlnorm", "qunif")) |> 
    mutate(min = ifelse(min <= 0, 1e-10, min))

  MC.parms <- parm.MC$Parameter
  MC.factors <- length(MC.parms)

  q <- parm.MC$Distribution

  # Define q.arg: q.arg needs to be a list of lists
  q.arg <- lapply(1:nrow(parm.MC), function(i) { # apply the function row-wise
    if (parm.MC$Distribution[i] == "qunif") { #if the distribution in uniform
      list(min = parm.MC$min[i], max = parm.MC$max[i]) #create a list containing min and max per row
    } else {
      list(meanlog = parm.MC$meanlog[i], sdlog = parm.MC$sdlog[i]) #if not then create a list containing mean and sdlog
    }
  })
  names(q.arg) <- MC.parms

  set.seed(10010)
  
  length(q) == length(q.arg)

  Experience <- fast99(
    model = NULL, 
    factors = MC.parms, # These are the parameters that will be varying in the model 
    n = 1000,# Integer giving the sample size, i.e. the length of the discretization of the s-space 
    q = q,
    q.arg = q.arg)

  
  dim(Experience$X) # dataframe of: n factors(number of parameters) * n n(1000) number of observations, of n factors(number of parameters) number of variables
  design <- Experience$X
  write.csv(design, file = here(OUTPUT, paste0(PFAS_i,"_Initial_MC_ExpDesign.csv")), row.names = FALSE)
  
  # check distributions
  design |> pivot_longer(everything(), 
                         names_to = "Parameter", 
                         values_to = "Value") |> 
    ggplot(aes(x=Value))+
    geom_histogram()+
    facet_wrap(~Parameter, scales = "free")+
    labs(title = "Monte Carlo parameter distributions") +
    theme_quick 

  # clean by removing infinite values
  Experience$X <- Experience$X |> mutate(across(everything(), ~na_if(., Inf))) |> drop_na()
  clean_design <- Experience$X
  write.csv(clean_design, file = here(OUTPUT, paste0(PFAS_i,"MC_ExpDesign.csv")), row.names = FALSE)
  Monte_Carlo_PFAS_design[[i]] <- clean_design
  
  # check distributions
  
  Design_plot <- clean_design |> 
    pivot_longer(everything(), 
                 names_to = "Parameter", 
                 values_to = "Value") |> 
    ggplot(aes(x=Value))+
    geom_histogram()+
    facet_wrap(~Parameter, scales = "free")+
    labs(title = "Monte Carlo parameter distributions",
         x = "",
         y = "") +
    theme_quick  
  Design_plot
  
  Design_plot <- Design_plot + CP_theme
  ggsave(plot = Design_plot,
         file = here(OUTPUT, paste0(PFAS_i,"_MC_ExpDesign.png")),
         dpi = 600,
         width = 18,
         height = 24,
         units = "cm")
  
  results <- apply(clean_design, 1, function(row) {
    run_Monte_Carlo(row, parm.phys = parm.phys)
  })
  
  y <- results
  savey <- mapply(cbind, y,
                  "Iteration" = seq(1,length(y), by = 1),
                  SIMPLIFY = FALSE)
  save(savey, file = here(OUTPUT, paste0(PFAS_i,"_Monte_Carlo_y.RData")))
  Monte_Carlo_PFAS_savey[[i]] <- savey
 }
names(Monte_Carlo_PFAS_savey) <- PFAS_list

library(pracma)


# Analyse results ----

calc_TK_results <- function(results) {
  
  data <-  results
  plasma <- data |> 
    select(time, CP)
  
  # AUC ----
  AUC <- trapz(data$time, data$CP) # ng * day/ml
  
  # Cmax and Tmax ----
  Cmax <- max(data$CP) # ng/ml
  Tmax <- filter(data, CP == max(data$CP)) %>% pull(time) # days
  
  # plasma half-life ----
  # elimination phase, i.e. after the tmax
  elimination_phase <- data |> 
    filter(time >= Tmax,
           CP >= 1e-9) |> 
    mutate(logCP = log(CP))
  
  # check that there is enough data points to run linear regression
  if (nrow(elimination_phase) < 3) {
    warning("not enough data-points to calculate half life")
    return(NA)
  }
  
  # linear regression 
  elimination_model <- lm(logCP ~ time, data = elimination_phase)
  
  # calculate elimination constant (i.e slope of the line)
  kel <- -elimination_model$coefficients[[2]]  
  
  # calculate the half-life
  half_life <- log(2) / kel
  
  # fraction absorbed ----
  
  # all amounts in tissues at Tmax
  A_Tmax <- data %>% filter(time == max(CP)) |> 
    select(starts_with("A")) |> 
    select(!c(AIL, AFe))
  
  # calculate the total amount absorbed at Tmax
  ATotal_Tmax <- rowSums(A_Tmax)
  
  # Absorbed fraction (as a percentage of total oral dose)
  OD <- 0.1
  Fa <- ATotal_Tmax/OD *100
  
  # fecal and urinary excretion ----
  # Calculate excretion (feces and urine)
  
  A_feces <- data |> select(AFe) |> 
    mutate(dailyFE = AFe - lag(AFe, n = 1, default = first(AFe))) |> 
    filter(dailyFE > 0)
  FfeEl <- mean(A_feces$AFe)/OD * 100
  
  A_urine <- data |> select(AUr) |> 
    mutate(dailyUr = AUr - lag(AUr, n = 1, default = first(AUr))) |> 
    filter(dailyUr > 0)
  
  FurEl <- last(A_urine$AUr)/OD * 100
  
  
  return(c(AUC = AUC,
           Cmax = Cmax,
           Tmax = Tmax,
           kel = kel,
           HL = half_life,
           Fa = Fa,
           FfeEl = FfeEl,
           FurEl = FurEl
  ))
}

## ADONA ----

load(here(OUTPUT, "ADONA_Monte_Carlo_y.RData"))
ADONA <- savey
ADONA_df <- do.call(rbind, ADONA) |> mutate(PFAS = "ADONA")
ADONA_CP <- ADONA_df |> mutate(PFAS = "ADONA") |> select(time, CP, PFAS, Iteration)
ADONA_CP_wide <- ADONA_CP |> select(!PFAS) |> 
  mutate(Iteration = paste0("n_", Iteration)) |> 
  pivot_wider(names_from = Iteration, values_from = CP) 
ADONA_CP_summary <- ADONA_CP_wide |> 
  select(!time) |> t() |> 
  as.data.frame() |> 
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75)))) |>
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value") |> 
  separate(Variable, 
           into = c("Time", "Variable")) |> 
  mutate(Time = as.numeric(str_replace(Time, "V", ""))) |> 
  mutate(Time = (Time*10)/365) |>  #unit correction, seq was per 10 days
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(ADONA_CP_summary_long, file = here(OUTPUT, "ADONA_CP_sumstat.csv"))



ADONA_TK_results <- sapply(ADONA, calc_TK_results)
ADONA_TK_results_all <- do.call(rbind, ADONA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", (Value*10)/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(ADONA_TK_results_all, file = here(OUTPUT, "ADONA_TK_sumstat.csv"))


ADONA_CP_plot <- ggplot(ADONA_CP_summary_long, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75), 
              fill = "grey70") +
  geom_path()+
  labs(title = "ADONA", 
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
ADONA_CP_plot <- ADONA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = ADONA_CP_plot,
       file = here(OUTPUT, "ADONA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")
## others ----


load(here(OUTPUT, "HFPO-DA_Monte_Carlo_y.RData"))
HFPO_DA <- savey
HFPO_DA_df <- do.call(rbind, HFPO_DA)|> mutate(PFAS = "HFPO_DA")

load(here(OUTPUT, "PFBA_Monte_Carlo_y.RData"))
PFBA <- savey
PFBA_df <- do.call(rbind, PFBA) |> mutate(PFAS = "PFBA")

load(here(OUTPUT, "PFBS_Monte_Carlo_y.RData"))
PFBS <- savey
PFBS_df <- do.call(rbind, PFBS) |> mutate(PFAS = "PFBS")

load(here(OUTPUT, "PFDA_Monte_Carlo_y.RData"))
PFDA <- savey
PFDA_df <- do.call(rbind, PFDA)|> mutate(PFAS = "PFDA")

load(here(OUTPUT, "PFDoDA_Monte_Carlo_y.RData"))
PFDoDA <- savey
PFDoDA_df <- do.call(rbind, PFDoDA) |> mutate(PFAS = "PFDoDA")

load(here(OUTPUT, "PFHpA_Monte_Carlo_y.RData"))
PFHpA <- savey
PFHpA_df <- do.call(rbind, PFHpA)|> mutate(PFAS = "PFHpA")

load(here(OUTPUT, "PFOA_Monte_Carlo_y.RData"))
PFOA <- savey
PFOA_df <- do.call(rbind, PFOA) |> mutate(PFAS = "PFOA")

load(here(OUTPUT, "PFOA_Monte_Carlo_y.RData"))
PFOA <- savey
PFOA_df <- do.call(rbind, PFOA) |> mutate(PFAS = "PFOA")


load(here(OUTPUT, "PFHxS_Monte_Carlo_y.RData"))
PFHxS <- savey
PFHxS_df <- do.call(rbind, PFHxS) |> mutate(PFAS = "PFHxS")

load(here(OUTPUT, "PFNA_Monte_Carlo_y.RData"))
PFNA <- savey
PFNA_df <- do.call(rbind, PFNA) |> mutate(PFAS = "PFNA")

load(here(OUTPUT, "PFOA_Monte_Carlo_y.RData"))
PFOA <- savey
PFOA_df <- do.call(rbind, PFOA) |> mutate(PFAS = "PFOA")

load(here(OUTPUT, "PFOS_Monte_Carlo_y.RData"))
PFOS <- savey
PFOS_df <- do.call(rbind, PFOS) |> mutate(PFAS = "PFOS")

load(here(OUTPUT, "PFPeA_Monte_Carlo_y.RData"))
PFPeA <- savey
PFPeA_df <- do.call(rbind, PFPeA) |> mutate(PFAS = "PFPeA")

load(here(OUTPUT, "PFPeS_Monte_Carlo_y.RData"))
PFPeS <- savey
PFPeS_df <- do.call(rbind, PFPeS) |> mutate(PFAS = "PFPeS")

load(here(OUTPUT, "PFTrDA_Monte_Carlo_y.RData"))
PFTrDA <- savey
PFTrDA_df <- do.call(rbind, PFTrDA) |> mutate(PFAS = "PFTrDA")

load(here(OUTPUT, "PFUnDA_Monte_Carlo_y.RData"))
PFUnDA <- savey
PFUnDA_df <- do.call(rbind, PFUnDA) |> mutate(PFAS = "PFUnDA")



HFPO_DA_TK_results <- sapply(HFPO_DA, calc_TK_results)
HFPO_DA_TK_results_all <- do.call(rbind, HFPO_DA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                         
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

HFPO_TK_plot <- ggplot(ADONA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = HFPO_TK_plot,
       file = here(OUTPUT, "HFPO_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")

PFBA_TK_results <- sapply(PFBA, calc_TK_results)
PFBA_TK_results_all <- do.call(rbind, PFBA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                         
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFBA_TK_plot <- ggplot(PFBA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFBA_TK_plot,
       file = here(OUTPUT, "PFBA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")


PFBS_TK_results <- sapply(PFBS, calc_TK_results)
PFBS_TK_results_all <- do.call(rbind, PFBS_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                          
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFBS_TK_plot <- ggplot(PFBS_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFBS_TK_plot,
       file = here(OUTPUT, "PFBS_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")


PFDA_TK_results <- sapply(PFDA, calc_TK_results)
PFDA_TK_results_all <- do.call(rbind, PFDA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                           
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFDA_TK_plot <- ggplot(PFDA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFDA_TK_plot,
       file = here(OUTPUT, "PFDA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")

PFDoDA_TK_results <- sapply(PFDoDA, calc_TK_results)
PFDoDA_TK_results_all <- do.call(rbind, PFDoDA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                          
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFDoDA_TK_plot <- ggplot(PFDoDA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFDoDA_TK_plot,
       file = here(OUTPUT, "PFDoDA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")

PFHpA_TK_results <- sapply(PFHpA, calc_TK_results)
PFHpA_TK_results_all <- do.call(rbind, PFHpA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                         
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFHpA_TK_plot <- ggplot(PFHpA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFHpA_TK_plot,
       file = here(OUTPUT, "PFHpA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")


PFHpS_TK_results <- sapply(PFHpS, calc_TK_results)
PFHpS_TK_results_all <- do.call(rbind, PFHpS_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                          
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFHpS_TK_plot <- ggplot(PFHpS_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFHpS_TK_plot,
       file = here(OUTPUT, "PFHpS_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")

PFHxA_TK_results <- sapply(PFHxA, calc_TK_results)
PFHxA_TK_results_all <- do.call(rbind, PFHxA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                          
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFHxA_TK_plot <- ggplot(PFHxA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFHxS_TK_plot,
       file = here(OUTPUT, "PFHxA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")

PFHxS_TK_results <- sapply(PFHxS, calc_TK_results)
PFHxS_TK_results_all <- do.call(rbind, PFHxS_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                          
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFHxS_TK_plot <- ggplot(PFHxS_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFHxS_TK_plot,
       file = here(OUTPUT, "PFHxS_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")

PFNA_TK_results <- sapply(PFNA, calc_TK_results)
PFNA_TK_results_all <- do.call(rbind, PFNA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                           
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFNA_TK_plot <- ggplot(PFNA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFNA_TK_plot,
       file = here(OUTPUT, "PFNA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")

PFOA_TK_results <- sapply(PFOA, calc_TK_results)
PFOA_TK_results_all <- do.call(rbind, PFOA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                          
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFOA_TK_plot <- ggplot(PFOA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFOA_TK_plot,
       file = here(OUTPUT, "PFOA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")



PFOS_TK_results <- sapply(PFOS, calc_TK_results)
PFOS_TK_results_all <- do.call(rbind, PFOS_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                           
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFOS_TK_plot <- ggplot(PFOS_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFOS_TK_plot,
       file = here(OUTPUT, "PFOS_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")




PFPeA_TK_results <- sapply(PFPeA, calc_TK_results)
PFPeA_TK_results_all <- do.call(rbind, PFPeA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                         
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFPeA_TK_plot <- ggplot(PFPeA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFPeA_TK_plot,
       file = here(OUTPUT, "PFPeA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")




PFPeS_TK_results <- sapply(PFPeS, calc_TK_results)
PFPeS_TK_results_all <- do.call(rbind, PFPeS_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                         
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFPeS_TK_plot <- ggplot(PFPeS_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFPeS_TK_plot,
       file = here(OUTPUT, "PFPeS_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")


PFTrDA_TK_results <- sapply(PFTrDA, calc_TK_results)
PFTrDA_TK_results_all <- do.call(rbind, PFTrDA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                         
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFTrDA_TK_plot <- ggplot(PFTrDA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFTrDA_TK_plot,
       file = here(OUTPUT, "PFTrDA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")

PFUnDA_TK_results <- sapply(PFUnDA, calc_TK_results)
PFUnDA_TK_results_all <- do.call(rbind, PFUnDA_TK_results) |> 
  na.omit() |> as.data.frame() |>  
  summarise(across(where(is.numeric), .fns = 
                     list(median = median,                       
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |> 
  pivot_longer(everything(), names_to = "Variable", 
               values_to = "Value")  |> 
  separate(Variable, 
           into = c("Parameter", "Variable"))|> 
  mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |> 
  pivot_wider(names_from = Variable, values_from = Value)

PFUnDA_TK_plot <- ggplot(PFUnDA_TK_results_all, aes(Parameter, mean)) +
  geom_pointrange(aes(ymin = q25,
                      ymax = q75)) +
  facet_wrap(~Parameter, scales = "free", nrow = 1) +
  labs(x = "", 
       y = "") +
  CP_theme +
  theme(axis.text.x = element_blank())
ggsave(plot = PFUnDA_TK_plot,
       file = here(OUTPUT, "PFUnDA_TK_plot.png"),
       dpi = 600,
       width = 14,
       height = 9,
       units = "cm")


Collect_all_tk <- rbind(ADONA_TK_results_ll,
                        HFPO_DA_TK_results_all,
                        PFBA_TK_results_all ,
                        PFBS_df,
                        PFDA_df,
                        PFDoDA_df,
                        PFHpA_df,
                        PFOA_df,
                        PFOA_df,
                        PFHxS_df,
                        PFNA_df,
                        PFOA_df,
                        PFOS_df,
                        PFPeA_df,
                        PFPeS_df,
                        PFTrDA_df,
                        PFUnDA_df)

ADONA_TK_results_all <-ADONA_TK_results_all |> mutate(PFAS = "ADONA")
HFPO_DA_TK_results_all <- HFPO_DA_TK_results_all |> mutate(PFAS = "HFPO_DA")
PFBA_TK_results_all <- PFBA_TK_results_all |> mutate(PFAS = "PFBA")
PFBS_TK_results_all <- PFBS_TK_results_all |> mutate(PFAS = "PFBS") 
PFDA_TK_results_all <- PFDA_TK_results_all |> mutate(PFAS = "PFDA")
PFDoDA_TK_results_all <- PFDoDA_TK_results_all |> mutate(PFAS = "PFDoDA") 
PFHpA_TK_results_all <- PFHpA_TK_results_all |> mutate(PFAS = "PFHpA") 
PFHpS_TK_results_all <- PFHpS_TK_results_all |> mutate(PFAS = "PFHpS")
PFHxA_TK_results_all <- PFHxA_TK_results_all |> mutate(PFAS = "PFHxA")
PFHxS_TK_results_all <- PFHxS_TK_results_all |> mutate(PFAS = "PFHxS")
PFNA_TK_results_all <- PFNA_TK_results_all |> mutate(PFAS = "PFNA")
PFOA_TK_results_all <- PFOA_TK_results_all |> mutate(PFAS = "PFOA")
PFOS_TK_results_all <- PFOS_TK_results_all |> mutate(PFAS = "PFOS")
PFPeA_TK_results_all <- PFPeA_TK_results_all |> mutate(PFAS = "PFPeA")
PFPeS_TK_results_all <- PFPeS_TK_results_all |> mutate(PFAS = "PFPeS")
PFTrDA_TK_results_all <- PFTrDA_TK_results_all |> mutate(PFAS = "PFTrDA")
PFUnDA_TK_results_all <- PFUnDA_TK_results_all |> mutate(PFAS = "PFUnDA")


Collect_all_tk <- rbind(ADONA_TK_results_all, 
               HFPO_DA_TK_results_all,
               PFBA_TK_results_all,
               PFBS_TK_results_all,
               PFDA_TK_results_all,
               PFDoDA_TK_results_all,
               PFHpA_TK_results_all,
               PFHpS_TK_results_all,
               PFHxA_TK_results_all,
               PFHxS_TK_results_all,
               PFNA_TK_results_all,
               PFOA_TK_results_all,
               PFOS_TK_results_all,
               PFPeA_TK_results_all,
               PFPeS_TK_results_all,
               PFTrDA_TK_results_all,
               PFUnDA_TK_results_all) 

write.csv(Collect_all_tk, file = here(OUTPUT,"Monte_Carlo_TK_results.csv"), row.names = FALSE)

HL_predicted <- Collect_all_tk |> 
  filter(Parameter == "HL", 
         PFAS != "PFTrDA")
Evaluation_data <- read_excel(here("Input","Evaluation_data.xlsx")) %>% 
  select(!c(PFAS_name, CAS_n)) %>% 
  rename(PFAS = Acronym)
Evaluation_data$PFAS = str_replace(Evaluation_data$PFAS, "HFPO-DA", "HFPO_DA")  

HL_measured <- Evaluation_data %>% 
  filter(Parameter == "Halflife",
         PFAS != "PFTrDA")

THL_measured <- HL_measured %>% 
  group_by(PFAS) %>% 
  summarise(mean = mean(Value), 
            median = median(Value)) %>% 
  mutate(HL_length = ifelse(mean < 0.01, "very_short",
                            ifelse(mean < 1, "short", 
                                   ifelse(mean < 5, "medium", "long")))) %>% 
  mutate(HL_length = factor(HL_length, levels = c("very_short", "short", "medium", "long"))) |> 
  select(PFAS, HL_length)

HL_measured <- full_join(HL_measured, THL_measured)

get_Category <- HL_measured |> select(PFAS, Category_OECD) |> distinct()
HL_predicted <- full_join(HL_predicted, THL_measured) 
HL_predicted <- full_join(HL_predicted, get_Category)

Plot_HL <- ggplot(HL_measured, aes(PFAS, Value)) +
  geom_jitter(color = "grey50", 
              width = 0.2, alpha = 0.4) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_errorbar(data = HL_predicted, aes(x = PFAS,
                      y = mean,
                      colour = Category_OECD,
                      ymin = q25,
                      ymax = q75),
                width = .1,
                linewidth = .5,
                alpha = 0.8) +
  geom_point(data = HL_predicted, 
             aes(x = PFAS,
                 y = mean,
                 colour = Category_OECD),
             size = 1.5,
             shape = 17) +
  facet_wrap(~ HL_length, scales = "free", nrow = 1 ) +
  labs(x = "",
       y = "Half life (years)") +
  theme_quick +
  theme(axis.text.x = element_text(angle = 45),
        strip.text = element_blank())


Plot_HL <- Plot_HL + CP_theme +
  theme(legend.position = "bottom")

ggsave(plot = Plot_HL,
       file = here(OUTPUT,"Plot_HL.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

HL_predicted_n <- HL_predicted |> select(!HL_length) |> 
  mutate(HL_length = ifelse(mean < 0.1, "very_short",
                            ifelse(mean < 1, "short", 
                                   ifelse(mean < 5, "medium", "long")))) %>% 
  mutate(HL_length = factor(HL_length, levels = c("very_short", "short", "medium", "long"))) 

Plot_predicted_HL <- ggplot(HL_predicted_n, aes(PFAS, mean,
                                                colour = Category_OECD)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_errorbar(aes(ymin = q25,
                    ymax = q75),
                width = .1,
                linewidth = .5,
                alpha = 0.8) +
  geom_point(size = 1.5,
             shape = 17) +
  facet_wrap(~ HL_length, scales = "free", nrow = 1 ) +
  labs(x = "",
       y = "Half life (years)") +
  theme_quick +
  theme(axis.text.x = element_text(angle = 45),
        strip.text = element_blank())


Plot_predicted_HL <- Plot_predicted_HL + CP_theme +
  theme(legend.position = "bottom")

ggsave(plot = Plot_predicted_HL,
       file = here(OUTPUT,"Plot_predicted_HL.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

ADONA_df <- ADONA_df |> mutate(PFAS = "ADONA") |> select(time, CP, PFAS)
HFPO_DA_df <- HFPO_DA_df |> mutate(PFAS = "HFPO_DA") |> select(time, CP, PFAS)
PFBA_df <- PFBA_df |> mutate(PFAS = "PFBA") |> select(time, CP, PFAS) 
PFBS_df <- PFBS_df |> mutate(PFAS = "PFBS") |> select(time, CP, PFAS)
PFDA_df <- PFDA_df |> mutate(PFAS = "PFDA") |> select(time, CP, PFAS)
PFDoDA_df <- PFDoDA_df |> mutate(PFAS = "PFDoDA") |> select(time, CP, PFAS)
PFHpA_df <- PFHpA_df |> mutate(PFAS = "PFHpA") |> select(time, CP, PFAS)
PFHpS_df <- PFHpS_df |> mutate(PFAS = "PFHpS")|> select(time, CP, PFAS)
PFHxA_df <- PFHxA_df |> mutate(PFAS = "PFHxA") |> select(time, CP, PFAS)
PFHxS_df <- PFHxS_df |> mutate(PFAS = "PFHxS") |> select(time, CP, PFAS)
PFNA_df <- PFNA_df |> mutate(PFAS = "PFNA") |> select(time, CP, PFAS)
PFOA_df <- PFOA_df |> mutate(PFAS = "PFOA") |> select(time, CP, PFAS)
PFOS_df <- PFOS_df |> mutate(PFAS = "PFOS") |> select(time, CP, PFAS)
PFPeA_df <- PFPeA_df |> mutate(PFAS = "PFPeA") |> select(time, CP, PFAS)
PFPeS_df <- PFPeS_df |> mutate(PFAS = "PFPeS") |> select(time, CP, PFAS)
PFTrDA_df <- PFTrDA_df |> mutate(PFAS = "PFTrDA") |> select(time, CP, PFAS)
PFUnDA_df <- PFUnDA_df |> mutate(PFAS = "PFUnDA") |> select(time, CP, PFAS)


Collect_all_CP <- rbind(ADONA_df, 
                        HFPO_DA_df,
                        PFBA_df,
                        PFBS_df,
                        PFDA_df,
                        PFDoDA_df,
                        PFHpA_df,
                        PFHpS_df,
                        PFHxA_df,
                        PFHxS_df,
                        PFNA_df,
                        PFOA_df,
                        PFOS_df,
                        PFPeA_df,
                        PFPeS_df,
                        PFTrDA_df,
                        PFUnDA_df) 

ggplot(ADONA_df, aes(time, CP)) +
  geom_path()+
  # facet_wrap(~PFAS) +
  theme_quick


