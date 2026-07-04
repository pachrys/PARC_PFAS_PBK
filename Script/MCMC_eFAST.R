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
OUTPUT <- here("Output", "MonteCarlo") #format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
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
  filter(Predicted == "yes") |> 
  select(PFAS, Parameter, Value, SD, Distribution)

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
  
  # print("Intestine - passive absorption")
  CL_IL <- (Papp_SI*SA_SI*1e-3)    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/d -> L/d /1000)
  
  CL_OATP2B1_I <- CL_OATP2B1c*SF_OATP2B1_I*VI # L/d
  # print("Intestine - active absorption")
  # print(QI)
  # print(CL_OATP2B1_I)
  
  # print("Liver uptake")
  CL_OATP1B1 <- CL_OATP1B1c*SF_OATP1B1*VL_ec   # L/d
  CL_OATP1B3 <- CL_OATP1B3c*SF_OATP1B3*VL_ec   # L/d
  CL_OATP2B1_L <- CL_OATP2B1c*SF_OATP2B1*VL_ec # L/d
  CL_NTCP <- CL_NTCPc*SF_NTCP*VL_ec            # L/d    
  # print(QI+QL)
  # print(CL_OATP1B1 + CL_OATP1B3 + CL_OATP2B1_L)
  
  # print("Liver excretion")
  CL_BSEP <- CL_BSEPc*SF_BSEP*VL_ic  # L/d
  # print(QI + QL)
  # print(CL_BSEP)
  
  # print("Kidney uptake")
  CL_OAT3 <- CL_OAT3c*SF_OAT3*VPT  # L/d 
  CL_OAT1 <- CL_OAT1c*SF_OAT1*VPT  # L/d 
  # print(QK)
  # print(CL_OAT1 + CL_OAT3)
  
  # print("Kidney out")
  CL_OAT4 <- CL_OAT4c*SF_OAT4*VPT  # L/d 
  # print(GFR)
  # print(CL_OAT4)
  
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
  
  PBK.out <- lsoda(y = A_init,
                   times = seq(0, 1825, by = 10),
                   func = PBK.model, 
                   parms = pars)
  
  PBK_OUTPUT <- as.data.frame(PBK.out) %>% 
    mutate(CP = ifelse(CP < 0-1000, 0, CP)) 
  
 return(PBK_OUTPUT)
  
}

Monte_Carlo_PFAS_design <- list()
Monte_Carlo_PFAS_RESULTS <- list()
Monte_Carlo_PFAS_y <- list()


# for (i in 1:nPFAS) {
  
  PFAS_i <- "PFOA"#PFAS_list[i]
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
           Distribution = ifelse(Distribution == "lognormal", "qlnorm", "qunif"))

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
    n = 100,# Integer giving the sample size, i.e. the length of the discretization of the s-space 
    q = q,
    q.arg = q.arg)

  
  dim(Experience$X) # dataframe of: n factors(number of parameters) * n n(1000) number of observations, of n factors(number of parameters) number of variables
  design <- Experience$X
  write.csv(design, file = here(OUTPUT, paste0(PFAS_i,"Initial_MC_ExpDesign.csv")), row.names = FALSE)
  
  # check distributions
  design |> pivot_longer(everything(), 
                         names_to = "Parameter", 
                         values_to = "Value") |> 
    ggplot(aes(x=Value))+
    geom_histogram()+
    facet_wrap(~Parameter, scales = "free")+
    labs(title = "Monte Carlo parameter distributions",
         x = "",
         y = "") +
    theme_quick 

  # clean by removing infinite values
  Experience$X <- Experience$X |> mutate(across(everything(), ~na_if(., Inf))) |> drop_na()
  clean_design <- Experience$X
  write.csv(clean_design, file = here(OUTPUT, paste0(PFAS_i,"MC_ExpDesign.csv")), row.names = FALSE)
  # Monte_Carlo_PFAS_design[[i]] <- clean_design
  
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
         file = here(OUTPUT, paste0(PFAS_i,"MC_ExpDesign.png")),
         dpi = 600,
         width = 12,
         height = 9,
         units = "cm")
  
  results <- apply(clean_design, 1, function(row) {
    run_Monte_Carlo(row, parm.phys = parm.phys)
  })
  
  y <- results
  # Monte_Carlo_PFAS_y[[i]] <- results
  save(y, file = here(OUTPUT, paste0(PFAS_i,"Monte_Carlo_y.RData")))
  
# }
  names(y) <- paste0("iter",(seq(1,length(y), by = 1)))
  
  Monte_Carlo_RESULTS <- do.call(rbind, y)
  # Morris_RESULTS <- cbind(rownames(Morris_RESULTS), 
  #                         data.frame(Morris_RESULTS, 
  #                                    row.names=NULL)) %>%
  #   separate("rownames(Morris_RESULTS)", 
  #            sep = "_", into = c("SF", "Parameter")) %>% 
  # select(!SF)
  
  # Cmax <- max(PBK_OUTPUT$CP)
  # message(glue("-Cmax {Cmax}---"))
  # 
  # conc <- PBK_OUTPUT$CP
  # 
  # time <- PBK_OUTPUT$time
  # Tmax <- filter(PBK_OUTPUT, 
  #                CP == max(PBK_OUTPUT$CP)) %>% pull(time)
  # message(glue("-Tmax {Tmax}---"))
  # 
  # tlast <- max(time[conc > 0])
  # 
  # half_life <- pk.calc.half.life(
  #   PBK_OUTPUT$CP,
  #   PBK_OUTPUT$time
  # )
  # half_life <- half_life$half.life
  # message(glue("-HL {half_life}---"))
  # 
  # print(plot(PBK_OUTPUT$time, PBK_OUTPUT$CP))
  # 
  # return(c(Cmax = Cmax, 
  #          half_life = half_life
  # ))
