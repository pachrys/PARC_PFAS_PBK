# --------------------------------------------------------------------------- #
# FUNCTIONS FOR CALCULATING PBK MODEL PARAMETERS
# By: Chrysanthi Pachoulide
# Date: 05-11-2025
# --------------------------------------------------------------------------- #

# !! max BW in Physio_params is 76 for Female

BASE_PARAMS <- function(expAGE = NULL, expBW = NULL, sex = NULL) { 
  
  ## Function conditions ####
  if (is.na(expAGE) && is.na(expBW)) {
    expBW <- 70  # Default body weight
  } 
  
  if (!is.na(expAGE) && expAGE > 80) {
    warning("Age above 80 which is the max age in the lifestage model. Physiological parameters assumed as those of an 80 year old")
    expAGE <- 80
  } 
  
  if (is.na(sex)) {
    sex <- "M"  # Default sex is male...
  }
  # 25ml plasma per cycle (29.2 days per cycle until; 12.5 cycles per year ); Ruark (12.89 +- 9.11 cycles per year)
  
  # Menstruation, reference Ruark et al 2016 http://dx.doi.org/10.1016/j.envint.2016.11.030
  if (sex == "F" && expAGE >= 12 && expAGE < 45){  
    CL_menses <-  25e-3/29.2 # L/day
  } else if (sex == "F" && expAGE >= 45 && expAGE <= 51){
    CL_menses <- 27.1e-3/35.1 # L/day
  } else {
    CL_menses <- 0 # L/day
  }
  
 
  ## Input
  Physio_params <- Physio.c
  
  suffix <- ifelse(sex == "F", "_F", "_M")
  
  if (!is.na(expAGE)) { # filter based on age
    Physio_params <- Physio_params %>% 
      mutate(age_diff = abs(age - expAGE)) %>% 
      filter(age_diff == min(age_diff)) %>% 
      slice(1)
  } else if (!is.na(expBW)) { # filter based on BW if age is not provided
    bw_col <- paste0("BW", suffix)
    Physio_params <- Physio_params %>%
      mutate(bw_diff = abs(!!sym(bw_col) - expBW)) %>% 
      filter(bw_diff == min(bw_diff)) %>% 
      slice(1)
  }
  
  Physio_params <- Physio_params %>% select(ends_with(suffix))
  
  
  ## Physiological Parameters ####
  Physio_params <- Physio_params %>% 
    mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
    mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
  
  BW <- if (!is.na(expBW)) {
    expBW # BW is the actual BW at the time of measurement
    } else {
    Physio_params[[paste0("BW", suffix)]]}
  # BW <- Physio_params[[paste0("BW", suffix)]]   # Body weight calculated only from the lifestage equations
  
  QC <- Physio_params[[paste0("CardOut", suffix)]]              # L/d, This is corrected for hematocrit already so it's plasma
  Hct <- Physio_params[[paste0("Hct", suffix)]]                 # Hematocrit
  
  BH <- Physio_params[[paste0("BH", suffix)]]                   # Height (cm)
  BSA <- exp(-3.75 + 0.42*log(BH)+0.52*log(BW))*1e4             # Body Surface area (cm2)
  
  GFR <- Physio_params[[paste0("GFR", suffix)]]                 # Age/creatinine dependent glomerular filtration rate (L/day)
  
  SAlb <- Physio_params[[paste0("SAlb", suffix)]]               # Serum albumin concentration (g_HSA /L_serum, or mg/ml)
  ### Organ volumes -------------------------
  
  
  # Intestine
  VIc <- Physio_params[[paste0("V_gutFraction", suffix)]] 
  
  # Intestinal lumen is taken as a compartment outside the intestine (weight of Intestinal lumen is outside of the bodyweight)
  L = 280                           # cm, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b
  R = 1.75                          # cm, (adult of 70kg mean value) Willmann2004 doi: 10.1021/jm030999b
  VIL = pi*L*(R^2)/1000             # L, cm3/1000 #Value is the same as Punt et al. 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  VILc = VIL/70                     # deriving the constant by dividing by BW, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b
  
  SA_SIc = 2*pi*R*L*25/70           # cm2, amplification factor of 25 for the microvilli in the intestinal lumen, Willmann2004 doi: 10.1021/jm030999b
  # SA_SIc = 70.1 * 1e3/70            # cm2 converted to a constant by dividing by BW, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b
  
  # Liver
  VLc <- Physio_params[[paste0("V_liverFraction", suffix)]]  # cm2, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b        
  
  # # Liver is divided in intracellular and extracellular compartments. 
  # # Extracellular compartment combines both the vascular and interstitial space
  # # Fractional volume of intracellular space was taken from Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  VL_icc = 0.573                   # Fractional volume of intracellular space in the liver
  VL_ecc = 0.161# 1- VL_icc #0.161             # Fractional volume of extracellular space in the liver
  
  
  # Kidney
  VKc <- Physio_params[[paste0("V_kidneyFraction", suffix)]]
  
  # Kidney is divided in to proximal tubule and rest of kidney. Each compartment is divided in to tissue and lumen (or filtrate). 
  # This separation is necessary because of the active reabsorption happening in the proximal tubule, but also due to the differences in composition and pH of primary and terminal urine.
  # Fractional volumes were recalculated from Pletz et al. 2020 https://doi.org/10.1016/j.comtox.2021.100172
  
  VPTc <- 0.398                   # Fractional volume proximal tubule 
  VPTBc <- 0.137                  # Fractional volume proximal tubule blood
  VPTCc <- 0.460                  # Fractional volume proximal tubule cell
  VPTTc <- VPTBc*(1-Hct) + VPTCc  # Fractional volume proximal tubule tissue
  VPTLc <- 0.403                  # Fractional volume proximal tubule lumen
  
  VRKLc <- 0.283                   # Fractional volume rest of kidney lumen
  
  
  # Adipose
  VAc <- Physio_params[[paste0("V_adiposeFraction", suffix)]]
  MAc <- Physio_params[[paste0("AdiposeMass", suffix)]]
  
  
  # Plasma
  VPc <- Physio_params[[paste0("V_bloodFraction", suffix)]]
  
  
  # Total body volume
  VTotc <- Physio_params$VolumesSum                  
  
  ### Organ blood flows -------------------------
  QTotc <- Physio_params$BloodFlowSum
  QIc <- Physio_params[[paste0("Q_gutFraction", suffix)]]
  QLc <- Physio_params[[paste0("Q_liverFraction", suffix)]]
  QKc <- Physio_params[[paste0("Q_kidneyFraction", suffix)]] 
  QAc <- Physio_params[[paste0("Q_adiposeFraction", suffix)]]
  
  
  ### Other physiological flows and constants -------------------------
  
  tco = 0.14*24               # /d, Bowel residence/transit time in the colon, (24*/h) or 7h willmann2004 doi: 10.1021/jm030999b
  
  QUrc = 0.022                # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
  GFRc = 0.18                 # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
  QTc = 0.36 # QT normalised to GFR, based on reference values of GFR 120 and QT of 43.2 from Scotcher et al. 2016 https://doi.org/10.1016/j.ejps.2016.03.018 and Pletz
  
  ### Physiological pHs in different matrices -------------------------
  
  pH_P <- 7.4     # plasma
  pH_IL <- 7      # intestinal lumen, average
  
  ### Albumin concentrations in different matrices -------------------------

  # From Akihiro Tojo and Satoshi Kinugasa 2012 doi:10.1155/2012/481520
  # In the same paper: the proximal tubule reabsorbes 71% of albumin, while LoH and DT 23% and the collecting duct 3%
  Calb_P <- SAlb      # g/L plasma (see above for values)
  Calb_PTL <- 14.4e-3 # g/L proximal tubule 
  Calb_IL <- 0.007e-3    # g/L intestinal lumen
  Calb_L_ec <- 0.0012e-3 # g/L Liver extracellular 

  # Based on the Poulin and Theil 2009, below Table 6
  # Albumin ratio
  R_PTL <- Calb_P/Calb_PTL # plasma to proximal tubule lumen albumin ratio
  R_L_ec <- Calb_P/Calb_L_ec # plasma to liver albumin ratio
  R_IL <- Calb_P/Calb_IL # plasma to intestinal lumen albumin ratio
  
  ## Chemical Specific ####
  

  MW <- PhysChemTK_params %>% 
    filter(Parameter == "MW") %>%
    pull(Value) %>% as.numeric()
  
  # Distribution coefficients to different matrices 
  logML <- PhysChemTK_params %>% 
    filter(Parameter == "logD_MembraneLipidWater") %>% pull(Value) %>% as.numeric()
  logSP <- PhysChemTK_params %>% 
    filter(Parameter == "logD_StructuralProteinWater") %>% pull(Value) %>% as.numeric()
  logALB <- PhysChemTK_params %>% 
    filter(Parameter == "logD_AlbuminWater") %>% pull(Value) %>% as.numeric()
  logSL <- PhysChemTK_params %>% 
    filter(Parameter == "logD_StorageLipidWater") %>% pull(Value) %>% as.numeric()
  logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954
  # !! need to change this !!
  logGLOB <- PhysChemTK_params %>% 
    filter(Parameter == "logD_GlobulinWater") %>% pull(Value) %>% as.numeric()
  
  k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
  k_SP <- 10^logSP       # structural protein:water partition coefficient
  k_ALB <- 10^logALB     # albumin:water partition coefficient
  k_SL <- 10^logSL       # structural lipids:water partition coefficient
  k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient
  k_GLOB <- 10^logGLOB   # globulin:water partition coefficient
  
  ### Fraction unbound in plasma -------------------------
  # Equation for calculating fraction unbound in plasma was taken from https://doi.org/10.1021/acs.est.3c07415
  
  DAlb <- 1.36 #kg/L density of albumin
  D_ALB <- k_ALB*DAlb # L_w/L_hsa; distribution coefficient
  FV_Alb <- SAlb*10^-3 # kg/kg of albumin in serum

  DGlob <- 1.20 #kg/L density of globulin
  D_Glob <- k_GLOB*DGlob # L_w/L_hsa; distribution coefficient
  FV_Glob <- 0.031 # kg/kg of globulin in serum
  
  FV_W <- 1 - FV_Alb - FV_Glob # serum mass balance
  
  fup <- 1/(1 + (D_ALB*(FV_Alb/FV_W)) + (D_Glob*(FV_Glob/FV_W)))

  
  ### Partition coefficients -------------------------
  # Code for calculating the partition coefficients is a slightly modified version from Utsey et al. 2020 https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_Schmitt.R

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
           (k_ALB*Physio.Tissues$f_ALB) +
           (k_FABP*Physio.Tissues$f_FABP)
  ) #* fup Removing fup from here, Kps will be recalculated directly in the model for correct sensitivity analysis
  
  Names <- Physio.Tissues$Tissue %>% substr(1,2)
  Names <- paste("Kp",Names ,sep="")
  Kp.df <- as.data.frame(list(Tissue = Names)) %>%
    mutate(Value = Kp) %>% 
    pivot_wider(names_from = Tissue, values_from = Value)
  
  
  # Calculate Plasma/Rest of the body partition coefficient
  # KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  KpRe <- (Kp.df$KpSk * Physio_params[[paste0("V_skinFraction", suffix)]] +
             Kp.df$KpLu * Physio_params[[paste0("V_lungFraction", suffix)]] +
             Kp.df$KpSt * Physio_params[[paste0("V_stomachFraction", suffix)]] +
             Kp.df$KpBr * Physio_params[[paste0("V_brainFraction", suffix)]] +
             Kp.df$KpHe * Physio_params[[paste0("V_heartFraction", suffix)]] + 
             Kp.df$KpMu * Physio_params[[paste0("V_muscleFraction", suffix)]] +
             Kp.df$KpSp * Physio_params[[paste0("V_spleenFraction", suffix)]] +
             Kp.df$KpGo * Physio_params[[paste0("V_reproFraction", suffix)]] +
             Kp.df$KpBo * Physio_params[[paste0("V_boneFraction", suffix)]]) / (
               Physio_params[[paste0("V_skinFraction", suffix)]] +
               Physio_params[[paste0("V_lungFraction", suffix)]] +
                 Physio_params[[paste0("V_stomachFraction", suffix)]] +
                 Physio_params[[paste0("V_brainFraction", suffix)]] +
                 Physio_params[[paste0("V_heartFraction", suffix)]] +
                 Physio_params[[paste0("V_muscleFraction", suffix)]] +
                 Physio_params[[paste0("V_spleenFraction", suffix)]] +
                 Physio_params[[paste0("V_reproFraction", suffix)]] +
                 Physio_params[[paste0("V_boneFraction", suffix)]])  
  Kp.df$KpRe <- KpRe
  
  PIc <- Kp.df$KpIn   # Intestine
  PLc <- Kp.df$KpLi   # Liver
  PKc <- Kp.df$KpKi   # Kidney
  PAc <- Kp.df$KpAd   # Adipose
  PRc <- Kp.df$KpRe   # Rest
  
  ### Uptake from the gastro-intestinal duct -------------------------
  
  # Input data, in vitro clearance
  Papp_SI <- PhysChemTK_params %>% 
    filter(Parameter == "Papp") %>% 
    pull(Value) %>% as.numeric() # cm/s
  
  CL_OATP2B1c <- PhysChemTK_params %>% 
      filter(Parameter == "OATP2B1") %>%
      pull(Value) %>% as.numeric()*1e-3 # L/min/mg protein

  # IVIVE scaling factors
  # OATP2B1_I = 0.5 # pmol/mg protein
  # REF_OATP2B1_I = OATP2B1_I/0.812                     # relative expression factor OATP2B1, calculated from: OATP2B1_vivo/OATP2B1_vitro (pmol/mg membrane protein / pmol/mg membrane protein); Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables7 and s6, ref23
  # mgProtein.mgI = 0.13                           # mg of protein per mg intestine; Utsey et al. 2020 http://dx.doi.org/10.1124/dmd.120.090498
  SF_OATP2B1_I = 79248.77 #mgProtein.mgI * REF_OATP2B1_I * 1e6       # (mg liver * 1e6 -> kg liver)
  
  
  ### Uptake to the liver -------------------------
  
  # Input data, in vitro clearance
  CL_OATP1B1c <- PhysChemTK_params %>% 
    filter(Parameter == "OATP1B1") %>%
    pull(Value) %>% as.numeric()*1e-3 # L/min/mg protein
  
  CL_OATP1B3c <- PhysChemTK_params %>% 
    filter(Parameter == "OATP1B3") %>%
    pull(Value) %>% as.numeric()*1e-3 # L/min/mg protein
  
  CL_NTCPc <-   PhysChemTK_params %>%
    filter(Parameter == "NTCP") %>%
    pull(Value) %>% as.numeric()*1e-3 # L/min/mg protein
  
  # IVIVE scaling factors fmol/ug protein
  
  # OATP1B1_L = 4.5 # pmol/mg membrane protein
  # OATP1B3_L = 1.75 # pmol/mg membrane protein
  # OATP2B1_L = 1.5 # pmol/mg membrane protein
  # NTCP_L = 2 # pmol/mg membrane protein
  # BSEP_L = 1.75 # pmol/mg membrane protein
  # 
  # REF_OATP1B1 = OATP1B1_L/0.120                #mean(c(19.6,2,9.7))             # relative expression factor OATP1B1, calculated from: OATP1B1_vivo / OATP1B1_vitro  (pmol/mg membrane protein/pmol/mg membrane protein); Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables7 and s6, ref23
  # REF_OATP1B3 = OATP1B3_L/0.719                #mean(c(2.7,1,6.3))             # relative expression factor OATP1B3, calculated from: OATP1B3_vivo / OATP1B3_vitro  (pmol/mg membrane protein/pmol/mg membrane protein); Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables7 and s6, ref23
  # REF_OATP2B1_L = OATP2B1_L/0.812                     # relative expression factor OATP2B1, calculated from: OATP2B1_vivo/OATP2B1_vitro (pmol/mg membrane protein / pmol/mg membrane protein); Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables7 and s6, ref23
  # REF_NTCP = NTCP_L/1.202  # relative expression factor NTCP, calculated from: OATP1B3_vivo / OATP1B3_vitro  (pmol/mg membrane protein/pmol/mg membrane protein); Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables7 and s6, ref23
  # 
  # mgProtein.mgL = 0.35 #Prasad et al. 2016 DOI:10.1002/cpt.409 # 0.18                                 # mg of protein per mg liver; Utsey et al. 2020 http://dx.doi.org/10.1124/dmd.120.090498
  SF_OATP1B1 = 11017310  #mgProtein.mgL * REF_OATP1B1 * 1e6       # (mg liver * 1e6 -> kg liver)
  SF_OATP1B3 = 770384.3 #mgProtein.mgL * REF_OATP1B3 * 1e6       # (mg liver * 1e6 -> kg liver)
  SF_OATP2B1 = 356421.6 #mgProtein.mgL * REF_OATP2B1_L * 1e6       # (mg liver * 1e6 -> kg liver)
  SF_NTCP = 392315.4 #mgProtein.mgL * REF_NTCP * 1e6             # (mg liver * 1e6 -> kg liver)
  

  
  ### Biliary clearance -------------------------
  # Input data, in vitro clearance
  Vmax_BSEPc <- mean(c(5.8, 7, 8.4))                            # umol/min/mg BSEP, Average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011 
  Km_BSEP <- mean(c(4.3, 25.9, 19))                            # umol/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  
  CL_BSEPc <- Vmax_BSEPc/Km_BSEP
  
  # MW_BSEP <- 146400 #g/mol ref: The Human Protein Atlas (https://www.proteinatlas.org/ENSG00000073734-ABCB11/structure)
  # pmolBSEP.mgprot <- 2.2/15/MW_BSEP * 1000000000
  # 
  # 
  # 
  # # IVIVE 
  # mgBSEP.HC <- 0.839 * 140000 * 1e-9         # mg BSEP/1e6 hepatocytes (0.839 pmole BSEP/1e6 hepatocytes (amound of BSEP per hepatocyte) * 140000 g/mole (MW BSEP) -> pg * 1e-9 -> mg) de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011]
  # HC.GL <- 99                                 # 1e6 hepatocytes in the liver ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  SF_BSEP <- 297.3746 #mgBSEP.HC * HC.GL * 1e3          # Scaling factor for BSEP mediated hepatic efflux for GCA and GCDC ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  
  
  ### Renal Clearance -------------------------
  
  ## Active transport
  
  CL_OAT4c <- PhysChemTK_params %>%
    filter(Parameter == "OAT4") %>%
    pull(Value) %>% as.numeric()*1e-3 # L/min/mg protein
    # 
  CL_OAT3c <- PhysChemTK_params %>% 
    filter(Parameter == "OAT3") %>%
    pull(Value) %>% as.numeric()*1e-3 # L/min/mg protein
    # 
  CL_OAT1c <- PhysChemTK_params %>% 
    filter(Parameter == "OAT1") %>%
    pull(Value) %>% as.numeric()*1e-3 # L/min/mg protein
  
  # IVIVE 
  # REF_OAT4 <- 0.29/0.063              # relative expression factor OAT, calculated from: OAT4_vivo / OAT4_vitro  (pmol/mg membrane protein/ pmol/mg membrane protein); in vivo: Al-Madjoub et al 2021 (https://doi.org/10.1002/cpt.2396) table S3 (prev 1.9; mean 3.2, sd 4.6, CV% 145% n = 9 #), in vitro:; Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables5
  # REF_OAT3 <- 3.50/0.711              # relative expression factor OAT, calculated from: OAT4_vivo / OAT4_vitro  (pmol/mg membrane protein/ pmol/mg membrane protein); in vivo: Al-Madjoub et al 2021 (https://doi.org/10.1002/cpt.2396) table S3 (mean 9.6, sd 9.8, CV 102%, n = 16), in vitro:; Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables5
  # REF_OAT1 <- 5.33/2.596              # relative expression factor OAT, calculated from: OAT4_vivo / OAT4_vitro  (pmol/mg membrane protein/ pmol/mg membrane protein); in vivo: Al-Madjoub et al 2021 (https://doi.org/10.1002/cpt.2396) table S3 (mean 10.2, sd 10.1, CV 99%, n = 15), in vitro:; Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables5
  # 
  # mgProtein.mgK <- 0.21                       # mg of protein per mg kidney; Utsey et al. 2020 http://dx.doi.org/10.1124/dmd.120.090498
  SF_OAT4 <- 21047620 #mgProtein.mgK * REF_OAT4 * 1e6   # (mg kidney * 1e6 -> kg kidney)
  SF_OAT3 <- 4638537 #mgProtein.mgK * REF_OAT3 * 1e6   # (mg kidney * 1e6 -> kg kidney)
  SF_OAT1 <- 1329353 #mgProtein.mgK * REF_OAT1 * 1e6   # (mg kidney * 1e6 -> kg kidney)

  
  ### Final parameter constants -------------------------
  
  parm.c <- list(
   BW = BW,
   QC = QC,
   Hct = Hct,
   VIc = VIc,
   VILc = VILc,
   SA_SIc = SA_SIc,
   VLc = VLc,
   VL_icc = VL_icc,
   VL_ecc = VL_ecc,
   VKc = VKc,
   VPTc = VPTc, 
   VPTTc = VPTTc, 
   VPTLc = VPTLc, 
   VRKLc = VRKLc,
   VAc = VAc,
   MAc = MAc,
   VPc = VPc,
   QIc = QIc,
   QLc = QLc,
   QKc = QKc,
   QAc = QAc,
   QUrc = QUrc,
   GFR = GFR,
   QTc = QTc,
   tco = tco,
   R_PTL = R_PTL,
   R_L_ec = R_L_ec,
   MW = MW,
   fup = fup,
   PIc = PIc,
   PLc = PLc,
   PKc = PKc,
   PAc = PAc,
   PRc = PRc,
   Papp_SI = Papp_SI,
   CL_OATP2B1c = CL_OATP2B1c,
   SF_OATP2B1_I = SF_OATP2B1_I,
   CL_OATP1B1c = CL_OATP1B1c,
   CL_OATP1B3c = CL_OATP1B3c,
   CL_NTCPc = CL_NTCPc,
   SF_OATP1B1 = SF_OATP1B1,
   SF_OATP1B3 = SF_OATP1B3,
   SF_OATP2B1 = SF_OATP2B1,
   SF_NTCP = SF_NTCP,
   CL_BSEPc = CL_BSEPc,
   SF_BSEP = SF_BSEP,
   CL_OAT4c = CL_OAT4c,
   CL_OAT3c = CL_OAT3c,
   CL_OAT1c = CL_OAT1c,
   SF_OAT4 = SF_OAT4,
   SF_OAT3 = SF_OAT3,
   SF_OAT1 = SF_OAT1,
   CL_menses = CL_menses
  )
  
  return(parm.c)
}
