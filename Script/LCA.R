  # --------------------------------------------------------------------------- #
  # Script for running the the Local Sensitivity Analysis
  # CP, 18-06-2026
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
  OUTPUT <- here("Output", "LSA") #, format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
  dir.create(OUTPUT, recursive = TRUE)

  source(here("Script", "PBK_model_LSA.R"), local = TRUE)  
  load(here("Input", "LSA_initial_params.RData"))
  # read.csv("Input", "LSA_parameters.csv")
  
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
  
  LSA_PFAS_RESULTS <- list(
    CALC_Parameters = vector("list", nPFAS), # list of length nPFAS
    Sensi_OUTPUT = vector("list", nPFAS) #
  )
  
  names(LSA_PFAS_RESULTS$CALC_Parameters) <- PFAS_list
  names(LSA_PFAS_RESULTS$Sensi_OUTPUT) <- PFAS_list
  
  FinalParam <- NULL
  
  LSA_parms <- c("SF_OATP2B1_I", "SF_OATP1B1", "SF_OATP1B3", "SF_OATP2B1", "SF_NTCP", "SF_BSEP", "SF_OAT4", "SF_OAT3", "SF_OAT1")
  
  # for (i in 1:nPFAS) {
    
    PFAS <- PFAS_list[7]
    message(glue("{PFAS}------------------------------------"))
    
    pars <- ModelParams[[PFAS]]
    
  # }
#   
  pars <- as.list(log10(unlist(pars)))
  

  solvePBK <- function(pars){
    
    
    derivs <- function(t, state, pars){
      with(as.list(c(state, pars)), {
        
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
        
        
        CL_OATP2B1_I <- CL_OATP2B1c*60*24*SF_OATP2B1_I*VI #L/d
        CL_OATP1B1 <- CL_OATP1B1c*60*24*SF_OATP1B1*VL_ec #L/d
        CL_OATP1B3 <- CL_OATP1B3c*60*24*SF_OATP1B3*VL_ec #L/d
        CL_OATP2B1_L <- CL_OATP2B1c*60*24*SF_OATP2B1*VL_ec #L/d
        CL_NTCP <- CL_NTCPc*60*24*SF_NTCP*VL_ec # L/d    
        CL_BSEP <- CL_BSEPc*60*24*SF_BSEP*VL_ic         # L/d
        CL_OAT4 <- CL_OAT4c*60*24*SF_OAT4*VPT           # L/d (min -> d; /mg protein 
        CL_OAT3 <- CL_OAT3c*60*24*SF_OAT3*VPT           # L/d (min -> d)
        CL_OAT1 <- CL_OAT1c*60*24*SF_OAT1*VPT           # L/d (min -> d)
        
        
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
        
        dOD <-  0 #- OD             # ug/d, Oral dose input OralD
        
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
        Atot <- OD +
          AIL + AI + AFe + 
          AL_ec + AL_ic +
          APTT + APTL + ARKT + ARKL + AUr +
          AA + 
          AR +
          AP +
          AMp
        
        dAin <- 0 #0.01 
        MB <- 0.1 - Atot + 1    
        
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
               dAin),
             c(CP = CP,
               MB = MB)
             )
        })
    }
    
    
    
    
    state <- c(OD = 0,
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
    
    tout    <- seq(0, 1825, by = 10)
    
    pars <- as.list(10^(unlist(pars)))
    
    PBK_OUTPUT <- as.data.frame(
      lsoda(y = state, 
            times = tout, 
            func = derivs,
            parms = pars
      ))
    
    PBK_OUTPUT <- PBK_OUTPUT %>% mutate(CP = ifelse(CP < 1e-100, 0, CP))
    # Cmax <- max(PBK_OUTPUT$CP)
    # conc <- PBK_OUTPUT$CP
    # 
    # time <- PBK_OUTPUT$time
    # Tmax <- filter(PBK_OUTPUT, CP == max(PBK_OUTPUT$CP)) %>% pull(time)
    # tlast <- max(time[conc > 0])
    # 
    # half_life <- pk.calc.half.life(
    #   conc,
    #   time,
    #   Tmax,
    #   tlast
    # )
    # half_life <- half_life$half.life
    # 
    return(PBK_OUTPUT)
    
    
  }
  out <- solvePBK(pars)
  
  # SENSI_all <- sensFun(
  #   func = SENSI_model,
  #   parms = parm.c, 
  #   tiny = 0.01)
  
  
  # lsa <- sensFun(func = SENSI_model,
  #                parms = parm.c,
  #                sensvar = "CP", 
  #                senspar = LSA_parms,
  #                tiny = 0.01 # change parameter by 1%
  #                  )
  # 
  # head(lsa)
  # summary(lsa)
  # 
  # LSArun <- function(parms){
  #     solve(model2, parms, inits=inits, events=ev, cores=4) # This is the solve function I use in the LSA?
  #   }
  
  SnsPBK <- sensFun(func = solvePBK, 
                    parms = pars,
                    senspar = LSA_parms,
                    sensvar = "CP", 
                    tiny = log10(99) #1.4
                    # ,10
                    # tiny = log10(c(19.80,
                    #          32.20,
                    #          102.00,
                    #          9600.00,
                    #          3200.00,
                    #          6830.05,
                    #          6.08,
                    #          495.00,
                    #          1506.53)
                    # )
                    )
  
  head(SnsPBK)
  plot(SnsPBK)
  plot(SnsPBK, type = "b", pch = 15:19, col = 2:6,
       main = "Sensitivity all vars")

  summary(SnsPBK)
  plot(summary(SnsPBK))
  
  
  
  Sns_summary <- summary(SnsPBK) %>% as.data.frame()
  Snsplot <- Sns_summary %>% select(Mean) %>%
    mutate(Parameter = rownames(Sns_summary)) %>%
    # filter(Mean >= 0.1 | Mean <= -0.1) %>%  
    mutate(Parameter = fct_reorder(Parameter, abs(Mean), .desc = TRUE)) %>%
    ggplot(aes(x = Mean, y = Parameter)) +
    geom_bar(stat = "identity") +
    coord_flip() + # Flip coordinates for better readability if there are many parameters
    labs(title = "Output: Plasma concentration", 
         x = "Mean sensitivity coefficient", 
         y = "Parameter") +
    CP_theme +
    theme( axis.text.x = element_text(angle = 45, hjust = 1)  )
  
  ggsave(plot = Snsplot,
         file = "Snsplot.png",
         dpi = 600,
         width = 12,
         height = 9,
         units = "cm")
  
 