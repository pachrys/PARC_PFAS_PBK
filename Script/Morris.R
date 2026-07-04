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
  library(readxl)
  
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
  OUTPUT <- here("Output", "Morris") #format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
  dir.create(OUTPUT, recursive = TRUE)
  
  # Model Code
  PBK.model <- function(t, state, pars){
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
  

  # SENSITIVITY FUNCTION ####
  # ---------------------------------------------------------------------------- #
  

  SENSI_model <- function(parm.Morris, pars){
    
    # Final parameters per test
    names(parm.Morris) <- ParametersMorris$Parameter
    parm.final <- pars
    print(parm.final)
    
    parm.final[names(parm.Morris)] <- parm.Morris
    print(parm.final)
    
    
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
    
    
    PBK.out <- lsoda(y = A_init,
                     times = seq(0, 1825, by = 10),
                     func = PBK.model, 
                     parms = parm.final)
    
    PBK_OUTPUT <- as.data.frame(PBK.out) %>% mutate(CP = ifelse(CP < 0-1000, 0, CP))

    Cmax <- max(PBK_OUTPUT$CP)
    message(glue("-Cmax {Cmax}---"))
    
    conc <- PBK_OUTPUT$CP
    
    time <- PBK_OUTPUT$time
    Tmax <- filter(PBK_OUTPUT, CP == max(PBK_OUTPUT$CP)) %>% pull(time)
    message(glue("-Tmax {Tmax}---"))
    
    tlast <- max(time[conc > 0])
    
    half_life <- pk.calc.half.life(
      PBK_OUTPUT$CP,
      PBK_OUTPUT$time,
      # Tmax,
      # tlast
    )
    half_life <- half_life$half.life
    message(glue("-HL {half_life}---"))
    
    print(plot(PBK_OUTPUT$time, PBK_OUTPUT$CP))
    
    
    return(c(Cmax = Cmax, 
             half_life = half_life
              ))
    
  }
  
  # Morris Parameters
  load(here("Input", "LSA_initial_params.RData"))
  
  InputMorrisParams <- read_excel(here("Input", "ProteinAbundance.xlsx"), 
                                 sheet = "Input_Morris")
  
  ParametersMorris <- InputMorrisParams %>% 
    select(!Tissue) %>% 
    rename("Parameter" = "Transporter",
           "value" = "Mean_SF",
           "binf" = "Minimum_SF", 
           "bsup" = "Maximum_SF")

  Morris_PFAS_design <- list()
  Morris_PFAS_RESULTS <- list()
  Morris_PFAS_y <- list()
  
  nPFAS <- length(PFAS_list)
  
  
  for (i in 1:nPFAS) {
    
    PFAS <- PFAS_list[i]
    message(glue("{PFAS}------------------------------------"))
    
    Parameters_df <- tibble(
      Parameter = names(ModelParams[[PFAS]]),
      Value = unlist(ModelParams[[PFAS]]))
    
    
    # SET MORRIS PARAMETERS ####
    # ---------------------------------------------------------------------------- #
    
    # Parameter upper and lower bounds
    
    P <- ParametersMorris %>% column_to_rownames("Parameter") # P is parameters to be tested in the morris test
    binf <- P$binf
    bsup <- P$bsup
    
    Morris.factors <- ParametersMorris$Parameter  
    N.par <- nrow(P)
    
    pars <- ModelParams[[PFAS]] # pars is all model parameters
    
    
    
    # SET MORRIS EXPERIMENT DESIGN ####
    # ---------------------------------------------------------------------------- #
    
    
    # Experiment design
    Morris.design <- list(type = "oat", levels = 6, grid.jump = 3)
    
    # Repetitions
    Morris.r <- 5
    
    # Morris test function
    Experience <- morris(model =  NULL, 
                         factors = Morris.factors, # number of parameters
                         r = Morris.r,  # number of repetitions
                         design = Morris.design,  
                         binf = binf,  
                         bsup = bsup, 
                         scale = TRUE)
    
    # Save Morris test design
    design <- Experience$X
    colnames(design) <- row.names(P)
    
    write.csv(design, here(OUTPUT, paste0(PFAS,"MorrisExpDesign.csv")), row.names = FALSE)
    Morris_PFAS_design[[i]] <- design
    
    
    pdf(here(OUTPUT, paste0(PFAS, "Distr_ParFull.pdf")))
    par(mfrow=c(4,3))
    for( i in 1:ncol(design)){ hist(design[,i], breaks=100, col="purple", main=colnames(design)[i] ) }
    dev.off()
    save(Experience, file = here(OUTPUT, paste0(PFAS, "ExperienceFull.RData")))
    
    # load(here(OUTPUT,paste0(PFAS,"ExperienceFull.RData")))
    
    
    # RUN MORRIS DESIGN ####
    # ---------------------------------------------------------------------------- #
    
    # Runs r (n of repetitions) * (param+1) simulations * (n) model outputs * (n) model outputs per time
    # results <- apply(design, 1,  SENSI_model) 
    results <- apply(design, 1, function(row) {
      SENSI_model(row, pars = pars)  
    })
    
    
    y <- results
    Morris_PFAS_y[[i]] <- results
    save(y, file = here(OUTPUT, paste0(PFAS,"y.RData")))
    
    y <- t(y) %>% na.omit() %>% t()
    sim.results.morris <- as.matrix(y)
    
    
    
    mu = mu.star = sigma = NULL
    
    # Tell morris y for each model output
    for (i in (1:length(sim.results.morris[,1])))
    {
      tell(Experience, sim.results.morris[i,])
      mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
      mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
      sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))
      
    }
    
    Experience$ee <- Experience$ee %>% na.omit() 
    print(PFAS)
    plot(Experience)
    print(Experience)

    rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(sim.results.morris)
    
    mu.star_scale = mu.star/apply(mu.star, 1, mean, na.rm=TRUE) 
    sigma_scale = sigma/apply(sigma, 1, mean, na.rm=TRUE)
    
    # Calculate global indices
    IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE)
    IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
    
    Morris_IGlobal.df <- data.frame(
      IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE),
      IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
    )
    
    IGlobal <- data.frame(IGlobal_sigma, IGlobal_mu.star) 
    IGlobal <- IGlobal %>% mutate(parm_names = rownames(IGlobal),
                                  PFAS = PFAS)
    Morris_PFAS_RESULTS[[i]] <- IGlobal
    
    Morris.indices <- list(
      mu = mu,
      mu_star = mu.star,
      sigma = sigma,
      mu_star_scale = mu.star_scale,
      sigma_scale = sigma_scale,
      IGlobal = IGlobal
    )
    print(IGlobal)
    
    save(Morris.indices, file = here(OUTPUT,paste0(PFAS, "Morris.indices.RData")))
    
  }

  for(i in 1:nPFAS){
    
    PFAS <- PFAS_list[i]
    load(here(OUTPUT,paste0(PFAS, "Morris.indices.RData")))
    
    Morris_PFAS_RESULTS[[i]] <- Morris.indices$IGlobal
    
  }


#   # Plotting ----

  Morris_RESULTS <- do.call(rbind, Morris_PFAS_RESULTS)
  # Morris_RESULTS <- cbind(rownames(Morris_RESULTS), 
  #                         data.frame(Morris_RESULTS, 
  #                                    row.names=NULL)) %>%
  #   separate("rownames(Morris_RESULTS)", 
  #            sep = "_", into = c("SF", "Parameter")) %>% 
    # select(!SF)
  
  library("viridis")
  
  morris_data <- Morris_RESULTS
  ALL <- morris_data %>% arrange(desc(IGlobal_mu.star)) %>%
    mutate(n = seq(from=1, to=153, by=1)) %>%
    filter(IGlobal_sigma > 0) %>% 
    separate(parm_names, sep = "SF_", into=c("SF", "parm_names")) %>%
    select(!c(Parameter, SF))

  Morris_plot <-  ALL %>%
    ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma,
               label = parm_names,
               color = parm_names)) +
    scale_color_viridis(discrete = TRUE, option = "B", end = 0.8) +
    geom_jitter(width = .5, alpha = 0.7, size = 1)  +
    geom_text_repel(aes(label = as.character(parm_names)),
                    size = 12,
                    max.overlaps = Inf
    ) +
    facet_wrap(~PFAS, scales = "free") +
    labs(
      x = expression(paste("Mean Effect (μ*)")),
      y = expression(paste("Standard Deviation (σ)"))
    )  +
    CP_theme +
    theme(legend.position = "none", 
          axis.text = element_text(size = 25),
          strip.text.x = element_text(size = 35)) 
  
  ggsave(plot = Morris_plot,
         file = here(OUTPUT, "Morris_plot.png"),
         dpi = 600,
         width = 18,
         height = 16,
         units = "cm")
