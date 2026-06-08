# --------------------------------------------------------------------------- #
# FUNCTIONS TO  CALCULATING MODEL RESULTS
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

# Without physiological changes with lifestage (depending on age) ####
RUNandOUT <- function(exposure_type, 
                      exp_Oral = NULL, exp_Dermal = NULL, exp, 
                       expSTOP, 
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {
  
  
  ## Add state ####
  Astate <- c(OD = 0,
              AIL = 0,
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
  
  
  # Calculate parameters
  base.parm.c = BASE_PARAMS(expAGE, expBW, sex) # common parameters
  parm.c <- c(base.parm.c, 
              expSTOP = expSTOP, 
              expOral = exp)
  # write.csv(parm.c, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
  
  
  
  ## Run the PBK model ####
  PBK_OUTPUT <- ORAL_PBK_RUN(y = Astate, 
                             parms = parm.c, 
                             times = seq(Tstart,Tstop,by=Dt))

  return(list(
    CALC_Parameters = parm.c,
    PBK_OUTPUT = PBK_OUTPUT
  ))
  
  
}


# With physiological changes with lifestage (depending on age) ####
RUNandOUT_lifestage <- function(exposure_type, 
                      exp_Oral = NULL, exp_Dermal = NULL, exp, 
                       expSTOP, 
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {

  ## Add state ####
  Astate <- c(OD = 0,
              AIL = 0,
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

  
  # Initial input for lifestage changes
  fullPBK_OUTPUT <- NULL
  newTstart <- Tstart
  newAstate <- Astate
  newAGE <- expAGE
  newBW <- expBW
  
  while (newTstart < Tstop) {
    
    yearlyTstop <- min(newTstart + 365, Tstop) 
    
    
    ## Calculare parameters, per newAGE ####
    
    base.parm.c = BASE_PARAMS(newAGE, newBW, sex) # common parameters
    parm.c <- c(base.parm.c, 
                expSTOP = expSTOP, 
                expOral = exp)
  
    # write.csv(parm.c, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
    
    
    
    ## Run the PBK model, per current age ####
    yearlyPBK_OUTPUT <- ORAL_PBK_RUN(y = newAstate, 
                                     parms = parm.c, 
                                     times = seq(newTstart, yearlyTstop, by=Dt))

    if(is.null(fullPBK_OUTPUT)){
      fullPBK_OUTPUT <- yearlyPBK_OUTPUT # this is for the first iteration
    } else {
      fullPBK_OUTPUT <- rbind(fullPBK_OUTPUT, yearlyPBK_OUTPUT[-1,]) # binding all iteration results 
    }
    
    
    # Update input for next iteration
    newAstate <- unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:17])   
    newTstart <- yearlyTstop
    newAGE <- expAGE + (newTstart/365)
    newBW <- NA
    
  }

  PBK_OUTPUT <- fullPBK_OUTPUT
  
  return(list(
    CALC_Parameters = parm.c,
    PBK_OUTPUT = PBK_OUTPUT # time is in days
  ))
  
  
}

# Post-run analysis and plots per PFAS ####

POST.Run <- function(exposure_type, 
                     exp_Oral = NULL, exp, 
                     expSTOP, 
                     expAGE = NULL, expBW = NULL, sex = "M", 
                     Tstart, Tstop, Dt,
                     Data){
  
  ## Input ####
  PBK_OUTPUT <- Data #$PBK_OUTPUT

  ## Plot Mass Balance/Error ####
  MB.df <- PBK_OUTPUT %>% select(time, Ain, Atot, MB)  # change days to years if needed
  MB.df$MB <- round(MB.df$MB, 10) 
  MB.df$ERROR <- (MB.df$Ain - MB.df$Atot) / MB.df$Atot * 100
  MB.df$ERROR <- round(MB.df$ERROR, 10) 
  MB_plot <- ggplot(data = MB.df)+
    geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
    geom_line(aes(x = time, y = MB, color = "MB")) +
    scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"), 
                       name = NULL) +
    CP_theme +
    ylab("MB / ERROR")
  MB_plot
  
  ## Plot Tissue concentrations ####
  
  # For liver and kidney
  BW <- RawData$CALC_Parameters$BW
  VKc <- RawData$CALC_Parameters$VKc
  VK <- VKc * BW                    # L, Volume of kidney
  VPTc <- RawData$CALC_Parameters$VPTc
  VPT <- VPTc * VK                  # L, Volume of proximal tubule
  VPTTc <- RawData$CALC_Parameters$VPTTc 
  VPTT <- VPT * VPTTc              # L, Volume of proximal tubule Tissue
  VPTLc <- RawData$CALC_Parameters$VPTLc  
  VPTL <- VPT * VPTLc
  VRKLc <- RawData$CALC_Parameters$VRKLc  
  VRKL <- (VK - VPT) *  VRKLc    # L, Volume of rest of kidney lumen
  cVRKL <- VK/VRKL
  VRKT <- VK - VPTT - VPTL - VRKL # L, Volume of rest of kidney Tissue, used in the model
  VRKTc <- VK/VRKT
  
  VLc <- RawData$CALC_Parameters$VLc
  VL <- VLc * BW
  VL_icc <- RawData$CALC_Parameters$VL_icc 
  VL_ic <- VL_icc * VL              # L, Volume liver intracellular                 
  VL_ecc <- RawData$CALC_Parameters$VL_ecc
  VL_ec <- VL_ecc * VL              # L, Volume liver extracellular 
  
  C_Tissues.df <- PBK_OUTPUT %>% 
    transmute(
      time = time,
      CI = CI,
      CK = CPTT, 
      CL = CL_ic,
      CA, CR, CP  
      ) %>% 
    pivot_longer(cols = CI:CP, 
                 names_to = "Tissue", 
                 values_to = "Concentration") %>% 
    mutate(time = time/365) #transforming time in years
  
  # Plot Tissue concentrations
  Plot_C_Tissues <- C_Tissues.df %>% 
    ggplot(aes(time, Concentration)) +
    geom_path(linewidth = 0.5) +
    facet_wrap(~Tissue) +
    labs(title = "Tissue concentrations",
         x = "Time (years)", # check that time is indeed in days and not years
         y = "Concentration (ng/ml)") +
    CP_theme 
  Plot_C_Tissues
  
  ## Calculate AUC and Half life ####
  AUC <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CP"])  # ug*day/L
  AUC <- round(AUC, 3)
  message(glue("AUC {AUC} ug*day/L"))
  
  
  time <- PBK_OUTPUT[ , "time"] # days
  conc <- PBK_OUTPUT[ , "CP"] # ug/L or ng/ml
  Cmax <- max(conc)
  Tmax <- time[which.max(conc)]
  tlast <- max(time[conc > 0])
  half_life <- pk.calc.half.life(
    conc,
    time,
    Tmax,
    tlast
  )
  
  Cmax <- round(Cmax, 3)
  message(glue("Cmax {Cmax} ug/L"))
  
  Tmax <- round(Tmax/365, 3)
  message(glue("Tmax {Tmax} years"))
  
  HalfLife <- half_life$half.life/365  # half-life in years
  HalfLife <- round(HalfLife,3)
  message(glue("Half-life {HalfLife} years"))
  
  return(list(
    ANALYSED_data = data.frame(exposure_type = paste(exposure_type),
                               exp = paste(round(exp, digits = 10),"ug/kg/day", sep = "_"),
                               expAGE = paste(expAGE, "years", sep = "_"),
                               Tstart = paste(round(Tstart/365, digits = 10), "years", sep = "_"),
                               expSTOP = paste(round(expSTOP/365, digits = 10), "years", sep = "_"),
                               Tstop = paste(round(Tstop/365, digits = 10), "years", sep = "_"),
                               expBW = paste(round(expBW, digits = 10), "kg", sep = "_"),
                               sex = paste(sex),
                               AUC = paste(round(AUC, digits = 3), "ug*day/L", sep = "_"),
                               Cmax = paste(round(Cmax, digits = 3), "ug/L", sep = "_"),
                               Tmax = paste(round(Tmax, digits = 3), "years", sep = "_"),
                               HalfLife = paste(round(HalfLife, digits = 3), "years", sep = "_")),
    OUT_Plots = list(MB_plot, Plot_C_Tissues) 
  ))
}




