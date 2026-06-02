  # --------------------------------------------------------------------------- #
  # SCRIPT FOR RUNNING THE PBK MODEL & OBTAINING RESULTS
  # By: Chrysanthi Pachoulide
  # Date: 05-11-2025
  # --------------------------------------------------------------------------- #
  
  rm(list=ls()) 
  
  # Packages
  library(readxl)
  library(here)
  library(tidyverse)
  library(deSolve)
  library(PKNCA)
  library(pracma)
  library(glue)
  library(patchwork)
  library(quarto)
  library(tinytex)
  library(webshot2)
  library(knitr)
  library(showtext)
  font_add(family = "Garamond", regular = "GARA.TTF")
  showtext_auto()

  
  CP_theme <- theme_minimal() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "bottom",
          plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm") 
    )
  
  # Set output storage directory
  OUTPUT <- here("Output") #, format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
  dir.create(OUTPUT, recursive = TRUE)
  
  # Manual input
  Lifestage = "Yes"   # Choose if physiology should change with age ("Yes" to include physiological changes due to age and "No" to assume the same physiology over time)
  Population = "No"   # Choose to include population or individual exposure ("Yes" to include population based and "No" to only run the model for one person)
  Test_study = ""     # Choose study from input file
  PFAS = "PFNA"       # Choose PFAS to run the model on

  # Load input files
  INPUT_dummy <- read.csv(here("Input", "INPUT_dummy.csv"))
  Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
  Tissue.c <- read_csv(here("Input", "TissueComposition.csv"))
  PhysChemTK.c <- read_xlsx(here("Input", "TK_parameters.xlsx")) 
  
  # Load R files 
  source(here("Script", "RUN_and_OUTPUT.R"))
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  source(here("Script", "PBK_model.R"), local = TRUE)  
  


  # Load input ----
  # Set exposure type, choose between "Oral", "Dermal", "Oral_Dermal" (if exposure is both Oral and Dermal), Inhalation"
  exposure_type = "Oral"
  
  # Add input information
  
  # Exposure-relevant information
  # Current input is that of the Abraham study; Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047 
  exposure_type = exposure_type # type of exposure
  exp_Oral = 0.00019 #0.0418 # ug/kg/day, (for Abraham: 3.96/BW of 82Kg),to be used only when both oral and dermal are used
  exp_Dermal = 0 # ug/kg/day, to be used only when both oral and dermal are used
  exp = exp_Oral + exp_Dermal # ug/kg/day 
  Tinput = 1 # for repeated exposure or so default = 1
  tinterval = 1 # for repeated exposure or so default = 1
  expSTOP = 5*365 #7300 #50*365 # time in days after which the exposure stopped
  
  # Subject-relevant information
  expAGE = 20 #65  # years, old age at exposure if not provided then age argument is not used physiology is based on BW
  expBW = NA # 82 # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
  sex = "F" # sex either "F" or "M" if none then default is "M"
  
  # Simulation relevant information
  Tstart = 0 # days, start of the simulation
  Tstop = 20*365 # days, stop of the simulation
  Dt = 1 #1/10000 # days, iteration steps (decrease/increase depending on run time)
  
  
  RawData <- list(
    CALC_Parameters = vector("list", 1), 
    PBK_OUTPUT = vector("list", 1) #,
    # ANALYSED_data = vector("list", 1),
    # OUT_Plots = vector("list", 1) 
  )
  
  if(Lifestage == "Yes"){
    
    expBW <- NA # forcing BW to be NA as BW is taken based on age in the lifestage model
    if(is.na(expAGE)){
      stop("expAGE required to run the lifestage option")
    }
    
  }

 
  MODEL_OUTPUT <- if(Lifestage == "Yes"){ # Individual and lifestage
    
    RUNandOUT_lifestage(
      
      # Exposure-relevant information
      exposure_type = exposure_type, 
      exp = exp, 
      exp_Oral = exp_Oral, 
      exp_Dermal = exp_Dermal, 
      Tinput = Tinput, 
      tinterval = tinterval, 
      expSTOP = expSTOP, 
      
      # Subject-relevant information
      expAGE = expAGE, 
      expBW = expBW, 
      sex = sex, 
      
      # Simulation relevant information
      Tstart = Tstart, 
      Tstop = Tstop, 
      Dt = Dt 
      
    )
    
  } else { # Individual, no lifestage
    
    RUNandOUT(
      
      # Exposure-relevant information
      exposure_type = exposure_type, 
      exp = exp, 
      exp_Oral = exp_Oral, 
      exp_Dermal = exp_Dermal, 
      Tinput = Tinput, 
      tinterval = tinterval, 
      expSTOP = expSTOP, 
      
      # Subject-relevant information
      expAGE = expAGE, 
      expBW = expBW, 
      sex = sex, 
      
      # Simulation relevant information
      Tstart = Tstart, 
      Tstop = Tstop, 
      Dt = Dt 
      
    )
    
  }
 
  # Save Results ----

  RawData <- MODEL_OUTPUT
  
  save(RawData, file = here(OUTPUT, "RawData_PBK.RData"))
  
  # Post-run Analysis ----

  Pers.POST.RUN_RESULTS <- POST.Run(exposure_type = exposure_type,
                                 exp = exp, 
                                 exp_Oral = exp_Oral, 
                                 exp_Dermal = exp_Dermal, 
                                 Tinput = Tinput, 
                                 tinterval = tinterval, 
                                 expSTOP = expSTOP, 
                                 
                                 # Subject-relevant information
                                 expAGE = expAGE, 
                                 expBW = expBW, 
                                 sex = sex, 
                                 
                                 # Simulation relevant information
                                 Tstart = Tstart, 
                                 Tstop = Tstop, 
                                 Dt = Dt, 
                                 RawData = RawData)
    
  save(Pers.POST.RUN_RESULTS, file = here(OUTPUT, "Pers.POST.RUN_RESULTS.RData"))
 

  
  # Create Report ----
  # quarto_render(input = (here("PBK_Results_Report.qmd")),
  #               output_format = "pdf",
  #               # output_file = (here(OUTPUT, "PBK_Results_Report.html")),
  #               execute_params = list(Lifestage = Lifestage,
  #                                     Population = Population,
  #                                     Test_study = Test_study)
  # )
  # 
  # quarto_render(input = (here("PBK_Results_Report.qmd")),
  #               output_format = "html",
  #               # output_file = (here(OUTPUT, "PBK_Results_Report.html")),
  #               execute_params = list(Lifestage = Lifestage,
  #                                     Population = Population,
  #                                     Test_study = Test_study)
  # )
  # browseURL(here("PBK_Results_Report.html"))
  # 
  # 
