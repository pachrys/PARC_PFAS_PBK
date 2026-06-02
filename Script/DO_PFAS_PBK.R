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
  
  # Load input files
  INPUT_dummy <- read.csv(here("Input", "INPUT_dummy.csv"))
  Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
  Tissue.c <- read_csv(here("Input", "TissueComposition.csv"))
  PhysChemTK.c <- read_xlsx(here("Input", "TK_parameters.xlsx")) 
  
  # Load R files 
  source(here("Script", "RUN_and_OUTPUT.R"))
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  source(here("Script", "PBK_model.R"), local = TRUE)  
  
  
  # Manual input
  Lifestage = "Yes"   # Choose if physiology should change with age ("Yes" to include physiological changes due to age and "No" to assume the same physiology over time)

  PFAS_list <- unique(PhysChemTK.c$Acronym)
  nPFAS <- length(PFAS_list)
  
  # Exposure ----
  exposure_type = "Oral"
  exp_Oral = 0.0418 # ug/kg/day, (for Abraham: 3.96/BW of 82Kg),to be used only when both oral and dermal are used
  exp = exp_Oral # ug/kg/day 
  expSTOP = 5*365 #7300 #50*365 # time in days after which the exposure stopped
  
  # Subject-relevant information
  sex = "F" # sex either "F" or "M" if none then default is "M"
  expAGE = 20 #65  # years, old age at exposure if not provided then age argument is not used physiology is based on BW
  expBW = NA # 82 # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one

  # Simulation relevant information
  Tstart = 0 # days, start of the simulation
  Tstop = 20*365 # days, stop of the simulation
  Dt = 1  # days, iteration steps (decrease/increase depending on run time)
  
  PFAS.RESULTS <- list(
    CALC_Parameters = vector("list", nPFAS), # list of length nPFAS
    PBK_OUTPUT = vector("list", nPFAS) #,
  )
  
  names(PFAS.RESULTS$CALC_Parameters) <- PFAS_list
  names(PFAS.RESULTS$PBK_OUTPUT) <- PFAS_list
 
  for (i in 1:nPFAS) {
    
    PFAS <- PFAS_list[i]
    message(glue("Simulating {PFAS} {i} of {nPFAS} PFAS"))
    PhysChemTK_params <- PhysChemTK.c %>% filter(Acronym == PFAS) 
  
  MODEL_OUTPUT <- if(Lifestage == "Yes"){ # Individual and lifestage
    
    RUNandOUT_lifestage(
      
      # Exposure-relevant information
      exposure_type = exposure_type, 
      exp = exp, 
      exp_Oral = exp_Oral, 
      exp_Dermal = exp_Dermal, 
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
  
  PFAS.RESULTS$CALC_Parameters[[i]] <- MODEL_OUTPUT$CALC_Parameters
  PFAS.RESULTS$PBK_OUTPUT[[i]] <- MODEL_OUTPUT$PBK_OUTPUT
  
  }
 
  # Save Results ----

  RawData <- PFAS.RESULTS
  
  save(RawData, file = here(OUTPUT, "RawData_PBK.RData"))
  
  # Post-run Analysis ----

  POST.RUN_RESULTS <- list()
  
  for (i in 1:nPFAS) {
    
    PFAS <- PFAS_list[i]
    message(glue("\
                 {PFAS}"))
    Data <- RawData$PBK_OUTPUT[[i]]
    
    POST.RUN_RESULTS[[i]] <- POST.Run(exposure_type,  
                                      exp,
                                      exp_Oral,
                                      expSTOP,   
                                      expAGE, 
                                      expBW, 
                                      sex,  
                                      Tstart, 
                                      Tstop,  
                                      Dt, 
                                      Data)
    
    
  }
  names(POST.RUN_RESULTS) <- PFAS_list
  
  
  save(POST.RUN_RESULTS, file = here(OUTPUT, "POST.RUN_RESULTS.RData"))
 

  
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
