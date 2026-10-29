# --------------------------------------------------------------------------- #
# Script for optimising the most sensitive parameters
# Chrysanthi Pachoulide, August 2026
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment


# Load packages
library(here)
library(ggplot2)
library(deSolve)
library(tidyverse)
library(glue)
library(readxl)
library(PKNCA)
library(parallel)
library(pracma)
library(viridis)
library(lhs)
library(data.table)
library(patchwork)

# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()


# Set output storage directory
# INPUT <- here("Output", "GoF_upTo3days") 
OUTPUT <- here("Output", "test")
dir.create(OUTPUT, recursive = TRUE)


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
    
    strip.text.x = element_text(size = 35)
  )


Patchwork_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 50, face = "bold", margin = margin(0, 0, b = 0.5, 0)),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    plot.title.position = "plot",
    panel.background = element_blank(), #rect(fill = 'white', color = NA)
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_line(color = "grey90",linewidth = 0.1),
    
    axis.title = element_text(size = 45, face = "bold"),
    axis.text.x = element_text(size = 35, face = "bold", margin = margin(0, 0, 0, 0)),
    axis.text.y = element_text(size = 35, face = "bold", margin = margin(0, 0, 0, 0)),
    axis.line = element_blank(), 
    axis.ticks = element_blank(),
    
    legend.position = "bottom",
    legend.title = element_blank(), 
    legend.text = element_text(size = 35),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.8, "lines"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
    
    # strip.text.x = element_blank()
    strip.text = element_text(hjust = 0, size = 45, face = "bold", margin = margin(0,0,0,0))
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


# Load data
Abraham_df <- read_excel(here("Input","Evaluation_data.xlsx")) %>% 
  filter(Reference == "Abraham_2024") |> 
  select(!c(PFAS_name, CAS_n)) %>% 
  rename(PFAS = Acronym) |> 
  mutate(Dose = round(as.numeric(Dose_ugPerkg_bw)*82, 3)) |>  #82 is the body weight of Abraham
  select(PFAS, Dose) |> 
  mutate(PFAS = str_replace_all(PFAS, "_", "-")) |> 
  distinct()

PFAS_list <- c(
  "PFTrDA",
  "PFDoDA",
  "PFUnDA",
  "PFDA",
  "PFOS",
  "PFNA",
  "PFHpS",
  "PFOA",
  "PFHxS",
  "ADONA",
  "PFHpA",
  "PFPeS",
  "HFPO-DA",
  "PFHxA",
  "PFBS",
  "PFPeA",
  "PFBA"
)

Dose_df <- tibble(
  PFAS = PFAS_list,
  Dose = 4) |> 
  mutate(Dose = ifelse(PFAS == "PFBS", 18.4,
                       ifelse(PFAS == "HFPO-DA" , 19.6,
                              Dose)))
TIMES <- seq(0, 20, by = 1) #1 day

nPFAS <- length(PFAS_list)


## Tissue constituents for partition coefficient calculation
Physio.data <- read_csv(here("Input", "TissueComposition.csv")) %>% 
  filter(Species == "Human") %>%
  filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
  select(- c(Comment, Species)) %>% 
  mutate(across(-Tissue, as.numeric))


## Physiological parameters
Physio_params <- read_xlsx(here("Input", "Physio_Parameters.xlsx")) 

parm.phys <- Physio_params$Value
names(parm.phys) <- Physio_params$Abbreviation

## MC parameters
Tune_TK_params <- read_xlsx(here("Input", "TK_Parameters_last.xlsx")) |> 
  filter(Optimize == "yes") 
Fix_TK_params <- read_xlsx(here("Input", "TK_Parameters.xlsx")) |> 
  filter(Optimize != "yes") 

Input_PhysicoChemical_Parameters <- c("logD_MembraneLipidWater",
                                      "logD_StructuralProteinWater",
                                      "logD_AlbuminWater",
                                      "logD_StorageLipidWater",
                                      "logD_GlobulinWater",
                                      "Papp_SI",
                                      "Clear_OATP2B1c",
                                      "Clear_OATP1B1c",
                                      "Clear_OATP1B3c",
                                      "Clear_NTCPc",
                                      "Clear_BSEPc",
                                      "Clear_OAT3c",
                                      "Clear_OAT1c",
                                      "Clear_OAT4c",
                                      "SF_OATP2B1_I",
                                      "SF_OATP1B1",
                                      "SF_OATP1B3",
                                      "SF_OATP2B1",
                                      "SF_NTCP",
                                      "SF_BSEP",
                                      "SF_OAT3",
                                      "SF_OAT1",
                                      "SF_OAT4")

Intermediate_TK_Parameters <- c("fup", 
                                "fu_PTL",     
                                "fu_L_ec",     
                                "fu_PT_ec",
                                "Clear_apparent",       
                                "Clear_OATP2B1_I",          
                                "Clear_OATP1B1",              
                                "Clear_OATP1B3",                
                                "Clear_OATP2B1_L",                
                                "Clear_NTCP",                  
                                "Clear_BSEP",            
                                "Clear_OAT3",              
                                "Clear_OAT1",                      
                                "Clear_OAT4")

Final_TK_Parameters <- c("PI",
                         "PL",
                         "PK",
                         "PA",
                         "PR",
                         "Clear_uptake_I",             
                         "Clear_uptake_L",
                         "Clear_excretion_L",
                         "GFR_free",
                         "Clear_reabsorption_K")

# load functions
source(here("Script", "functions_NG_PBK_model.R"))

# Run UA ----
# Create Monte Carlo function
run_Monte_Carlo <- function(parm.chem, parm.phys, Dose_ug_day){
  
  A_init <- c(AIL = Dose_ug_day,
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
              InMB = Dose_ug_day)
  
  
  pars <- Calc_final_params(parm.chem = parm.chem, 
                            parm.phys = parm.fix)
  
  # save parameters to check when analysing results
  run_pars <- as.data.frame(pars)
  
  PBK.out <- lsoda(y = A_init,
                   times = TIMES, 
                   func = PBK.model,
                   parms = pars,
                   rtol = 1e-4,        
                   atol = 1e-6,
                   hmax = 0.5,         
                   maxsteps = 100000)
  
  PBK_OUTPUT <- as.data.frame(PBK.out) 
  PBK_OUTPUT[,-1] <- pmax(PBK_OUTPUT[,-1], 0)
  
  PBK_OUTPUT <- cbind(PBK_OUTPUT, run_pars)
  
  return(PBK_OUTPUT)
  
}

Monte_Carlo_PFAS_savey <- list()

for (i in 1:nPFAS) {
  
  
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}-n{i}-----------------------------------"))
  
  # physiological parameters; they are kept fix
  Dose_ug_day <- filter(Dose_df, PFAS == PFAS_i) |> 
    select(Dose) |> pull()
  
  parm.phys <- parm.phys
  
  names_parm.chem.fix <- Fix_TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(Parameter) |> pull()
    
  parm.chem.fix <- Fix_TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(Value_fix) |> 
    pull()
  
  names(parm.chem.fix) <- names_parm.chem.fix
  
  parm.fix <- c(parm.phys, parm.chem.fix)
  
  parm.chem <- Tune_TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(!PFAS) |> 
    mutate(Value = ifelse(PBK_mechanism == "Transport", log(Value), Value)) |> 
    mutate(meanlog = ifelse(Distribution == "lognormal", Value, NA),
           sdlog = ifelse(Distribution == "lognormal", RSE, NA),
           min = ifelse(PBK_mechanism == "Scaling", Value-RSE, NA), 
           max = ifelse(PBK_mechanism == "Scaling", Value+RSE,NA),
           Distribution = ifelse(Distribution == "lognormal", "qlnorm", "qunif"))
  
  MC.parms <- parm.chem$Parameter
  MC.factors <- length(MC.parms) 
  
  ## Latin Hypercube Sampling ---- 
  # this optimises the Monte Carlo sampling, also used by McNally & Loizou
  
  set.seed(10010) # set seeds for reproducibility
  n = 10#00 # number of samples per parameter
  X.LHS.MC.design <- randomLHS(n, MC.factors, preserveDraw = TRUE) # create lhs grid
  Y.LHS.MC.design <- matrix(0, nrow = n, ncol = MC.factors) # empty matrix to populate the parameters with
  
  colnames(X.LHS.MC.design) = colnames(Y.LHS.MC.design) = MC.parms
  
  for (i in 1:MC.factors) {
    
    parm <-   MC.parms[i]
    parm.data <- filter(parm.chem, Parameter == parm) 
    pdf <- parm.data |> pull(Distribution)
    
    Y.LHS.MC.design[,i] <- if(pdf == "qlnorm"){
      
      qlnorm(X.LHS.MC.design[,i], 
             meanlog = parm.data$meanlog, 
             sdlog = parm.data$sdlog)
      
    } else{
      qunif(X.LHS.MC.design[,i], 
            min = parm.data$min, 
            max = parm.data$max)
    }
    
  }
  
  # create the final MC design and remove infinite values if any
  MC.design <- Y.LHS.MC.design |> 
    as.data.frame() |> 
    mutate(across(everything(), ~na_if(., Inf))) |> 
    drop_na()
  
  # visualise and safe design
  Design_plot <- MC.design |>  
    pivot_longer(everything(), 
                 names_to = "Parameter", 
                 values_to = "Value") |> 
    ggplot(aes(x=Value))+
    geom_histogram()+
    facet_wrap(~ factor(Parameter, Input_PhysicoChemical_Parameters), 
               scales = "free") +    coord_flip()+
    labs(title = paste0(PFAS_i," Monte Carlo parameter distributions"),
         y = "Frequency",
         x = "Value (different unit per parameter)") +
    theme_quick 
  Design_plot
  
  Design_plot <- Design_plot + Patchwork_theme +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(plot = Design_plot,
         file = here(OUTPUT, paste0(PFAS_i,"_LHS_MC_ExpDesign.png")),
         dpi = 600,
         width = 24,
         height = 16,
         units = "cm")
  
  write.csv(MC.design, 
            file = here(OUTPUT, paste0(PFAS_i,"_LHS_MC_ExpDesign.csv")), 
            row.names = FALSE)
  
  results <- apply(MC.design, 1, function(row) {
    run_Monte_Carlo(row,
                    parm.phys = parm.phys,
                    Dose_ug_day = Dose_ug_day)
  })
  
  y <- results
  savey <- mapply(cbind, y,
                  "Iteration" = seq(1,length(y), by = 1),
                  "PFAS" = PFAS_i,
                  SIMPLIFY = FALSE)
  save(savey, file = here(OUTPUT, paste0(PFAS_i,"_LHS_MC_y.RData")))
  
  Monte_Carlo_PFAS_savey[[i]] <- savey
}

# names(Monte_Carlo_PFAS_savey) <- PFAS_list
# 
# Save_MC_PBK_results <- rbindlist(Monte_Carlo_PFAS_savey, idcol = "PFAS")


# Produce results ----


## functions ----

source(here("Script", "functions_calculate_results.R"))


## produce result files and plots for all PFAS ----
INPUT <- OUTPUT #here("Output", "Second_try_UA")
# OUTPUT <- here("Output","Second_try_UA", "Processed_results") 
# dir.create(OUTPUT, recursive = TRUE)
# 
# PFAS_list2 <- c(
#   "PFHxS",
#   "ADONA",
#   "PFHpA",
#   "PFPeS",
#   "HFPO-DA",
#   "PFHxA",
#   "PFBS",
#   "PFPeA"  
# )
nPFAS <- length(PFAS_list)
for (i in 1:nPFAS) {
  
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}-n{i}-----------------------------------"))
  
  # Inputs ----
  
  # load results (called savey)
  load(here(INPUT, paste0(PFAS_i, "_LHS_MC_y.RData" )
  ))
  # load oral dose (different per PFAS if comparing results to Abraham)
  OD <- filter(Dose_df, PFAS == PFAS_i) |>
    select(Dose) |> pull()
  
  # TK results ----
  # calculate AUC, half-life, elimination rates, fraction absorbed etc

  # To see if there are any errors

  TK_results <- sapply(savey, check_error, OD)

  # if(is.list(TK_results)){
  #   errored_iterations <- which(sapply(TK_results, function(x) length(x) < 3)) #
  #   message(glue("Removing iterations n{errored_iterations} that give errors. \\"))
  #   write.csv(errored_iterations, here(OUTPUT, paste0(PFAS_i, "_errored_iterations.csv")), row.names = FALSE)
  #   savey <- savey[-errored_iterations]
  #   TK_results <- sapply(savey, calc_TK_results, OD)
  # } else if(!is_empty(which(is.na(TK_results[2,]) | TK_results[2,]<=2))){
  #   errored_iterations <- which(is.na(TK_results[2,]) | TK_results[2,]<=2)
  #   write.csv(errored_iterations, here(OUTPUT, paste0(PFAS_i, "_errored_iterations.csv")), row.names = FALSE)
  #   message(glue("Removing iterations n{errored_iterations} that give errors. \\"))
  #   savey <- savey[-errored_iterations]
  #   TK_results <- sapply(savey, calc_TK_results, OD)
  # } else {
    TK_results <- sapply(savey, calc_TK_results, OD)
  # }


  # to clean TK_results
  clean_TK_results <- TK_results |>
    t() |>
    as.data.frame() |>
    filter(CLfe < 1e5,
           CLur < 1e5) |>
    t()

  # find which iterations went wrong
  errored_TK_results <- TK_results |>
    t() |>
    as.data.frame() |>
    rownames_to_column() |>
    rename(Iteration = rowname) |>
    filter(CLfe > 1e5,
           CLur > 1e5)

  errored_TK_result_iterations <- errored_TK_results$Iteration

  if(is.numeric(errored_TK_results$Iteration)){
    message(glue("Removing iterations n{errored_TK_result_iterations} that give unrealistic clearance values (CL > 1e5)"))
    write.csv(errored_TK_results, file = here(OUTPUT, paste0(PFAS_i, "errored_TK_results.csv")))
  }

  TK_summary_results <- calc_TK_summary_results(clean_TK_results) |>
    as.data.frame() |> mutate(PFAS = PFAS_i)
  write.csv(TK_summary_results, file = here(OUTPUT, paste0(PFAS_i, "_TK_sumstat.csv")))


  Results_df <- do.call(rbind, savey) 

  # Actual parameters ----
  
  # Parameters_df <- Results_df |>
  #   select(VI_L:Clear_OAT4) |>
  #   distinct()
  # 
  # write.csv(Parameters_df, file = here(OUTPUT, paste0(PFAS_i, "_final_model_parameters.csv")))
  # 
  # 
  # Input_PhysicoChemical_Parameters_plot <- Parameters_df |> 
  #   select(logD_MembraneLipidWater:SF_OAT4) |> 
  #   pivot_longer(everything(),
  #                names_to = "Parameter",
  #                values_to = "Value") |> 
  #   ggplot() +
  #   geom_histogram(aes(x = Value)) +
  #   facet_wrap(~ factor(Parameter,Input_PhysicoChemical_Parameters), 
  #              scales = "free") +
  #   labs(title = paste0(PFAS_i, " Input physico-chemical parameters"),
  #        x = "Value",
  #        y = "Frequency") +
  #   coord_flip()+
  #   theme_quick
  # 
  # Intermediate_TK_Parameters_plot <- Parameters_df |> 
  #   select(fup:fu_PT_ec,Clear_apparent:Clear_OAT4) |> 
  #   pivot_longer(everything(),
  #                names_to = "Parameter",
  #                values_to = "Value") |> 
  #   ggplot() +
  #   geom_histogram(aes(x = Value)) +
  #   facet_wrap(~ factor(Parameter, Intermediate_TK_Parameters), scales = "free") +
  #   labs(title = paste0(PFAS_i, " intermediate TK parameters"),
  #        x = "Value",
  #        y = "Frequency") +
  #   coord_flip()+
  #   theme_quick
  # 
  # Final_TK_Parameters_plot <- Parameters_df |> 
  #   select(PI:Clear_reabsorption_K) |> 
  #   pivot_longer(everything(),
  #                names_to = "Parameter",
  #                values_to = "Value") |> 
  #   ggplot() +
  #   geom_histogram(aes(x = Value)) +
  #   facet_wrap(~ factor(Parameter, Final_TK_Parameters), 
  #              scales = "free") +
  #   labs(title = paste0(PFAS_i, " final TK parameters"),
  #        x = "Value",
  #        y = "Frequency") +
  #   coord_flip()+
  #   theme_quick
  # 
  # Final_Physio_Parameters_plot <- Parameters_df |> 
  #   select(VI_L:Calb_P, BW:QC) |> 
  #   pivot_longer(everything(),
  #                names_to = "Parameter",
  #                values_to = "Value") |> 
  #   ggplot() +
  #   geom_histogram(aes(x = Value)) +
  #   facet_wrap(~ Parameter, scales = "free") +
  #   labs(title = paste0(PFAS_i, " physiological parameters"),
  #        x = "Value",
  #        y = "Frequency") +
  #   coord_flip()+
  #   theme_quick
  # 
  # 
  # Input_PhysicoChemical_Parameters_plot = Input_PhysicoChemical_Parameters_plot +
  #   Patchwork_theme + theme(axis.text.x = element_text(angle = 0))
  # 
  # Intermediate_TK_Parameters_plot  = Intermediate_TK_Parameters_plot +
  #   Patchwork_theme + theme(axis.text.x = element_text(angle = 0))
  # 
  # Final_TK_Parameters_plot = Final_TK_Parameters_plot +
  #   Patchwork_theme + theme(axis.text.x = element_text(angle = 0))
  # 
  # Final_Physio_Parameters_plot = Final_Physio_Parameters_plot +
  #   Patchwork_theme + theme(axis.text.x = element_text(angle = 0))
  # 
  # ggsave(plot = Input_PhysicoChemical_Parameters_plot,
  #        file = here(OUTPUT, paste0(PFAS_i, "_Input_PhysicoChemical_Parameters_plot.png")),
  #        dpi = 600,
  #        width = 24,
  #        height = 18,
  #        units = "cm")
  # 
  # ggsave(plot = Intermediate_TK_Parameters_plot,
  #        file = here(OUTPUT, paste0(PFAS_i, "_Intermediate_TK_Parameters_plot.png")),
  #        dpi = 600,
  #        width = 18,
  #        height = 18,
  #        units = "cm")
  # 
  # ggsave(plot = Final_TK_Parameters_plot,
  #        file = here(OUTPUT, paste0(PFAS_i, "_Final_TK_Parameters_plot.png")),
  #        dpi = 600,
  #        width = 16,
  #        height = 16,
  #        units = "cm")
  # 
  # ggsave(plot = Final_Physio_Parameters_plot,
  #        file = here(OUTPUT, paste0(PFAS_i, "_Final_Physio_Parameters_plot.png")),
  #        dpi = 600,
  #        width = 24,
  #        height = 18,
  #        units = "cm")
  
  # Predicted tissue concentrations ----
  
  VL <- parm.phys[["VLc"]]*parm.phys[["BW"]]
  VL_ic <- parm.phys[["VL_icc"]]*VL
  VL_ec <- parm.phys[["VL_ecc"]]*VL
  
  VK <- parm.phys[["VKc"]]*parm.phys[["BW"]]
  VPT <- parm.phys[["VPTc"]]*VK
  VPTT <- VPT*parm.phys[["VPTTc"]]
  VPTL <- VPT*parm.phys[["VPTLc"]]
  VRKL <- (VK - VPT)*parm.phys[["VRKLc"]]
  VRKT <- VK - VPTT - VPTL - VRKL
  
  Initial_Tissues_list <- c("CI", "CL_ec", "CL_ic", "CPTL", "CPTT", "CRKL", "CRKT", "CA", "CR", "CP")
  
  
  C_all_df <- Results_df |>
    select(time, CI, CL_ec, CL_ic, CPTL, CPTT, CRKL, CRKT, CA, CR, CP, Iteration) |>
    mutate(PFAS = PFAS_i) |>
    mutate(CL = (CL_ec*VL_ec + CL_ic*VL_ic)/VL,
           CK = (CPTL*VPTL + CPTT*VPTT + CRKL*VRKL + CRKT*VRKT)/VK,
           time = time) |> # time to years
    select(PFAS, Iteration, time, CA, CR, CI, CL, CK, CP) |>
    pivot_longer(CA:CP,
                 names_to = "Tissue",
                 values_to = "Concentration")
  
  Tissues_list <- c("CA", "CR", "CI", "CL", "CK", "CP")
  
  Tissues_names <- c("Adipose", "Rest","Intestine", "Liver", "Kidney", "Plasma")
  
  # plasma concentration
  nTissues <- length(Tissues_list)
  
  C_Tissues_list <- NULL
  
  for (j in 1:nTissues) {
    
    Tissue_j <- Tissues_list[j]
    Tissue_j_name <- Tissues_names[j]
    
    C_summary <- calc_C_summary(C_all_df, Tissue_j) |>
      mutate(Tissue = Tissue_j_name,
             PFAS = PFAS_i)
    
    C_Tissues_list[[j]] <- C_summary
    
  }
  
  C_Tissues <- do.call(rbind, C_Tissues_list) |> 
    na.omit() |> 
    mutate(Time = rep(TIMES, nTissues),
           # Time = rep(TIMES, nTissues),
           Time = Time/365)
  
  C_Tissues_plot <- plot_all_tissues(C_Tissues)
  
  C_Tissues_plot <- C_Tissues_plot + Patchwork_theme +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(plot = C_Tissues_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_TissueC_plot.png")),
         dpi = 600,
         width = 16,
         height = 11,
         units = "cm")
  
  write.csv(C_Tissues, file = here(OUTPUT, paste0(PFAS_i, "_TissueC_sumstat.csv")))
  
  
}

# run summary results ----
## load input files ----
nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

TK_df <-
  list.files(path = here(OUTPUT),
             pattern = "\\_TK_sumstat.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .)))|> select(!"...1")
write.csv(TK_df, here(OUTPUT, "TK_results_all_PFAS.csv"), row.names = FALSE)

C_df <-
  list.files(path = here(OUTPUT),
             pattern = "\\_TissueC_sumstat.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .)))

Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |>
  rename(PFAS = Acronym) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))


## plot results ----
# ### Rest ----
# C_Rest_predicted_df <- C_df |>
#   filter(Tissue == "Rest") 
# 
# C_Rest_predicted_df <- full_join(C_Rest_predicted_df, nC_MW_df) 
# C_Rest_predicted_df <- full_join(C_Rest_predicted_df, Dose_df) 
# 
# CR_plot <- ggplot() +
#   geom_ribbon(data = C_Rest_predicted_df,
#               aes(x = Time,
#                   ymin = q025/Dose,
#                   ymax = q975/Dose,
#                   fill = Category_OECD),
#               alpha = 0.4) +
#   geom_ribbon(data = C_Rest_predicted_df,
#               aes(x = Time,
#                   ymin = q25/Dose,
#                   ymax = q75/Dose,
#                   fill = Category_OECD),
#               alpha = 0.3) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   geom_path(data = C_Rest_predicted_df,
#             aes(Time, median/Dose,
#                 colour = Category_OECD)) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA")
#   ),
#   nrow = 2) +
#   labs(x = "Time (years)",
#        y = "Concentration / dose",
#        title = "Rest") +
#   theme_quick
# CR_plot
# # CR_plot <- CR_plot + Patchwork_theme +
# #   theme(legend.position = "bottom")
# # 
# # ggsave(plot = CR_plot,
# #        file = here(OUTPUT, "CR_plot_UA.png"),
# #        dpi = 600,
# #        width = 24,
# #        height = 11,
# #        units = "cm")
# 
# 
# 
# ### Adipose ----
# C_Adipose_predicted_df <- C_df |>
#   filter(Tissue == "Adipose") 
# 
# C_Adipose_predicted_df <- full_join(C_Adipose_predicted_df, nC_MW_df) 
# C_Adipose_predicted_df <- full_join(C_Adipose_predicted_df, Dose_df) 
# 
# CA_plot <- ggplot() +
#   geom_ribbon(data = C_Adipose_predicted_df,
#               aes(x = Time,
#                   ymin = q025/Dose,
#                   ymax = q975/Dose,
#                   fill = Category_OECD),
#               alpha = 0.4) +
#   geom_ribbon(data = C_Adipose_predicted_df,
#               aes(x = Time,
#                   ymin = q25/Dose,
#                   ymax = q75/Dose,
#                   fill = Category_OECD),
#               alpha = 0.3) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   geom_path(data = C_Adipose_predicted_df,
#             aes(Time, median/Dose,
#                 colour = Category_OECD)) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA")
#   ),
#   nrow = 2) +
#   labs(x = "Time (years)",
#        y = "Concentration / dose",
#        title = "Adipose") +
#   theme_quick
# CA_plot
# # CA_plot <- CA_plot + Patchwork_theme +
# #   theme(legend.position = "bottom")
# # 
# # ggsave(plot = CA_plot,
# #        file = here(OUTPUT, "CA_plot_UA.png"),
# #        dpi = 600,
# #        width = 24,
# #        height = 11,
# #        units = "cm")
# 
# 
# ### Kidney ----
# C_Kidney_predicted_df <- C_df |>
#   filter(Tissue == "Kidney") 
# 
# C_Kidney_predicted_df <- full_join(C_Kidney_predicted_df, nC_MW_df) 
# C_Kidney_predicted_df <- full_join(C_Kidney_predicted_df, Dose_df) 
# 
# CK_plot <- ggplot() +
#   geom_ribbon(data = C_Kidney_predicted_df,
#               aes(x = Time,
#                   ymin = q025/Dose,
#                   ymax = q975/Dose,
#                   fill = Category_OECD),
#               alpha = 0.4) +
#   geom_ribbon(data = C_Kidney_predicted_df,
#               aes(x = Time,
#                   ymin = q25/Dose,
#                   ymax = q75/Dose,
#                   fill = Category_OECD),
#               alpha = 0.3) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   geom_path(data = C_Kidney_predicted_df,
#             aes(Time, median/Dose,
#                 colour = Category_OECD)) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA")
#   ),
#   nrow = 2) +
#   labs(x = "Time (years)",
#        y = "Concentration / dose",
#        title = "Kidney") +
#   theme_quick
# CK_plot
# # CK_plot <- CK_plot + Patchwork_theme +
# #   theme(legend.position = "bottom")
# # 
# # ggsave(plot = CK_plot,
# #        file = here(OUTPUT, "CK_plot_UA.png"),
# #        dpi = 600,
# #        width = 24,
# #        height = 11,
# #        units = "cm")
# 
# ### Liver ----
# C_liver_predicted_df <- C_df |>
#   filter(Tissue == "Liver") 
# 
# C_liver_predicted_df <- full_join(C_liver_predicted_df, nC_MW_df) 
# C_liver_predicted_df <- full_join(C_liver_predicted_df, Dose_df) 
# 
# CL_plot <- ggplot() +
#   geom_ribbon(data = C_liver_predicted_df,
#               aes(x = Time,
#                   ymin = q025/Dose,
#                   ymax = q975/Dose,
#                   fill = Category_OECD),
#               alpha = 0.4) +
#   geom_ribbon(data = C_liver_predicted_df,
#               aes(x = Time,
#                   ymin = q25/Dose,
#                   ymax = q75/Dose,
#                   fill = Category_OECD),
#               alpha = 0.3) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   geom_path(data = C_liver_predicted_df,
#             aes(Time, median/Dose,
#                 colour = Category_OECD)) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA")
#   ),
#   nrow = 2) +
#   labs(x = "Time (years)",
#        y = "Concentration / dose",
#        title = "Liver") +
#   theme_quick
# CL_plot
# # CL_plot <- CL_plot + Patchwork_theme +
# #   theme(legend.position = "bottom")
# # 
# # ggsave(plot = CL_plot,
# #        file = here(OUTPUT, "CL_plot_UA.png"),
# #        dpi = 600,
# #        width = 24,
# #        height = 11,
# #        units = "cm")
# 
# 
# 
# 
# 
# ### Intestine ----
# C_Intestine_predicted_df <- C_df |>
#   filter(Tissue == "Intestine") 
# 
# C_Intestine_predicted_df <- full_join(C_Intestine_predicted_df, nC_MW_df) 
# C_Intestine_predicted_df <- full_join(C_Intestine_predicted_df, Dose_df) 
# 
# CI_plot <- ggplot() +
#   geom_ribbon(data = C_Intestine_predicted_df,
#               aes(x = Time,
#                   ymin = q025/Dose,
#                   ymax = q975/Dose,
#                   fill = Category_OECD),
#               alpha = 0.4) +
#   geom_ribbon(data = C_Intestine_predicted_df,
#               aes(x = Time,
#                   ymin = q25/Dose,
#                   ymax = q75/Dose,
#                   fill = Category_OECD),
#               alpha = 0.3) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   geom_path(data = C_Intestine_predicted_df,
#             aes(Time, median/Dose,
#                 colour = Category_OECD)) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA")
#   ),
#   nrow = 2) +
#   labs(x = "Time (years)",
#        y = "Concentration / dose",
#        title = "Intestine") +
#   theme_quick
# CI_plot
# # CL_plot <- CL_plot + Patchwork_theme +
# #   theme(legend.position = "bottom")
# # 
# # ggsave(plot = CL_plot,
# #        file = here(OUTPUT, "CL_plot_UA.png"),
# #        dpi = 600,
# #        width = 24,
# #        height = 11,
# #        units = "cm")
# 
# 
# 
# 
# 

### Plasma ----
C_plasma_predicted_df <- C_df |>
  filter(Tissue == "Plasma") 

C_plasma_predicted_df <- full_join(C_plasma_predicted_df, nC_MW_df) 
C_plasma_predicted_df <- full_join(C_plasma_predicted_df, Dose_df) 

# CP_plot <- ggplot() +
#   geom_ribbon(data = C_plasma_predicted_df,
#               aes(x = Time,
#                   ymin = q025/Dose,
#                   ymax = q975/Dose,
#                   fill = Category_OECD),
#               alpha = 0.4) +
#   geom_ribbon(data = C_plasma_predicted_df,
#               aes(x = Time,
#                   ymin = q25/Dose,
#                   ymax = q75/Dose,
#                   fill = Category_OECD),
#               alpha = 0.3) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   geom_path(data = C_plasma_predicted_df,
#             aes(Time, median/Dose,
#                 colour = Category_OECD)) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA")
#   ),
#   nrow = 2) +
#   labs(x = "Time (years)",
#        y = "Concentration / dose",
#        title = "Plasma") +
#   theme_quick
# CP_plot
# CP_plot <- CP_plot + Patchwork_theme +
#   theme(legend.position = "top")
# 
# ggsave(plot = CP_plot,
#        file = here(OUTPUT, "CP_plot_UA.png"),
#        dpi = 600,
#        width = 24,
#        height = 11,
#        units = "cm")

### Gather all compartment Cs ----
# 
# OneC <- CP_plot
# TwoC <- CL_plot
# ThreeC <- CK_plot
# FourC <- CI_plot
# FiveC <- CA_plot
# SixC <- CR_plot
# 
# Plot_Cs <- OneC/TwoC/ThreeC/FourC/FiveC/SixC +
#   plot_layout(nrow = 6, guides = "collect") &
#   Patchwork_theme +
#   theme(legend.position = "bottom")
# 
# ggsave(plot = Plot_Cs,
#        file = here(OUTPUT, "Plot_Cs_UA.png"),
#        dpi = 600,
#        width = 18,
#        height = 28,
#        units = "cm")

### Half life ----
# 
# HL_predicted_df <- TK_df |>
#   filter(Parameter == "HL")
# 
# HL_predicted_df <- full_join(HL_predicted_df, nC_MW_df)
# 
# HL_predicted_plot <- ggplot() +
#   geom_boxplot(data = HL_predicted_df,
#                aes(x = reorder(PFAS, MW), y = log10(median),
#                    ymin = log10(q025), lower = log10(q25), 
#                    middle = log10(median), 
#                    upper = log10(q75), ymax = log10(q975),
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   guides(color = guide_legend(override.aes = list(fill = NA))) +
#   theme_quick
# HL_predicted_plot
# 
# HL_predicted_plot <- HL_predicted_plot + Patchwork_theme +
#   theme(legend.position = "top",
#         axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0))
#         )
# 
# ggsave(plot = HL_predicted_plot,
#        file = here(OUTPUT,"HL_predicted_plot.png"),
#        dpi = 600,
#        width = 12,
#        height = 11,
#        units = "cm")



### Fraction absorbed and eliminated in urine and feces ----
# 
# FelFa_predicted_df <- TK_df |> 
#   filter(Parameter %in% c("Fa", 
#                           # "FelTotal",
#                           "Ffeces",
#                           "FbilEl",
#                           "Furine")) |> 
#   pivot_wider(names_from = Parameter, values_from = c(median, mean, q25, q75)) |> 
#   mutate(median_Fpreintfeces = median_Ffeces-median_FbilEl,
#          mean_Fpreintfeces = mean_Ffeces-mean_FbilEl,
#          q25_Fpreintfeces = q25_Ffeces-q25_FbilEl,
#          q75_Fpreintfeces = q75_Ffeces-q75_FbilEl,
#          median_cumFpreintfeces = median_Fpreintfeces + median_Furine,
#          mean_cumFpreintfeces = mean_Fpreintfeces + mean_Furine,
#          q25_cumFpreintfeces = q25_Fpreintfeces + q25_Furine,
#          q75_cumFpreintfeces = q75_Fpreintfeces + q75_Furine,
#          median_cumFbile = median_Fpreintfeces + median_Furine + median_FbilEl,
#          mean_cumFbile = mean_Fpreintfeces + mean_Furine + mean_FbilEl,
#          q25_cumFbile = q25_Fpreintfeces + q25_Furine + q25_FbilEl,
#          q75_cumFbile = q75_Fpreintfeces + q75_Furine + q75_FbilEl,
#   ) |>  
#   pivot_longer(median_Fa:q75_cumFbile, names_to = c("What", "Parameter"), values_to = "Value", names_sep = "_") |> 
#   pivot_wider(names_from = "What", values_from = "Value") |> 
#   filter(Parameter != "Ffeces") |> 
#   mutate(Direction = ifelse(Parameter == "Fa", "In",
#                             ifelse(Parameter == "Furine", "Urine",
#                                    "Feces"))) 
# 
# FelFa_predicted_df <- full_join(FelFa_predicted_df, nC_MW_df)
# 
# Clean_FelFa <- FelFa_predicted_df |> 
#   filter(Parameter %in% c("FbilEl", "Fpreintfeces", "Furine"))
# Cum_FelFa <- FelFa_predicted_df |> 
#   filter(Parameter %in% c("cumFbile", "cumFpreintfeces", "Furine")) |> 
#   mutate(Parameter = str_replace_all(Parameter, "cum", ""),
#          Parameter = str_replace_all(Parameter, "Fbile", "FbilEl"))
# 
# Fa_predicted_df <- TK_df |> 
#   filter(Parameter %in% c("Fa")) 
# 
# Fa_predicted_df <- full_join(Fa_predicted_df, nC_MW_df)
# 
# Fa_plot <- ggplot()+
#   geom_boxplot(data = Fa_predicted_df,
#                aes(x = reorder(PFAS, MW), y = median,
#                    ymin = q025, lower = q25, 
#                    middle = median, 
#                    upper = q75, ymax = q975,
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   theme_quick
# 
# Fe_plot <- ggplot() +
#   geom_bar(data = Clean_FelFa,
#            aes(x = reorder(PFAS, MW), y = median, 
#                fill = Parameter),
#            stat = "identity",
#            position = "stack",
#            width = 0.5) +
#   scale_fill_manual(name = "Parameter",
#                     values = c("FbilEl" = "#B198B9",
#                                "Fpreintfeces" = "#8C6897",
#                                "Furine" = "#89D499"),
#                     labels = c("FbilEl" = "Feces: bile",
#                                "Fpreintfeces" = "Feces: pre-systemic",
#                                "Furine" = "Urine")) +
#   geom_errorbar(data = Cum_FelFa,
#                 aes(x = reorder(PFAS, MW),
#                     ymin = q25, ymax = q75,
#                     colour = Parameter),
#                 position = "dodge",
#                 width = .3,
#                 linewidth = .4
#   ) +
#   scale_color_manual(name = "Parameter",
#                      values = c("FbilEl" = "#683975",
#                                 "Fpreintfeces" = "#440154",
#                                 "Furine" = "#2A7741"),
#                      labels = c("FbilEl" = "Feces: bile",
#                                 "Fpreintfeces" = "Feces: pre-systemic",
#                                 "Furine" = "Urine")) +
#   labs(x = "", 
#        y = "Cumulative fraction excreted") +
#   theme_quick +
#   theme(legend.position = "right")
# Fe_plot
# 
# Fe_plot <- Fe_plot + Patchwork_theme +
#   theme(legend.position = "top",
#         legend.text = element_text(size = 35),
#         legend.key.size = unit(0.8, "lines"),
#         legend.spacing.y = unit(0.1, "cm"),
#         legend.spacing.x = unit(0.1, "cm"),
#         legend.margin = margin(0,0,0,0),
#         axis.text.x = element_text(size = 35, face = "bold", margin = margin(0, 0, 0, 0)),
#         axis.text.y = element_text(size = 35, margin = margin(0, 0, 0, 0)),   
#         axis.title = element_text(size = 45) #, face = "bold"
#   )
# 
# 
# ggsave(plot = Fe_plot,
#        file = here(OUTPUT,"Fe_plot.png"),
#        dpi = 600,
#        width = 16,
#        height = 12,
#        units = "cm")

### Urinary and fecal clearances ----
# 
# CLurCLfe_predicted_df <- TK_df |> 
#   filter(Parameter %in% c("CLfe", 
#                           "CLur")) 
# 
# CLurCLfe_predicted_df <- full_join(CLurCLfe_predicted_df, nC_MW_df)
# 
# RurRfe_predicted_df <- TK_df |> 
#   filter(Parameter %in% c("RfeEl", 
#                           "RurEl"))
# RurRfe_predicted_df <- full_join(RurRfe_predicted_df, nC_MW_df)
# 
# CL_plot <- ggplot() +
#   geom_boxplot(data = CLurCLfe_predicted_df,
#                aes(x = reorder(PFAS, MW), y = log10(median),
#                    ymin = log10(q025), lower = log10(q25), 
#                    middle = log10(median), 
#                    upper = log10(q75), ymax = log10(q975),
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ Parameter, nrow = 2, 
#              labeller = as_labeller(c("CLfe" = "Biliary/ fecal",
#                                       "CLur" = "Renal")),
#              scales = "free_y") +
#   theme_quick
# CL_plot
# CL_plot <- CL_plot + Patchwork_theme +
#    theme(axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0)),
#          legend.position = "top"
#    )
#  
#  
#  ggsave(plot = CL_plot,
#         file = here(OUTPUT,"CL_plot.png"),
#         dpi = 600,
#         width = 12,
#         height = 11,
#         units = "cm")

### Final plot Uncertainty Analysis ----
# UA_one <- CP_plot + 
#   labs(title = "A.",
#        x = "Time (years)",
#        y = "Plasma concentration (μg/L)") +
#   guides(fill = FALSE) + guides(colour = FALSE) + 
#   Patchwork_theme + 
#   theme(legend.position = "none")
# 
# UA_two <- Fa_plot + 
#   labs(title = "B." ,
#        x = "",
#        y = "Fraction absorbed") +
#   guides(colour = FALSE) + guides(linetype = FALSE) + 
#   Patchwork_theme + 
#   theme(legend.position = "none", 
#         axis.text.x = element_text(angle = 90))
# 
# UA_three <- CL_plot + 
#   labs(title = "C." ,
#        x = "",
#        y = "Clearance (L/day) (log10 scale)") +
#   guides(fill = FALSE) + guides(colour = FALSE) + 
#   Patchwork_theme + 
#   theme(legend.position = "none", 
#         axis.text.x = element_text(angle = 90))
# 
# UA_four <- HL_predicted_plot + 
#   labs(title = "D.",
#        x = "",
#        y = "Half-life (years) (log10 scale)") +
#   Patchwork_theme + 
#   theme(legend.position = "bottom", 
#         axis.text.x = element_text(angle = 90))
# 
# UA_plot <- UA_one / (UA_two + UA_three + UA_four) +
#   plot_layout(heights = c(0.7, 1),
#               guides = "collect")  &
#   theme(plot.margin = margin(0.1, 0.1, b = 0, 0.1, "cm"),
#         legend.box.margin = margin(t = -0.5, 0, b = 0.1, 0, "cm"),
#         legend.position='bottom')
# 
# ggsave(plot = UA_plot,
#        file = here(OUTPUT,"UA_plot_log.png"),
#        dpi = 600,
#        width = 24,
#        height = 16,
#        units = "cm")

## Post process parameters ----

# P_df <- list.files(path = here(OUTPUT),
#                    pattern = "\\_final_model_parameters.csv$") |>
#   map_df(~ read_csv(file = here(OUTPUT, .),
#                     id = "name"))|>
#   select(!"...2") |>
#   # rename(Iteration = "...2") |>
#   mutate(PFAS = str_replace_all(name,
#                                 "C:/Users/Admin/Desktop/PhD_Chrysa/CP_L_R/PARC_PFAS_PBK/Output/GoF_upTo3days/",
#                                 "")) |>
#   mutate(PFAS = str_replace_all(PFAS,
#                                 "_final_model_parameters.csv","")) |>
#   select(!name) #|> unite("PFAS_Iteration", c(PFAS,Iteration), sep = "_")
# 
# 
# Visualise_P <- P_df |>
#   # select(!Clear_net_reabsorption_unscaled) |>
#   rename_with(~str_remove(.x, "_")) |>
#   rename_with(~str_remove(.x, "_")) |>
#   group_by(PFAS) |>
#   summarise(across(where(is.numeric),
#                    .fns =
#                      list(median = ~median(.,na.rm = TRUE),
#                           q25 = ~quantile(., 0.25, na.rm = TRUE),
#                           q75 = ~quantile(.,0.75, na.rm = TRUE),
#                           q025 = ~quantile(., 0.025, na.rm = TRUE),
#                           q975 = ~quantile(.,0.975, na.rm = TRUE)))
#   ) |>
#   pivot_longer(!PFAS, names_to = c("Parameter", "stat"),
#                names_sep = "_",values_to = "Value") |>
#   pivot_wider(names_from = stat, values_from = Value)
# 
# Visualise_P <- full_join(Visualise_P, nC_MW_df)
# write.csv(Visualise_P, here(OUTPUT,"All_PFAS_finalParameters.csv"), row.names = FALSE)
# 
# Input_PhysicoChemical_Parameters <- str_remove_all(Input_PhysicoChemical_Parameters, "_")
# log_Input_PhysicoChemical_Parameters <- c("logDMembraneLipidWater",
#                                           "logDStructuralProteinWater",
#                                           "logDAlbuminWater",
#                                           "logDStorageLipidWater",
#                                           "logDGlobulinWater",
#                                           "PappSI")
# Visualise_P <- Visualise_P |> 
#   mutate(median=ifelse(Parameter %in% log_Input_PhysicoChemical_Parameters, log(median), median),
#          q025=ifelse(Parameter %in% log_Input_PhysicoChemical_Parameters, log(q025), q025),
#          q25=ifelse(Parameter %in% log_Input_PhysicoChemical_Parameters, log(q25), q25),
#          q75=ifelse(Parameter %in% log_Input_PhysicoChemical_Parameters, log(q75), q75),
#          q975=ifelse(Parameter %in% log_Input_PhysicoChemical_Parameters, log(q975), q975))
# 
# 
# 
# Input_PhysicoChemical_Parameters_plot <- Visualise_P |>
#   filter(Parameter %in% Input_PhysicoChemical_Parameters) |>
#   ggplot()+
#   geom_boxplot(aes(x = reorder(PFAS, MW),
#                    y = median,
#                    ymin = q025, lower = q25,
#                    middle = median,
#                    upper = q75, ymax = q975,
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   guides(color = guide_legend(override.aes = list(fill = NA))) +
#   # facet_wrap(~Parameter, scales = "free_y")+
#   facet_wrap(~factor(Parameter,Input_PhysicoChemical_Parameters),
#              scales = "free_y") +
#   labs(x = "",
#        y = "") +
#   theme_quick +
#   theme(axis.text.x = element_text(angle = 90))
# 
# Intermediate_TK_Parameters <- str_remove_all(Intermediate_TK_Parameters, "_")
# 
# Intermediate_TK_Parameters_plot <- Visualise_P |>
#   filter(Parameter %in% Intermediate_TK_Parameters) |>
#   ggplot()+
#   geom_boxplot(aes(x = reorder(PFAS, MW),
#                    y = median,
#                    ymin = q025, lower = q25,
#                    middle = median,
#                    upper = q75, ymax = q975,
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   guides(color = guide_legend(override.aes = list(fill = NA))) +
#   # facet_wrap(~Parameter, scales = "free_y")+
#   facet_wrap(~factor(Parameter,Intermediate_TK_Parameters),
#              scales = "free_y") +
#   labs(x = "",
#        y = "") +
#   theme_quick +
#   theme(axis.text.x = element_text(angle = 90))
# Intermediate_TK_Parameters_plot
# 
# Final_TK_Parameters <- str_remove_all(Final_TK_Parameters, "_")
# Final_TK_Parameters_plot <- Visualise_P |>
#   filter(Parameter %in% Final_TK_Parameters) |>
#   ggplot()+
#   geom_boxplot(aes(x = reorder(PFAS, MW),
#                    y = median,
#                    ymin = q025, lower = q25,
#                    middle = median,
#                    upper = q75, ymax = q975,
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   guides(color = guide_legend(override.aes = list(fill = NA))) +
#   # facet_wrap(~Parameter, scales = "free_y")+
#   facet_wrap(~factor(Parameter,Final_TK_Parameters),
#              scales = "free_y") +
#   labs(x = "",
#        y = "") +
#   theme_quick +
#   theme(axis.text.x = element_text(angle = 90))
# Final_TK_Parameters_plot
# 
# 
# Final_Physio_Parameters_plot <- Visualise_P |>
#   filter_out(Parameter %in% Final_TK_Parameters) |>
#   filter_out(Parameter %in% Input_PhysicoChemical_Parameters) |>
#   filter_out(Parameter %in% Intermediate_TK_Parameters) |>
#   ggplot()+
#   geom_boxplot(aes(x = reorder(PFAS, MW),
#                    y = median,
#                    ymin = q025, lower = q25,
#                    middle = median,
#                    upper = q75, ymax = q975,
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   guides(color = guide_legend(override.aes = list(fill = NA))) +
#   # facet_wrap(~Parameter, scales = "free_y")+
#   facet_wrap(~Parameter,
#              scales = "free_y") +
#   labs(x = "",
#        y = "") +
#   theme_quick +
#   theme(axis.text.x = element_text(angle = 90))
# Final_Physio_Parameters_plot
# 
# 
# Input_PhysicoChemical_Parameters_plot = Input_PhysicoChemical_Parameters_plot +
#   Patchwork_theme + theme(axis.text.x = element_text(angle = 90))
# 
# Intermediate_TK_Parameters_plot  = Intermediate_TK_Parameters_plot +
#   Patchwork_theme + theme(axis.text.x = element_text(angle = 90))
# 
# Final_TK_Parameters_plot = Final_TK_Parameters_plot +
#   Patchwork_theme + theme(axis.text.x = element_text(angle = 90))
# 
# Final_Physio_Parameters_plot = Final_Physio_Parameters_plot +
#   Patchwork_theme + theme(axis.text.x = element_text(angle = 90))
# 
# ggsave(plot = Input_PhysicoChemical_Parameters_plot,
#        file = here(OUTPUT, paste0("all_PFAS", "_Input_PhysicoChemical_Parameters_plot.png")),
#        dpi = 600,
#        width = 24,
#        height = 18,
#        units = "cm")
# 
# ggsave(plot = Intermediate_TK_Parameters_plot,
#        file = here(OUTPUT, paste0("all_PFAS", "_Intermediate_TK_Parameters_plot.png")),
#        dpi = 600,
#        width = 18,
#        height = 18,
#        units = "cm")
# 
# ggsave(plot = Final_TK_Parameters_plot,
#        file = here(OUTPUT, paste0("all_PFAS", "_Final_TK_Parameters_plot.png")),
#        dpi = 600,
#        width = 16,
#        height = 16,
#        units = "cm")
# 
# ggsave(plot = Final_Physio_Parameters_plot,
#        file = here(OUTPUT, paste0("all_PFAS", "_Final_Physio_Parameters_plot.png")),
#        dpi = 600,
#        width = 24,
#        height = 18,
#        units = "cm")
# 

# GoF ----
Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |>
  rename(PFAS = Acronym) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

## Against Fa ----

# Fa_measured_df <- Measured_df |>
#   filter(Parameter %in% "Fa") 
# 
# Fa_measured_df <- full_join(Fa_measured_df, nC_MW_df) |>
#   select(PFAS, Category_OECD,
#          Parameter,Value, Unit,
#          MW, n_C,
#          Reference) |>
#   # na.omit() |> 
#   mutate(PFAS = factor(PFAS, rev(PFAS_list)))
# Fa_predicted_df <- Fa_predicted_df |> 
#   mutate(PFAS = factor(PFAS, rev(PFAS_list)))
# 
# Fa_MesPred_plot <- ggplot() +
#   geom_point(data = Fa_measured_df,
#              aes(x = reorder(PFAS, MW), y = Value,
#                  shape = Reference), #,             color = Category_OECD),
#              color = "grey60") +
#   geom_boxplot(data = Fa_predicted_df,
#                aes(x =  reorder(PFAS, MW), y = median,
#                    ymin = q025, lower = q25,
#                    middle = median,
#                    upper = q75, ymax = q975,
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   labs(x = "",y = "Fraction absorbed")+
#   theme_quick
# 
# Fa_MesPred_plot
# 
# 
# 
# Fa_MesPred_plot <- ggplot() +
#   geom_point(data = Fa_measured_df,
#              aes(x = reorder(PFAS, MW), y = Value,
#                  shape = Reference), #,             color = Category_OECD), 
#              color = "grey60") +
#   geom_boxplot(data = Fa_predicted_df,
#                aes(x = reorder(PFAS, MW), y = median,
#                    ymin = q025, lower = q25, 
#                    middle = median, 
#                    upper = q75, ymax = q975,
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   # facet_wrap(~ Parameter, nrow = 2,
#   #            labeller = as_labeller(c("CLfe" = "Biliary/ fecal",
#   #                                     "CLur" = "Renal")),
#   #            scales = "free_y") +
#   labs(x = "",y = "Fraction absorbed")+
#   theme_quick
# Fa_MesPred_plot
# # Fa_MesPred_plot <- Fa_MesPred_plot + Patchwork_theme +
# #   theme(legend.position = "right",
# #         axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0))
# #   )
# # 
# # ggsave(plot = Fa_MesPred_plot,
# #        file = here(OUTPUT,"Fa_MesPred_plot.png"),
# #        dpi = 600,
# #        width = 16,
# #        height = 11,
# #        units = "cm")
# 
# 
# 
# save_Fa_measured_df <- Fa_measured_df |> 
#   mutate(Fa = Value, 
#          Origin = "Measured") |> 
#   select(PFAS, Fa, Origin)
# 
# save_Fa_predicted_df <- Fa_predicted_df |> 
#   select(PFAS, median) |> 
#   mutate(Origin = "Predicted") |> 
#   rename(Fa = median)
# 
# save_Fa <- rbind(save_Fa_measured_df, save_Fa_predicted_df) |> 
#   pivot_wider(names_from = Origin, values_from = Fa) |> 
#   mutate(Predicted_q025 = Fa_predicted_df$q025,
#          Predicted_q975 = Fa_predicted_df$q975)
# 
# write.csv(save_Fa, 
#           file = here(OUTPUT, "Fa.csv"), 
#           row.names = FALSE)
# write_csv(save_Fa, 
#           file = here(OUTPUT, "Fa.csv"))
# 
# write.csv(MC.design, 
#           file = here(OUTPUT, paste0(PFAS_i,"_LHS_MC_ExpDesign.csv")), 
#           row.names = FALSE)
# 

## Against HL measured ----
# 
# raw_HL_measured <- Measured_df |>
#   filter(Parameter == "Halflife",
#          PFAS != "PFTrDA") |>
#   select(PFAS, Category_OECD, Reference,
#          Value, Range_CV_Std_min, Range_CV_Std_max,
#          n, sex, Age_years) |>
#   mutate(n = round(n,0))
# 
# set.seed(1234)
# 
# HL_measured <- raw_HL_measured %>%
#   mutate(n = as.numeric(n),
#          Value = as.numeric(Value),
#          min = as.numeric(Range_CV_Std_min),
#          max = as.numeric(Range_CV_Std_max)) |>
#   mutate(meandistance = (Value-min + max-Value)/2,
#          sd = ifelse(meandistance > Value, meandistance/4, meandistance), #Value/4, meandistance/2)
#          sd = sd/2 #,
#          # sd = ifelse(Reference == "Yao_2023" & PFAS == "PFBA", sd*2, sd),
#          # sd = ifelse(Reference == "Yao_2023" & PFAS == "PFPeA", sd*2, sd),
#          # sd = ifelse(Reference == "Xu_2020" & PFAS == "PFBS", sd*1.1, sd)
#   ) |>
#   rowwise() |>
#   mutate(r_sex = case_when( n != 1 & sex == "F" ~ list(rep("F", n)),
#                             n != 1 & sex == "M" ~ list(rep("M", n)),
#                             TRUE ~ list(sex)
#   ),
#   r_Value = ifelse(n != 1,
#                    list(rnorm(n, Value, sd)),
#                    list(Value))
#   ) |>
#   unnest(cols = r_sex:r_Value,
#          keep_empty = TRUE) |>
#   mutate(Value = round(Value, 3),
#          r_Value = round(r_Value, 3),
#          min = ifelse(min <= 0, 0.00001, min),
#          r_Value = ifelse(n != 1 & r_Value < min, min, r_Value) ,
#          r_Value = ifelse(n != 1 & r_Value > max, max, r_Value),
#          r_Value = round(r_Value, 3)
#   ) |>
#   select(!sex) |>
#   rename(mean = Value,
#          Value = r_Value,
#          sex = r_sex)
# 
# 
# HL_measured_df <- full_join(HL_measured, nC_MW_df)
# 
# HL_MesPred_plot <- 
#   ggplot() +
#   geom_jitter(data = HL_measured_df,
#               aes(reorder(PFAS, MW),
#                   log10(Value),
#                   shape = Category_OECD),
#               color = "grey60",
#               width = 0.2,
#               alpha = 0.9,
#               size = .3) +
#   geom_boxplot(data = HL_predicted_df,
#                aes(x = reorder(PFAS, MW), y = log10(median),
#                    ymin = log10(q025), lower = log10(q25), 
#                    middle = log10(median), 
#                    upper = log10(q75), ymax = log10(q975),
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   guides(color = guide_legend(override.aes = list(fill = NA))) +
#   labs(x = "",
#        y = "Half life (years) (log10 scale)") +
#   theme_quick
# 
# HL_MesPred_plot
# 
# HL_MesPred_plot <- HL_MesPred_plot + Patchwork_theme +
#   theme(legend.position = "right",
#         axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0))
#   )
# 
# ggsave(plot = HL_MesPred_plot,
#        file = here(OUTPUT,"HL_MesPred_plot.png"),
#        dpi = 600,
#        width = 16,
#        height = 11,
#        units = "cm")

## Against measured clearances ----
# 
# CL_UrFe_measured_df <- Measured_df |>
#   filter(Parameter %in% c("CL_urinary", 
#                           "CL_fecal")) |> 
#   # filter_out(Reference == "Abraham_2024") |>
#   mutate(Parameter = str_replace_all(Parameter, "CL_fecal", "CLfe"),
#          Parameter = str_replace_all(Parameter, "CL_urinary", "CLur")
#   )
# 
# CL_UrFe_measured_df <- full_join(CL_UrFe_measured_df, nC_MW_df) |>
#   select(PFAS, Category_OECD,
#          Parameter,Value, Unit,
#          MW, n_C,
#          Reference) |>
#   na.omit()
# 
# CL_MesPred_plot <- ggplot() +
#   geom_point(data = CL_UrFe_measured_df,
#              aes(x = reorder(PFAS, MW), y = log10(Value),
#                  shape = Reference), #,             color = Category_OECD), 
#              color = "grey60") +
#   geom_boxplot(data = CLurCLfe_predicted_df,
#                aes(x = reorder(PFAS, MW), y = log10(median),
#                    ymin = log10(q025), lower = log10(q25), 
#                    middle = log10(median), 
#                    upper = log10(q75), ymax = log10(q975),
#                    fill = Category_OECD, colour = Category_OECD),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ Parameter, nrow = 2,
#              labeller = as_labeller(c("CLfe" = "Biliary/ fecal",
#                                       "CLur" = "Renal")),
#              scales = "free_y") +
#   labs(x = "",y = "Clearance (L/day) (log10 scale)")+
#   theme_quick
# 
# CL_MesPred_plot

# CL_MesPred_plot <- CL_MesPred_plot + Patchwork_theme +
#   theme(legend.position = "right",
#         axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0))
#   )
# 
# ggsave(plot = CL_MesPred_plot,
#        file = here(OUTPUT,"CL_MesPred_plot.png"),
#        dpi = 600,
#        width = 16,
#        height = 11,
#        units = "cm")


## Against plasma measured ----

Measured_Abraham_df <- Measured_df |>
  filter(Reference == "Abraham_2024")

C_plasma_measured_df <- Measured_df |>
  filter(Parameter == "C_plasma") |>
  mutate(Time = as.numeric(Time_days)/365) |>
  select(PFAS, Time, Value)
C_plasma_measured_df <- full_join(C_plasma_measured_df, nC_MW_df)
C_plasma_measured_df <- full_join(C_plasma_measured_df, Dose_df)

C_plasma_measured_df <- C_plasma_measured_df |> 
  filter(Time < 0.006)

C_plasma_predicted_df <- C_df |>
  filter(Tissue == "Plasma")

C_plasma_predicted_df <- full_join(C_plasma_predicted_df, nC_MW_df)
C_plasma_predicted_df <- full_join(C_plasma_predicted_df, Dose_df)

C_plasma_predicted_df <- C_plasma_predicted_df |> 
  filter_out(PFAS %in% c("PFPeS", "PFHpS", "PFTrDA"))
C_plasma_measured_df <- C_plasma_measured_df |> 
  filter_out(PFAS %in% c("PFPeS", "PFHpS", "PFTrDA"))

write.csv(C_plasma_predicted_df, here(OUTPUT, "Plasma_predicted.csv"), row.names = FALSE)

ObsVsPred_Plasma_plot <- ggplot() +
  geom_ribbon(data = C_plasma_predicted_df,
              aes(x = Time*365*24,
                  ymin = q025,
                  ymax = q975,
                  fill = Category_OECD),
              alpha = 0.4) +
  geom_ribbon(data = C_plasma_predicted_df,
              aes(x = Time*365*24,
                  ymin = q25,
                  ymax = q75,
                  fill = Category_OECD),
              alpha = 0.3) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
  geom_path(data = C_plasma_predicted_df,
            aes(Time, median,
                colour = Category_OECD)) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  geom_point(data = C_plasma_measured_df, #filter(C_plasma_measured_df, Category_OECD == "PFCAs"),
             aes(Time*365*24, Value),
             colour = "grey10",
             size = 0.05) +
  labs(x = "Time (hours)",
       y = "Concentration (ng/ml)") +
  facet_wrap(~ factor(PFAS, rev(PFAS_list)),
             nrow = 3,
             scales = "free_y") +
  theme_quick
ObsVsPred_Plasma_plot

ObsVsPred_Plasma_plot <- ObsVsPred_Plasma_plot + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 35, face = "bold", margin = margin(0, 0, b = 0.5, 0)),
    plot.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
    plot.title.position = "plot",
    panel.background = element_blank(), #rect(fill = 'white', color = NA)
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_line(color = "grey90",linewidth = 0.1),
    
    axis.title = element_text(size = 35, face = "bold"),
    axis.text.x = element_text(size = 25, face = "bold", margin = margin(0, 0, 0, 0)),
    axis.text.y = element_text(size = 25, face = "bold", margin = margin(0, 0, 0, 0)),
    axis.line = element_blank(), 
    axis.ticks = element_blank(),
    
    legend.position = "bottom",
    legend.title = element_blank(), 
    legend.text = element_text(size = 25),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.spacing.x = unit(0.05, "cm"),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
    
    # strip.text.x = element_blank()
    strip.text = element_text(hjust = 0, size = 35, face = "bold", margin = margin(0,0,0,0))
  )
  

ggsave(plot = ObsVsPred_Plasma_plot,
       file = here(OUTPUT, "CP_plot_UA_FreeY.png"),
       dpi = 600,
       width = 14,
       height = 8,
       units = "cm")

### obs vs pred lines ----
ObsPlasma <- C_plasma_measured_df |> 
  mutate(Time = signif(as.numeric(Time),3),
         q025 = Value, 
         q975 = Value) |> 
  rename(Observed = Value)

Obstimepoints <- signif(unique(ObsPlasma$Time),3)


PredPlasma <- C_plasma_predicted_df |> 
  mutate(Time = signif(Time, 3)) |> 
  filter(Time %in% Obstimepoints) |> 
  select(PFAS, Time, median, q025, q975, Category_OECD, MW, n_C) |> 
  rename(Predicted = median)

Predtimepoints <- signif(unique(PredPlasma$Time),3)

ObsPlasma <- ObsPlasma |> 
  filter(Time %in% Predtimepoints)

# MedianPredPlasma <- PredPlasma |> select(!c(q025, q975)) 
# MedianObsPlasma <- ObsPlasma |> select(!c(q025, q975))

ObsPredPlasma <- left_join(ObsPlasma, PredPlasma, by = c("PFAS", "Time", "Category_OECD", "MW", "n_C")) 

s_reg <- NULL
s_r.sq <- NULL
s_RSE <- NULL
s_RMSE <- NULL
PFAS_listreg <- unique(ObsPlasma$PFAS)
nPFASreg <- length(PFAS_listreg)
for (i in 1:nPFASreg) {
  PFAS_i <- PFAS_listreg[i]
  
  Forreg <- filter(ObsPredPlasma, PFAS == PFAS_i) |> 
    select(Predicted, Observed) 
  
  reg <- lm(Predicted ~ Observed, data = Forreg)
  r.sq <- summary(reg)$r.squared
  RSE <- summary(reg)$sigma
  RMSE <- sqrt(mean((Forreg$Observed - Forreg$Predicted)^2))
  
  s_reg[i] <- reg
  s_r.sq[i] <- r.sq
  s_RSE[i] <- RSE
  s_RMSE[i] <- RMSE
}

names(s_r.sq) = names(s_RSE) = names(s_RMSE) = PFAS_listreg

RegResults <- tibble(
  PFAS = PFAS_listreg, 
  R_squared = s_r.sq,
  RSE = s_RSE,
  RMSE = s_RMSE
)
write.csv(RegResults, here(OUTPUT, "RegResults.csv"))
write.csv(ObsPredPlasma, here(OUTPUT, "RawObsPred.csv"))


ObsPredPlasma_regr <- ggplot(ObsPredPlasma,
                             aes(#x = Predicted, y = Observed, 
                               x = log10(Observed), y = log10(Predicted),
                               colour = Category_OECD)) +
  # geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.5) + #, color = "grey50") +
  # geom_abline(intercept = 0, slope = 2, linetype = "dashed", linewidth = 0.3) +#, color = "grey50") +
  # geom_abline(intercept = 0, slope = 1/2, linetype = "dashed", linewidth = 0.3) + #, color = "grey50") +
  # geom_abline(intercept = 0, slope = 10, linetype = "dotted", linewidth = 0.3) + #, color = "grey50") +
  # geom_abline(intercept = 0, slope = 1/10, linetype = "dotted", linewidth = 0.3) + #, color = "grey50") +
  # 
  geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.5) + #, color = "grey50") +
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", linewidth = 0.3) +#, color = "grey50") +
  geom_abline(intercept = log10(0.5), slope = 1, linetype = "dashed", linewidth = 0.3) + #, color = "grey50") +
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotted", linewidth = 0.3) + #, color = "grey50") +
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotted", linewidth = 0.3) + #, color = "grey50") +
  # geom_errorbar(aes(ymin = log10(Cmin), ymax = log10(Cmax)),
  #               alpha = 0.3)+
  facet_wrap(~ factor(PFAS, rev(PFAS_list)),
             nrow = 2) +
  geom_point(size = 1, alpha = 0.9) +
  scale_color_viridis(option = "D", discrete = TRUE, end=0.7) +
  coord_cartesian(xlim = c(-3, 0.5), ylim = c(-3, 0.5)) +
  # coord_cartesian(xlim = c(0, 0.60), ylim = c(0, 0.60)) +
  labs(#title = "B. ZEN Blood",
    x = expression(paste("Observed concentration (log"[10] ~ "(ng/ml))")),
    y = expression(paste("Predicted concentration (log"[10] ~ "(ng/ml))")),
    color = "Category OECD"
  ) +
  theme_quick 
ObsPredPlasma_regr
ObsPredPlasma_regr <- ObsPredPlasma_regr + Patchwork_theme +
  theme(legend.position = "bottom")

ggsave(plot = ObsPredPlasma_regr,
       file = here(OUTPUT, "ObsPredPlasma_regr.png"),
       dpi = 600,
       width = 24,
       height = 8,
       units = "cm")





## Ratios ----

# Measured_ratios <- Measured_df |>
#   filter(Parameter %in% c("KPR", "LPR")) |>
#   select(PFAS, Parameter, Value, Reference) |>
#   mutate(Value = as.numeric(Value))
# 
# Measured_ratios <- full_join(Measured_ratios, nC_MW_df) |> 
#   na.omit()
# 
# Ratios_PFAS <- unique(Measured_ratios$PFAS)
# 
# Ratios <- C_df |>
#   select(!"...1") |>
#   filter(Tissue %in% c("Plasma", "Kidney", "Liver"),
#          PFAS %in% Ratios_PFAS) |>
#   filter(Time > 0.2) |>
#   pivot_longer(median:q975, names_to = "stat", values_to = "Value") |>
#   pivot_wider(names_from = Tissue, values_from = Value) |>
#   mutate("LPR" = Liver/Plasma,
#          "KPR" = Kidney/Plasma) |>
#   select(!c(Kidney, Liver, Plasma)) |>
#   pivot_longer(LPR:KPR, names_to = "Tissue", values_to = "Value") |>
#   pivot_wider(names_from = stat, values_from = Value)
# 
# KPR <- Ratios |>
#   filter(Tissue == "KPR") |>
#   group_by(PFAS) |>
#   summarise(across(median:q975, .fns =
#                      list(median = ~median(.,na.rm = TRUE),
#                           q25 = ~quantile(., 0.25, na.rm = TRUE),
#                           q75 = ~quantile(., 0.75, na.rm = TRUE),
#                           q025 = ~quantile(., 0.025, na.rm = TRUE),
#                           q975 = ~quantile(.,0.975, na.rm = TRUE)))) |>
#   pivot_longer(!PFAS,
#                names_to = "Variable",
#                values_to = "Value") |>
#   filter(Variable %in% c("median_median",
#                          "q25_q25",
#                          "q75_q75",
#                          "q025_q025",
#                          "q975_q975")) |>
#   separate(Variable, into = c("Keep", "not")) |>
#   select(!not) |>
#   pivot_wider(names_from = Keep, values_from = Value) |>
#   mutate(Parameter = "KPR")
# 
# LPR <- Ratios |>
#   filter(Tissue == "LPR") |>
#   group_by(PFAS) |>
#   summarise(across(median:q975, .fns =
#                      list(median = ~median(.,na.rm = TRUE),
#                           q25 = ~quantile(., 0.25, na.rm = TRUE),
#                           q75 = ~quantile(., 0.75, na.rm = TRUE),
#                           q025 = ~quantile(., 0.025, na.rm = TRUE),
#                           q975 = ~quantile(.,0.975, na.rm = TRUE)
#                      ))) |>
#   pivot_longer(!PFAS,
#                names_to = "Variable",
#                values_to = "Value") |>
#   filter(Variable %in% c("median_median",
#                          "q25_q25",
#                          "q75_q75",
#                          "q025_q025",
#                          "q975_q975")) |>
#   separate(Variable, into = c("Keep", "not")) |>
#   select(!not) |>
#   pivot_wider(names_from = Keep, values_from = Value) |>
#   mutate(Parameter = "LPR")
# 
# Ratios <- rbind(KPR, LPR)
# Ratios <- full_join(nC_MW_df, Ratios) |> na.omit()
# 
# ObsVsPred_Ratios_plot <- ggplot() +
#   geom_point(data = Measured_ratios,
#              aes(x = PFAS,
#                  y = Value,
#                  shape = Reference),
#              colour = "grey60",
#              size = 1) +
#   geom_boxplot(data = Ratios,
#                aes(x = factor(PFAS, c("PFUnDA",
#                                       "PFDA",
#                                       "PFOS",
#                                       "PFNA",
#                                       "PFOA",
#                                       "PFHxS",
#                                       "PFHpA"
#                )), 
#                y = median,
#                ymin = q025, lower = q25, 
#                middle = median, 
#                upper = q75, ymax = q975,
#                fill = Category_OECD,
#                colour = Category_OECD
#                ),
#                stat = "identity",
#                linewidth = 0.3,
#                width = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.5) +
#   scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~Parameter,
#              labeller = as_labeller(c("KPR" = "Kidney/Plasma ratio",
#                                       "LPR" = "Liver/Plasma ratio")),
#              scales = "free_y"
#   ) +
#   labs(x = "",
#        y = "Tissue / Plasma ratio")+
#   theme_quick
# ObsVsPred_Ratios_plot
# 
# 
# ObsVsPred_Ratios_plot <- ObsVsPred_Ratios_plot + CP_theme
# 
# ggsave(plot = ObsVsPred_Ratios_plot,
#        file = here(OUTPUT, "ObsVsPred_Ratios_plot.png"),
#        dpi = 600,
#        width = 12,
#        height = 9,
#        units = "cm")
# 
# 
# save_Measured_ratios <- Measured_ratios |> 
#   select(PFAS, Parameter, Value, Reference) |> 
#   rename(Ratio = Value)  
# # mutate(#Reference = "Measured",
# #        q025 = NA,
# #        q975 = NA)
# 
# save_Predicted_ratios <- Ratios |> 
#   select(PFAS, Parameter, median) |>  #, q025, q975) |> 
#   rename(Ratio  = median) |> 
#   mutate(Reference = "Predicted")
# 
# save_ratios <- rbind(save_Measured_ratios, save_Predicted_ratios) |> 
#   pivot_wider(names_from = Reference, values_from = c(Ratio)) |> 
#   mutate(Predicted_q025 = Ratios$q025,
#          Predicted_q975 = Ratios$q975)
# 
# write.csv(save_ratios, here(OUTPUT, "Ratios.csv"), row.names = FALSE)


## Final plot GoF ----
# 
# GoF_one <- ObsVsPred_Plasma_plot + 
#   labs(title = "A.",
#        x = "Time (years)",
#        y = "Plasma concentration (μg/L)") +
#   guides(fill = FALSE) + guides(colour = FALSE) + 
#   Patchwork_theme + 
#   theme(legend.position = "bottom")
# 
# GoF_two <- Fa_MesPred_plot + 
#   labs(title = "B." ,
#        x = "",
#        y = "Fraction absorbed") +
#   guides(colour = FALSE) + guides(linetype = FALSE) + 
#   Patchwork_theme + 
#   theme(legend.position = "bottom", 
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# GoF_three <- CL_MesPred_plot + 
#   labs(title = "C." ,
#        x = "",
#        y = "Clearance (L/day) (log10 scale)") +
#   guides(fill = FALSE) + guides(colour = FALSE) + 
#   Patchwork_theme + 
#   theme(legend.position = "bottom", 
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# GoF_four <- HL_MesPred_plot + 
#   labs(title = "D.",
#        x = "",
#        y = "Half-life (years) (log10 scale)") +
#   Patchwork_theme + 
#   theme(legend.position = "bottom", 
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# GoF_plot <- GoF_one / (GoF_two + GoF_three + GoF_four) +
#   plot_layout(heights = c(0.7, 1))  &
#   theme(plot.margin = margin(0.1, 0.1, b = 0, 0.1, "cm"),
#         legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
#         legend.box.margin = margin(t = -0.5, 0, b = 0.1, 0, "cm"), 
#         legend.position = "bottom")
# 
# ggsave(plot = GoF_plot,
#        file = here(OUTPUT,"GoF_plot_log3.png"),
#        dpi = 600,
#        width = 28,
#        height = 18,
#        units = "cm")
# 
# 
# 
# 









