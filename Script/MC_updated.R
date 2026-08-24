# --------------------------------------------------------------------------- #
# Script for running the Uncertainty Analysis using LHS
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

# Set output storage directory
OUTPUT <- here("Output", "UA_GoF_22Aug") 
dir.create(OUTPUT, recursive = TRUE)

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
# 
# PFAS_list <- c(
#   "PFBA",
#   "PFBS"
# )

Dose_df <- tibble(
  PFAS = PFAS_list,
  Dose = 4) |> 
  mutate(Dose = ifelse(PFAS == "PFBS", 18.4,
                       ifelse(PFAS == "HFPO-DA" , 19.6,
                              Dose)))

# PFAS_list <- Dose_df$PFAS

nPFAS <- length(PFAS_list)


## Tissue constituents for partition coefficient calculation
Physio.data <- read_csv(here("Input", "TissueComposition.csv")) %>% 
  filter(Species == "Human") %>%
  filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
  select(- c(Comment, Species)) %>% 
  mutate(across(-Tissue, as.numeric))


## All PBK parameters
parm.phys <- c(
  BW = 64.41978,
  QC = 5486.301,
  Hct = 0.3645031,
  VIc = 0.016,
  VILc = 0.03848451,
  SA_SIc = 1099.557,
  VLc = 0.02458126,
  VL_icc = 0.835,
  VL_ecc = 0.161,
  VKc = 0.004630086,
  VPTc = 0.398,
  VPTTc = 0.5470631,
  VPTLc = 0.403,
  VRKLc = 0.283,
  VAc = 0.3426787,
  MAc = 1.487552,
  VPc = 0.065,
  QIc = 0.165,
  QLc = 0.06857434,
  QKc = 0.1761446, 
  QAc = 0.09413654,
  QUrc = 0.022, 
  GFR = 193.2765, # = 0.2*0.1761446*5486.301, #156.6414, 
  QTc = 0.36, 
  tco = 3.36,
  Qmenses = 0.0008561644)

## MC parameters
MC_TK_params <- read_xlsx(here("Input", "UA_TK_Parameters.xlsx")) |> 
  mutate(Predicted = ifelse(Specific == "calculated", "yes", Predicted))|>
  filter(Predicted == "yes") 

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
                            parm.phys = parm.phys)
  
  # save parameters to check when analysing results
  run_pars <- as.data.frame(pars)
  
  PBK.out <- lsoda(y = A_init,
                   times = seq(0, 365*3, by = 0.1), 
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

Monte_Carlo_PFAS_design <- list()
Monte_Carlo_PFAS_RESULTS <- list()
Monte_Carlo_PFAS_y <- list()
Monte_Carlo_PFAS_savey <- list()#

for (i in 1:nPFAS) {
  
  
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}-n{i}-----------------------------------"))
  
  # physiological parameters; they are kept fix
  parm.phys <- parm.phys
  Dose_ug_day <- filter(Dose_df, PFAS == PFAS_i) |> 
    select(Dose) |> pull()
  
  PhysChemTK_params <- MC_TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(!PFAS) 
  
  parm.chem <- PhysChemTK_params |> 
    mutate(Value = ifelse(PBK_mechanism == "Transport", log(Value), Value)) |> 
    mutate(cValue = ifelse(PBK_mechanism == "Transport", log(cValue), cValue)) 
    
  
  parm.chem <- parm.chem |> 
    mutate(meanlog = ifelse(Distribution == "lognormal", cValue, NA),
           sdlog = ifelse(Distribution == "lognormal", RSE, NA),
           # Value = ifelse(PBK_mechanism == "Scaling", Value*1e-6, Value),
           # RSE = ifelse(PBK_mechanism == "Scaling", RSE*1e-6, RSE),
           min = ifelse(PBK_mechanism == "Scaling", 1, NA),#1, NA),
           # min = ifelse(PBK_mechanism == "Scaling", 
           #              ifelse((Value - RSE)>0, Value - RSE, 0.1), NA),#1, NA),
           max = ifelse(PBK_mechanism == "Scaling", cValue,NA),#Value, NA),
           # max = ifelse(PBK_mechanism == "Scaling", Value + RSE,NA),#Value, NA),
           Distribution = ifelse(Distribution == "lognormal", "qlnorm", "qunif"))

  MC.parms <- parm.chem$Parameter
  MC.factors <- length(MC.parms) 
  
  ## Latin Hypercube Sampling ---- 
  # this optimises the Monte Carlo sampling, also used by McNally & Loizou
  
  set.seed(10010) # set seeds for reproducibility
  n = 1000 # number of samples per parameter
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
    facet_wrap(~Parameter, scales = "free")+
    coord_flip()+
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
  
names(Monte_Carlo_PFAS_savey) <- PFAS_list

Save_MC_PBK_results <- rbindlist(Monte_Carlo_PFAS_savey, idcol = "PFAS")


# Produce results ----


## functions ----

source(here("Script", "functions_calculate_results.R"))


## produce result files and plots for all PFAS ----
INPUT <- OUTPUT #here("Output", "Second_try_UA")
# OUTPUT <- here("Output","Second_try_UA", "Processed_results") 
# dir.create(OUTPUT, recursive = TRUE)


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
  
  if(is.list(TK_results)){
    errored_iterations <- which(sapply(TK_results, function(x) length(x) < 3)) #
    message(glue("Removing iterations n{errored_iterations} that give errors. \\"))
    write.csv(errored_iterations, here(OUTPUT, paste0(PFAS_i, "_errored_iterations.csv")), row.names = FALSE)
    savey <- savey[-errored_iterations] 
    TK_results <- sapply(savey, calc_TK_results, OD) 
  } else if(!is_empty(which(is.na(TK_results[2,]) | TK_results[2,]<=2))){
    errored_iterations <- which(is.na(TK_results[2,]) | TK_results[2,]<=2)
    write.csv(errored_iterations, here(OUTPUT, paste0(PFAS_i, "_errored_iterations.csv")), row.names = FALSE)
    message(glue("Removing iterations n{errored_iterations} that give errors. \\"))
    savey <- savey[-errored_iterations] 
    TK_results <- sapply(savey, calc_TK_results, OD) 
  } else {
    TK_results <- sapply(savey, calc_TK_results, OD) 
  }
  
  
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
  
  Parameters_df <- Results_df |>
    select(GFR:Clear_reabsorption) |>
    select(!starts_with("CF_")) |> 
    distinct() |> 
    mutate(fup_GFR = GFR) |> 
    select(!GFR) 

  write.csv(Parameters_df, file = here(OUTPUT, paste0(PFAS_i, "_final_model_parameters.csv")))
  
  P_long_df <- Parameters_df |>
    pivot_longer(everything(),
                 names_to = "Parameter",
                 values_to = "Value")
  
  Parameters_plot <- ggplot(P_long_df, aes(x = Value)) +
    geom_histogram() +
    facet_wrap(~ factor(Parameter, 
                        c("PI",   
                          "PL",        
                          "PK",          
                          "PA",         
                          "PR",  
                          "fup", 
                          "fu_PTL",     
                          "fu_L_ec",     
                          "fu_PT_ec",
                          "SF_OATP2B1_I",          
                          "SF_OATP2B1",              
                          "SF_OATP1B3",                     
                          "SF_OATP1B1",             
                          "SF_NTCP",               
                          "SF_BSEP",                        
                          "SF_OAT1",             
                          "SF_OAT3",                
                          "SF_OAT4",    
                          "tco",
                          "Clear_apparent",       
                          "Clear_OATP2B1_I",          
                          "Clear_uptake_I",                 
                          "Clear_OATP1B1",              
                          "Clear_OATP1B3",                
                          "Clear_OATP2B1_L",                
                          "Clear_NTCP",                  
                          "Clear_BSEP",                 
                          "Clear_uptake_L",                 
                          "Clear_excretion_L",     
                          "fup_GFR",
                          "Clear_OAT3",              
                          "Clear_OAT1",                      
                          "Clear_OAT4",              
                          "Clear_reabsorption",
                          "QT"
                        )), scales = "free") +
    labs(title = paste0(PFAS_i, "final model parameters"),
         x = "Value",
         y = "Frequency") +
    coord_flip()+
    theme_quick
  
  Parameters_plot = Parameters_plot + Patchwork_theme +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(plot = Parameters_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_final_model_parameters_plot.png")),
         dpi = 600,
         width = 24,
         height = 18,
         units = "cm")
  
  
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
  
  
  C_df <- Results_df |>
    select(time, CI, CL_ec, CL_ic, CPTL, CPTT, CRKL, CRKT, CA, CR, CP, Iteration) |>
    mutate(PFAS = PFAS_i) |>
    mutate(CL = (CL_ec*VL_ec + CL_ic*VL_ic)/VL,
           CK = (CPTL*VPTL + CPTT*VPTT + CRKL*VRKL + CRKT*VRKT)/VK,
           time = time/365) |> # time to years
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
    
    C_summary <- calc_C_summary(C_df, Tissue_j) |>
      mutate(Tissue = Tissue_j_name,
             PFAS = PFAS_i)
    
    C_Tissues_list[[j]] <- C_summary
    
  }
  
  C_Tissues <- do.call(rbind, C_Tissues_list)
  
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
# TK_df <- read_csv("Output/CleanInitialParams_UncertaintyAnalysis/TK_results_all_PFAS.csv",
#                                 trim_ws = FALSE)
C_df <-
  list.files(path = here(OUTPUT),
             pattern = "\\_TissueC_sumstat.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .)))

Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |>
  rename(PFAS = Acronym) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))


## plot results ----
### Plasma ----
C_plasma_predicted_df <- C_df |>
  filter(Tissue == "Plasma") 

C_plasma_predicted_df <- full_join(C_plasma_predicted_df, nC_MW_df) 



CP_plot <- ggplot() +
  geom_ribbon(data = C_plasma_predicted_df,
              aes(x = Time/365,
                  ymin = q25,
                  ymax = q75,
                  fill = Category_OECD),
              alpha = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
  geom_path(data = C_plasma_predicted_df,
            aes(Time/365, median,
                colour = Category_OECD)) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  facet_wrap(~ factor(PFAS, c( #factor(n_C),
    "PFBA",
    "PFPeA",
    "PFBS",
    "PFHxA",
    "HFPO-DA",
    "PFPeS",
    "PFHpA",
    "ADONA",
    "PFHxS",
    "PFOA",
    "PFHpS",
    "PFNA",
    "PFOS",
    "PFDA",
    "PFUnDA",
    "PFDoDA",
    "PFTrDA")
  ),
  nrow = 2) +
  theme_quick

# CP_plot <- CP_plot + Patchwork_theme +
#   theme(legend.position = "top")
# 
# ggsave(plot = CP_plot,
#        file = here(OUTPUT, "CP_plot_UA.png"),
#        dpi = 600,
#        width = 24,
#        height = 11,
#        units = "cm")

### Half life ----

HL_predicted_df <- TK_df |>
  filter(Parameter == "HL")

HL_predicted_df <- full_join(HL_predicted_df, nC_MW_df)

HL_predicted_plot <- ggplot() +
  geom_errorbar(data = HL_predicted_df,
                aes(x = reorder(PFAS, MW),
                    y = log10(median),
                    colour = Category_OECD,
                    ymin = log10(q25),
                    ymax = log10(q75)),
                width = .01,
                linewidth = .4) +
  geom_point(data = HL_predicted_df,
             aes(x = reorder(PFAS, MW),
                 y = log10(median),
                 shape = Category_OECD,
                 colour = Category_OECD),
             size = 1.5) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_quick
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



### Fraction eliminated in urine and feces ----

FelFa_predicted_df <- TK_df |> 
  filter(Parameter %in% c("Fa", 
                          # "FelTotal",
                          "Ffeces",
                          "FbilEl",
                          "Furine")) |> 
  pivot_wider(names_from = Parameter, values_from = c(median, mean, q25, q75)) |> 
  mutate(median_Fpreintfeces = median_Ffeces-median_FbilEl,
         mean_Fpreintfeces = mean_Ffeces-mean_FbilEl,
         q25_Fpreintfeces = q25_Ffeces-q25_FbilEl,
         q75_Fpreintfeces = q75_Ffeces-q75_FbilEl,
         median_cumFpreintfeces = median_Fpreintfeces + median_Furine,
         mean_cumFpreintfeces = mean_Fpreintfeces + mean_Furine,
         q25_cumFpreintfeces = q25_Fpreintfeces + q25_Furine,
         q75_cumFpreintfeces = q75_Fpreintfeces + q75_Furine,
         median_cumFbile = median_Fpreintfeces + median_Furine + median_FbilEl,
         mean_cumFbile = mean_Fpreintfeces + mean_Furine + mean_FbilEl,
         q25_cumFbile = q25_Fpreintfeces + q25_Furine + q25_FbilEl,
         q75_cumFbile = q75_Fpreintfeces + q75_Furine + q75_FbilEl,
  ) |>  
  pivot_longer(median_Fa:q75_cumFbile, names_to = c("What", "Parameter"), values_to = "Value", names_sep = "_") |> 
  pivot_wider(names_from = "What", values_from = "Value") |> 
  filter(Parameter != "Ffeces") |> 
  mutate(Direction = ifelse(Parameter == "Fa", "In",
                            ifelse(Parameter == "Furine", "Urine",
                                   "Feces"))) 

FelFa_predicted_df <- full_join(FelFa_predicted_df, nC_MW_df)

Clean_FelFa <- FelFa_predicted_df |> 
  filter(Parameter %in% c("FbilEl", "Fpreintfeces", "Furine"))
Cum_FelFa <- FelFa_predicted_df |> 
  filter(Parameter %in% c("cumFbile", "cumFpreintfeces", "Furine")) |> 
  mutate(Parameter = str_replace_all(Parameter, "cum", ""),
         Parameter = str_replace_all(Parameter, "Fbile", "FbilEl"))

Fa_predicted_df <- TK_df |> 
  filter(Parameter %in% c("Fa")) 

Fa_predicted_df <- full_join(Fa_predicted_df, nC_MW_df)

Fa_plot <- ggplot()+
  geom_boxplot(data = Fa_predicted_df,
               aes(x = reorder(PFAS, MW), y = log10(median),
                   ymin = log10(q25), lower = log10(q25), 
                   middle = log10(median), 
                   upper = log10(q75), ymax = log10(q75),
                   fill = Category_OECD, colour = Category_OECD),
               stat = "identity",
               linewidth = 0.3,
               width = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.7) +
  scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
  theme_quick



Fe_plot <- ggplot() +
  geom_bar(data = Clean_FelFa,
           aes(x = reorder(PFAS, MW), y = median, 
               fill = Parameter),
           stat = "identity",
           position = "stack",
           width = 0.5) +
  scale_fill_manual(name = "Parameter",
                    values = c("FbilEl" = "#B198B9",
                               "Fpreintfeces" = "#8C6897",
                               "Furine" = "#89D499"),
                    labels = c("FbilEl" = "Feces: bile",
                               "Fpreintfeces" = "Feces: pre-systemic",
                               "Furine" = "Urine")) +
  geom_errorbar(data = Cum_FelFa,
                aes(x = reorder(PFAS, MW),
                    ymin = q25, ymax = q75,
                    colour = Parameter),
                position = "dodge",
                width = .3,
                linewidth = .4
  ) +
  scale_color_manual(name = "Parameter",
                     values = c("FbilEl" = "#683975",
                                "Fpreintfeces" = "#440154",
                                "Furine" = "#2A7741"),
                     labels = c("FbilEl" = "Feces: bile",
                                "Fpreintfeces" = "Feces: pre-systemic",
                                "Furine" = "Urine")) +
  labs(x = "", 
       y = "Cumulative fraction excreted") +
  theme_quick +
  theme(legend.position = "right")


Fe_plot <- Fe_plot + Patchwork_theme +
  theme(legend.position = "top",
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.8, "lines"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        axis.text.x = element_text(size = 35, face = "bold", margin = margin(0, 0, 0, 0)),
        axis.text.y = element_text(size = 35, margin = margin(0, 0, 0, 0)),   
        axis.title = element_text(size = 45) #, face = "bold"
  )


ggsave(plot = Fe_plot,
       file = here(OUTPUT,"Fe_plot.png"),
       dpi = 600,
       width = 16,
       height = 12,
       units = "cm")

### Urinary and fecal clearances ----

CLurCLfe_predicted_df <- TK_df |> 
  filter(Parameter %in% c("CLfe", 
                          "CLur")) 

CLurCLfe_predicted_df <- full_join(CLurCLfe_predicted_df, nC_MW_df)

RurRfe_predicted_df <- TK_df |> 
  filter(Parameter %in% c("RfeEl", 
                          "RurEl"))
RurRfe_predicted_df <- full_join(RurRfe_predicted_df, nC_MW_df)

CL_plot <- ggplot() +
  geom_boxplot(data = CLurCLfe_predicted_df,
               aes(x = reorder(PFAS, MW), y = log10(median),
                   ymin = log10(q25), lower = log10(q25), 
                   middle = log10(median), 
                   upper = log10(q75), ymax = log10(q75),
                   fill = Category_OECD, colour = Category_OECD),
               stat = "identity",
               linewidth = 0.3,
               width = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.7) +
  scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
  facet_wrap(~ Parameter, nrow = 2, 
             labeller = as_labeller(c("CLfe" = "Biliary/ fecal",
                                      "CLur" = "Renal")),
             scales = "free_y") +
  theme_quick
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
UA_one <- CP_plot + 
  labs(title = "A.",
       x = "Time (years)",
       y = "Plasma concentration (μg/L)") +
  guides(fill = FALSE) + guides(colour = FALSE) + 
  Patchwork_theme + 
  theme(legend.position = "none")

UA_two <- Fa_plot + 
  labs(title = "B." ,
       x = "",
       y = "Fraction absorbed (log10 scale)") +
  guides(colour = FALSE) + guides(linetype = FALSE) + 
  Patchwork_theme + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90))

UA_three <- CL_plot + 
  labs(title = "C." ,
       x = "",
       y = "Clearance (L/day) (log10 scale)") +
  guides(fill = FALSE) + guides(colour = FALSE) + 
  Patchwork_theme + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90))

UA_four <- HL_predicted_plot + 
  labs(title = "D.",
       x = "",
       y = "Half-life (years) (log10 scale)") +
  Patchwork_theme + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))

UA_plot <- UA_one / (UA_two + UA_three + UA_four) +
  plot_layout(heights = c(0.7, 1),
              guides = "collect")  &
  theme(plot.margin = margin(0.1, 0.1, b = 0, 0.1, "cm"),
        legend.box.margin = margin(t = -0.5, 0, b = 0.1, 0, "cm"),
        legend.position='bottom')

ggsave(plot = UA_plot,
       file = here(OUTPUT,"UA_plot_log.png"),
       dpi = 600,
       width = 24,
       height = 16,
       units = "cm")

### Against HL measured ----


Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |>
  rename(PFAS = Acronym) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

raw_HL_measured <- Measured_df |>
  filter(Parameter == "Halflife",
         PFAS != "PFTrDA") |>
  select(PFAS, Category_OECD, Reference,
         Value, Range_CV_Std_min, Range_CV_Std_max,
         n, sex, Age_years) |>
  mutate(n = round(n,0))

set.seed(1234)

HL_measured <- raw_HL_measured %>%
  mutate(n = as.numeric(n),
         Value = as.numeric(Value),
         min = as.numeric(Range_CV_Std_min),
         max = as.numeric(Range_CV_Std_max)) |>
  mutate(meandistance = (Value-min + max-Value)/2,
         sd = ifelse(meandistance > Value, meandistance/4, meandistance), #Value/4, meandistance/2)
         sd = sd/2 #,
         # sd = ifelse(Reference == "Yao_2023" & PFAS == "PFBA", sd*2, sd),
         # sd = ifelse(Reference == "Yao_2023" & PFAS == "PFPeA", sd*2, sd),
         # sd = ifelse(Reference == "Xu_2020" & PFAS == "PFBS", sd*1.1, sd)
  ) |>
  rowwise() |>
  mutate(r_sex = case_when( n != 1 & sex == "F" ~ list(rep("F", n)),
                            n != 1 & sex == "M" ~ list(rep("M", n)),
                            TRUE ~ list(sex)
  ),
  r_Value = ifelse(n != 1,
                   list(rnorm(n, Value, sd)),
                   list(Value))
  ) |>
  unnest(cols = r_sex:r_Value,
         keep_empty = TRUE) |>
  mutate(Value = round(Value, 3),
         r_Value = round(r_Value, 3),
         min = ifelse(min <= 0, 0.00001, min),
         r_Value = ifelse(n != 1 & r_Value < min, min, r_Value) ,
         r_Value = ifelse(n != 1 & r_Value > max, max, r_Value),
         r_Value = round(r_Value, 3)
  ) |>
  select(!sex) |>
  rename(mean = Value,
         Value = r_Value,
         sex = r_sex)


HL_measured_df <- full_join(HL_measured, nC_MW_df)

HL_MesPred_plot <- 
  ggplot() +
  geom_jitter(data = HL_measured_df,
              aes(reorder(PFAS, MW),
                  log10(Value),
                  shape = Category_OECD),
              color = "grey60",
              width = 0.2,
              alpha = 0.9,
              size = .3) +
  geom_errorbar(data = HL_predicted_df,
                aes(x = reorder(PFAS, MW),
                    y = log10(median),
                    colour = Category_OECD,
                    ymin = log10(q25),
                    ymax = log10(q75)),
                width = .01,
                linewidth = .4) +
  geom_point(data = HL_predicted_df,
             aes(x = reorder(PFAS, MW),
                 y = log10(median),
                 shape = Category_OECD,
                 colour = Category_OECD),
             size = 1.5) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  labs(x = "",
       y = "Half life (years) (log10 scale)") +
  theme_quick

HL_MesPred_plot

HL_MesPred_plot <- HL_MesPred_plot + Patchwork_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0))
  )

ggsave(plot = HL_MesPred_plot,
       file = here(OUTPUT,"HL_MesPred_plot.png"),
       dpi = 600,
       width = 16,
       height = 11,
       units = "cm")

# Post process parameters ----

P_df <- list.files(path = here(OUTPUT),
                   pattern = "\\_final_model_parameters.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .),
                    id = "name"))|> 
  select(!"...2") |> 
  # rename(Iteration = "...2") |> 
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/Admin/Desktop/PhD_Chrysa/CP_L_R/PARC_PFAS_PBK/Output/UA_GoF_22Aug/",
                                "")) |> 
  mutate(PFAS = str_replace_all(PFAS,
                                "_final_model_parameters.csv","")) |>
  select(!name) #|> unite("PFAS_Iteration", c(PFAS,Iteration), sep = "_") 


Visualise_P <- P_df |>  
  # select(!Clear_net_reabsorption_unscaled) |> 
  rename_with(~str_remove(.x, "_")) |> 
  rename_with(~str_remove(.x, "_")) |> 
  group_by(PFAS) |> 
  summarise(across(where(is.numeric), 
                   .fns =
                     list(median = ~median(.,na.rm = TRUE),
                          q25 = ~quantile(., 0.25, na.rm = TRUE),
                          q75 = ~quantile(.,0.75, na.rm = TRUE)))
  ) |>  
  pivot_longer(!PFAS, names_to = c("Parameter", "stat"), 
               names_sep = "_",values_to = "Value") |> 
  pivot_wider(names_from = stat, values_from = Value)

Visualise_P <- full_join(Visualise_P, nC_MW_df)
write.csv(Visualise_P, here(OUTPUT,"All_PFAS_finalParameters.csv"), row.names = FALSE)



SA_SI <-1099.557*64.41978
VI_L <- 0.03848451 * 64.41978      
VL_ic <- 0.03848451 * 64.41978 * 0.835               
VPT <- 0.004630086 * 64.41978 * 0.398


Selected_final_P <- Visualise_P |> 
  filter(Parameter %in% c(
    "fuPTec",
    "fuLec",
    "fuPTL",
    "SFOATP2B1I",
    "SFOATP2B1",
    "SFOATP1B3",
    "SFOATP1B1", 
    "SFNTCP",    
    "SFBSEP",
    "SFOAT1",           
    "SFOAT3",
    "SFOAT4",
    "Clearapparent",
    "ClearOATP2B1I",
    "ClearuptakeI",
    "ClearOATP1B1",
    "ClearOATP1B3",
    "ClearOATP2B1L",
    "ClearNTCP",
    "ClearBSEP",
    "ClearOAT3",
    "ClearOAT1",
    "ClearOAT4"  
  )) |> 
  pivot_wider(names_from = Parameter, values_from = c(median, q25, q75)) |> 
  mutate(
    median__Papp_SI = median_Clearapparent / (SA_SI*1e-3),    
    q25__Papp_SI = q25_Clearapparent / (SA_SI*1e-3),    
    q75__Papp_SI = q75_Clearapparent / (SA_SI*1e-3),    
    
    median__Clear_OATP2B1c = median_ClearOATP2B1I / (median_SFOATP2B1I*VI_L),  
    q25__Clear_OATP2B1c = q25_ClearOATP2B1I / (q25_SFOATP2B1I*VI_L),  
    q75__Clear_OATP2B1c = q75_ClearOATP2B1I / (q75_SFOATP2B1I*VI_L),  
    
    median__Clear_OATP1B1c = median_ClearOATP1B1 / (median_SFOATP1B1*VL_ic*median_fuLec),
    q25__Clear_OATP1B1c = q25_ClearOATP1B1 / (q25_SFOATP1B1*VL_ic*q25_fuLec),
    q75__Clear_OATP1B1c = q75_ClearOATP1B1 / (q75_SFOATP1B1*VL_ic*q75_fuLec),
    
    median__Clear_OATP1B3c = median_ClearOATP1B3 / (median_SFOATP1B3*VL_ic*median_fuLec),
    q25__Clear_OATP1B3c = q25_ClearOATP1B3 / (q25_SFOATP1B3*VL_ic*q25_fuLec),
    q75__Clear_OATP1B3c = q75_ClearOATP1B3 / (q75_SFOATP1B3*VL_ic*q75_fuLec),
    
    median__Clear_NTCPc = median_ClearNTCP / (median_SFNTCP*VL_ic*median_fuLec),   
    q25__Clear_NTCPc = q25_ClearNTCP / (q25_SFNTCP*VL_ic*q25_fuLec),    
    q75__Clear_NTCPc = q75_ClearNTCP / (q75_SFNTCP*VL_ic*q75_fuLec),    
    
    median__Clear_BSEPc = median_ClearBSEP / (median_SFBSEP*VL_ic),
    q25__Clear_BSEPc = q25_ClearBSEP / (q25_SFBSEP*VL_ic),
    q75__Clear_BSEPc = q75_ClearBSEP / (q75_SFBSEP*VL_ic),
    
    median__Clear_OAT3c = median_ClearOAT3 / (median_SFOAT3*VPT*median_fuPTec),
    q25__Clear_OAT3c = q25_ClearOAT3 / (q25_SFOAT3*VPT*q25_fuPTec),
    q75__Clear_OAT3c = q75_ClearOAT3 / (q75_SFOAT3*VPT*q75_fuPTec),
    
    median__Clear_OAT1c = median_ClearOAT1 / (median_SFOAT1*VPT*median_fuPTec),
    q25__Clear_OAT1c = q25_ClearOAT1 / (q25_SFOAT1*VPT*q25_fuPTec),
    q75__Clear_OAT1c = q75_ClearOAT1 / (q75_SFOAT1*VPT*q75_fuPTec),
    
    median__Clear_OAT4c = median_ClearOAT4 / (median_SFOAT4*VPT*median_fuPTL),
    q25__Clear_OAT4c = q25_ClearOAT4 / (q25_SFOAT4*VPT*q25_fuPTL),
    q75__Clear_OAT4c = q75_ClearOAT4 / (q75_SFOAT4*VPT*q75_fuPTL)
  ) |> 
  pivot_longer(!c(PFAS, Category_OECD, MW, n_C),
               names_to = "Parameter", 
               values_to = "Value") |> 
  separate(Parameter, into = c("stat", "Parameter"), sep = "__") |> 
  separate(stat, into = c("stat1", "Parameter2"), sep = "_") |> 
  mutate(Parameter = ifelse(!is.na(Parameter), Parameter, Parameter2)) |> 
  select(!Parameter2) |> 
  pivot_wider(names_from = stat1, values_from = Value) |> 
  mutate(Parameter = str_replace_all(Parameter,
    c("SFOATP2B1I"="SF_OATP2B1_I",
      "SFOATP2B1"="SF_OATP2B1" ,
      "SFOATP1B3"= "SF_OATP1B3",
      "SFOATP1B1"="SF_OATP1B1",
      "SFNTCP"="SF_NTCP",    
      "SFBSEP"="SF_BSEP",
      "SFOAT1"="SF_OAT1",           
      "SFOAT3"="SF_OAT3",
      "SFOAT4"="SF_OAT4",
      "fuPTL"="fu_PTL",
      "fuLec"="fu_L_ec",
      "fuPTec"="fu_PT_ec"   
    ))) |> 
  filter_out(Parameter %in% c(
    "Clearapparent", "ClearOATP2B1I",
    "ClearuptakeI",  "ClearOATP1B1",  
    "ClearOATP1B3",  "ClearOATP2B1L", 
    "ClearNTCP",     "ClearBSEP",     
    "ClearOAT3",        "ClearOAT1",     "ClearOAT4")
    )
  
  write.csv(Selected_final_P, here(OUTPUT, "Calibrated_P.csv"), row.names = FALSE)

Parameters_plot <- ggplot() +
  geom_pointrange(data = Visualise_P,
                  aes(x = reorder(PFAS, MW),
                      y = median,
                      colour = Category_OECD,
                      ymin = q25,
                      ymax = q75),
                  size = 0.1,
                  # width = .01,
                  linewidth = .4) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  # facet_wrap(~Parameter, scales = "free_y")+
  facet_wrap(~factor(Parameter,
                     c( "SFOATP2B1I",    "SFOATP2B1",       "SFOATP1B3",     "SFOATP1B1",
                        "SFNTCP",        "SFBSEP",          "SFOAT1",        "SFOAT3",
                        "SFOAT4",
                        "PI",            "PL",              "PK",            "PA",        "PR",
                        "tco",           "Clearapparent",   "ClearOATP2B1I",   "ClearuptakeI",
                        "fuLec", "ClearOATP1B1",  "ClearOATP1B3",    "ClearOATP2B1L", "ClearNTCP", "ClearuptakeL",
                        "ClearBSEP",     "ClearexcretionL",
                        "fup", "GFR",           "fupGFR",   "fuPTec", "ClearOAT3",     "ClearOAT1",  "ClearOAT4",
                        "fuPTL", "Clearreabsorption", "QT")
  ), scales = "free_y") +
  labs(x = "", 
       y = "") +
  theme_quick +
  theme(axis.text.x = element_text(angle = 90))

Parameters_plot = Parameters_plot + Patchwork_theme +
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = Parameters_plot,
       file = here(OUTPUT, "All_PFAS_calculated_model_parameters_plot.png"),
       dpi = 600,
       width = 24,
       height = 18,
       units = "cm")

SelectedP_df <- P_df |> 
  select(PFAS,Clear_uptake_I, Clear_uptake_L, Clear_excretion_L, fup_GFR, Clear_reabsorption) |> 
  rename(UI = Clear_uptake_I, UL = Clear_uptake_L, EL = Clear_excretion_L, fGFR = fup_GFR, RK = Clear_reabsorption) |> 
  group_by(PFAS) |> 
  summarise(across(where(is.numeric), .fns =
                     list(median = ~median(.,na.rm = TRUE),
                          mean = ~mean(.,na.rm = TRUE),
                          q25 = ~quantile(., 0.25, na.rm = TRUE),
                          q75 = ~quantile(.,0.75, na.rm = TRUE)))) |> 
  pivot_longer(!PFAS, names_to = c("Parameter", "stat"), 
               names_sep = "_",values_to = "Value") |> 
  mutate(Parameter = str_replace_all(Parameter, c("UI"= "Clear_uptake_I", 
                                                  "UL"= "Clear_uptake_L", 
                                                  "EL"= "Clear_excretion_L", 
                                                  "fGFR"= "fup_GFR", 
                                                  "RK"= "Clear_reabsorption"
  ))) |> 
  pivot_wider(names_from = stat, values_from = Value)

SelectedP_df <- full_join(SelectedP_df, nC_MW_df)

Key_params_plot <- ggplot() +
  geom_pointrange(data = SelectedP_df,
                  aes(x = reorder(PFAS, MW),
                      y = median,
                      colour = Category_OECD,
                      ymin = q25,
                      ymax = q75),
                  size = 0.1,
                  # width = .01,
                  linewidth = .4) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  facet_wrap(~Parameter, scales = "free_y") +
  labs(x = "", 
       y = "Clearance (L/day)") +
  theme_quick +
  theme(axis.text.x = element_text(angle = 90))


Key_params_plot <- Key_params_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = Key_params_plot,
       file = here(OUTPUT, "All_PFAS_key_params_plot.png"),
       dpi = 600,
       width = 16,
       height = 12,
       units = "cm")


### Against measured clearances ----

CL_UrFe_measured_df <- Measured_df |>
  filter(Parameter %in% c("CL_urinary", 
                          "CL_fecal")) |> 
  # filter_out(Reference == "Abraham_2024") |>
  mutate(Parameter = str_replace_all(Parameter, "CL_fecal", "CLfe"),
         Parameter = str_replace_all(Parameter, "CL_urinary", "CLur")
         )

CL_UrFe_measured_df <- full_join(CL_UrFe_measured_df, nC_MW_df) |>
  select(PFAS, Category_OECD,
         Parameter,Value, Unit,
         MW, n_C,
         Reference) |>
  na.omit()

CL_MesPred_plot <- ggplot() +
  geom_point(data = CL_UrFe_measured_df,
             aes(x = reorder(PFAS, MW), y = log10(Value),
                 shape = Reference), #,             color = Category_OECD), 
             color = "grey60") +
  geom_pointrange(data = CLurCLfe_predicted_df,
               aes(x = reorder(PFAS, MW), y = log10(median),
                   ymin = log10(q25), #lower = log10(q25),
                   # middle = log10(median), upper = log10(q75), 
                   ymax = log10(q75),
                   fill = Category_OECD, colour = Category_OECD),
               # stat = "identity",
               linewidth = 0.3,
               width = .01,
               size = .3) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.7) +
  scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
  facet_wrap(~ Parameter, nrow = 2,
             labeller = as_labeller(c("CLfe" = "Biliary/ fecal",
                                      "CLur" = "Renal")),
             scales = "free_y") +
  labs(x = "",y = "Clearance (L/day) (log10 scale)")+
  theme_quick


CL_MesPred_plot

CL_MesPred_plot <- CL_MesPred_plot + Patchwork_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0))
  )

ggsave(plot = CL_MesPred_plot,
       file = here(OUTPUT,"CL_MesPred_plot.png"),
       dpi = 600,
       width = 16,
       height = 11,
       units = "cm")


### Against plasma measured ----

Measured_Abraham_df <- Measured_df |>
  filter(Reference == "Abraham_2024")

C_plasma_measured_df <- Measured_df |>
  filter(Parameter == "C_plasma") |>
  mutate(Time = as.numeric(Time_days)) |>
  select(PFAS, Time, Value)
C_plasma_measured_df <- full_join(C_plasma_measured_df, nC_MW_df)

C_plasma_predicted_df <- C_df |>
  filter(Tissue == "Plasma",
         Time < 450)

C_plasma_predicted_df <- full_join(C_plasma_predicted_df, nC_MW_df)


ObsVsPred_Plasma_plot <- ggplot() +
  geom_ribbon(data = C_plasma_predicted_df,
              aes(x = Time/365,
                  ymin = q25,
                  ymax = q75,
                  fill = Category_OECD),
              alpha = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
  geom_path(data = C_plasma_predicted_df,
            aes(Time/365, median,
                colour = Category_OECD)) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  geom_point(data = C_plasma_measured_df, #filter(C_plasma_measured_df, Category_OECD == "PFCAs"),
             aes(Time/365, Value),
             colour = "black",
             size = 0.1) +
  labs(x = "Time (years)",
       y = "Concentration (μg/ml)") +
  facet_wrap(~ factor(PFAS, c( #factor(n_C),
    "PFBA",
    "PFPeA",
    "PFBS",
    "PFHxA",
    "HFPO-DA",
    "PFPeS",
    "PFHpA",
    "ADONA",
    "PFHxS",
    "PFOA",
    "PFHpS",
    "PFNA",
    "PFOS",
    "PFDA",
    "PFUnDA",
    "PFDoDA",
    "PFTrDA")
  ),
  nrow = 2) +
  theme_quick


ObsVsPred_Plasma_plot <- ObsVsPred_Plasma_plot + Patchwork_theme +
  theme(legend.position = "bottom")

ggsave(plot = ObsVsPred_Plasma_plot,
       file = here(OUTPUT, "CP_plot_UA.png"),
       dpi = 600,
       width = 24,
       height = 11,
       units = "cm")



























