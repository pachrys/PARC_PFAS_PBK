# --------------------------------------------------------------------------- #
# Script for running the MCMC using eFAST
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
    
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 30, margin = margin(0, 0, 0, 0)),
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
    
    strip.text.x = element_text(size = 30, face = "bold", margin = margin(0,0,0,0))
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
OUTPUT <- here("Output", "SA_Morris") 
dir.create(OUTPUT, recursive = TRUE)

# Load data

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
  Dose = 0.1) 

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
SA_TK_params <- read_xlsx(here("Input", "SA_TK_Parameters.xlsx")) |>
  mutate(Predicted = ifelse(Specific == "calculated", "yes", Predicted))|>
  filter(Predicted == "yes") |>
  select(PFAS, Parameter, Value, SD, Distribution)

# load functions
source(here("Script", "functions_NG_PBK_model.R"))


# Function to run the morris design
run_Morris <- function(parms, Dose_ug_day){
  
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
  
  
  pars <- Morris_Calc_final_params(parms = parms)
  
  # save parameters to check when analysing results
  run_pars <- as.data.frame(pars)
  
  PBK.out <- lsoda(y = A_init,
                   times = seq(0, 3650, by = 1), 
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

Morris_PFAS_design <- list()
Morris_PFAS_RESULTS <- list()
Morris_PFAS_y <- list()
Morris_PFAS_savey <- list()#

for (i in 1:nPFAS) {
  
  
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}------------------------------------"))
  
  # physiological parameters; they are kept fix
  PhysTK_params <- tibble(
    Parameter = names(parm.phys),
    Value = parm.phys,
    binf = Value - (0.1*Value),
    bsup = Value + (0.1*Value)
  )
  
  Dose_ug_day <- filter(Dose_df, PFAS == PFAS_i) |> 
    select(Dose) |> pull()
  
  PhysChemTK_params <- SA_TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(!PFAS) |> 
    transmute(Parameter = Parameter,
              Value = Value,
              binf = ifelse(Distribution == "lognormal", 
                            exp(Value) - (0.1*exp(Value)), 
                            Value - SD),
              bsup = ifelse(Distribution == "lognormal", 
                            exp(Value) + (0.1*exp(Value)), 
                            Value + SD)) |> 
    mutate(binf = ifelse(binf <= 0, 1e-10, binf))
  
  Model_params <- rbind(PhysTK_params, PhysChemTK_params)
  
  Morris_parms <- Model_params$Parameter
  binf <- Model_params$binf
  bsup <- Model_params$bsup
  Morris_factors <- length(Morris_parms) 
  
  ## Morris experiment design ---- 

  set.seed(10010) # set seeds for reproducibility
  n = 1000 # number of samples per parameter
  
  # Experiment design
  Morris_design <- list(type = "oat", levels = 6, grid.jump = 3)
  
  # Repetitions
  Morris_r <- 5
  
  # Morris test function
  Experience <- morris(model =  NULL, 
                       factors = Morris_factors, # number of parameters
                       r = Morris_r,  # number of repetitions
                       design = Morris_design,  
                       binf = binf,  
                       bsup = bsup, 
                       scale = TRUE)
  

  # create the final MC design and remove infinite values if any
  Morris_design <- Experience$X |> 
    as.data.frame() |> 
    mutate(across(everything(), ~na_if(., Inf))) |> 
    drop_na()
  colnames(Morris_design) <- Morris_parms

  # visualise and save design
  Design_plot <- Morris_design |>  
    pivot_longer(everything(), 
                         names_to = "Parameter", 
                         values_to = "Value") |> 
    ggplot(aes(x=Value))+
    geom_histogram()+
    facet_wrap(~Parameter, scales = "free_y")+
    coord_flip()+
    labs(title = paste0(PFAS_i," Morris parameter distributions"),
         y = "Frequency",
         x = "Value (different unit per parameter)") +
    theme_quick 
  Design_plot
  
  Design_plot <- Design_plot + CP_theme +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(plot = Design_plot,
         file = here(OUTPUT, paste0(PFAS_i,"_Morris_ExpDesign.png")),
         dpi = 600,
         width = 22,
         height = 18,
         units = "cm")
  
  write.csv(Morris_design, 
            file = here(OUTPUT, paste0(PFAS_i,"_Morris_ExpDesign.csv")), 
            row.names = FALSE)
  
  results <- apply(Morris_design, 1, function(row) {
    run_Morris(row,
               Dose_ug_day = Dose_ug_day)
  })

  y <- results
  savey <- mapply(cbind, y,
                  "Iteration" = seq(1,length(y), by = 1),
                  "PFAS" = PFAS_i,
                  SIMPLIFY = FALSE)
  save(savey, file = here(OUTPUT, paste0(PFAS_i,"_Morris_y.RData")))

  Morris_PFAS_savey[[i]] <- savey
  
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
  
  Morris_indices <- list(
    mu = mu,
    mu_star = mu.star,
    sigma = sigma,
    mu_star_scale = mu.star_scale,
    sigma_scale = sigma_scale,
    IGlobal = IGlobal
  )
  print(IGlobal)
  
  write.csv(Morris_indices, here(OUTPUT,paste0(PFAS, "_Morris_indices.csv")), row.names = FALSE)
  
  }
  
names(Morris_PFAS_savey) <= names(Morris_PFAS_RESULTS) <- PFAS_list

Save_Morris_PBK_results <- rbindlist(Morris_PFAS_savey, 
                                     idcol = "PFAS")
Save_Morris_MorrisSA_results <- rbindlist(Morris_PFAS_RESULTS, 
                                     idcol = "PFAS")


























