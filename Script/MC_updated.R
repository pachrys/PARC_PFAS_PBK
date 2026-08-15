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
OUTPUT <- here("Output", "trial_MC") 
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

Dose_df <- tibble(
  PFAS = PFAS_list,
  Dose = 4) |>
  mutate(Dose = ifelse(PFAS == "PFBS", 18.4, 
                       ifelse(PFAS == "HFPO-DA" , 19.6,
                              Dose)))

# Dose_df <- Abraham_df
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
  GFR = 156.6414, 
  QTc = 0.36, 
  tco = 3.36,
  Qmenses = 0.0008561644)


## MC parameters
MC_TK_params <- read_xlsx(here("Input", "MC_TK_Parameters.xlsx")) |>
  mutate(Predicted = ifelse(Specific == "calculated", "yes", Predicted))|>
  filter(Predicted == "yes") |>
  select(PFAS, Parameter, Value, SD, Distribution)
# mutate(SD = ifelse(Distribution == "lognormal", SD-0.2, SD)) #SD = ifelse((Distribution == "uniform") & (Value - SD <=0), Value - (Value/2), SD),

## MC parameters
MC_TK_params <- read_xlsx(here("Input", "MC_TK_Parameters.xlsx")) |>
  mutate(Predicted = ifelse(Specific == "calculated", "yes", Predicted))|>
  filter(Predicted == "yes") |>
  select(PFAS, Parameter, Value, SD, Distribution)

# load functions
source(here("Script", "functions_NG_PBK_model.R"))


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
                   times = seq(0, 3650, by = 0.4), # in days 0.04
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
  message(glue("{PFAS_i}------------------------------------"))
  
  # physiological parameters; they are kept fix
  parm.phys <- parm.phys
  Dose_ug_day <- filter(Dose_df, PFAS == PFAS_i) |> 
    select(Dose) |> pull()
  
  PhysChemTK_params <- MC_TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(!PFAS) 
  
  parm.chem <- PhysChemTK_params |> 
    mutate(meanlog = ifelse(Distribution == "lognormal", Value, NA),
           sdlog = ifelse(Distribution == "lognormal", SD, NA),
           min = ifelse(Distribution == "uniform", Value - SD, NA),
           max = ifelse(Distribution == "uniform", Value + SD, NA),
           Distribution = ifelse(Distribution == "lognormal", "qlnorm", "qunif")) |> 
    mutate(min = ifelse(min <= 0, 1e-10, min))
  
  MC.parms <- parm.chem$Parameter
  MC.factors <- length(MC.parms) 
  
  ## Latin Hypercube Sampling ---- 
  # this optimises the Monte Carlo sampling, also used by McNally & Loizou
  
  set.seed(10010) # set seeds for reproducibility
  n = 10 #1000 # number of samples per parameter
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
    labs(title = paste0(PFAS_i,"Monte Carlo parameter distributions"),
         y = "Frequency",
         x = "Value (different unit per parameter)") +
    theme_quick 
  Design_plot
  
  Design_plot <- Design_plot + CP_theme +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(plot = Design_plot,
         file = here(OUTPUT, paste0(PFAS_i,"_LHS_MC_ExpDesign.png")),
         dpi = 600,
         width = 18,
         height = 24,
         units = "cm")
  
  write.csv(MC.design, 
            file = here(OUTPUT, paste0(PFAS_i,"MC_ExpDesign.csv")), 
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
  save(savey, file = here(OUTPUT, paste0(PFAS_i,"_Monte_Carlo_y.RData")))

  Monte_Carlo_PFAS_savey[[i]] <- savey
  }
  
names(Monte_Carlo_PFAS_savey) <- PFAS_list

Save_MC_PBK_results <- rbindlist(Monte_Carlo_PFAS_savey, idcol = "PFAS")



























