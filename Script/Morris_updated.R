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
library(sensitivity)
library(ggrepel)


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

Patchwork_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 50, face = "bold", margin = margin(0, 0, b = 0.5, 0)),
    plot.subtitle = element_text(hjust = 0, size = 35, margin = margin(0, 0, b = 0.5, 0)),
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


# Set output storage directory
OUTPUT <- here("Output", "26August_Run_Morris") 
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
SA_Physio_params <- read_xlsx(here("Input", "Physio_Parameters.xlsx")) 
# 
# parm.phys <- c(
#   BW = 64.41978,
#   QC = 5486.301,
#   Hct = 0.3645031,
#   VIc = 0.016,
#   VILc = 0.03848451,
#   SA_SIc = 1099.557,
#   VLc = 0.02458126,
#   VL_icc = 0.835,
#   VL_ecc = 0.161,
#   VKc = 0.004630086,
#   VPTc = 0.398,
#   VPTTc = 0.5470631,
#   VPTLc = 0.403,
#   VRKLc = 0.283,
#   VAc = 0.3426787,
#   MAc = 1.487552,
#   VPc = 0.065,
#   QIc = 0.165,
#   QLc = 0.06857434,
#   QKc = 0.1761446,
#   QAc = 0.09413654,
#   QUrc = 0.022,
#   GFR = 193.2765, # 144 - 288 L/day (=0.2*0.1761446*5486.301)
#   QTc = 0.36,
#   tco = 3.36,
#   Qmenses = 0.0008561644,
#   Calb_P = 40 # 35-50 g/L plasma
#   )




## MC parameters
SA_TK_params <- read_xlsx(here("Input", "TK_Parameters.xlsx")) |>
  filter(Predicted == "yes") 

# load functions
source(here("Script", "functions_NG_PBK_model.R"))


# Function to run the morris design
run_Morris <- function(parms){
  
  Dose_ug_day = 0.1
  
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
  
  
  pars <- GSA_Calc_final_params(parms = parms)
  
  # save parameters to check when analysing results
  run_pars <- as.data.frame(pars)
  
  PBK.out <- lsoda(y = A_init,
                   times = seq(0, 365*3, by = 1), 
                   func = PBK.model,
                   parms = pars,
                   rtol = 1e-4,        
                   atol = 1e-6,
                   hmax = 0.5,         
                   maxsteps = 100000)
  
  PBK_OUTPUT <- as.data.frame(PBK.out) 
  PBK_OUTPUT[,-1] <- pmax(PBK_OUTPUT[,-1], 0)

  data <- PBK_OUTPUT |> drop_na() 
  
  # AUC ----
  AUC <- trapz(data$time, data$CP) # ng * day/ml
  
  # Cmax, Tmax, Tlast ----
  Cmax <- max(data$CP) # ng/ml
  Tmax <- filter(data, CP == max(data$CP)) %>% pull(time) # days
  Tlast <- filter(data, CP > 1e-9) |> select(time) |> last() |> pull() # 1e-9 ng/ml is already a very small number
  
  # calculate half life ----
  elimination_phase <- data |>
    filter(time >= Tmax & time <= Tlast) |>
    mutate(logCP = log(CP))
  
  # linear regression to find the slope
  elimination_model <- lm(logCP ~ time, data = elimination_phase)
  
  # calculate elimination constant (i.e. slope of the line)
  kel <- -elimination_model$coefficients[[2]]
  
  # calculate the half-life
  half_life <- (log(2) / kel) / 365 # unit correction to years
  
  return( c(Cmax=Cmax, 
            half_life = half_life, 
            AUC = AUC))
}



for (i in 1:nPFAS) {
  
  
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}------------------------------------"))
  
  # physiological parameters; they are kept fix
  Physio_params <- tibble(
    Parameter = SA_Physio_params$Abbreviation,
    Value = SA_Physio_params$Value,
    binf = SA_Physio_params$Binf,
    bsup = SA_Physio_params$Bsup
  )

  Dose_ug_day <- filter(Dose_df, PFAS == PFAS_i) |>
    select(Dose) |> pull()

  PhysChemTK_params <- SA_TK_params |>
    filter(PFAS == PFAS_i) |>
    mutate(Value = ifelse(PBK_mechanism == "Distribution", exp(Value),
                          ifelse(PBK_mechanism == "Absorption", exp(Value),
                                 Value)
                          ),
           binf = min,
           bsup = max
           ) |>
    select(Parameter, Value, binf, bsup)

  Model_params <- rbind(Physio_params, PhysChemTK_params)

  Morris_parms <- Model_params$Parameter
  binf <- Model_params$binf
  bsup <- Model_params$bsup
  Morris_factors <- length(Morris_parms)

  ## Morris experiment design ----

  set.seed(10010) # set seeds for reproducibility

  # Experiment design
  Morris_design <- list(type = "oat", levels = 6, grid.jump = 3)

  # Repetitions
  Morris_r <- 100
    
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
    labs(title = paste0(PFAS_i," Morris test parameter distributions"),
         y = "Frequency",
         x = "Value (different unit per parameter)") +
    theme_quick
  Design_plot

  Design_plot <- Design_plot + Patchwork_theme +
    theme(axis.text.x = element_text(angle = 0))

  ggsave(plot = Design_plot,
         file = here(OUTPUT, paste0(PFAS_i,"_Morris_ExpDesign.png")),
         dpi = 600,
         width = 34,
         height = 24,
         units = "cm")

  write.csv(Morris_design,
            file = here(OUTPUT, paste0(PFAS_i,"_Morris_ExpDesign.csv")),
            row.names = FALSE)

  results <- apply(Morris_design, 1,  run_Morris)

  y <- results
  savey = y
  # savey <- mapply(cbind, y,
  #                 "Iteration" = seq(1,length(y), by = 1),
  #                 "PFAS" = PFAS_i,
  #                 SIMPLIFY = FALSE)
  save(savey, file = here(OUTPUT, paste0(PFAS_i,"_Morris_y.RData")))

  # load(here(OUTPUT, paste0(PFAS_i,"_Morris_y.RData")))
  # y <- savey
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
  plot(Experience, main = paste0(PFAS_i))
  print(PFAS_i)
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
  rownames(IGlobal) <- Morris_parms
  IGlobal <- IGlobal %>% mutate(parm_names = rownames(IGlobal),
                                PFAS = PFAS_i)
  
  write.csv(IGlobal, here(OUTPUT,paste0(PFAS_i, "_Morris_indices.csv")), row.names = FALSE)
  
  }


for (i in 1:nPFAS) {
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}------------------------------------"))
  
  data <- read.csv(here(OUTPUT,paste0(PFAS_i, "_Morris_indices.csv")))
  # load(here(OUTPUT,paste0(PFAS_i, "_Morris_indices.RData")))
  morris_data <- data |> #Morris_indices$IGlobal |> 
    arrange(desc(IGlobal_mu.star)) %>% 
    mutate(n = seq(from=1, to=50, by=1)) 
  write.csv(morris_data, here(OUTPUT,paste0(PFAS_i, "_arranged_Morris_indices.csv")), row.names = FALSE)
  
  # 
  # morris_plot <- morris_data %>% filter(IGlobal_mu.star > 1) %>% 
  #   ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
  #              label = parm_names)) +
  #   geom_jitter(width = .5, alpha = 0.7, size = 1)  +
  #   geom_text_repel(aes(label = as.character(parm_names)), 
  #                   size=15, 
  #                   max.overlaps = Inf,
  #                   # min.segment.length = 0.05,
  #                   segment.size = 0.1
  #   ) +
  #   labs(
  #     title = paste0(PFAS_i, " Morris elementary effects sensitivity analysis"),
  #     subtitle = expression(paste("Parameters with μ* >= 1")),
  #     x = expression(paste("Mean Effect (μ*)")),
  #     y = expression(paste("Standard Deviation (σ)"))
  #   ) +
  #   theme_quick 
  # morris_plot
  # 
  # morris_plot = morris_plot +
  #   Patchwork_theme + 
  #   theme(axis.text.x = element_text(angle = 0))
  # 
  # ggsave(plot = morris_plot,
  #        file = here(OUTPUT, paste0(PFAS_i, "_morris_plot.png")),
  #        dpi = 600,
  #        width = 12,
  #        height = 12,
  #        units = "cm")
  
}


Indices_df <-
  list.files(path = here(OUTPUT),
             pattern = "\\_arranged_Morris_indices.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .)))


Wide_Indices_n_df <- Indices_df |> 
  select(!c(IGlobal_sigma, IGlobal_mu.star)) |> 
  mutate(parm_names = str_replace_all(parm_names,"_", "--")) |> 
  pivot_wider(names_from = parm_names, values_from = n) |> 
  summarise(across(where(is.numeric), .fns =
                     list(mean = ~mean(.,na.rm = TRUE),
                          min = ~min(.,na.rm = TRUE),
                          max = ~max(.,na.rm = TRUE)
                     ))) |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value") |> 
  separate(Variable,
           into = c("parm_names", "Variable"), 
           sep = "_") |> 
  mutate(parm_names = str_replace_all(parm_names, "--", "_")) |> 
  pivot_wider(names_from = Variable, values_from = Value) |> 
  rename(n_mean = mean,
         n_min = min,
         n_max = max)


Wide_Indices_IGlobal_mu.star_df <- Indices_df |> 
  select(!c(IGlobal_sigma, n)) |> 
  mutate(parm_names = str_replace_all(parm_names,"_", "--")) |> 
  pivot_wider(names_from = parm_names, values_from = IGlobal_mu.star) |> 
  summarise(across(where(is.numeric), .fns =
                     list(mean = ~mean(.,na.rm = TRUE),
                          min = ~min(.,na.rm = TRUE),
                          max = ~max(.,na.rm = TRUE)
                     ))) |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value") |> 
  separate(Variable,
           into = c("parm_names", "Variable"), 
           sep = "_") |> 
  mutate(parm_names = str_replace_all(parm_names, "--", "_")) |> 
  pivot_wider(names_from = Variable, values_from = Value) |> 
  rename(mu_mean = mean,
         mu_min = min,
         mu_max = max)


Wide_Indices_IGlobal_sigma_df <- Indices_df |> 
  select(!c(IGlobal_mu.star, n)) |> 
  mutate(parm_names = str_replace_all(parm_names,"_", "--")) |> 
  pivot_wider(names_from = parm_names, values_from = IGlobal_sigma) |> 
  summarise(across(where(is.numeric), .fns =
                     list(mean = ~mean(.,na.rm = TRUE),
                          min = ~min(.,na.rm = TRUE),
                          max = ~max(.,na.rm = TRUE)
                     ))) |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value") |> 
  separate(Variable,
           into = c("parm_names", "Variable"), 
           sep = "_") |> 
  mutate(parm_names = str_replace_all(parm_names, "--", "_")) |> 
  pivot_wider(names_from = Variable, values_from = Value) |> 
  rename(sd_mean = mean,
         sd_min = min,
         sd_max = max)


Wide_Indices_df <- full_join(Wide_Indices_n_df, Wide_Indices_IGlobal_mu.star_df)
Wide_Indices_df <- full_join(Wide_Indices_df, Wide_Indices_IGlobal_sigma_df)
write.csv(Wide_Indices_df, here(OUTPUT,"Processed_indices.csv"), row.names = FALSE)
