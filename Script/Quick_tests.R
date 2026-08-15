# --------------------------------------------------------------------------- #
# Script for running the PFAS-PBK model
# Chrysanthi Pachoulide, June 2026
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment


# Load packages
library(here)
library(ggplot2)
library(deSolve)
library(tidyverse)
library(sensitivity)
library(glue)
library(readxl)
library(PKNCA)
library(parallel)
library(pracma)
library(viridis)
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
    
    # strip.text.x = element_blank()
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
WhatAmITesting <- "FINAL_BSEPoo1"
OUTPUT <- here("Output", WhatAmITesting) 
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

## PFAS

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

# load functions
source(here("Script", "functions_NG_PBK_model.R"))

# Create to collect all inputs and functions needed to run the model
run_PBK <- function(parm.chem, parm.phys, Dose_ug_day){
  
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
                   times = seq(0, 365*10, by = 0.4), # in days 0.04
                   func = PBK.model,
                   parms = pars,
                   rtol = 1e-4,        # Looser tolerance
                   atol = 1e-6,
                   hmax = 0.5,         # Smaller max step
                   maxsteps = 100000)

  PBK_OUTPUT <- as.data.frame(PBK.out) 
  PBK_OUTPUT[,-1] <- pmax(PBK_OUTPUT[,-1], 0)
  # %>%
  #   mutate(CP = ifelse(CP < 1e-9, 0, CP))

  PBK_OUTPUT <- cbind(PBK_OUTPUT, run_pars)
  
  return(PBK_OUTPUT)
  # return(run_pars) just to check parameters

}

PFAS_PBK_model_results <- list()

for (i in 1:nPFAS) {
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}------------------------------------"))
  
  # physiological parameters; they are kept fix
  parm.phys <- parm.phys 
  Dose_ug_day <- filter(Dose_df, PFAS == PFAS_i) |> 
    select(Dose) |> pull()
  
  PhysChemTK_params <- MC_TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(!PFAS) |> 
    mutate(Value = ifelse(Distribution == "lognormal", exp(Value), Value)) 
  
  # if(PFAS_i %in% c("PFBA", "PFPeA", "PFBS", "PFHxA", "HFPO-DA")){
    # PhysChemTK_params <- PhysChemTK_params |>
    # mutate(Value = ifelse(Parameter == "Papp_SI", Value+exp(SD), Value))
# }
  
  parm.chem <- PhysChemTK_params$Value
  names(parm.chem) <- PhysChemTK_params$Parameter
  
  output <- run_PBK(parm.chem = parm.chem, 
                     parm.phys = parm.phys,
                     Dose_ug_day = Dose_ug_day)
  
  results <- output
  save(results, file = here(OUTPUT, paste0(WhatAmITesting,"_PBK_results.RData")))
  
  PFAS_PBK_model_results[[i]] <- results
}

names(PFAS_PBK_model_results) <- PFAS_list

Save_PBK_results <- rbindlist(PFAS_PBK_model_results, idcol = "PFAS")
save(Save_PBK_results, file = here(OUTPUT,paste0(WhatAmITesting,"_PBK_results.RData")))

TK_results_list <- NULL

# produce result files and plots for all PFAS ----
source(here("Script", "functions_calculate_results.R"))

for (i in 1:length(PFAS_list)) {


  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}------------------------------------"))
  # load(here(INPUT, paste0(PFAS_i, "_model_results.RData")))
  
  # Inputs ----

  # Select results
  PBK_OUTPUT <- PFAS_PBK_model_results[[PFAS_i]] #results

  # load oral dose (different per PFAS if comparing results to Abraham)
  OD <- filter(Dose_df, PFAS == PFAS_i) |>
    select(Dose) |> pull()

  # Check actual model parameters ----
  Parameters_df <- PBK_OUTPUT |>
    select(fup, fu_L_ec, fu_PTL, fu_PT_ec,
           QI, Clear_apparent, Clear_OATP2B1_I, Clear_uptake_I,
           QL, Clear_OATP1B1, Clear_OATP1B3, Clear_OATP2B1_L, Clear_NTCP,  Clear_uptake_L,
           Clear_BSEP, Clear_excretion_L,
           GFR, QK, Clear_OAT3, Clear_OAT1, 
           Clear_OAT4, Clear_net_reabsorption) |> 
    distinct() |> 
    transmute(Q_I = QI,
              # Clear_apparent_I_in = Clear_apparent,
              # Clear_OATP2B1_I_in = Clear_OATP2B1_I,
              Clear_uptake_I = Clear_uptake_I,
              Q_L = QL+QI,
              # Clear_OATP1B1_L_in = Clear_OATP1B1,
              # Clear_OATP1B3_L_in = Clear_OATP1B3,
              # Clear_OATP2B1_L_in =  Clear_OATP2B1_L,
              # Clear_NTCP_L_in = Clear_NTCP,
              Clear_uptake_L = Clear_uptake_L,
              # Clear_BSEP_L_out = Clear_BSEP,
              Clear_excretion_L = Clear_excretion_L,
              Q_K = QK, # - (fup*GFR),
              # Clear_OAT3_K_out = Clear_OAT3,
              # Clear_OAT1_K_out = Clear_OAT1,
              # Clear_elimination_K = Clear_OAT1 +  Clear_OAT3,
              # GFR = fup*GFR,
              # Clear_OAT4_K_in = Clear_OAT4,
              Clear_net_reabsorption_K = Clear_net_reabsorption
    )  
  
  P_long_df <- Parameters_df |>
    pivot_longer(everything(),
                 values_to = "Rate",
                 names_to = "Parameter") |> 
    mutate(
      Tissue = case_when(
        str_detect(Parameter,"_I") ~ "Intestine",
        
        str_detect(Parameter, "_K") ~ "Kidney",
        # str_detect(Parameter, "GFR") ~ "Kidney",
        
        TRUE ~ "Liver"),
      Direction = case_when(
        str_detect(Parameter,"_uptake") ~ "Uptake",
        str_detect(Parameter, "_excretion") ~ "Elimination",
        str_detect(Parameter, "_reabsorption") ~ "Reabsorption",
        
        TRUE ~ "Blood-flow")) |> 
    # filter(Parameter %in% c("Q_I", "Q_L", "Q_K", "GFR")) |> 
    mutate(Parameter = ifelse(Parameter %in% c("Q_I", "Q_L", "Q_K"),
                              "Blood-flow", Parameter),
           Parameter = str_replace_all(Parameter, "_uptake", ""),
           Parameter = str_replace_all(Parameter, "_excretion", ""),
           Parameter = str_replace_all(Parameter, "_reabsorption", ""),
           Parameter = str_replace_all(Parameter, "_I", ""), 
           Parameter = str_replace_all(Parameter, "_K", ""), 
           Parameter = str_replace_all(Parameter, "_L", "") 
    )
  
  Parameters_plot <- ggplot() +
    geom_bar(data = filter(P_long_df, Parameter != "Blood-flow"), #, "GFR")),
             # Direction %in% c("Uptake")),
             aes(x = Direction, y = Rate, fill = Parameter),
             stat = "identity",
             position = "dodge",
             width = 0.5) +
    scale_fill_viridis(option = "D", discrete = TRUE) +
    geom_bar(data = filter(P_long_df, Parameter == "Blood-flow"), #, "GFR")),
             aes(x = Direction, y = Rate,
                 colour = Parameter),
             # colour = "grey",
             fill = "grey",
             # alpha = 0.8,
             stat = "identity",
             position = "dodge",
             width = 0.5) +
    scale_color_manual(values = c(
      "Blood-flow" = "grey",
      "GFR" = "grey9"
    )) +
    coord_flip() +
    facet_wrap(~ factor(Tissue, c("Intestine", "Liver", "Kidney")), 
               scales = "free_x") +
    labs(title = paste0(PFAS_i),
         x = "",
         y = "Rate (L/day)") +
    theme_quick
  
  Parameters_plot = Parameters_plot + CP_theme +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "none")
  
  ggsave(plot = Parameters_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_Parameters_plot.png")),
         dpi = 600,
         width = 16,
         height = 10,
         units = "cm")
  
  write.csv(Parameters_df, file = here(OUTPUT, paste0(PFAS_i, "_Parameters.csv")))
  # Predicted tissue concentrations ----
  
  Tissues_list <- c("CI", "CL_ec", "CL_ic", "CPTL", "CPTT", "CRKL", "CRKT", "CA", "CR", "CP")
  
  VL <- parm.phys[["VLc"]]*parm.phys[["BW"]]
  VL_ic <- parm.phys[["VL_icc"]]*VL
  VL_ec <- parm.phys[["VL_ecc"]]*VL
  
  VK <- parm.phys[["VKc"]]*parm.phys[["BW"]]
  VPT <- parm.phys[["VPTc"]]*VK
  VPTT <- VPT*parm.phys[["VPTTc"]]
  VPTL <- VPT*parm.phys[["VPTLc"]]
  VRKL <- (VK - VPT)*parm.phys[["VRKLc"]]
  VRKT <- VK - VPTT - VPTL - VRKL
  
  C_Tissues <- PBK_OUTPUT |> select(time, Tissues_list) |>
    mutate(CL = (CL_ec*VL_ec + CL_ic*VL_ic)/VL,
           CK = (CPTL*VPTL + CPTT*VPTT + CRKL*VRKL + CRKT*VRKT)/VK,
           Time = time/365) |> # time to years
    select(!c(time,CL_ec, CL_ic, CPTL, CPTT, CRKL, CRKT))
  
  C_Tissues_long <- C_Tissues |>
    pivot_longer(!Time,
                 names_to = "Tissue",
                 values_to = "Concentration"
    )
  
  Tissue_labeller <- c(
    "CP" = "Plasma",
    "CI" = "Intestine",
    "CL" = "Liver",
    "CK" = "Kidney",
    "CA" = "Adipose",
    "CR" = "Rest")
  
  C_Tissues_plot <- ggplot(C_Tissues_long, aes(Time, Concentration)) +
    geom_path() +
    facet_wrap(~ Tissue,
               scales = "free_y",
               nrow = 2,
               labeller = as_labeller(Tissue_labeller)
    ) +
    labs(title = paste0(PFAS_i),
         x = "Time (years)",
         y = "Concentration (ng/ml)") +
    theme_quick
  
  C_Tissues_plot <- C_Tissues_plot + CP_theme +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(plot = C_Tissues_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_C_Tissues_plot.png")),
         dpi = 600,
         width = 16,
         height = 11,
         units = "cm")
  
  write.csv(C_Tissues, file = here(OUTPUT, paste0(PFAS_i, "_C_Tissues.csv")))
  
  # TK results ----
  # calculate AUC, half-life, elimination rates, fraction absorbed etc
  TK_results <- calc_TK_results(PBK_OUTPUT, OD)
  TK_df <- tibble(Parameter = names(TK_results),
                  Value = TK_results)
  TK_results_list[[i]] <- TK_results

  TK_results_plot <- ggplot(TK_df, aes(Parameter, Value)) +
    geom_point() +
    labs(title = paste0(PFAS_i),
         x = "",
         y = "") +
    coord_flip() +
    theme_quick

  TK_results_plot <- TK_results_plot + CP_theme +
    theme(axis.text.x = element_text(angle = 0))

  ggsave(plot = TK_results_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_TK_results_plot.png")),
         dpi = 600,
         width = 10,
         height = 16,
         units = "cm")

}

names(TK_results_list) <- PFAS_list
Save_TK_results <- as.data.frame(do.call(rbind, TK_results_list)) |>
  mutate(Test = WhatAmITesting)

write.csv(Save_TK_results, here(OUTPUT,paste0(WhatAmITesting,"_TK_results.csv")))


# Evaluate against measured ####

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

  nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>
    mutate(PFAS = str_replace_all(PFAS, "_", "-"))
  
  HL_measured_df <- full_join(HL_measured, nC_MW_df)

  HL_predicted_df <- Save_TK_results |> 
    rownames_to_column() |>
    select(rowname, HL) |> 
    rename("PFAS" = rowname, 
           "Value" = HL)

  HL_predicted_df <- full_join(HL_predicted_df, nC_MW_df)


  # HL_predicted_plot <- ggplot() +
  #   geom_point(data = HL_predicted_df,
  #                 aes(x = reorder(PFAS, -MW),
  #                     y = log10(Value)),
  #                 width = .01,
  #                 linewidth = .4) +
  #   scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
  #   labs(x = "",
  #        y = "Half life (years)") +
  #   coord_flip() +
  #   theme_quick
  # HL_predicted_plot

  # HL_predicted_plot <- HL_predicted_plot + CP_theme +
  #   theme(legend.position = "right")
  # 
  # ggsave(plot = HL_predicted_plot,
  #        file = here(OUTPUT,"HL_predicted_plott.png"),
  #        dpi = 600,
  #        width = 12,
  #        height = 12,
  #        units = "cm")


  HL_MeasPred_plot <- ggplot(HL_measured_df,
                             aes(reorder(PFAS, -MW),
                                 log10(Value))) +
    geom_jitter(color = "grey60",
                width = 0.2,
                alpha = 0.7,
                size = .1) +
    geom_point(data = HL_predicted_df,
               aes(x = reorder(PFAS, -MW),
                   y = log10(Value)),
                   colour = "purple",
               width = .01,
               linewidth = .4) +
    scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
    labs(x = "",
         y = "Half life (years)") +
    coord_flip() +
    theme_quick
  HL_MeasPred_plot

  HL_MeasPred_plot <- HL_MeasPred_plot + CP_theme +
    theme(axis.text.x = element_text(angle = 0))
  

  ggsave(plot = HL_MeasPred_plot,
         file = here(OUTPUT,"HL_MeasPred_plot.png"),
         dpi = 600,
         width = 12,
         height = 12,
         units = "cm")



  ## Plasma ----
  Measured_Abraham_df <- Measured_df |>
    filter(Reference == "Abraham_2024")

  C_plasma_measured_df <- Measured_df |>
    filter(Parameter == "C_plasma") |>
    mutate(Time = as.numeric(Time_days)) |> # time in days
    select(PFAS, Time, Value)

  C_plasma_predicted_df <- Save_PBK_results |>
    select(PFAS, time, CP) |> 
    filter(time <= 450)

  ObsPred_plasma <- 
    ggplot() +
    geom_path(data = C_plasma_predicted_df,
              aes(time, CP)) +
    geom_point(data = C_plasma_measured_df,
               aes(Time, Value),
               color = "purple",
               size = 0.5) +
    facet_wrap(~ PFAS, scales = "free_y") +
    labs(x = "Time (days)",
         y = "Plasma concentration (μg/L)") +
    theme_quick
  ObsPred_plasma
  ObsPred_plasma <- ObsPred_plasma + CP_theme
  ggsave(plot = ObsPred_plasma,
                file = here(OUTPUT, "ObsPred_plasma_plot.png"),
                dpi = 600,
                width = 18,
                height = 12,
                units = "cm")

  
  
  
  
  
  
  
  
  