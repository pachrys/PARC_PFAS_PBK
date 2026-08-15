# --------------------------------------------------------------------------- #
# Script for analysing Quick_test results
# Chrysanthi Pachoulide, August 2026
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
    # axis.text.x = element_text(size = 35, face = "bold", angle = 30),
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
OUTPUT <- here("Output", "AnalyseMC") 
dir.create(OUTPUT, recursive = TRUE)

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

nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

# # TK ----
# 
# INPUT <- here("Output", "AnalyseParamInfluence") 
# 
# TK_df <- list.files(path = here(INPUT),
#                            pattern = "\\_TK_results.csv$") |>
#   map_df(~ read_csv(file = here(INPUT, .))) |> 
#   rename(PFAS = ...1)
# 
# Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |>
#   rename(PFAS = Acronym) |>
#   mutate(PFAS = str_replace_all(PFAS, "_", "-"))
# 
# ## Half life ----
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
# HL_measured_df <- full_join(HL_measured, nC_MW_df)
# 
# HL_Measured_plot <- ggplot(HL_measured_df,
#          aes(reorder(PFAS, -MW),
#              log10(Value),
#              color = Reference)) +
#   geom_jitter(width = 0.2,
#               alpha = 0.7,
#               size = .5) +
#   scale_color_viridis(option = "H", discrete = TRUE, end  = 0.9) +
#   labs(x = "",
#        y = "Half life log10(years)") +
#   coord_flip() +
#   theme_quick
# HL_Measured_plot
# 
# HL_predicted_df <- TK_df |>
#   pivot_longer(AUC:FbilEl, names_to = "Parameter", values_to = "Value") |> 
#   filter(Parameter == "HL")
# 
# HL_predicted_df <- full_join(HL_predicted_df, nC_MW_df)
# 
# 
# HL_predicted_plot <- ggplot() +
#   geom_point(data = HL_predicted_df,
#                 aes(x = reorder(PFAS, -MW),
#                     y = log10(Value),
#                     colour = Test),
#                 width = .01,
#                 linewidth = .4) +
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
#   labs(x = "",
#        y = "Half life (years)") +
#   coord_flip() +
#   theme_quick
# HL_predicted_plot
# 
# # HL_predicted_plot <- HL_predicted_plot + CP_theme +
# #   theme(legend.position = "right")
# # 
# # ggsave(plot = HL_predicted_plot,
# #        file = here(OUTPUT,"HL_predicted_plott.png"),
# #        dpi = 600,
# #        width = 12,
# #        height = 12,
# #        units = "cm")
# 
# 
# HL_MeasPred_plot <- ggplot(HL_measured_df,
#                            aes(reorder(PFAS, -MW),
#                                log10(Value))) +
#   geom_jitter(color = "grey60",
#               width = 0.2,
#               alpha = 0.7,
#               size = .1) +
#   geom_point(data = HL_predicted_df,
#              aes(x = reorder(PFAS, -MW),
#                  y = log10(Value),
#                  colour = Test),
#              width = .01,
#              linewidth = .4) +
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
#   labs(x = "",
#        y = "Half life (years)") +
#   coord_flip() +
#   theme_quick
# HL_MeasPred_plot
# 
# HL_MeasPred_plot <- HL_MeasPred_plot + CP_theme
# 
# ggsave(plot = HL_MeasPred_plot,
#        file = here(OUTPUT,"HL_MeasPred_plot.png"),
#        dpi = 600,
#        width = 12,
#        height = 12,
#        units = "cm")
# 
# 
# 









# Concentrations ----

INPUT <- here("Output", "BSEP0_OAT4toQK") # BSEPs, BSEPsOAT4onQK, OAT4onQK, Inital

BSEP0_OAT4toQK <- list.files(path = here(INPUT),
                                pattern = "\\_C_Tissues.csv$") |>
  map_df(~ read_csv(file = here(INPUT, .), id = "name"))  |>
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBK/Output/BSEP0_OAT4toQK/",
                                "")) |>
  mutate(PFAS = str_replace_all(PFAS,
                                "_C_Tissues.csv",""),
         Test = "BSEP0_OAT4toQK") |>
  select(CI:Test)

INPUT <- here("Output", "BSEPdiv10OAT4toQK") # BSEPs, BSEPsOAT4onQK, OAT4onQK, Inital

BSEPdiv10OAT4toQK <- list.files(path = here(INPUT),
                    pattern = "\\_C_Tissues.csv$") |>
  map_df(~ read_csv(file = here(INPUT, .), id = "name"))  |>
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBK/Output/BSEPdiv10OAT4toQK/",
                                "")) |>
  mutate(PFAS = str_replace_all(PFAS,
                                "_C_Tissues.csv",""),
         Test = "BSEPdiv10OAT4toQK") |>
  select(CI:Test)

# INPUT <- here("Output", "BSEPs") # BSEPs, BSEPsOAT4onQK, OAT4onQK, Inital
# 
# BSEPs <- list.files(path = here(INPUT),
#                            pattern = "\\_C_Tissues.csv$") |>
#   map_df(~ read_csv(file = here(INPUT, .), id = "name"))  |>
#   mutate(PFAS = str_replace_all(name,
#                                "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBK/Output/BSEPs/",
#                                 "")) |>
#   mutate(PFAS = str_replace_all(PFAS,
#                                 "_C_Tissues.csv",""),
#          Test = "BSEPs") |>
#   select(CI:Test)

INPUT <- here("Output", "BSEPsOAT4onQK") # BSEPs, BSEPsOAT4onQK, OAT4onQK, Inital

BSEPsOAT4onQK <- list.files(path = here(INPUT),
                    pattern = "\\_C_Tissues.csv$") |>
  map_df(~ read_csv(file = here(INPUT, .), id = "name"))  |>
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBK/Output/BSEPsOAT4onQK/",
                                "")) |>
  mutate(PFAS = str_replace_all(PFAS,
                                "_C_Tissues.csv",""),
         Test = "BSEPsOAT4onQK") |>
  select(CI:Test)

INPUT <- here("Output", "OAT4onQK") # BSEPs, BSEPsOAT4onQK, OAT4onQK, Inital

OAT4onQK <- list.files(path = here(INPUT),
                    pattern = "\\_C_Tissues.csv$") |>
  map_df(~ read_csv(file = here(INPUT, .), id = "name"))  |>
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBK/Output/OAT4onQK/",
                                "")) |>
  mutate(PFAS = str_replace_all(PFAS,
                                "_C_Tissues.csv",""),
         Test = "OAT4onQK") |>
  select(CI:Test)

INPUT <- here("Output", "InitalParams") # BSEPs, BSEPsOAT4onQK, OAT4onQK, Inital

InitalParams <- list.files(path = here(INPUT),
                    pattern = "\\_C_Tissues.csv$") |>
  map_df(~ read_csv(file = here(INPUT, .), id = "name"))  |>
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBK/Output/InitalParams/",
                                "")) |>
  mutate(PFAS = str_replace_all(PFAS,
                                "_C_Tissues.csv",""),
         Test = "Inital") |>
  select(CI:Test)

C_df <- rbind(InitalParams, BSEPsOAT4onQK, OAT4onQK, BSEPdiv10OAT4toQK, BSEP0_OAT4toQK) |> #BSEPs
  pivot_longer(CI:CK, names_to = "Tissue", values_to = "Concentration")


for (i in 1:length(PFAS_list)) {

  PFAS_i <- PFAS_list[i]
  PFAS_df <- C_df |> filter(PFAS == PFAS_i)
  save(PFAS_df, file = here(OUTPUT, paste0(PFAS_i,"_Tissue_concentrations.RData")))

  Plot_i <- ggplot(PFAS_df,aes(Time, Concentration,
                               colour = Test,
                               linetype = Test)) +
    geom_path() +
    scale_color_viridis(option = "B",  discrete = TRUE,
                        begin = 0.2,  end = 0.8, alpha = 0.9) +
    scale_linetype_manual(values = c("Inital" = "solid",
                                     "BSEPdiv1000_OAT4toQK" = "dotdash", #BSEPs
                                     "BSEPsOAT4onQK" = "longdash",
                                     "OAT4onQK" = "dotted",
                                     "BSEPdiv10OAT4toQK" = "dashed"))+
    facet_wrap(~Tissue, scales = "free_y") +
    labs(title = paste0(PFAS_i),
         x = "Time (years)",
         y = "Concentration (ng/ml)") +
  theme_quick

  Plot_i <- Plot_i + CP_theme +
    theme(legend.position = "bottom")

  ggsave(plot = Plot_i,
         file = here(OUTPUT, paste0(PFAS_i, "_Tissue_concentrations.png")),
         dpi = 600,
         width = 12,
         height = 9,
         units = "cm")

}





