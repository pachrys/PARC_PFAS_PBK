# # MCMC Results Analysis 
# # -----------------------------
# 
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
    axis.text.x = element_text(size = 35, face = "bold"), #, angle = 30),
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
    legend.position = "none", # legend.position = c(0.9, 0.9),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
    legend.margin = margin(1, 1, 1, 1),

    # strip.text.x = element_blank()
    strip.text.x = element_text(size = 35),
    strip.text.y = element_text(size = 35)

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

    legend.position = "none"
  )

# Set output storage directory
INPUT <- here("Output", "MonteCarlo_1000_Abraham")
OUTPUT <- here("Output","MonteCarlo_1000_Abraham", "SumResults") #format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

Dose_df <- read_excel(here("Input","Evaluation_data.xlsx")) %>% 
  filter(Reference == "Abraham_2024") |> 
  select(!c(PFAS_name, CAS_n)) %>% 
  rename(PFAS = Acronym) |> 
  mutate(Dose = round(Dose_ugPerkg_bw*82, 3)) |>  #82 is the body weight of Abraham
  select(PFAS, Dose) |> 
  mutate(PFAS = str_replace_all(PFAS, "_", "-")) |> 
  distinct()


# Analyse results ----

calc_TK_results <- function(results, OD) {

  data <-  results 
  plasma <- data |>
    select(time, CP)

  # AUC ----
  AUC <- trapz(data$time, data$CP) # ng * day/ml

  # Cmax and Tmax ----
  Cmax <- max(data$CP) # ng/ml
  Tmax <- filter(data, CP == max(data$CP)) %>% pull(time) # days

  # plasma half-life ----
  # elimination phase, i.e. after the tmax
  elimination_phase <- data |>
    filter(time >= Tmax,
           CP >= 1e-9) |>
    mutate(logCP = log(CP))

  # check that there is enough data points to run linear regression
  if (nrow(elimination_phase) < 3) {
    warning("not enough data-points to calculate half life")
    return(NA)
  }

  # linear regression
  elimination_model <- lm(logCP ~ time, data = elimination_phase)

  # calculate elimination constant (i.e slope of the line)
  kel <- -elimination_model$coefficients[[2]]

  # calculate the half-life
  half_life <- (log(2) / kel) / 365 # unit correction to years

  # fraction absorbed ----

  # all amounts in tissues at Tmax
  A_Tmax <- data %>% filter(time == Tmax) |>
    select(starts_with("A")) |>
    select(!c(AIL, AFe))

  # calculate the total amount absorbed at Tmax
  ATotal_Tmax <- rowSums(A_Tmax)

  # Absorbed fraction (as a percentage of total oral dose)
  OD <- OD
  Fa <- ATotal_Tmax/OD *100

  # fecal and urinary excretion ----
  # Calculate excretion (feces and urine)

  A_feces <- data |> filter(time > Tmax) |> select(AFe) |>
    mutate(dailyFE = (AFe - lag(AFe, n = 1, default = first(AFe)))*0.1) |>
    filter(dailyFE > 0)
  FfeEl <- last(A_feces$AFe)/OD * 100
  First10daysFE <- data |> filter(time == Tmax) |> select(AFe) |> pull()
  AvDailyFE <- mean(A_feces$dailyFE) # average fecal excretion per day
  
  A_urine <- data |> filter(time > Tmax) |> select(AUr) |>
    mutate(dailyFE = (AUr - lag(AUr, n = 1, default = first(AUr)))*0.1) |>
    filter(dailyFE > 0)
  FurEl <- last(A_urine$AUr)/OD * 100
  First10daysUR <- data |> filter(time == Tmax) |> select(AUr) |> pull()
  AvDailyUR <- mean(A_urine$dailyFE) # average fecal excretion per day
  

  return(c(AUC = AUC,
           Cmax = Cmax,
           Tmax = Tmax,
           kel = kel,
           HL = half_life,
           Fa = Fa,
           FfeEl = FfeEl,
           First10daysFE = First10daysFE,
           AvDailyFE = AvDailyFE,
           FurEl = FurEl,
           First10daysUR = First10daysUR,
           AvDailyUR = AvDailyUR
  ))
}

# ADONA ----
# 
# load(here(INPUT, "ADONA_Monte_Carlo_y.RData"))
# ADONA <- savey
# rm(savey)
# OD_ADONA <- filter(Dose_df, PFAS == "ADONA") |> 
#   select(Dose) |> pull()
# 
# # Plasma concentration
# ADONA_df <- do.call(rbind, ADONA) |> mutate(PFAS = "ADONA")
# ADONA_CP <- ADONA_df |> mutate(PFAS = "ADONA") |> select(time, CP, PFAS, Iteration)
# ADONA_CP_wide <- ADONA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# ADONA_CP_summary <- ADONA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(ADONA_CP_summary, file = here(OUTPUT, "ADONA_CP_sumstat.csv"))
# 

# TK parameters
# ADONA_TK_results <- sapply(ADONA, calc_TK_results, OD = OD_ADONA)
ADONA_TK_results_all <- do.call(rbind, ADONA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value") |> 
  separate(Variable,
           into = c("Parameter", "Variable")) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(ADONA_TK_results_all, file = here(OUTPUT, "ADONA_TK_sumstat.csv"))

ADONA_CP_plot <- ggplot(ADONA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "ADONA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
ADONA_CP_plot <- ADONA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = ADONA_CP_plot,
       file = here(OUTPUT, "ADONA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")


# HFPO_DA ----
# 
# load(here(INPUT, "HFPO-DA_Monte_Carlo_y.RData"))
# HFPO_DA <- savey
# rm(savey)
# OD_HFPO_DA <- filter(Dose_df, PFAS == "HFPO_DA") |> 
#   select(Dose) |> pull()
# 
# # Plasma concentration
# HFPO_DA_df <- do.call(rbind, HFPO_DA) |> mutate(PFAS = "HFPO_DA")
# HFPO_DA_CP <- HFPO_DA_df |> mutate(PFAS = "HFPO_DA") |> select(time, CP, PFAS, Iteration)
# HFPO_DA_CP_wide <- HFPO_DA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# HFPO_DA_CP_summary <- HFPO_DA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(HFPO_DA_CP_summary, file = here(OUTPUT, "HFPO_DA_CP_sumstat.csv"))


# TK parameters
# HFPO_DA_TK_results <- sapply(HFPO_DA, calc_TK_results, OD = OD_HFPO_DA)
HFPO_DA_TK_results_all <- do.call(rbind, HFPO_DA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(HFPO_DA_TK_results_all, file = here(OUTPUT, "HFPO_DA_TK_sumstat.csv"))

HFPO_DA_CP_plot <- ggplot(HFPO_DA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "HFPO_DA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
HFPO_DA_CP_plot <- HFPO_DA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = HFPO_DA_CP_plot,
       file = here(OUTPUT, "HFPO_DA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")




# PFBA ----

load(here(INPUT, "PFBA_Monte_Carlo_y.RData"))
PFBA <- savey
rm(savey)
OD_PFBA <- filter(Dose_df, PFAS == "PFBA") |>
  select(Dose) |> pull()
# # TK parameters
# PFBA <- PFBA
# PFBA_TK_results <- sapply(PFBA, calc_TK_results, OD = OD_PFBA)
PFBA_TK_results_all <- do.call(rbind, PFBA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFBA_TK_results_all, file = here(OUTPUT, "PFBA_TK_sumstat.csv"))



# Plasma concentration
PFBA_df <- do.call(rbind, PFBA) |> mutate(PFAS = "PFBA")
PFBA_CP <- PFBA_df |> mutate(PFAS = "PFBA") |> select(time, CP, PFAS, Iteration)

PFBA_summary <- PFBA_df |> filter(CP > 0.00025) |>
  select(!starts_with("A")) |>
  select(!time) |>
  select(!CIL:CP) |>
  distinct()
write.csv(PFBA_summary, file = here(OUTPUT, "PFBA_summary.csv"))

PFBA_CP <- PFBA_CP |> 
  filter(CP < 0.00025) |> 
  mutate(CP = ifelse(CP <= 1e-9, 0, CP))

PFBA_CP_wide <- PFBA_CP |> select(!PFAS) |>
  mutate(Iteration = paste0("n_", Iteration)) |>
  pivot_wider(names_from = Iteration, values_from = CP)
PFBA_CP_summary <- PFBA_CP_wide |>
  select(!time) |> t() |>
  as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = ~median(.,na.rm = TRUE),
                          mean = ~mean(.,na.rm = TRUE),
                          q25 = ~quantile(., 0.25, na.rm = TRUE),
                          q75 = ~quantile(.,0.75, na.rm = TRUE)))) |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value") |>
  separate(Variable,
           into = c("Time", "Variable")) |>
  mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
  mutate(Time = ((Time-1)*10) /365) |>  #unit correction, to years
  pivot_wider(names_from = Variable, values_from = Value)


PFBA_CP_plot <- 
  
  ggplot(PFBA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFBA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick

PFBA_CP_plot

PFBA_CP_plot <- PFBA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFBA_CP_plot,
       file = here(OUTPUT, "PFBA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")
write.csv(PFBA_CP_summary, file = here(OUTPUT, "PFBA_CP_sumstat.csv"))



#
#
# PFBS ----
#
# 
# load(here(INPUT, "PFBS_Monte_Carlo_y.RData"))
# PFBS <- savey
# rm(savey)
# 
# # Plasma concentration
# PFBS_df <- do.call(rbind, PFBS) |> mutate(PFAS = "PFBS")
# 
# PFBS_summary <- PFBS_df |> filter(CP > 0.004) |>
#   select(!starts_with("A")) |>
#   select(!time) |>
#   select(!CIL:CP) |>
#   distinct()
# write.csv(PFBS_summary, file = here(OUTPUT, "PFBS_summary.csv"))
# 
# PFBS_CP <- PFBS_df |> mutate(PFAS = "PFBS") |> select(time, CP, PFAS, Iteration)
# PFBS_CP <- PFBS_CP |> filter(CP < 0.004) |>
#   mutate(CP = ifelse(CP <= 1e-09, 0, CP))
# 
# PFBS_CP_wide <- PFBS_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFBS_CP_summary <- PFBS_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = ~median(.,na.rm = TRUE),
#                           mean = ~mean(.,na.rm = TRUE),
#                           q25 = ~quantile(., 0.25, na.rm = TRUE),
#                           q75 = ~quantile(.,0.75, na.rm = TRUE)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFBS_CP_summary, file = here(OUTPUT, "PFBS_CP_sumstat.csv"))
# 

# TK parameters
# OD_PFBS <- filter(Dose_df, PFAS == "PFBS") |> 
#   select(Dose) |> pull()
# PFBS <- PFBS
# PFBS_TK_results <- sapply(PFBS, calc_TK_results, OD = OD_PFBS)
PFBS_TK_results_all <- do.call(rbind, PFBS_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFBS_TK_results_all, file = here(OUTPUT, "PFBS_TK_sumstat.csv"))

PFBS_CP_plot <- ggplot(PFBS_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFBS",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFBS_CP_plot

PFBS_CP_plot <- PFBS_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFBS_CP_plot,
       file = here(OUTPUT, "PFBS_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



#
# PFDA ----
#
# 
# load(here(INPUT, "PFDA_Monte_Carlo_y.RData"))
# PFDA <- savey
# rm(savey)
# 
# # Plasma concentration
# PFDA_df <- do.call(rbind, PFDA) |> mutate(PFAS = "PFDA")
# PFDA_CP <- PFDA_df |> mutate(PFAS = "PFDA") |> select(time, CP, PFAS, Iteration)
# PFDA_CP_wide <- PFDA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFDA_CP_summary <- PFDA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFDA_CP_summary, file = here(OUTPUT, "PFDA_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFDA <- filter(Dose_df, PFAS == "PFDA") |> 
#   select(Dose) |> pull()
# PFDA_TK_results <- sapply(PFDA, calc_TK_results, OD = OD_PFDA)
PFDA_TK_results_all <- do.call(rbind, PFDA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFDA_TK_results_all, file = here(OUTPUT, "PFDA_TK_sumstat.csv"))

PFDA_CP_plot <- ggplot(PFDA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFDA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFDA_CP_plot <- PFDA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFDA_CP_plot,
       file = here(OUTPUT, "PFDA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



#
# PFDoDA ----
#
# 
# load(here(INPUT, "PFDoDA_Monte_Carlo_y.RData"))
# PFDoDA <- savey
# rm(savey)
# 
# # Plasma concentration
# PFDoDA_df <- do.call(rbind, PFDoDA) |> mutate(PFAS = "PFDoDA")
# PFDoDA_CP <- PFDoDA_df |> mutate(PFAS = "PFDoDA") |> select(time, CP, PFAS, Iteration)
# PFDoDA_CP_wide <- PFDoDA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFDoDA_CP_summary <- PFDoDA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFDoDA_CP_summary, file = here(OUTPUT, "PFDoDA_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFDoDA <- filter(Dose_df, PFAS == "PFDoDA") |> 
# select(Dose) |> pull()
# PFDoDA_TK_results <- sapply(PFDoDA, calc_TK_results, OD = OD_PFDoDA)
PFDoDA_TK_results_all <- do.call(rbind, PFDoDA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFDoDA_TK_results_all, file = here(OUTPUT, "PFDoDA_TK_sumstat.csv"))

PFDoDA_CP_plot <- ggplot(PFDoDA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFDoDA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFDoDA_CP_plot <- PFDoDA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFDoDA_CP_plot,
       file = here(OUTPUT, "PFDoDA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



#
# PFHpA ----
# load(here(INPUT, "PFHpA_Monte_Carlo_y.RData"))
# PFHpA <- savey
# rm(savey)
# 
# # Plasma concentration
# PFHpA_df <- do.call(rbind, PFHpA) |> mutate(PFAS = "PFHpA")
# PFHpA_CP <- PFHpA_df |> mutate(PFAS = "PFHpA") |> select(time, CP, PFAS, Iteration)
# PFHpA_CP_wide <- PFHpA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFHpA_CP_summary <- PFHpA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFHpA_CP_summary, file = here(OUTPUT, "PFHpA_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFHpA <- filter(Dose_df, PFAS == "PFHpA") |> 
#   select(Dose) |> pull()
# PFHpA_TK_results <- sapply(PFHpA, calc_TK_results, OD = OD_PFHpA)
PFHpA_TK_results_all <- do.call(rbind, PFHpA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFHpA_TK_results_all, file = here(OUTPUT, "PFHpA_TK_sumstat.csv"))

PFHpA_CP_plot <- ggplot(PFHpA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFHpA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFHpA_CP_plot <- PFHpA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFHpA_CP_plot,
       file = here(OUTPUT, "PFHpA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



#
# PFHxA ----
#
# 
# load(here(INPUT, "PFHxA_Monte_Carlo_y.RData"))
# PFHxA <- savey
# rm(savey)
# 
# # Plasma concentration
# PFHxA_df <- do.call(rbind, PFHxA) |> mutate(PFAS = "PFHxA")
# PFHxA_CP <- PFHxA_df |> mutate(PFAS = "PFHxA") |> select(time, CP, PFAS, Iteration)
# 
# PFHxA_summary <- PFHxA_df |> filter(CP > 0.020) |>
#   select(!starts_with("A")) |>
#   select(!time) |>
#   select(!CIL:CP) |>
#   distinct()
# write.csv(PFHxA_summary, file = here(OUTPUT, "PFHxA_summary.csv"))
# 
# PFHxA_CP <- PFHxA_df |> mutate(PFAS = "PFHxA") |> select(time, CP, PFAS, Iteration)
# 
# PFHxA_CP <- PFHxA_CP |> filter(CP < 0.020) |>
#   mutate(CP = ifelse(CP <= 1e-09, 0, CP))
# 
# PFHxA_CP_wide <- PFHxA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFHxA_CP_summary <- PFHxA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = ~median(.,na.rm = TRUE),
#                           mean = ~mean(.,na.rm = TRUE),
#                           q25 = ~quantile(., 0.25, na.rm = TRUE),
#                           q75 = ~quantile(.,0.75, na.rm = TRUE)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# 
# write.csv(PFHxA_CP_summary, file = here(OUTPUT, "PFHxA_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFHxA <- filter(Dose_df, PFAS == "PFHxA") |> 
#   select(Dose) |> pull()
# PFHxA <- PFHxA #[-It_PFHxA]
# PFHxA_TK_results <- sapply(PFHxA, calc_TK_results, OD_PFHxA)
PFHxA_TK_results_all <- do.call(rbind, PFHxA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFHxA_TK_results_all, file = here(OUTPUT, "PFHxA_TK_sumstat.csv"))

PFHxA_CP_plot <- ggplot(PFHxA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFHxA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFHxA_CP_plot

PFHxA_CP_plot <- PFHxA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFHxA_CP_plot,
       file = here(OUTPUT, "PFHxA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")

#
# PFHxS ----
# 
# load(here(INPUT, "PFHxS_Monte_Carlo_y.RData"))
# PFHxS <- savey
# rm(savey)
# 
# # Plasma concentration
# PFHxS_df <- do.call(rbind, PFHxS) |> mutate(PFAS = "PFHxS")
# PFHxS_CP <- PFHxS_df |> mutate(PFAS = "PFHxS") |> select(time, CP, PFAS, Iteration)
# PFHxS_CP_wide <- PFHxS_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFHxS_CP_summary <- PFHxS_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFHxS_CP_summary, file = here(OUTPUT, "PFHxS_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFHxS <- filter(Dose_df, PFAS == "PFHxS") |> 
#   select(Dose) |> pull()
# PFHxS_TK_results <- sapply(PFHxS, calc_TK_results, OD = OD_PFHxS)
PFHxS_TK_results_all <- do.call(rbind, PFHxS_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFHxS_TK_results_all, file = here(OUTPUT, "PFHxS_TK_sumstat.csv"))

PFHxS_CP_plot <- ggplot(PFHxS_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFHxS",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFHxS_CP_plot <- PFHxS_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFHxS_CP_plot,
       file = here(OUTPUT, "PFHxS_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



# PFNA ----
# 
# load(here(INPUT, "PFNA_Monte_Carlo_y.RData"))
# PFNA <- savey
# rm(savey)
# 
# # Plasma concentration
# PFNA_df <- do.call(rbind, PFNA) |> mutate(PFAS = "PFNA")
# PFNA_CP <- PFNA_df |> mutate(PFAS = "PFNA") |> select(time, CP, PFAS, Iteration)
# PFNA_CP_wide <- PFNA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFNA_CP_summary <- PFNA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFNA_CP_summary, file = here(OUTPUT, "PFNA_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFNA <- filter(Dose_df, PFAS == "PFNA") |> 
#   select(Dose) |> pull()
# PFNA_TK_results <- sapply(PFNA, calc_TK_results, OD = OD_PFNA)
PFNA_TK_results_all <- do.call(rbind, PFNA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFNA_TK_results_all, file = here(OUTPUT, "PFNA_TK_sumstat.csv"))

PFNA_CP_plot <- ggplot(PFNA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFNA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFNA_CP_plot <- PFNA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFNA_CP_plot,
       file = here(OUTPUT, "PFNA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



# PFOA ----
# 
# load(here(INPUT, "PFOA_Monte_Carlo_y.RData"))
# PFOA <- savey
# rm(savey)
# 
# # Plasma concentration
# PFOA_df <- do.call(rbind, PFOA) |> mutate(PFAS = "PFOA")
# PFOA_CP <- PFOA_df |> mutate(PFAS = "PFOA") |> select(time, CP, PFAS, Iteration)
# PFOA_CP_wide <- PFOA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFOA_CP_summary <- PFOA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFOA_CP_summary, file = here(OUTPUT, "PFOA_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFOA <- filter(Dose_df, PFAS == "PFOA") |> 
#   select(Dose) |> pull()
# PFOA_TK_results <- sapply(PFOA, calc_TK_results, OD = OD_PFOA)
PFOA_TK_results_all <- do.call(rbind, PFOA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFOA_TK_results_all, file = here(OUTPUT, "PFOA_TK_sumstat.csv"))

PFOA_CP_plot <- ggplot(PFOA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFOA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFOA_CP_plot <- PFOA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFOA_CP_plot,
       file = here(OUTPUT, "PFOA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



# PFOS ----
# 
# load(here(INPUT, "PFOS_Monte_Carlo_y.RData"))
# PFOS <- savey
# rm(savey)
# 
# # Plasma concentration
# PFOS_df <- do.call(rbind, PFOS) |> mutate(PFAS = "PFOS")
# PFOS_CP <- PFOS_df |> mutate(PFAS = "PFOS") |> select(time, CP, PFAS, Iteration)
# PFOS_CP_wide <- PFOS_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFOS_CP_summary <- PFOS_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFOS_CP_summary, file = here(OUTPUT, "PFOS_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFOS <- filter(Dose_df, PFAS == "PFOS") |> 
#   select(Dose) |> pull()

# PFOS_TK_results <- sapply(PFOS, calc_TK_results, OD = OD_PFOS)
PFOS_TK_results_all <- do.call(rbind, PFOS_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFOS_TK_results_all, file = here(OUTPUT, "PFOS_TK_sumstat.csv"))

PFOS_CP_plot <- ggplot(PFOS_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFOS",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFOS_CP_plot <- PFOS_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFOS_CP_plot,
       file = here(OUTPUT, "PFOS_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



# PFPeA ----
# 
# load(here(INPUT, "PFPeA_Monte_Carlo_y.RData"))
# PFPeA <- savey
# rm(savey)
# 
# # Plasma concentration
# PFPeA_df <- do.call(rbind, PFPeA) |> mutate(PFAS = "PFPeA")
# PFPeA_CP <- PFPeA_df |> mutate(PFAS = "PFPeA") |> select(time, CP, PFAS, Iteration)
# 
# PFPeA_summary <- PFPeA_df |> filter(CP > 0.02) |>
#   select(!starts_with("A")) |>
#   select(!time) |>
#   select(!CIL:CP) |>
#   distinct()
# write.csv(PFPeA_summary, file = here(OUTPUT, "PFPeA_summary.csv"))
# 
# PFPeA_CP <- PFPeA_df |> mutate(PFAS = "PFPeA") |> select(time, CP, PFAS, Iteration)
# 
# PFPeA_CP <- PFPeA_CP |> filter(CP < 0.02) |>
#   mutate(CP = ifelse(CP <= 1e-09, 0, CP))
# 
# 
# 
# PFPeA_CP_wide <- PFPeA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFPeA_CP_summary <- PFPeA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = ~median(.,na.rm = TRUE),
#                           mean = ~mean(.,na.rm = TRUE),
#                           q25 = ~quantile(., 0.25, na.rm = TRUE),
#                           q75 = ~quantile(.,0.75, na.rm = TRUE)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFPeA_CP_summary, file = here(OUTPUT, "PFPeA_CP_sumstat.csv"))
# 
# 
# # TK parameters
# PFPeA <- PFPeA
# OD_PFPeA <- filter(Dose_df, PFAS == "PFPeA") |> 
#   select(Dose) |> pull()

# PFPeA_TK_results <- sapply(PFPeA, calc_TK_results, OD = OD_PFPeA)
PFPeA_TK_results_all <- do.call(rbind, PFPeA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFPeA_TK_results_all, file = here(OUTPUT, "PFPeA_TK_sumstat.csv"))

PFPeA_CP_plot <- ggplot(PFPeA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFPeA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFPeA_CP_plot

PFPeA_CP_plot <- PFPeA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )


ggsave(plot = PFPeA_CP_plot,
       file = here(OUTPUT, "PFPeA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")



# PFUnDA ----
# 
# load(here(INPUT, "PFUnDA_Monte_Carlo_y.RData"))
# PFUnDA <- savey
# rm(savey)
# 
# # Plasma concentration
# PFUnDA_df <- do.call(rbind, PFUnDA) |> mutate(PFAS = "PFUnDA")
# PFUnDA_CP <- PFUnDA_df |> mutate(PFAS = "PFUnDA") |> select(time, CP, PFAS, Iteration)
# PFUnDA_CP_wide <- PFUnDA_CP |> select(!PFAS) |>
#   mutate(Iteration = paste0("n_", Iteration)) |>
#   pivot_wider(names_from = Iteration, values_from = CP)
# PFUnDA_CP_summary <- PFUnDA_CP_wide |>
#   select(!time) |> t() |>
#   as.data.frame() |>
#   summarise(across(where(is.numeric), .fns =
#                      list(median = median,
#                           mean = mean,
#                           q25 = ~quantile(., 0.25),
#                           q75 = ~quantile(.,0.75)))) |>
#   pivot_longer(everything(), names_to = "Variable",
#                values_to = "Value") |>
#   separate(Variable,
#            into = c("Time", "Variable")) |>
#   mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
#   mutate(Time = Time/365) |>  #unit correction, to years
#   pivot_wider(names_from = Variable, values_from = Value)
# write.csv(PFUnDA_CP_summary, file = here(OUTPUT, "PFUnDA_CP_sumstat.csv"))
# 
# 
# # TK parameters
# OD_PFUnDA <- filter(Dose_df, PFAS == "PFUnDA") |> 
#   select(Dose) |> pull()

# PFUnDA_TK_results <- sapply(PFUnDA, calc_TK_results, OD_PFUnDA)
PFUnDA_TK_results_all <- do.call(rbind, PFUnDA_TK_results) |>
  na.omit() |> as.data.frame() |>
  summarise(across(where(is.numeric), .fns =
                     list(median = median,
                          mean = mean,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(.,0.75))))   |>
  pivot_longer(everything(), names_to = "Variable",
               values_to = "Value")  |>
  separate(Variable,
           into = c("Parameter", "Variable"))|>
  # mutate(Value = ifelse(Parameter == "HL" | Parameter == "Tmax", Value/365, Value)) |>
  pivot_wider(names_from = Variable, values_from = Value)
write.csv(PFUnDA_TK_results_all, file = here(OUTPUT, "PFUnDA_TK_sumstat.csv"))

PFUnDA_CP_plot <- ggplot(PFUnDA_CP_summary, aes(Time, mean)) +
  geom_ribbon(aes(ymin = q25,
                  ymax = q75),
              fill = "grey70") +
  geom_path()+
  labs(title = "PFUnDA",
       x = "Time (years)",
       y = "Plasma concentration (ng/ml)") +
  theme_quick
PFUnDA_CP_plot <- PFUnDA_CP_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 0)
  )
ggsave(plot = PFUnDA_CP_plot,
       file = here(OUTPUT, "PFUnDA_CP_plot.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")

# Iterations <- read.csv(here(OUTPUT, "It_dirty.csv"))

# 
# design_PFBA <- read.csv(here(INPUT, "PFBAMC_ExpDesign.csv"))
# design_PFBA <- design_PFBA |>
#   mutate(Iteration = rownames(design_PFBA))
# select_design_PFBA <- design_PFBA |>
#   filter(Iteration %in% It_PFBA)
# 
# design_PFBS <- read.csv(here(INPUT, "PFBSMC_ExpDesign.csv"))
# design_PFBS <- design_PFBS |>
#   mutate(Iteration = rownames(design_PFBS))
# select_design_PFBS <- design_PFBS |>
#   filter(Iteration %in% It_PFBS)
# 
# design_PFHxA <- read.csv(here(INPUT, "PFHxAMC_ExpDesign.csv"))
# design_PFHxA <- design_PFHxA |>
#   mutate(Iteration = rownames(design_PFHxA))
# select_design_PFHxA <- design_PFHxA |>
#   filter(Iteration %in% It_PFHxA)
# 
# design_PFPeA <- read.csv(here(INPUT, "PFPeAMC_ExpDesign.csv"))
# design_PFPeA <- design_PFPeA |>
#   mutate(Iteration = rownames(design_PFPeA))
# select_design_PFPeA <- design_PFPeA |>
#   filter(Iteration %in% It_PFPeA)


# summary results ----
ADONA_TK_results_all <-ADONA_TK_results_all |> mutate(PFAS = "ADONA")
HFPO_DA_TK_results_all <- HFPO_DA_TK_results_all |> mutate(PFAS = "HFPO_DA")
PFBA_TK_results_all <- PFBA_TK_results_all |> mutate(PFAS = "PFBA")
PFBS_TK_results_all <- PFBS_TK_results_all |> mutate(PFAS = "PFBS")
PFDA_TK_results_all <- PFDA_TK_results_all |> mutate(PFAS = "PFDA")
PFDoDA_TK_results_all <- PFDoDA_TK_results_all |> mutate(PFAS = "PFDoDA")
PFHpA_TK_results_all <- PFHpA_TK_results_all |> mutate(PFAS = "PFHpA")
PFHxA_TK_results_all <- PFHxA_TK_results_all |> mutate(PFAS = "PFHxA")
PFHxS_TK_results_all <- PFHxS_TK_results_all |> mutate(PFAS = "PFHxS")
PFNA_TK_results_all <- PFNA_TK_results_all |> mutate(PFAS = "PFNA")
PFOA_TK_results_all <- PFOA_TK_results_all |> mutate(PFAS = "PFOA")
PFOS_TK_results_all <- PFOS_TK_results_all |> mutate(PFAS = "PFOS")
PFPeA_TK_results_all <- PFPeA_TK_results_all |> mutate(PFAS = "PFPeA")
PFUnDA_TK_results_all <- PFUnDA_TK_results_all |> mutate(PFAS = "PFUnDA")

Collect_all_tk <- rbind(ADONA_TK_results_all,
                        HFPO_DA_TK_results_all,
                        PFBA_TK_results_all,
                        PFBS_TK_results_all,
                        PFDA_TK_results_all,
                        PFDoDA_TK_results_all,
                        PFHpA_TK_results_all,
                        PFHxA_TK_results_all,
                        PFHxS_TK_results_all,
                        PFNA_TK_results_all,
                        PFOA_TK_results_all,
                        PFOS_TK_results_all,
                        PFPeA_TK_results_all,
                        PFUnDA_TK_results_all)
nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv"))

nC_MW_df <- nC_MW_df |> mutate(PFAS = str_replace_all(PFAS, "HFPO-DA", "HFPO_DA") )
Collect_all_tk <- full_join(Collect_all_tk, nC_MW_df)

write.csv(Collect_all_tk, file = here(OUTPUT,"Monte_Carlo_TK_results.csv"), row.names = FALSE)
# 
# Collect_all_tk <- read.csv(here(OUTPUT,"Monte_Carlo_TK_results.csv"))
HL_predicted <- Collect_all_tk |>
  filter(Parameter == "HL",
         PFAS != "PFTrDA")

Evaluation_data <- read_excel(here("Input","Evaluation_data.xlsx")) %>%
  select(!c(PFAS_name, CAS_n)) %>%
  rename(PFAS = Acronym)

raw_HL_measured <- Evaluation_data %>%
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

THL_measured <- HL_measured %>%
  group_by(PFAS) %>%
  summarise(mean = mean(Value),
            median = median(Value)) %>%
  mutate(HL_length = ifelse(median < 0.01, "very_short",
                            ifelse(median < 1, "short",
                                   ifelse(median < 5, "medium", "long")))) %>%
  mutate(HL_length = factor(HL_length, levels = c("very_short", "short", "medium", "long"))) |>
  select(PFAS, HL_length)

nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv"))

HL_measured <- full_join(HL_measured, nC_MW_df)
HL_measured <- full_join(HL_measured, THL_measured) |>
  filter(PFAS != "PFTrDA") |>
  mutate(PFAS = str_replace(PFAS, "_", "-"))

HL_predicted <- full_join(HL_predicted, THL_measured) |>
  filter(PFAS != "PFTrDA") |>
  mutate(PFAS = str_replace(PFAS, "_", "-"))


# plot measured:
HL_measured_plot <- ggplot(HL_measured, aes(reorder(PFAS, -MW), log(Value),
                                            color = Reference)) +
  geom_jitter(#color = "grey70",
              width = 0.2,
              alpha = 0.9,
              size = .1) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
  labs(x = "",
       y = "Half life (years)") +
  coord_flip() +
  theme_quick


# plot measured & predicted:
HL_MeasPred_plot <- ggplot(HL_measured, aes(reorder(PFAS, -MW), log(Value))) +
  geom_jitter(color = "grey60",
    width = 0.2,
    alpha = 0.9,
    size = .1) +
  geom_errorbar(data = HL_predicted,
                aes(x = reorder(PFAS, -MW),
                    y = log(mean),
                    colour = Category_OECD,
                    ymin = log(q25),
                    ymax = log(q75)),
                width = .01,
                linewidth = .4) +
  geom_point(data = HL_predicted,
             aes(x = reorder(PFAS, -MW),
                 y = log(mean),
                 colour = Category_OECD),
             size = 1.5,
             shape = 17) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
  labs(x = "",
       y = "Half life (years)") +
  coord_flip() +
  theme_quick

HL_measured_plot <- HL_measured_plot + CP_theme +
  theme(legend.position = "right")
HL_MeasPred_plot <- HL_MeasPred_plot + CP_theme

ggsave(plot = HL_measured_plot,
       file = here(OUTPUT,"HL_measured_plot.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")

ggsave(plot = HL_MeasPred_plot,
       file = here(OUTPUT,"HL_MeasPred_plot.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")