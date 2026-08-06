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
    axis.text.x = element_text(size = 35), #, face = "bold"), #, angle = 30),
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

# Set input and output directory
INPUT <- here("Output", "test_MC")
OUTPUT <- here("Output","test_MC", "Processed_results") 
dir.create(OUTPUT, recursive = TRUE)

Dose_df <- read_excel(here("Input","Evaluation_data.xlsx")) %>%
  filter(Reference == "Abraham_2024") |>
  select(!c(PFAS_name, CAS_n)) %>%
  rename(PFAS = Acronym) |>
  mutate(Dose = round(Dose_ugPerkg_bw*82, 3)) |>  # 82 is the body weight of Abraham
  select(PFAS, Dose) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-")) |>
  distinct()


# functions ----

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
  AvDailyUR <- mean(A_urine$dailyFE) # average urinary excretion per day


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
  )
  )
}
calc_TK_summary_results <- function(raw_TK_results){

  Results_all <- do.call(rbind, raw_TK_results) |>
    na.omit() |> as.data.frame() |>
    summarise(across(where(is.numeric), .fns =
                       list(median = ~median(.,na.rm = TRUE),
                            mean = ~mean(.,na.rm = TRUE),
                            q25 = ~quantile(., 0.25, na.rm = TRUE),
                            q75 = ~quantile(.,0.75, na.rm = TRUE)))) |>
    pivot_longer(everything(), names_to = "Variable",
                 values_to = "Value")  |>
    separate(Variable,
             into = c("Parameter", "Variable"))|>
    pivot_wider(names_from = Variable, values_from = Value)

  return(Results_all)
}
calc_C_summary <- function(C_df, Tissue){

  C_df <- C_df  |> select(time, Tissue, PFAS, Iteration)
  C_wide_df <- C_df |> select(!PFAS) |>
    mutate(Iteration = paste0("n_", Iteration)) |>
    pivot_wider(names_from = Iteration, values_from = Tissue)
  C_summary <- C_wide_df |>
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
    mutate(Time = (Time*10)/365) |>  #unit correction, to years
    pivot_wider(names_from = Variable, values_from = Value)

  return(C_summary)
}
plot_per_tissue <- function(C_sum){

  Tissue_i <- unique(C_sum$Tissue)

  C_plot <- ggplot(C_sum, aes(Time, mean)) +
    geom_ribbon(aes(ymin = q25,
                    ymax = q75),
                fill = "grey70") +
    geom_path()+
    labs(title = paste0(PFAS_i),
         x = "Time (years)",
         y = paste0(Tissue_i, " Concentration (ng/ml)")) +
    theme_quick

  return(C_plot)

}
plot_all_tissues <- function(C_ts){

  C_ts_plot <- ggplot(C_ts, aes(Time, mean)) +
    geom_ribbon(aes(ymin = q25,
                    ymax = q75),
                fill = "grey70") +
    geom_path() +
    facet_wrap(~ Tissue, scales = "free", nrow = 2) +
    labs(title = paste0(PFAS_i),
         x = "Time (years)",
         y = "Concentration (ng/ml)") +
    theme_quick

  return(C_ts_plot)

}

PFAS_list <- c("ADONA",
               "HFPO-DA",
               "PFBA",
               "PFBS",
                "PFDA",
                "PFDoDA",
                "PFHpA",
               "PFHpS",
               "PFHxA",
               "PFHxS",
               "PFNA",
               "PFOA",
               "PFOS",
               "PFPeA",
               "PFPeS",
               "PFTrDA",
               "PFUnDA"
               )
nPFAS <- length(PFAS_list)

# produce result files and plots for all PFAS ----
for (i in 1:nPFAS) {


  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}------------------------------------"))

  # Inputs ----

  # load results (called savey)
  load(here(INPUT, paste0(PFAS_i, "_Monte_Carlo_y.RData" )
            ))
  # load oral dose (different per PFAS if comparing results to Abraham)
  OD <- filter(Dose_df, PFAS == PFAS_i) |>
    select(Dose) |> pull()

  # TK results ----
  # calculate AUC, half-life, elimination rates, fraction absorbed etc
  TK_results <- sapply(savey, calc_TK_results, OD)
  TK_summary_results <- calc_TK_summary_results(TK_results) |>
    as.data.frame() |> mutate(PFAS = PFAS_i)
  write.csv(TK_summary_results, file = here(OUTPUT, paste0(PFAS_i, "_TK_sumstat.csv")))

  # Predicted tissue concentrations ----
  Results_df <- do.call(rbind, savey) |> mutate(PFAS = PFAS_i)
  if(PFAS_i == "PFBA"){
    Results_df <- Results_df |> filter(CP < 0.0025) #if dose 0.1 -> 0.00025
  } else if(PFAS_i == "PFBS"){
    Results_df <- Results_df |> filter(CP < 0.4) #if dose 0.1 -> 0.004
  } else if(PFAS_i == "PFHxA"){
    Results_df <- Results_df |> filter(CP < 0.20) # if dose 0.1 -> 0.02
  } else if(PFAS_i == "PFPeA"){
    Results_df <- Results_df |> filter(CP < 0.2) # if dose 0.1 -> 0.02
  } else {
    Results_df <- Results_df
  }

  Tissues_list <- c("CI", "CL_ic", "CPTT", "CA", "CR", "CP")
  Tissues_names <- c("Intestine", "Liver IC", "Proximal tubule", "Adipose", "Rest", "Plasma")
  # plasma concentration
  nTissues <- length(Tissues_list)

  C_Tissues_list <- NULL

  for (j in 1:nTissues) {

    Tissue_j <- Tissues_list[j]#"PFOA"
    Tissue_j_name <- Tissues_names[j]

    C_summary <- calc_C_summary(Results_df, Tissue_j) |>
      mutate(Tissue = Tissue_j_name,
             PFAS = PFAS_i)

    # C_plot <- plot_per_tissue(C_summary)
    # C_plot <- C_plot + CP_theme +
    #   theme(axis.text.x = element_text(angle = 0))
    #
    # ggsave(plot = C_plot,
    #        file = here(OUTPUT, paste0(PFAS_i, Tissue_j, "_C_plot.png")),
    #        dpi = 600,
    #        width = 11,
    #        height = 9,
    #        units = "cm")

    C_Tissues_list[[j]] <- C_summary

  }

  C_Tissues <- do.call(rbind, C_Tissues_list)

  C_Tissues_plot <- plot_all_tissues(C_Tissues)

  C_Tissues_plot <- C_Tissues_plot + CP_theme +
    theme(axis.text.x = element_text(angle = 0))

  ggsave(plot = C_Tissues_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_C_Tissues_plot.png")),
         dpi = 600,
         width = 16,
         height = 11,
         units = "cm")

  write.csv(C_Tissues, file = here(OUTPUT, paste0(PFAS_i, "_C_Tissues_sumstat.csv")))

  # Actual parameters ----

  Parameters_df <- Results_df |>
    select(CF_I:PFAS) |> select(!CL_menses) |>
    distinct()

  P_long_df <- Parameters_df |> select(!c(PFAS, Iteration)) |>
    pivot_longer(everything(),
                 names_to = "Parameter",
                 values_to = "Value") |>
    mutate(Parameter = factor(Parameter,
                              levels = c(
                                "CL_IL",
                                "CF_I",
                                "CL_OATP2B1_I",
                                "CF_L_in",
                                "CL_OATP1B1",
                                "CL_OATP1B3",
                                "CL_OATP2B1_L",
                                "CL_NTCP",
                                "CF_L_out",
                                "CL_BSEP",
                                "CF_K_in",
                                "CL_OAT1",
                                "CL_OAT3",
                                "CF_K_out",
                                "CL_OAT4")))

  Parameters_plot <- ggplot(P_long_df, aes(x = Value)) +
    geom_histogram() +
    facet_wrap(~ Parameter, scales = "free", nrow = 3) +
    labs(x = "Parameter value",
         y = "Frequency")
    theme_quick

  Parameters_plot = Parameters_plot + CP_theme +
    theme(axis.text.x = element_text(angle = 0))

  ggsave(plot = Parameters_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_C_Parameters_plot.png")),
         dpi = 600,
         width = 16,
         height = 11,
         units = "cm")

  write.csv(Parameters_df, file = here(OUTPUT, paste0(PFAS_i, "_Parameters.csv")))


}

# run summary results ----

## load input files ----
# nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>  
#   mutate(PFAS = str_replace_all(PFAS, "_", "-"))
# 
# TK_df <-
#   list.files(path = here(INPUT),
#              pattern = "\\_TK_sumstat.csv$") |>
#   map_df(~ read_csv(file = here(INPUT, .))) 
# 
# # C_df <-
# #   list.files(path = here(INPUT),
# #              pattern = "\\_C_Tissues_sumstat.csv$") |>
# #   map_df(~ read_csv(file = here(INPUT, .))) 
# 
# Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |> 
#   rename(PFAS = Acronym) |> 
#   mutate(PFAS = str_replace_all(PFAS, "_", "-"))

## Plasma ----
# Measured_Abraham_df <- Measured_df |> 
#   filter(Reference == "Abraham_2024")
# 
# C_plasma_measured_df <- Measured_df |> 
#   filter(Parameter == "C_plasma") |> 
#   mutate(Time = as.numeric(Time_days)/367) |> 
#   select(PFAS, Time, Value) 
# 
# C_plasma_predicted_df <- C_df |> 
#   filter(Tissue == "Plasma", 
#          Time < 1.5)
# 
# ObsPred_plasma <- ggplot() +
#   geom_ribbon(data = C_plasma_predicted_df, 
#               aes(x = Time*365,
#                   ymin = q25,
#                   ymax = q75),
#               fill = "grey80") +
#   geom_path(data = C_plasma_predicted_df,
#             aes(Time*365, median)) +
#   geom_point(data = C_plasma_measured_df,
#              aes(Time*365, Value),
#              color = "purple",
#              size = 0.5) +
#   facet_wrap(~PFAS, scales = "free_y") +
#   labs(x = "Time (days)",
#        y = "Plasma concentration (μg/L)") +
#   theme_quick
# 
# ObsPred_plasma <- ObsPred_plasma + CP_theme
# ggsave(plot = ObsPred_plasma,
#               file = here(OUTPUT, "ObsPred_plasma_plot.png"),
#               dpi = 600,
#               width = 18,
#               height = 12,
#               units = "cm")
# 

## Urinary and fecal clearances ----
# 
# CL_UrFe_measured_df <- Measured_df |> 
#   filter(Parameter %in% c("CL_urinary", # in ml/min
#                           "CL_fecal"))|> 
#   mutate(bValue = ifelse(Unit == "ml/day/kg",
#                          (Value/1000) * 82,#L/day
#                          ifelse(Unit == "ml/min",
#                                 (Value/1000)*60*24,
#                                 ifelse(Unit == "L/day/kg",
#                                        Value * 82,
#                                        Value)
#                          )),
#          bUnit = "L/day")|> 
#   select(!c(Value, Unit)) |> 
#   rename(Value = bValue, 
#          Unit = bUnit)
# 
# CL_UrFe_measured_df <- full_join(CL_UrFe_measured_df, nC_MW_df) |> 
#   select(PFAS, Category_OECD, 
#          Parameter,Value, Unit, 
#          MW, n_C, 
#          Reference) |> 
#   na.omit()
# 
# CL_UrFe_predicted_df <- TK_df |> 
#   filter(Parameter %in% c(
#     "First10daysFE",
#     "AvDailyFE",
#     "First10daysUR",
#     "AvDailyUR"
#   )) |> 
#   mutate(mean = ifelse(Parameter == "First10daysFE" | 
#                          Parameter == "First10daysUR",
#                        mean*0.1, mean),
#          q25 = ifelse(Parameter == "First10daysFE" | 
#                         Parameter == "First10daysUR",
#                       q25*0.1, q25),
#          q75 = ifelse(Parameter == "First10daysFE" | 
#                         Parameter == "First10daysUR",
#                       q75*0.1, q75),
#     When = ifelse(Parameter == "First10daysFE" | 
#                          Parameter == "First10daysUR",
#                        "First_10", "Average")) |> 
#   mutate(Parameter = ifelse(Parameter == "First10daysFE" | 
#                               Parameter ==  "AvDailyFE",
#                             "CL_fecal",
#                             "CL_urinary"
#                               ))
# 
# CL_UrFe_predicted_df <- full_join(CL_UrFe_predicted_df, nC_MW_df) |> 
#   na.omit()
# 
# CL_measured_plot <- ggplot() +
#   geom_point(data = CL_UrFe_measured_df,
#            aes(x = reorder(PFAS, -MW), y = log10(Value), 
#                shape = Parameter,
#                color = Reference)) +
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
#   coord_flip() +
#   labs(y = "Clearance log10(L/day)", 
#        x = "") +
#   theme_quick
# CL_measured_plot
# 
# CL_measured_plot <- CL_measured_plot + CP_theme +
#   theme(legend.position = "right")
# 
# ggsave(plot = CL_measured_plot,
#        file = here(OUTPUT,"CL_measured_log10_plot.png"),
#        dpi = 600,
#        width = 12,
#        height = 12,
#        units = "cm")
# 
# CL_Pred_plot <- 
#   ggplot() +
#   geom_pointrange(data = CL_UrFe_predicted_df,
#                   aes(x = reorder(PFAS, -MW),
#                       y = log10(mean), 
#                       ymin = log10(q25),
#                       ymax = log10(q75),
#                       shape = Parameter,
#                       colour = Category_OECD),
#                   linewidth = .4,
#                   size = .5
#   ) +
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
#   facet_wrap(~When, scales = "free_x") +
#   coord_flip() +
#   labs(y = "Clearance log10(L/day)",
#        x = "") +
#   theme_quick
# CL_Pred_plot
# 
# 
# CL_Pred_plot <- CL_Pred_plot + CP_theme +
#   theme(legend.position = "right")
# 
# ggsave(plot = CL_Pred_plot,
#        file = here(OUTPUT,"CL_Pred_log10_plot.png"),
#        dpi = 600,
#        width = 12,
#        height = 12,
#        units = "cm")
# 
# 
# CL_UrFe_measured_df <- CL_UrFe_measured_df |> 
#   mutate(When = ifelse(Reference == "Abraham_2024",
#                        "First_10", "Average")) 
# 
# CL_MeasPred_plot <- ggplot() +
#   geom_jitter(data = CL_UrFe_measured_df,
#              aes(x = reorder(PFAS, -MW), y = log10(Value), 
#                  shape = Parameter),
#              color = "grey60",
#              width = 0.2,
#              alpha = 0.9,
#              size = 1) +
#   geom_pointrange(data = CL_UrFe_predicted_df,
#                   aes(x = reorder(PFAS, -MW),
#                       y = log10(mean), 
#                       ymin = log10(q25),
#                       ymax = log10(q75),
#                       shape = Parameter,
#                       colour = Category_OECD),
#                   linewidth = .4,
#                   size = .3
#   ) +
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
#   # facet_wrap(~When, scales = "free_x") +
#   coord_flip() +
#   labs(y = "Clearance log10(L/day)", 
#        x = "") +
#   theme_quick
# CL_MeasPred_plot
# 
# CL_MeasPred_plot <- CL_MeasPred_plot + CP_theme +
#   theme(legend.position = "right")
# 
# ggsave(plot = CL_MeasPred_plot,
#        file = here(OUTPUT,"CL_MeasPred_plot_log10_plot.png"),
#        dpi = 600,
#        width = 12,
#        height = 12,
#        units = "cm")

## Half life ----

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
# HL_Measured_plot <- HL_Measured_plot + CP_theme +
#   theme(legend.position = "right")
# 
# ggsave(plot = HL_Measured_plot,
#        file = here(OUTPUT,"HL_measured_plot.png"),
#        dpi = 600,
#        width = 12,
#        height = 12,
#        units = "cm")
# 
# 
# HL_predicted_df <- TK_df |> 
#   filter(Parameter == "HL")
# 
# HL_predicted_df <- full_join(HL_predicted_df, nC_MW_df)
# 
# 
# HL_predicted_plot <- ggplot() +
#   geom_errorbar(data = HL_predicted_df,
#                 aes(x = reorder(PFAS, -MW),
#                     y = log10(mean),
#                     colour = Category_OECD,
#                     ymin = log10(q25),
#                     ymax = log10(q75)),
#                 width = .01,
#                 linewidth = .4) +
#   geom_point(data = HL_predicted_df,
#              aes(x = reorder(PFAS, -MW),
#                  y = log10(mean),
#                  colour = Category_OECD),
#              size = 1.5,
#              shape = 17) +
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
#   labs(x = "",
#        y = "Half life (years)") +
#   coord_flip() +
#   theme_quick
# HL_predicted_plot
# 
# HL_predicted_plot <- HL_predicted_plot + CP_theme +
#   theme(legend.position = "right")
# 
# ggsave(plot = HL_predicted_plot,
#        file = here(OUTPUT,"HL_predicted_plott.png"),
#        dpi = 600,
#        width = 12,
#        height = 12,
#        units = "cm")
# 
# 
# HL_MeasPred_plot <- ggplot(HL_measured_df, 
#                            aes(reorder(PFAS, -MW),
#                                log10(Value))) +
#   geom_jitter(color = "grey60",
#               width = 0.2,
#               alpha = 0.7,
#               size = .1) +
#   geom_errorbar(data = HL_predicted_df,
#                 aes(x = reorder(PFAS, -MW),
#                     y = log10(mean),
#                     colour = Category_OECD,
#                     ymin = log10(q25),
#                     ymax = log10(q75)),
#                 width = .01,
#                 linewidth = .4) +
#   geom_point(data = HL_predicted_df,
#              aes(x = reorder(PFAS, -MW),
#                  y = log10(mean),
#                  colour = Category_OECD),
#              size = 1.5,
#              shape = 17) +
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
# nC_MW_plot <- ggplot(nC_MW_df, aes(n_C, MW, color = Category_OECD))+
#   geom_point()+
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.9) +
#   theme_quick +
#   labs(x = "npf-C",
#        y = "MW (g/mol)")
# 
# nC_MW_plot <- nC_MW_plot + CP_theme +
#   theme(legend.position = "right")
# 
# ggsave(plot = nC_MW_plot,
#        file = here("Input","nC_MW_plot.png"),
#        dpi = 600,
#        width = 9,
#        height = 9,
#        units = "cm")
# 
# 
# 







