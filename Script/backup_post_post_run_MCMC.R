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







