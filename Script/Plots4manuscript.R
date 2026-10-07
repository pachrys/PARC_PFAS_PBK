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
INPUT <- here("Output", "GoF_0709") 
OUTPUT <- here("Output", "Plots")
dir.create(OUTPUT, recursive = TRUE)

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
    
    legend.position = "right",
    legend.title = element_blank(), 
    legend.text = element_text(size = 35),
    # legend.margin = margin(0, 0, 0, 0, "cm"),
    legend.margin = margin(0.05, r = 0.05, l = -0.01, 0.05, "cm"),
    legend.key.spacing = unit(0.05, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.08, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
    legend.box.margin = margin(t = 0, 0, b = 0, l = -0.3, "cm"),
    
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

# Inputs ----
nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))


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



TK_df <- read.csv(here(INPUT, "TK_results_all_PFAS.csv"))

C_df <-
  list.files(path = here(INPUT),
             pattern = "\\_TissueC_sumstat.csv$") |>
  map_df(~ read_csv(file = here(INPUT, .)))


Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |>
  rename(PFAS = Acronym) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

# Plot 1: HL and CL ----
HL_predicted_df <- TK_df |>
  filter(Parameter == "HL")

HL_predicted_df <- full_join(HL_predicted_df, nC_MW_df) |> 
  filter (PFAS != "PFTrDA")

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

HL_measured_df <- full_join(HL_measured, nC_MW_df)|> 
  filter(PFAS != "PFTrDA")



CLurCLfe_predicted_df <- TK_df |> 
  filter(Parameter %in% c("CLfe", 
                          "CLur")) 

CLurCLfe_predicted_df <- full_join(CLurCLfe_predicted_df, nC_MW_df)

RurRfe_predicted_df <- TK_df |> 
  filter(Parameter %in% c("RfeEl", 
                          "RurEl"))
RurRfe_predicted_df <- full_join(RurRfe_predicted_df, nC_MW_df)



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


## Plots ----

shapes <- c(
 "Abraham. 2024" = 19,
 "Chiu. 2022"  = 21,
 "Delaere. 2025" = 13,
 "Li. 2022"   = 7,
 "Lyu. 2026"   = 24,
 "Nilsson. 2022" = 17,
 "Olsen. 2007"  = 4,
 "Olsen. 2009"  = 5,
 "Rosato. 2024" = 25,
 "Russell. 2013" = 3,
 "Worley. 2017" = 14,
 "Xu. 2020"   = 12,
 "Yao. 2023"   = 8,
 "Zhang. 2013"  = 12,
 "Andersson. 2025" = 9,
 "Fujii. 2015"   = 15
)

HL_measured_df <- HL_measured_df |> 
  mutate(Reference = str_replace_all(Reference, "_", ". "))
CL_UrFe_measured_df <- CL_UrFe_measured_df |> 
  mutate(Reference = str_replace_all(Reference, "_", ". "))

HL_MesPred_plot <- 
  ggplot() +
  geom_jitter(data = HL_measured_df,
              aes(x = reorder(PFAS, MW),
                  y = Value,
                  shape = Reference
                  ),
              color = "grey70",
              stroke = 0.2,
              width = 0.2,
              size = .8) +
  geom_boxplot(data = HL_predicted_df,
               aes(x = reorder(PFAS, MW), 
                   y = median,
                   ymin = q025, lower = q25, 
                   middle = median, 
                   upper = q75, ymax = q975,
                   fill = Category_OECD, 
                   colour = Category_OECD),
               stat = "identity",
               linewidth = 0.3,
               width = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.6, alpha = 0.8) +
  scale_colour_viridis(option = "D", discrete = TRUE, end = 0.6) +
  scale_shape_manual(values = shapes)+
  scale_y_log10() +
  labs(x = "",
       y = "Half life (years) ( scale)") +
  theme_quick
HL_MesPred_plot


CL_MesPred_plot <- ggplot() +
  geom_point(data = CL_UrFe_measured_df,
             aes(x = reorder(PFAS, MW), 
                 y = Value,
                 shape = Reference),
             color = "grey70",
             stroke = 0.2,
             size = .8) +
  geom_boxplot(data = CLurCLfe_predicted_df,
               aes(x = reorder(PFAS, MW), 
                   y = median,
                   ymin = q025, lower = q25, 
                   middle = median, 
                   upper = q75, ymax = q975,
                   fill = Category_OECD, 
                   colour = Category_OECD),
               stat = "identity",
               linewidth = 0.3,
               width = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.6, alpha = 0.8) +
  scale_colour_viridis(option = "D", discrete = TRUE, end = 0.6) +
  scale_shape_manual(values = shapes) +
  facet_wrap(~ Parameter, nrow = 2,
             labeller = as_labeller(c("CLfe" = "Biliary/ fecal",
                                      "CLur" = "Renal")),
             scales = "free_y") +
  scale_y_log10() +
  labs(x = "",
       y = "Clearance (L/day) ( scale)")+
  theme_quick
  CL_MesPred_plot

## Print ---- 
 
  GoF_one <- HL_MesPred_plot + 
    labs(title = "A. Half-life",
         x = "",
         y = "Half-life (years) (log10 scale)") 
  # +
  #   guides(colour = FALSE,
  #          fill = FALSE) 
  # 
  GoF_two <- CL_MesPred_plot + 
    labs(title = "B. Clearances" ,
         x = "",
         y = "Clearance (L/day) (log10 scale)") #+
    # guides(shape = FALSE) 
  
  # GoF_plot <- GoF_one + GoF_two +
  #   plot_layout(guides = "collect") &
  #   Patchwork_theme &
  #   theme(plot.margin = margin(0.1, 0.1, b = 0, 0.1, "cm"),
  #         legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
  #         legend.box.margin = margin(t = -0.5, 0, b = 0.1, 0, "cm"), 
  #         legend.position = "right",
  #         axis.text.x = element_text(angle = 45, hjust = 1)
  #         )
  # 
  # ggsave(plot = GoF_plot,
  #        file = here(OUTPUT,"Elimination.png"),
  #        dpi = 600,
  #        width = 18,
  #        height = 9,
  #        units = "cm")
  
  
  
# Plot 2: plasma C ----  
  
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
  
  Measured_Abraham_df <- Measured_df |>
    filter(Reference == "Abraham_2024")
  
  C_plasma_measured_df <- Measured_df |>
    filter(Parameter == "C_plasma") |>
    mutate(Time = as.numeric(Time_days)/365) |>
    select(PFAS, Time, Value)
  C_plasma_measured_df <- full_join(C_plasma_measured_df, nC_MW_df)
  C_plasma_measured_df <- full_join(C_plasma_measured_df, Dose_df) |> 
    filter_out(PFAS %in% c("PFTrDA", "PFPeS", "PFHpS"))
  
  
  C_plasma_predicted_df <- C_df |>
    filter(Tissue == "Plasma",
           Time < 1.25)
  
  C_plasma_predicted_df <- full_join(C_plasma_predicted_df, nC_MW_df)
  C_plasma_predicted_df <- full_join(C_plasma_predicted_df, Dose_df)|> 
    filter_out(PFAS %in% c("PFTrDA", "PFPeS", "PFHpS"))
  
  ## Plot ----
  
  ObsVsPred_Plasma_plot <- ggplot() +
    geom_ribbon(data = C_plasma_predicted_df,
                aes(x = Time,
                    ymin = q025,
                    ymax = q975,
                    fill = Category_OECD),
                alpha = 0.4) +
    geom_ribbon(data = C_plasma_predicted_df,
                aes(x = Time,
                    ymin = q25,
                    ymax = q75,
                    fill = Category_OECD),
                alpha = 0.3) +
    scale_fill_viridis(option = "D", discrete = TRUE, end = 0.6) +
    geom_path(data = C_plasma_predicted_df,
              aes(Time, median,
                  colour = Category_OECD)) +
    scale_color_viridis(option = "D", discrete = TRUE, end = 0.6) +
    geom_point(data = C_plasma_measured_df, #filter(C_plasma_measured_df, Category_OECD == "PFCAs"),
               aes(Time, Value),
               colour = "grey50",
               size = .5,
               stroke = .2) +
    labs(x = "Time (years)",
         y = "Concentration (ng/ml)") +
    facet_wrap(~ factor(PFAS, rev(PFAS_list)),
               nrow = 3,
               scales = "free_y") +
    theme_quick
  ObsVsPred_Plasma_plot
  
  ## Print ----
  
  GoF_three <- ObsVsPred_Plasma_plot + 
    labs(title = "C. Plasma concentration over time curves") +
    Patchwork_theme &
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
          legend.box.margin = margin(t = -0.5, 0, b = 0.1, 0, "cm"), 
          legend.position = "right"#,
          # axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # ggsave(plot = GoF_three,
  #        file = here(OUTPUT,"Plasma_long.png"),
  #        dpi = 600,
  #        width = 14,
  #        height = 11,
  #        units = "cm")
  
  # Print all ---
  
  GoF_cl <- GoF_one + GoF_two &
    Patchwork_theme &
    theme(axis.text.x = element_text(angle = 45, hjust = 1)
          )
  
  GoF_three <- ObsVsPred_Plasma_plot + 
    labs(title = "C. Plasma concentration over time curves") +
    Patchwork_theme 
  
  GoF <- GoF_cl / GoF_three 
    # plot_layout(guides = "collect") +
    # plot_annotation(title = "NG-PBK predicted and measured PFAS toxicokinetics align",
    #                 theme = theme(plot.title = element_text(size = 52, 
    #                                                         face = "bold",
    #                                                         hjust = 0.5)
    #                               )
    #                 ) 
  
  ggsave(plot = GoF,
         file = here(OUTPUT,"GoF2.png"),
         dpi = 600,
         width = 18,
         height = 18,
         units = "cm") 
  
  
  # Plot 3: fractions ----

  FelFa_predicted_df <- TK_df |>
    filter(Parameter %in% c("Fa",
                            # "FelTotal",
                            "Ffeces",
                            "FbilEl",
                            "Furine")) |> 
    select(Parameter, median, PFAS) |>
    pivot_wider(names_from = Parameter, 
                values_from = median) |> 
    mutate(Fpreintfeces = Ffeces - FbilEl,
           cumFpreintfeces = Fpreintfeces + Furine,
           cumFbile = Fpreintfeces + Furine + FbilEl) |> 
    select(!c(Fa, Ffeces)) |> 
    pivot_longer(FbilEl:cumFbile, 
                 names_to = "Parameter", 
                 values_to = "Value")
  
  FelFa_predicted_df <- full_join(FelFa_predicted_df, nC_MW_df)
  
  Clean_FelFa <- FelFa_predicted_df |>
    filter(Parameter %in% c("FbilEl", "Fpreintfeces", "Furine"))
  Cum_FelFa <- FelFa_predicted_df |>
    filter(Parameter %in% c("cumFbile", "cumFpreintfeces", "Furine")) |>
    mutate(Parameter = str_replace_all(Parameter, "cum", ""),
           Parameter = str_replace_all(Parameter, "Fbile", "FbilEl"))

  Fe_plot <- ggplot() +
    geom_bar(data = Clean_FelFa,
             aes(x = reorder(PFAS, MW), y = Value,
                 fill = Parameter),
             stat = "identity",
             position = "stack",
             width = 0.5) +
    scale_fill_viridis(option = "D", direction = 1, discrete = "TRUE", 
                       end = 0.8, 
                       labels = c("FbilEl" = "Feces: bile",
                                  "Fpreintfeces" = "Feces: pre-systemic",
                                  "Furine" = "Urine")) +
    labs(x = "",
         y = "Cumulative fraction excreted") +
    theme_quick +
    theme(legend.position = "right")
  Fe_plot

  # Plot 4: distributions ----
  
  INPUT <- here("Output", "0.00019")
  
  Cdistr_df <-
    list.files(path = here(INPUT),
               pattern = "\\Tissue_Distribution.csv$") |>
    map_df(~ read_csv(file = here(INPUT, .))) 
  
  Cdistr_df <- Cdistr_df |> 
    filter(time == max(Cdistr_df$time),
           Tissue != "Plasma_total") |> 
    mutate(SepTissue = Tissue) |> 
    separate(SepTissue, into = c("Tissue", "Component")) |> 
    mutate(Component = str_replace_all(
      Component,
      c("ALB" = "albumin",
        "GLOB" = "globulin",
        "ML" = "membrane lipids",
        "SL" = "storage lipids",
        "SP" = "structural proteins")
    )) |> 
    mutate(Component = factor(Component, 
                              levels = c("albumin",
                                         "globulin",
                                         "structural proteins",
                                         "storage lipids", 
                                         "membrane lipids",
                                         "free")))

  Cdistr_df <- full_join(Cdistr_df, nC_MW_df)
  
  
  ### liver----
  
  Liver_distri <- filter(Cdistr_df, Tissue == "Liver") |> 
    ggplot(aes(reorder(PFAS, MW), Concentration, fill = Component)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_viridis(option = "D", discrete = TRUE,
                       end = 0.9) +
    scale_y_continuous(labels = scales::percent) +
    labs(#title = "Distribution in liver",
      x = "",
      y = "Mass balance", 
      fill = "Tissue component") +
    theme_quick
  Liver_distri
  
  
  ### kidney ----
  Kidney_distri <- filter(Cdistr_df, Tissue == "Kidney") |> 
    ggplot(aes(reorder(PFAS, MW), as.numeric(Concentration), fill = Component)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_viridis(option = "D", discrete = TRUE,
                       end = 0.9) +
    # facet_wrap(~factor(PFAS, levels = rev(order))) +
    scale_y_continuous(labels = scales::percent) +
    labs(#title = "Distribution in kidney",
      x = "",
      y = "Mass balance", 
      fill = "Tissue component")+
    theme_quick
  Kidney_distri
  
  
  ### plasma ----
  Plasma_distri <- filter(Cdistr_df, Tissue == "Plasma") |> 
    ggplot(aes(reorder(PFAS, MW), Concentration, fill = Component)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_viridis(option = "D", discrete = TRUE,
                       end = 0.9) +
    # facet_wrap(~factor(PFAS, levels = rev(order))) +
    scale_y_continuous(labels = scales::percent) +
    labs(#title = "Distribution in plasma",
      x = "",
      y = "Mass balance", 
      fill = "Tissue component")+
    theme_quick
  Plasma_distri
  
  # All tissues ----
  
  Tdistr_df <-
    list.files(path = here(INPUT),
               pattern = "\\Tissue_Cs.csv$") |>
    map_df(~ read_csv(file = here(INPUT, .))) 
  
  Tdistr_df <- Tdistr_df |> 
    mutate(Tissue = str_replace_all(
      Tissue,
      c("CA" = "adipose",
        "CI" = "intestine",
        "CL" = "liver",
        "CK" = "kidney",
        "CR" = "rest",
        "CP" = "plasma")
    )) |> 
    mutate(Tissue = factor(Tissue, 
                           levels = c("liver",
                                      "intestine", 
                                      "kidney", 
                                      "rest",
                                      "plasma",
                                      "adipose"
                                      )))
  Tdistr_df <- full_join(Tdistr_df, nC_MW_df)
  
  ### Plasma ----
  
  # install.packages("geomtextpath")
  library(geomtextpath)
  library(ggrepel)
  data4ss_plot = filter(Tdistr_df,
                        Tissue == "plasma") 
  
  data4ss_plot_label <- data4ss_plot |> 
    mutate(label = ifelse(time == max(data4ss_plot$time), PFAS, NA),
           nudge_x = 1)|> 
    mutate(PFAS = factor(PFAS, levels = rev(PFAS_list)),
           time = time/365)
  
  Plasma_perCat_Cs <- ggplot(data = data4ss_plot_label, 
                             aes(time, Concentration, 
                                 colour = PFAS)) + #factor(PFAS, levels = rev(PFAS_list)))) +
    geom_path(linewidth = 0.1, alpha = 0.8) +
    geom_text_repel(aes(label = label),
                    min.segment.length = 0.1,
                    show.legend = FALSE,
                    size = unit(9, "mm"),
                    segment.size = 0.1,
                    nudge_x = 1.5,
                    direction = "y",
                    hjust = 0,
                    na.rm = TRUE) +
    scale_colour_viridis(option = "D", discrete = TRUE,
                         direction = -1,
                         end = 0.9) +
    # facet_wrap(~Category_OECD, #factor(PFAS, levels = rev(PFAS_list))
    #            # scales = "free_y",
    #            nrow = 1) +
    labs(x = "Time (years)",
         y = "Concentration (μmol/L)") +
    theme_quick +
    theme(legend.position = "none")
  Plasma_perCat_Cs

  
  
  ### Tissues ----
  
  lastTdistr_df <- Tdistr_df |> 
    filter(time == max(Tdistr_df$time))
  
  Tdistr_plot <- lastTdistr_df |> 
    ggplot(aes(reorder(PFAS, MW), Concentration, fill = Tissue)) + #factor(PFAS, levels = rev(order))
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_viridis(option = "D", discrete = TRUE
    ) +
    # facet_wrap(~factor(PFAS, levels = rev(order))) +
    scale_y_continuous(labels = scales::percent) +
    labs(#title = "Distribution in tissues",
         x = "",
         y = "Mass balance", 
         fill = "Tissue")+
    theme_quick
  Tdistr_plot 
  
  ### Print ----
  
  #### left A ----
  
  LA_l <- Plasma_perCat_Cs + Patchwork_theme +
    theme(legend.position = "null")
  LA_r <- Plasma_distri + Patchwork_theme +
    theme(plot.margin = margin(0.05, 0.05, b = 0, 0.05, "cm"),
          legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
          legend.box.margin = margin(t = -0.5, r = -0.5, b = 0.1, l = -0.5, "cm"),
          legend.position = c(0.89, 0.915),
          axis.text.x = element_text(angle = 45, hjust = 1)
          )

  LA <- LA_l + LA_r +
    plot_layout(widths = c(1.3, 0.7))
  
  ggsave(plot = LA,
         file = here("ManuscriptPlots","LA.png"),
         dpi = 600,
         width = 18,
         height = 11,
         units = "cm") 
  
  #### right B ----
  
  RB <- Tdistr_plot + Patchwork_theme +
    theme(plot.margin = margin(0.05, 0.05, b = 0, 0.05, "cm"),
          legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
          legend.box.margin = margin(t = -0.5, r = -0.5, b = 0.1, l = -0.5, "cm"),
          legend.position = c(0.962, 0.88),
          axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(plot = RB,
         file = here("ManuscriptPlots","RB.png"),
         dpi = 600,
         width = 18,
         height = 11,
         units = "cm") 
  
  ### right C ----
  
  RC <- Liver_distri + Kidney_distri & 
    Patchwork_theme &
    theme(plot.margin = margin(0.05, 0.05, b = 0, 0.05, "cm"),
          legend.text = element_text(size = 28),
          legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
          legend.box.margin = margin(t = -0.5, r = -0.5, b = 0.1, l = -0.5, "cm"),
          legend.position = c(0.11, 0.9),
          legend.key.spacing = unit(0.02, "cm"),
          legend.key.size = unit(0.4, "lines"),
          legend.spacing.y = unit(0.06, "cm"),
          legend.spacing.x = unit(0.06, "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1)
    )
    
  ggsave(plot = RC,
         file = here("ManuscriptPlots","RC.png"),
         dpi = 600,
         width = 18,
         height = 11,
         units = "cm")   
    
  ### left D ----
  
  LD <- Fe_plot +  Patchwork_theme +
    theme(plot.margin = margin(0.05, 0.05, b = 0, 0.05, "cm"),
          legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
          legend.box.margin = margin(t = -0.5, r = -0.5, b = 0.1, l = -0.5, "cm"),
          legend.position = c(0.934, 0.92),
          axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(plot = LD,
         file = here("ManuscriptPlots","LD.png"),
         dpi = 600,
         width = 18,
         height = 11,
         units = "cm") 
    
    
    
    
    
    
    
    
    