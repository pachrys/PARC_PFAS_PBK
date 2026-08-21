# Produce parameter plots for all PFAS together ---
# Find correlations between HL and clearance and parameters ---


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

INPUT <- here("Output", "Removing_fupfromGFR")
OUTPUT <- INPUT #here(INPUT, "Correlations")

nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

P_df <- list.files(path = here(INPUT),
                               pattern = "\\_final_model_parameters.csv$") |>
  map_df(~ read_csv(file = here(INPUT, .),
                    id = "name"))|> 
  select(!"...2") |> 
  # rename(Iteration = "...2") |> 
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/Admin/Desktop/PhD_Chrysa/CP_L_R/PARC_PFAS_PBK/Output/Removing_fupfromGFR/",
                                "")) |> 
  mutate(PFAS = str_replace_all(PFAS,
                                "_final_model_parameters.csv","")) |>
  select(!name) #|> unite("PFAS_Iteration", c(PFAS,Iteration), sep = "_") 


# TK_df <- list.files(path = here(INPUT),
#                     pattern = "\\_TK_sumstat.csv$") |>
#   map_df(~ read_csv(file = here(INPUT, .),
#                     id = "name"))|> 
#   select(!"...2") |> 
#   # rename(Iteration = "...2") |> 
#   mutate(PFAS = str_replace_all(name,
#                                 "C:/Users/Admin/Desktop/PhD_Chrysa/CP_L_R/PARC_PFAS_PBK/Output/Re_analyse_UA/",
#                                 "")) |> 
#   mutate(PFAS = str_replace_all(PFAS,
#                                 "_TK_sumstat.csv","")) |>
#   select(!name) #|> 
#   # unite("PFAS_Iteration", c(PFAS,Iteration), sep = "_") 

# TK_df <- full_join(nC_MW_df, TK_df)

# plot the actual model parameters ----
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


Key_params_plot = Key_params_plot + CP_theme +
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = Key_params_plot,
       file = here(OUTPUT, "All_PFAS_key_params_plot.png"),
       dpi = 600,
       width = 16,
       height = 12,
       units = "cm")



# check all parameters ----
Visualise_P <- P_df |>  
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
