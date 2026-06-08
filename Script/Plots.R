
rm(list=ls())

# Packages
library(here)
library(tidyverse)
library(showtext)
# library(plotly)
library(readxl)
library(patchwork)
library(ggrepel)
library(scales)
library(viridis)
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
    
    axis.title = element_text(size = 35), #, face = "bold"
    axis.text = element_text(size = 32),
    # axis.text.y = element_blank(),
    axis.text.y = element_text(size = 25, margin = margin(0,-10,0,-10)),
    axis.text.x = element_text(size = 32, angle = 45, hjust = 1),
    axis.line = element_blank(), #element_line(linewidth = 0.5),
    axis.ticks = element_line(color = "grey90", linewidth = 0.1),
    axis.ticks.length = unit(0.05, "cm"),
    
    legend.title = element_blank(), 
    legend.text = element_text(size = 22), #30
    # legend.margin = margin(0, 0, 0, 0, "cm"),
    # legend.key.spacing = unit(0.1, "cm"),
    # legend.key.size = unit(0.01, "lines"),
    # legend.spacing.y = unit(0.1, "cm"),
    # legend.spacing.x = unit(0.1, "cm"),
    legend.position = "top",
    legend.direction = "vertical",
    legend.box = "horizontal", 
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "cm"), 
    legend.spacing.x = unit(0.2, "cm"),
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
    legend.margin = margin(0.5, 0.5, 0.5, 0.5),
    strip.text.x = element_text(size = 35)
  )

load(here("Output","laIVIVE_PredictedParams_POST.RUN_RESULTS.RData")) #Predicted
PFAS <- names(POST.RUN_RESULTS)
# list2env(POST.RUN_RESULTS, envir = .GlobalEnv) to separate the result list
PredParam_POST.RUN_RESULTS <- POST.RUN_RESULTS
Results_PredictedParams <- tibble(
  PFAS = PFAS, 
  AUC = lapply(names(PredParam_POST.RUN_RESULTS), function(x) {
    PredParam_POST.RUN_RESULTS[[x]]$ANALYSED_data$AUC
  }),
  Cmax = lapply(names(PredParam_POST.RUN_RESULTS), function(x) {
    PredParam_POST.RUN_RESULTS[[x]]$ANALYSED_data$Cmax
  }),
  Tmax = lapply(names(PredParam_POST.RUN_RESULTS), function(x) {
    PredParam_POST.RUN_RESULTS[[x]]$ANALYSED_data$Tmax
  }),
  Half_life = lapply(names(PredParam_POST.RUN_RESULTS), function(x) {
    PredParam_POST.RUN_RESULTS[[x]]$ANALYSED_data$HalfLife
  })
) %>% 
  pivot_longer(AUC:Half_life, names_to = "TK", values_to = "Value") %>% 
  separate(col = Value, into = c("Value", "Unit"), sep = "_") %>% 
  mutate(Parameter_source = "Predicted")

load(here("Output","laIVIVE_ObservedParams_POST.RUN_RESULTS.RData"))
ObsParam_POST.RUN_RESULTS <- POST.RUN_RESULTS
Results_ObservedParams <- tibble(
  PFAS = PFAS, 
  AUC = lapply(names(ObsParam_POST.RUN_RESULTS), function(x) {
    ObsParam_POST.RUN_RESULTS[[x]]$ANALYSED_data$AUC
  }),
  Cmax = lapply(names(ObsParam_POST.RUN_RESULTS), function(x) {
    ObsParam_POST.RUN_RESULTS[[x]]$ANALYSED_data$Cmax
  }),
  Tmax = lapply(names(ObsParam_POST.RUN_RESULTS), function(x) {
    ObsParam_POST.RUN_RESULTS[[x]]$ANALYSED_data$Tmax
  }),
  Half_life = lapply(names(ObsParam_POST.RUN_RESULTS), function(x) {
    ObsParam_POST.RUN_RESULTS[[x]]$ANALYSED_data$HalfLife
  })
) %>% 
  pivot_longer(AUC:Half_life, names_to = "TK", values_to = "Value") %>% 
  separate(col = Value, into = c("Value", "Unit"), sep = "_") %>% 
  mutate(Parameter_source = "Observed")

Results_df <- rbind(Results_PredictedParams, Results_ObservedParams) %>% 
  mutate(Value = as.numeric(Value)) %>% 
  mutate(Value = ifelse(TK == "Cmax", Value/0.04829268, Value)) #normalise Cmax by dose

ModelParameters_df <- read_xlsx(here("ModelParameters.xlsx")) %>% 
  mutate(Value = as.numeric(Value)) %>% 
  filter(Parameter != "MW")

ValidationData_df <- read_xlsx(here("NotOnGit","PFAS_validation_data.xlsx"))

MeasuredHL <- ValidationData_df %>% filter(Parameter == "Halflife")
HL <- MeasuredHL %>%
  # filter(Acronym == "PFOA") %>%
  select(c(Acronym, Value, Unit, "n Subjects", Category_OECD)) %>%
  rename(n = "n Subjects") %>% 
  mutate(n = ifelse(is.na(n), 1, n),
         Category_OECD = factor(Category_OECD, 
                                levels = c("PFCAs", "PFSAs", "PFECAs"))) %>% 
  arrange(desc(Category_OECD))

MeasuredCmax <- ValidationData_df %>% filter(Parameter == "Cmax_plasma")

Cmax <- MeasuredCmax %>%
  # filter(Acronym == "PFOA") %>%
  select(c(Acronym, Value, Unit, "n Subjects", Category_OECD, Reference, "Dose (μg/kg bw)")) %>%
  rename(n = "n Subjects",
         Dose = "Dose (μg/kg bw)") %>% 
  mutate(n = ifelse(is.na(n), 1, n),
         Category_OECD = factor(Category_OECD, 
                                levels = c("PFCAs", "PFSAs", "PFECAs")))

Cmax$Dose <- as.numeric(Cmax$Dose)

Cmax <- Cmax %>% mutate(
  nValue = Value/Dose) %>% 
  arrange(desc(Category_OECD))

Select_Cmax <- Cmax %>% select(Acronym, nValue) %>% 
  mutate(Value = nValue, 
         TK = "Cmax") %>% 
  select(!nValue)

Select_HL <- HL %>% select(Acronym, Value) %>% 
  mutate(TK = "Half_life")

Measured <- rbind(Select_Cmax, Select_HL) %>% 
  mutate(Parameter_source = "Measured",
         Source = "Measured") %>% 
  rename(PFAS = Acronym)

Results_df <- Results_df %>% mutate(Source = "Predicted")

Results_all <- rbind(Measured, select(Results_df, !Unit))



PFAS_list <- unique(Results_df$PFAS)
nPFAS <- length(PFAS_list)
Plot_list <- NULL

for(i in 1:nPFAS){
  
  PFAS_t <- PFAS_list[i]
  Results <- Results_all %>% filter(PFAS == PFAS_t) 
  ModelParameters <- ModelParameters_df %>% filter(Acronym == PFAS_t) 

  Results_plot <- ggplot(filter(Results, Source == "Predicted"), 
                         aes(x = TK, 
                             y = Value, 
                             color = Parameter_source, 
                             shape = Source)) +
    geom_point(data = filter(Results, Source == "Measured"),
               aes(x = TK,
                   y = Value,
                   shape = Source),
               alpha = 0.5, size = 3, color = "grey50") +
    geom_point(alpha = 0.9, position = position_dodge(width = 0.3), size = 2) +
    scale_color_viridis_d(option = "D",
                          labels = c("Observed" = "Predicted, measured parameters where available", 
                                     "Predicted" = "Predicted, only estimated parameters")) +
    scale_shape_manual(values = c("Measured" = 18, 
                                  "Predicted" = 17),
                       labels = c("Observed" = "Measured",
                                  "Predicted" = "Estimated")) +
    guides(color = guide_legend(override.aes = list(size = 0.5)),
           shape = guide_legend(override.aes = list(size = 0.5))
           ) +
    facet_wrap(~TK, scales = "free", nrow = 1) +
    scale_x_discrete(expand = FALSE) +
    labs(title = paste0(PFAS_t, ": predicted and measured TK data"),
         x = "",
         y = "") +
    CP_theme +
    theme(axis.text.x = element_blank())
  
  
  Params <- ggplot(ModelParameters, 
                   aes(x = Parameter, 
                       y = Value, 
                       colour = Name,
                       shape = Source)) +
    # geom_col(position = "dodge") +
    geom_point(alpha = 0.9, position = position_dodge(width = 0.3), size = 2) +
    scale_color_viridis_d(option = "D",
                          labels = c("Observed" = "Predicted, measured parameters where available", 
                                     "Predicted" = "Predicted, only estimated parameters"))+
    scale_shape_manual(values = c("Observed" = 19, 
                                  "Predicted" = 17),
                       labels = c("Observed" = "Measured parameter",
                                  "Predicted" = "Estimated parameter")) +
    guides(color = guide_legend(override.aes = list(size = 0.5)),
           shape = guide_legend(override.aes = list(size = 0.5))
    ) +
    facet_wrap(~TK_process, scales = "free", nrow = 1) +
    labs(title = paste0(PFAS_t, ": input parameters"),
         x = "",
         y = "") +
    CP_theme
  
  Plot_I <- Results_plot + Params + 
    plot_layout(widths = c(1,1))
  
  # ,
  #               guides = "collect") &
  #   theme(legend.position = 'right')

  Plot_list[[i]] <- Plot_I
  
}

for(i in 1:nPFAS){
  PFAS_t <- PFAS_list[i]
    ggsave(plot = Plot_list[[i]], file = paste("laIVIVE_Predictions_Params",PFAS_t,".png",sep=""),
     dpi = 600,
     width = 18,
     height = 9,
     units = "cm")
}



# HL_exp <- HL %>% 
#   uncount(n) # point for each count

# HL_plot <- ggplot() +
#   geom_point(data = HL,
#              aes(x = Acronym, 
#                  y = Value, colour = Category_OECD),
#              size = 0.5,
#              position = position_jitter(width = 0.1)) +
#   scale_colour_manual(values = c("PFCAs" = "#4ea9cf",
#                                  "PFSAs" = "#b45c73",
#                                  "PFECAs"  = "#55b483")) +
#   CP_theme +
#   labs(y = "Half life (years)",
#        x = "") +
#   theme(legend.title = element_blank(), 
#         legend.position = "bottom", 
#         axis.text.x = element_text(size = 32, angle = 45, hjust = 1))
# HL_plot
# ggsave(plot = HL_plot, file = "MeasuredHL.png",
#        dpi = 600,
#        width = 12,
#        height = 9,
#        units = "cm")

 
# Cmax_plot <- ggplot() +
#   geom_point(data = Cmax,
#              aes(x = Acronym, 
#                  y = nValue, colour = Category_OECD),
#              size = 2,
#              position = position_jitter(width = 0.1)) +
#   scale_colour_manual(values = c("PFCAs" = "#4ea9cf",
#                                  "PFSAs" = "#b45c73",
#                                  "PFECAs"  = "#55b483")) +
#   CP_theme +
#   labs(y = "dose-normalised Cmax (μg/L /μg/kg bw)",
#        x = "") +
#   theme(legend.title = element_blank(), 
#         legend.position = "bottom", 
#         axis.text.x = element_text(size = 32, angle = 45, hjust = 1))
# Cmax_plot
# ggsave(plot = Cmax_plot, 
#        file = "MeasuredCmax.png",
#        dpi = 600,
#        width = 12,
#        height = 9,
#        units = "cm")
