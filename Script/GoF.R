# --------------------------------------------------------------------------- #
# Script for running the Morris Test
# CP, 30-06-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Load packages
library(here)
library(ggplot2)
library(deSolve)
library(tidyverse)
library(sensitivity)
library(ggrepel)
library(scales)
library(patchwork)
library(pracma)
library(PKNCA)
library(glue)
library(readxl)

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
    
    axis.title = element_text(size = 35), #, face = "bold"
    axis.text = element_text(size = 35),
    # axis.text.y = element_blank(),
    # axis.text.x = element_text(size = 64, face = "bold"),
    axis.line = element_blank(), #element_line(linewidth = 0.5),
    axis.ticks = element_line(color = "grey90", linewidth = 0.1),
    axis.ticks.length = unit(0.05, "cm"),
    
    legend.title = element_text(size = 35), 
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


# Set output storage directory
OUTPUT <- here("Output", "Testing clearance values and model stability") 
# dir.create(OUTPUT, recursive = TRUE)

Evaluation_data <- read_excel(here("Input","Evaluation_data.xlsx")) %>% 
  select(!c(PFAS_name, CAS_n, Category_OECD)) %>% 
  rename(PFAS = Acronym)
Evaluation_data$PFAS = str_replace(Evaluation_data$PFAS, "HFPO-DA", "HFPO.DA")  

load(here(OUTPUT, "Results_mean_values.RData"))

HL_measured <- Evaluation_data %>% filter(Parameter == "Halflife") 

THL_measured <- HL_measured %>% 
  group_by(PFAS) %>% 
  summarise(mean = mean(Value), 
            median = median(Value)) %>% 
  mutate(HL_length = ifelse(mean < 0.01, "very_short",
                            ifelse(mean < 1, "short", 
                            ifelse(mean < 5, "medium", "long")))) %>% 
  select(PFAS, HL_length)

HL_measured <- full_join(HL_measured, THL_measured)

HL_predicted <- as.data.frame(Results$HL) %>% 
  pivot_longer(everything(),names_to = "PFAS", values_to = "Value")
HL_predicted <- full_join(HL_predicted, THL_measured)


Results <- ggplot(HL_measured, aes(PFAS, Value)) +
  geom_jitter(aes(colour = Reference), width = 0.2, size = 1) +
  geom_point(data = HL_predicted, shape = 3) + 
  facet_wrap(~ HL_length, scales = "free") +
  CP_theme +
  theme(axis.text.x = element_text(angle = 45))

Results


  ggsave(plot = Results, 
         file = here(OUTPUT,"Results.png"),
         dpi = 600,
         width = 16,
         height = 12,
         units = "cm")

