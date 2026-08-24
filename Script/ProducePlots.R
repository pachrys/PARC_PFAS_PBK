# To compare initial and scaled parameters ----

rm(list=ls()) # to clear out the global environment


# Load packages
library(here)
library(ggplot2)
library(tidyverse)
library(glue)
library(readxl)
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

# Set output storage directory
OUTPUT <- here("Output", "RunTheUA") 
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

nPFAS <- length(PFAS_list)

Final_P <- read_csv("FinalModelParameters_allPFAS.csv") |> select(!"...1")
Initial_P <- read_excel("Input/MC_TK_Parameters_run_inThefirstUA.xlsx")


SA_SI <-1099.557*64.41978
VI_L <- 0.03848451 * 64.41978      
VL_ic <- 0.03848451 * 64.41978 * 0.835               
VPT <- 0.004630086 * 64.41978 * 0.398


Selected_final_P <- Final_P |> 
  filter(Parameter %in% c(
    "fuPTec",
    "fuLec",
    "fuPTL",
    "SFOATP2B1I",
    "SFOATP2B1",
    "SFOATP1B3",
    "SFOATP1B1", 
    "SFNTCP",    
    "SFBSEP",
    "SFOAT1",           
    "SFOAT3",
    "SFOAT4",
    "Clearapparent",
    "ClearOATP2B1I",
    "ClearuptakeI",
    "ClearOATP1B1",
    "ClearOATP1B3",
    "ClearOATP2B1L",
    "ClearNTCP",
    "ClearBSEP",
    "ClearOAT3",
    "ClearOAT1",
    "ClearOAT4"  
    )) |> 
  select(!c(q25, q75)) |> 
  pivot_wider(names_from = Parameter, values_from = median) |> 
  mutate(
    Papp_SI = Clearapparent / (SA_SI*1e-3),    
    Clear_OATP2B1c = ClearOATP2B1I / (SFOATP2B1I*VI_L),  
    Clear_OATP1B1c = ClearOATP1B1 / (SFOATP1B1*VL_ic*fuLec),
    Clear_OATP1B3c = ClearOATP1B3 / (SFOATP1B3*VL_ic*fuLec),
    Clear_NTCPc = ClearNTCP / (SFNTCP*VL_ic*fuLec),    
    Clear_BSEPc = ClearBSEP / (SFBSEP*VL_ic),
    Clear_OAT3c = ClearOAT3 / (SFOAT3*VPT*fuPTec),
    Clear_OAT1c = ClearOAT1 / (SFOAT1*VPT*fuPTec),
    Clear_OAT4c = ClearOAT4 / (SFOAT4*VPT*fuPTL) 
  ) |> 
  rename(SF_OATP2B1_I = SFOATP2B1I,
         SF_OATP2B1 = SFOATP2B1,
         SF_OATP1B3 = SFOATP1B3,
         SF_OATP1B1 = SFOATP1B1,
         SF_NTCP = SFNTCP,    
         SF_BSEP = SFBSEP,
         SF_OAT1 = SFOAT1,           
         SF_OAT3 = SFOAT3,
         SF_OAT4 = SFOAT4
         ) |> 
  mutate(SF_OAT1 = 809372.8,
         SF_OAT3 = 2782974,
         SF_OAT4 = 10462711) |> 
  pivot_longer(fuPTL:Clear_OAT4c,
               names_to = "Parameter", 
               values_to = "Value")

Selected_initial_P <- Initial_P |> 
  filter_out(Parameter %in% c("logD_AlbuminWater",
                               "logD_GlobulinWater",
                               "logD_MembraneLipidWater",
                               "logD_StorageLipidWater",
                               "logD_StructuralProteinWater")
         )

Selected_final_P <- Selected_final_P |> 
  filter(Parameter %in% unique(Selected_initial_P$Parameter))

Selected_initial_P <- Selected_initial_P |> 
  mutate(Predicted = ifelse(Specific == "calculated", "yes",
                            Predicted),
         Predicted = ifelse(Predicted == "yes", "QSPR-model",
                            "Experimental")) |> 
  select(!c(Specific, Reference, 
            SD_comment, Unit, PBK_mechanism, 
            Distribution, Distribution_comment)) |> 
  rename(Source = Predicted)

Selected_final_P <- Selected_final_P |> 
  mutate(Source = "Calibrated",
         SD = 0) |> 
  select(!c(q25, q75)) |> 
  rename(Value = median)

nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

Selected_initial_P <- full_join(nC_MW_df, Selected_initial_P)
Compare_P <- rbind(Selected_final_P, Selected_initial_P)



labels <- c(
  "Papp_SI" = "Papp log(cm/day)",     
  "Clear_OATP2B1c" = "CLint OATP2B1 (L/day/mg)",
  "Clear_OATP1B1c" = "CLint OATP1B1 (L/day/mg)",               
  "Clear_OATP1B3c" = "CLint OATP1B3 (L/day/mg)",
  "Clear_NTCPc" = "CLint NTCP (L/day/mg)",
  "Clear_BSEPc" = "CLint BSEP (L/day/mg)",
  "Clear_OAT1c" = "CLint OAT1 (L/day/mg)",            
  "Clear_OAT3c" = "CLint OAT3 (L/day/mg)",
  "Clear_OAT4c" = "CLint OAT4 (L/day/mg)",
  "SF_OATP2B1_I" = "SF_OATP2B1_I", 
  "SF_OATP2B1" = "SF_OATP2B1",   
  "SF_OATP1B3" = "SF_OATP1B3",   
  "SF_OATP1B1" = "SF_OATP1B1",   
  "SF_NTCP" = "SF_NTCP",     
  "SF_BSEP" = "SF_BSEP",     
  "SF_OAT1" = "SF_OAT1",
  "SF_OAT3" = "SF_OAT3",
  "SF_OAT4" = "SF_OAT4"
  )


Compare_P_plot <- ggplot() +
  geom_pointrange(data = Compare_P,
                  aes(x = reorder(PFAS, MW),
                      y = Value,
                      ymin = Value-SD,
                      ymax = Value+SD,
                      color =factor(Source,
                                    c("Experimental", "QSPR-model", "Calibrated")),
                      shape =factor(Source,
                                      c("Experimental", "QSPR-model", "Calibrated"))),
                  size = 0.3,
                  alpha = 0.8,
                  linewidth = .2) +  
  geom_pointrange(data = Selected_final_P,
                  aes(x = reorder(PFAS, MW),
                      y = Value,
                      ymin = q25,
                      ymax = q75),
                  colour = "grey50",
                  size = 0.3,
                  alpha = 0.8,
                  linewidth = .2) +  
  scale_color_viridis(option = "D", discrete = TRUE, direction = -1 , end = 0.8) +
  facet_wrap(~factor(Parameter,
                     c("Papp_SI",
                       "Clear_OATP2B1c",
                       "Clear_OATP1B1c",
                       "Clear_OATP1B3c",
                       "Clear_NTCPc",
                       "Clear_BSEPc",
                       "Clear_OAT3c",
                       "Clear_OAT1c",
                       "Clear_OAT4c",
                       "SF_OATP2B1_I",
                       "SF_OATP2B1",
                       "SF_OATP1B3",
                       "SF_OATP1B1",
                       "SF_NTCP",
                       "SF_BSEP",
                       "SF_OAT1",
                       "SF_OAT3",
                       "SF_OAT4"
                     )),
             scales = "free_y",
             labeller = as_labeller(labels)) +
  labs(x = "",
       y = "") +
  theme_quick +
  theme(axis.text.x = element_text(angle = 90)
        )


Compare_P_plot = Compare_P_plot + Patchwork_theme +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "right")

ggsave(plot = Compare_P_plot,
       file = here(OUTPUT,"ComparePplot.png"),
       dpi = 600,
       width = 24,
       height = 18,
       units = "cm")


