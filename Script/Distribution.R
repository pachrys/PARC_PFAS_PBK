# --------------------------------------------------------------------------- #
# Script for running QIVIVE 4 Nynke
# Chrysanthi Pachoulide, September 2026
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

# Set output storage directory
OUTPUT <- here("Output", "0.00019")
dir.create(OUTPUT, recursive = TRUE)

# Load data

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

PFAS_i <- "PFOA"
OD <- 0.00019 #14450*2203.137

nPFAS <- length(PFAS_list)

## Tissue constituents for partition coefficient calculation
Physio.data <- read_csv(here("Input", "TissueComposition.csv")) %>% 
  filter(Species == "Human") %>%
  filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
  select(- c(Comment, Species)) %>% 
  mutate(across(-Tissue, as.numeric))


## Input parameters
Physio_params <- read_xlsx(here("Input", "Physio_Parameters.xlsx")) 
TK_params <- read_csv(here("Input", "Final_TK_Parameters.csv")) 
ncMW <- read_csv(here("Input", "nC_MW_df.csv"))
Tissue_component_distributions <- read_xlsx(here("Input", "TK_Parameters_last.xlsx"))


# load functions
source(here("Script", "functions_NG_PBK_model.R"))


for (i in 1:nPFAS) {
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}-n{i}-----------------------------------"))
  
  # parameters
  
  parm.phys <- Physio_params$Value
  names(parm.phys) <- Physio_params$Abbreviation
  
  names_parm.TK <- TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(Parameter) |> pull()
  
  parm.TK <- TK_params |> 
    filter(PFAS == PFAS_i) |> 
    select(median) |> pull()
   
  names(parm.TK) <- names_parm.TK
  
  pars <- c(parm.phys, parm.TK)
  
  
  Dose_ug_day <- OD*parm.phys[["BW"]]
  
  TIMES <- seq(0, 365*40, by = 1)
  
  ODevents <- data.frame(var = "AIL",
                         time = 0:365*40, # daily exposure
                         value = Dose_ug_day,
                         method = "add") 
  
  A_init <- c(AIL = 0,
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
              InMB = 0)
  
  PBK.out <- lsoda(y = A_init,
                   times = TIMES, 
                   func = PBK.model.f,
                   parms = pars,
                   events = list(data = ODevents),
                   rtol = 1e-4,        
                   atol = 1e-6,
                   hmax = 0.5,         
                   maxsteps = 100000)
  
  PBK_OUTPUT <- as.data.frame(PBK.out) 
  PBK_OUTPUT[,-1] <- pmax(PBK_OUTPUT[,-1], 0)
  
  MW <- ncMW |> 
    filter(PFAS == PFAS_i) |> 
    select(MW) |> pull()
  
  Ds <- Tissue_component_distributions |> 
    filter(PFAS == PFAS_i,
           PBK_mechanism == "Distribution"#,
           # Parameter != "logD_GlobulinWater"
           ) |> 
    mutate(Value = exp(Value),
           D = case_when(
             str_detect(Parameter, "logD_AlbuminWater") ~ "D_ALB",
             str_detect(Parameter, "logD_MembraneLipidWater") ~ "D_ML",
             str_detect(Parameter, "logD_StructuralProteinWater") ~ "D_SP",
             str_detect(Parameter, "logD_StorageLipidWater") ~ "D_SL",
             str_detect(Parameter, "logD_GlobulinWater") ~ "D_GLOB"
             )
           ) |> 
    select(D, Value) |> 
    pivot_wider(names_from = D, values_from = Value)
  
  LiverComp <- Physio.data |> 
    filter(Tissue == "Liver")
  KidneyComp <- Physio.data |> 
    filter(Tissue == "Kidney")
  PlasmaComp <- Physio.data |> 
    filter(Tissue == "Plasma")
  
  P_ALB <- Ds$D_ALB*PlasmaComp$f_ALB
  f_Glob <- 0.031
  f_W <- 1 - PlasmaComp$f_ALB - f_Glob
  P_Glob <- Ds$D_GLOB*f_Glob
  P_Total <- P_ALB + P_Glob + f_W
  P_fra_ALB <- P_ALB/P_Total
  P_fra_Glob <- P_Glob/P_Total
  P_fra_W <- f_W/P_Total
  
  L_W <- LiverComp$f_W
  L_ALB <- Ds$D_ALB*LiverComp$f_ALB
  L_ML <- Ds$D_ML*LiverComp$f_ML
  L_SL <- Ds$D_SL*LiverComp$f_SL
  L_SP <- Ds$D_SP*LiverComp$f_SP
  L_Total <- L_ALB + L_ML + L_SL + L_SP + L_W
  L_fra_ML <- L_ML/L_Total
  L_fra_SL <- L_SL/L_Total
  L_fra_ALB <- L_ALB/L_Total
  L_fra_SP <- L_SP/L_Total
  L_fra_W <- L_W/L_Total
  
  K_W <- KidneyComp$f_W
  K_ALB <- Ds$D_ALB*KidneyComp$f_ALB
  K_ML <- Ds$D_ML*KidneyComp$f_ML
  K_SL <- Ds$D_SL*KidneyComp$f_SL
  K_SP <- Ds$D_SP*KidneyComp$f_SP
  K_Total <- K_ALB + K_ML + K_SL + K_SP + K_W
  K_fra_ML <- K_ML/K_Total * 10
  K_fra_SL <- K_SL/K_Total * 10
  K_fra_ALB <- K_ALB/K_Total * 10
  K_fra_SP <- K_SP/K_Total * 10
  K_fra_W <- K_W/K_Total * 10
  
  VL <- parm.phys[["VLc"]]*parm.phys[["BW"]]
  VL_ic <- parm.phys[["VL_icc"]]*VL
  VL_ec <- parm.phys[["VL_ecc"]]*VL
  
  VK <- parm.phys[["VKc"]]*parm.phys[["BW"]]
  VPT <- parm.phys[["VPTc"]]*VK
  VPTT <- VPT*parm.phys[["VPTTc"]]
  VPTL <- VPT*parm.phys[["VPTLc"]]
  VRKL <- (VK - VPT)*parm.phys[["VRKLc"]]
  VRKT <- VK - VPTT - VPTL - VRKL
  
  Working <- PBK_OUTPUT |> 
    mutate(CL = (CL_ec*VL_ec + CL_ic*VL_ic)/VL,
           CK = (CPTL*VPTL + CPTT*VPTT + CRKL*VRKL + CRKT*VRKT)/VK
           ) |> 
    select(time, CP, CL, CK) |> 
    mutate(Plasma_ALB = CP*P_fra_ALB,
           Plasma_GLOB = CP*P_fra_Glob,
           Plasma_free = CP*P_fra_W,
           Liver_ALB = CL*L_fra_ALB,
           Liver_SP = CL*L_fra_SP,
           Liver_ML = CL*L_fra_ML,
           Liver_SL = CL*L_fra_SL,
           Liver_free = CL*L_fra_W,
           Kidney_ALB = CK*K_fra_ALB,
           Kidney_SP = CK*K_fra_SP,
           Kidney_ML = CK*K_fra_ML,
           Kidney_SL = CK*K_fra_SL,
           Kidney_free = CK*K_fra_W
           ) |> 
    select(!c(CP, CL, CK)) |>
    pivot_longer(Plasma_ALB:Kidney_free, 
                 names_to = "Tissue", 
                 values_to = "Concentration") |> 
    mutate(Concentration = Concentration /MW, # to umol/L
           time = time/365)
  
  Last <- Working |> 
    filter(time == max(Working$time),
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
  
  Distribution_plot <- Last |> 
    ggplot(aes(Tissue, Concentration, fill = as.factor(Component))) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_viridis(option = "D", discrete = TRUE,
                       end = 0.9) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "",
         y = "Mass balance", 
         fill = "Tissue component")+
    theme_quick
  Distribution_plot 
  
  Distribution_plot <- Distribution_plot + Patchwork_theme +
    theme(legend.position = "right")
  
  ggsave(plot = Distribution_plot,
         file = here(OUTPUT, paste0(PFAS_i,"Distribution_plot.png")),
         dpi = 600,
         width = 16,
         height = 12,
         units = "cm")
  
  Working <- Working |> 
    mutate(PFAS = PFAS_i)
  write.csv(Working, here(OUTPUT, paste0(PFAS_i,"Tissue_Distribution.csv")), row.names = FALSE)

  Cs <- PBK_OUTPUT |> 
    select(time, CI:CP) |> 
    mutate(CL = (CL_ec*VL_ec + CL_ic*VL_ic)/VL,
           CK = (CPTL*VPTL + CPTT*VPTT + CRKL*VRKL + CRKT*VRKT)/VK) |> 
    select(!c(CL_ec, CL_ic, CPTL, CPTT, CRKL, CRKT, CVI, CVL_ec, 
              CVRKT, CVA, CVR)) |> 
    pivot_longer(CI:CK, names_to = "Tissue", values_to = "Concentration") |> 
    mutate(PFAS = PFAS_i)
  write.csv(Cs, here(OUTPUT, paste0(PFAS_i,"Tissue_Cs.csv")), row.names = FALSE)
  
}

# summarise ----

Cdistr_df <-
  list.files(path = here(OUTPUT),
             pattern = "\\Tissue_Distribution.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .))) 

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
# 
# test <- Cdistr_df |> 
#   select(!time) |> 
#   pivot_wider(names_from = Tissue, values_from = Concentration)

Cdistr_df <- full_join(Cdistr_df, ncMW)


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








Liver_distri <- Liver_distri + Patchwork_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 30,
                                   face = "bold", 
                                   margin = margin(t = 2, 0, 0, 0)
        ))

ggsave(plot = Liver_distri,
       file = here(OUTPUT, "Liver_distri.png"),
       dpi = 600,
       width = 16,
       height = 12,
       units = "cm")



Kidney_distri <- Kidney_distri + Patchwork_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 30,
                                   face = "bold", 
                                   margin = margin(t = 2, 0, 0, 0)
        ))

ggsave(plot = Kidney_distri,
       file = here(OUTPUT, "Kidney_distri.png"),
       dpi = 600,
       width = 16,
       height = 12,
       units = "cm")


Plasma_distri <- Plasma_distri + Patchwork_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 30,
                                   face = "bold", 
                                   margin = margin(t = 2, 0, 0, 0)
        ))

ggsave(plot = Plasma_distri,
       file = here(OUTPUT, "Plasma_distri.png"),
       dpi = 600,
       width = 16,
       height = 12,
       units = "cm")


# All tissues ----

Tdistr_df <-
  list.files(path = here(OUTPUT),
             pattern = "\\Tissue_Cs.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .))) 

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
                                       "adipose",
                                       "plasma")))
Tdistr_df <- full_join(Tdistr_df, ncMW)

## Check steady state plasma ----

### wide plot ----

install.packages("geomtextpath")
library(geomtextpath)
library(ggrepel)
data4ss_plot = filter(Tdistr_df,
                      Tissue == "plasma") 



data4ss_plot_label <- data4ss_plot |> 
  mutate(label = ifelse(time == max(data4ss_plot$time), PFAS, NA),
         # nudge_y = case_when(PFAS == "HFPO-DA" ~ 0.001,
         #                     PFAS == "ADONA"~ 0.001,
         #                     PFAS == "PFBA" ~0.001, 
         #                     PFAS == "PFPeA" ~ 0.001,
         #                     PFAS == "PFBS" ~ 0.001,
         #                     PFAS == "PFHxA" ~ 0.001),
         nudge_x = 1)|> 
  mutate(PFAS = factor(PFAS, levels = rev(PFAS_list)),
         time = time/365)

install.packages("ggforce")
library(ggforce)

Plasma_perCat_Cs <- ggplot(data = data4ss_plot_label, 
                          aes(time, Concentration, 
                              colour = PFAS)) + #factor(PFAS, levels = rev(PFAS_list)))) +
  geom_path(linewidth = 0.1, alpha = 0.8) +
  geom_text_repel(aes(label = label),
                  min.segment.length = 0.1,
                  show.legend = FALSE,
                  size = unit(8, "mm"),
                  # label.size = unit(0.25, "mm"),
                  # label.padding = unit(0.8, "mm"),
                  # label.r = unit(0.5, "mm"),
                  segment.size = 0.1,
                  nudge_x = 2.5,
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

Plasma_perCat_Cs <- Plasma_perCat_Cs + Patchwork_theme +
  theme(legend.position = "none")


ggsave(plot = Plasma_perCat_Cs,
       file = here(OUTPUT, "Plasma_Cons.png"),
       dpi = 600,
       width = 11,
       height = 9,
       units = "cm")  

### long plot ----
Plasma_Cs <- ggplot(data = filter(Tdistr_df,
                                  Tissue == "plasma"), 
                    aes(time/365, Concentration, 
                        colour = Category_OECD)) + #factor(PFAS, levels = rev(PFAS_list)))) +
  geom_path(linewidth = 0.25, alpha = 0.9) +
  scale_colour_viridis(option = "D", discrete = TRUE,
                       end = 0.7) +
  facet_wrap(~factor(PFAS, levels = rev(PFAS_list)), #factor(PFAS, levels = rev(PFAS_list))
             # scales = "free_y",
             nrow = 3) +
  labs(x = "Time (years)",
       y = "Concentration (μmol/L)") +
  theme_quick
Plasma_Cs

Plasma_Cs <- Plasma_Cs + Patchwork_theme +
  theme(legend.position = c(0.955, 0.05),
        strip.text = element_text(size = 30),
        axis.text.y = element_text(size = 25),
        axis.text.x = element_text(size = 25))

ggsave(plot = Plasma_Cs,
       file = here(OUTPUT, "Plasma_Cwide.png"),
       dpi = 600,
       width = 14,
       height = 6,
       units = "cm")  


## tissue mass balance ----
Tdistr_df <- Tdistr_df |> 
    filter(time == max(Tdistr_df$time))

Tdistr_plot <- Tdistr_df |> 
  ggplot(aes(reorder(PFAS, MW), Concentration, fill = Tissue)) + #factor(PFAS, levels = rev(order))
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_viridis(option = "D", discrete = TRUE
                     ) +
  # facet_wrap(~factor(PFAS, levels = rev(order))) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Distribution in tissues",
       x = "",
       y = "Mass balance", 
       fill = "Tissue")+
  theme_quick
Tdistr_plot 

Tdistr_plot <- Tdistr_plot + Patchwork_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 30,
                                   face = "bold", 
                                   margin = margin(t = 2, 0, 0, 0)
        ))

ggsave(plot = Tdistr_plot,
       file = here(OUTPUT, paste0(PFAS_i,"Tdistr_plot.png")),
       dpi = 600,
       width = 16,
       height = 12,
       units = "cm")


Tissue <- Plasma_distri + Liver_distri + Kidney_distri +
  plot_layout(nrow = 1, 
              guides = "collect") &
  Patchwork_theme &
  theme(plot.margin = margin(t = -0.5, 0.1, b = -0.5, 0.1, "cm"))

Tissues <- Tdistr_plot + Patchwork_theme +
  theme(plot.margin = margin(0.1, 0.1, b = -0.5, 0.1, "cm"))

P <- Tissues / Tissue

ggsave(plot = P,
       file = here(OUTPUT, "FinalTissueN.png"),
       dpi = 600,
       width = 16,
       height = 18,
       units = "cm")

# 


# Distribution_plot <- Cdistr_df |> 
#   ggplot(aes(Tissue, Concentration, fill = Component)) +
#   geom_bar(stat = "identity", position = "fill") +
#   scale_fill_viridis(option = "D", discrete = TRUE,
#                      end = 0.9) +
#   facet_wrap(~factor(PFAS, levels = rev(order))) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "",
#        y = "Mass balance", 
#        fill = "Tissue component")+
#   theme_quick
# Distribution_plot 
# 
# Distribution_plot <- Distribution_plot + Patchwork_theme +
#   theme(legend.position = "right",
#         axis.text.x = element_text(angle = 30,
#                                    face = "bold", 
#                                    margin = margin(t = 2, 0, 0, 0)
#         ))
# 
# ggsave(plot = Distribution_plot,
#        file = here(OUTPUT, "Distribution_plot.png"),
#        dpi = 600,
#        width = 16,
#        height = 12,
#        units = "cm")