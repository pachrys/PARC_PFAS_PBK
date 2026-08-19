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
library(patchwork)

# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

theme_quick <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'),
    plot.title = element_text(hjust = 0, face = "bold"),

    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_line(color = "grey90"),
    
    axis.line = element_blank(),
    axis.ticks=element_blank(),
    
    legend.position = "none"
  )


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
    axis.text = element_text(size = 30),
    # axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30), #, face = "bold"), #, angle = 30),
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
    legend.margin = margin(.1, .1, .1, .1),
    
    # strip.text.x = element_blank()
    strip.text.x = element_text(size = 35),
    strip.text.y = element_text(size = 35)
    
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
    

# Set input and output directory
INPUT <- here("Output", "Second_try_UA")
OUTPUT <- here("Output","Second_try_UA", "Processed_results") 
dir.create(OUTPUT, recursive = TRUE)

# Dose_df <- read_excel(here("Input","Evaluation_data.xlsx")) %>%
#   filter(Reference == "Abraham_2024") |>
#   select(!c(PFAS_name, CAS_n)) %>%
#   rename(PFAS = Acronym) |>
#   mutate(Dose = round(Dose_ugPerkg_bw*82, 3)) |>  # 82 is the body weight of Abraham
#   select(PFAS, Dose) |>
#   mutate(PFAS = str_replace_all(PFAS, "_", "-")) |>
#   distinct()


# functions ----

source(here("Script", "functions_calculate_results.R"))

# PFAS ----


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
  Dose = 0.1) 

# Dose_df <- tibble(
#   PFAS = PFAS_list,
#   Dose = 4) |>
#   mutate(Dose = ifelse(PFAS == "PFBS", 18.4, 
#                        ifelse(PFAS == "HFPO-DA" , 19.6,
#                               Dose)))

# Dose_df <- Abraham_df
# PFAS_list <- Dose_df$PFAS

nPFAS <- length(PFAS_list)


# produce result files and plots for all PFAS ----
for (i in 1:nPFAS) {


  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}-n{i}-----------------------------------"))

  # Inputs ----

  # load results (called savey)
  load(here(INPUT, paste0(PFAS_i, "_Monte_Carlo_y.RData" )
            ))
  # load oral dose (different per PFAS if comparing results to Abraham)
  OD <- filter(Dose_df, PFAS == PFAS_i) |>
    select(Dose) |> pull()

  # TK results ----
  # calculate AUC, half-life, elimination rates, fraction absorbed etc
  # 
  # To see if there are any errors
  # clean_savey <- savey[sapply(savey, function(x) nrow(x) >= 3)]
  # # paste0("Iteration ", which(sapply(savey, function(x) nrow(x) < 3)), " gave error")
  #  TK_results <- sapply(savey, check_error, OD)
  # TK_results <- sapply(savey, check_error, OD)
  # proper_savey <- savey[-c(3)] #,147,466,549,660,969)]
  
  TK_results <- sapply(savey, calc_TK_results, OD) #proper_savey

  # to clean TK_results
  clean_TK_results <- TK_results |> 
    t() |> 
    as.data.frame() |>
    filter(CLfe < 1e5, 
           CLur < 1e5) |> 
    t()
  
  # find which iterations went wrong
  errored_TK_results <- TK_results |> 
    t() |> 
    as.data.frame() |> 
    rownames_to_column() |> 
    rename(Iteration = rowname) |>
    filter(CLfe > 1e5, 
           CLur > 1e5) 
  
  print(errored_TK_results$Iteration)
  if(is.numeric(errored_TK_results$Iteration)){
    write.csv(errored_TK_results, file = here(OUTPUT, paste0(PFAS_i, "errored_TK_results.csv")))
  }
  
  
  TK_summary_results <- calc_TK_summary_results(clean_TK_results) |>
    as.data.frame() |> mutate(PFAS = PFAS_i)
  write.csv(TK_summary_results, file = here(OUTPUT, paste0(PFAS_i, "_TK_sumstat.csv")))

  # Predicted tissue concentrations ----

  parm.phys <- c(
    BW = 64.41978,
    QC = 5486.301,
    Hct = 0.3645031,
    VIc = 0.016,
    VILc = 0.03848451,
    SA_SIc = 1099.557,
    VLc = 0.02458126,
    VL_icc = 0.835,
    VL_ecc = 0.161,
    VKc = 0.004630086,
    VPTc = 0.398,
    VPTTc = 0.5470631,
    VPTLc = 0.403,
    VRKLc = 0.283,
    VAc = 0.3426787,
    MAc = 1.487552,
    VPc = 0.065,
    QIc = 0.165,
    QLc = 0.06857434,
    QKc = 0.1761446,
    QAc = 0.09413654,
    QUrc = 0.022,
    GFR = 193.2765, # = 0.2*0.1761446*5486.301, #156.6414,
    QTc = 0.36,
    tco = 3.36,
    Qmenses = 0.0008561644)

  VL <- parm.phys[["VLc"]]*parm.phys[["BW"]]
  VL_ic <- parm.phys[["VL_icc"]]*VL
  VL_ec <- parm.phys[["VL_ecc"]]*VL

  VK <- parm.phys[["VKc"]]*parm.phys[["BW"]]
  VPT <- parm.phys[["VPTc"]]*VK
  VPTT <- VPT*parm.phys[["VPTTc"]]
  VPTL <- VPT*parm.phys[["VPTLc"]]
  VRKL <- (VK - VPT)*parm.phys[["VRKLc"]]
  VRKT <- VK - VPTT - VPTL - VRKL

  Initial_Tissues_list <- c("CI", "CL_ec", "CL_ic", "CPTL", "CPTT", "CRKL", "CRKT", "CA", "CR", "CP")

  Results_df <- do.call(rbind, savey)

  C_df <- Results_df |>
    select(time, CI, CL_ec, CL_ic, CPTL, CPTT, CRKL, CRKT, CA, CR, CP, Iteration) |>
    mutate(PFAS = PFAS_i) |>
    mutate(CL = (CL_ec*VL_ec + CL_ic*VL_ic)/VL,
           CK = (CPTL*VPTL + CPTT*VPTT + CRKL*VRKL + CRKT*VRKT)/VK,
           time = time/365) |> # time to years
    select(PFAS, Iteration, time, CA, CR, CI, CL, CK, CP) |>
    pivot_longer(CA:CP,
                 names_to = "Tissue",
                 values_to = "Concentration")

  Tissues_list <- c("CA", "CR", "CI", "CL", "CK", "CP")

  Tissues_names <- c("Adipose", "Rest","Intestine", "Liver", "Kidney", "Plasma")

  # plasma concentration
  nTissues <- length(Tissues_list)

  C_Tissues_list <- NULL

  for (j in 1:nTissues) {

    Tissue_j <- Tissues_list[j]
    Tissue_j_name <- Tissues_names[j]

    C_summary <- calc_C_summary(C_df, Tissue_j) |>
      mutate(Tissue = Tissue_j_name,
             PFAS = PFAS_i)

    C_Tissues_list[[j]] <- C_summary

  }

  C_Tissues <- do.call(rbind, C_Tissues_list)

  C_Tissues_plot <- plot_all_tissues(C_Tissues)

  C_Tissues_plot <- C_Tissues_plot + Patchwork_theme +
    theme(axis.text.x = element_text(angle = 0))

  ggsave(plot = C_Tissues_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_C_Tissues_plot.png")),
         dpi = 600,
         width = 16,
         height = 11,
         units = "cm")

  write.csv(C_Tissues, file = here(OUTPUT, paste0(PFAS_i, "_C_Tissues_sumstat.csv")))

  # # Actual parameters ----

  Parameters_df <- Results_df |>
    select(fup:CF_K_reabsorption) |>
    select(!starts_with("CF_")) |> 
    distinct()

  write.csv(Parameters_df, file = here(OUTPUT, paste0(PFAS_i, "_Parameters.csv")))

  P_long_df <- Parameters_df |>
    rename(Clear_reabsorption = Clear_net_reabsorption_unscaled) |> 
    pivot_longer(everything(),
                 names_to = "Parameter",
                 values_to = "Value")

  Parameters_plot <- ggplot(P_long_df, aes(x = Value)) +
    geom_histogram() +
    facet_wrap(~ factor(Parameter, 
                        c("PI",   
                          "PL",        
                          "PK",          
                          "PA",         
                          "PR",  
                          "fup", 
                          "fu_PTL",     
                          "fu_L_ec",     
                          "fu_PT_ec",
                          "SF_OATP2B1_I",          
                          "SF_OATP2B1",              
                          "SF_OATP1B3",                     
                          "SF_OATP1B1",             
                          "SF_NTCP",               
                          "SF_BSEP",                        
                          "SF_OAT1",             
                          "SF_OAT3",                
                          "SF_OAT4",                        
                          "Clear_apparent",       
                          "Clear_OATP2B1_I",          
                          "Clear_uptake_I",                 
                          "Clear_OATP1B1",              
                          "Clear_OATP1B3",                
                          "Clear_OATP2B1_L",                
                          "Clear_NTCP",                  
                          "Clear_BSEP",                 
                          "Clear_uptake_L",                 
                          "Clear_excretion_L",         
                          "Clear_OAT3",              
                          "Clear_OAT1",                      
                          "Clear_OAT4",              
                          "Clear_reabsorption"
                         )), scales = "free") +
    labs(x = "Value",
         y = "Frequency") +
    coord_flip()+
    theme_quick

  Parameters_plot = Parameters_plot + Patchwork_theme +
    theme(axis.text.x = element_text(angle = 0))

  ggsave(plot = Parameters_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_C_Parameters_plot.png")),
         dpi = 600,
         width = 24,
         height = 18,
         units = "cm")

}

# run summary results ----
## load input files ----
nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

TK_df <-
  list.files(path = here(OUTPUT),
             pattern = "\\_TK_sumstat.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .)))|> select(!"...1")
write.csv(TK_df, here(OUTPUT, "TK_results_all_PFAS.csv"), row.names = FALSE)

C_df <-
  list.files(path = here(OUTPUT),
             pattern = "\\_C_Tissues_sumstat.csv$") |>
  map_df(~ read_csv(file = here(OUTPUT, .)))

Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |>
  rename(PFAS = Acronym) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))


# plot results ----
## Plasma ----
C_plasma_predicted_df <- C_df |>
  filter(Tissue == "Plasma") 

C_plasma_predicted_df <- full_join(C_plasma_predicted_df, nC_MW_df) 



CP_plot <- ggplot() +
  geom_ribbon(data = C_plasma_predicted_df,
              aes(x = Time/365,
                  ymin = q25,
                  ymax = q75,
                  fill = Category_OECD),
              alpha = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
  geom_path(data = C_plasma_predicted_df,
            aes(Time/365, median,
                colour = Category_OECD)) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  facet_wrap(~ factor(PFAS, c( #factor(n_C),
    "PFBA",
    "PFPeA",
    "PFBS",
    "PFHxA",
    "HFPO-DA",
    "PFPeS",
    "PFHpA",
    "ADONA",
    "PFHxS",
    "PFOA",
    "PFHpS",
    "PFNA",
    "PFOS",
    "PFDA",
    "PFUnDA",
    "PFDoDA",
    "PFTrDA")
  ),
  nrow = 2) +
  theme_quick

CP_plot <- CP_plot + Patchwork_theme +
  theme(legend.position = "top")

ggsave(plot = CP_plot,
       file = here(OUTPUT, "CP_plot_UA.png"),
       dpi = 600,
       width = 24,
       height = 11,
       units = "cm")

## Half life ----

HL_predicted_df <- TK_df |>
  filter(Parameter == "HL")

HL_predicted_df <- full_join(HL_predicted_df, nC_MW_df)

HL_predicted_plot <- ggplot() +
  geom_errorbar(data = HL_predicted_df,
                aes(x = reorder(PFAS, MW),
                    y = log10(median),
                    colour = Category_OECD,
                    ymin = log10(q25),
                    ymax = log10(q75)),
                width = .01,
                linewidth = .4) +
  geom_point(data = HL_predicted_df,
             aes(x = reorder(PFAS, MW),
                 y = log10(median),
                 shape = Category_OECD,
                 colour = Category_OECD),
             size = 1.5) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_quick
# HL_predicted_plot
# 
# HL_predicted_plot <- HL_predicted_plot + Patchwork_theme +
#   theme(legend.position = "top",
#         axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0))
#         )
# 
# ggsave(plot = HL_predicted_plot,
#        file = here(OUTPUT,"HL_predicted_plot.png"),
#        dpi = 600,
#        width = 12,
#        height = 11,
#        units = "cm")



## Fraction eliminated in urine and feces ----

FelFa_predicted_df <- TK_df |> 
  filter(Parameter %in% c("Fa", 
                          # "FelTotal",
                          "Ffeces",
                          "FbilEl",
                          "Furine")) |> 
  pivot_wider(names_from = Parameter, values_from = c(median, mean, q25, q75)) |> 
  mutate(median_Fpreintfeces = median_Ffeces-median_FbilEl,
         mean_Fpreintfeces = mean_Ffeces-mean_FbilEl,
         q25_Fpreintfeces = q25_Ffeces-q25_FbilEl,
         q75_Fpreintfeces = q75_Ffeces-q75_FbilEl,
         median_cumFpreintfeces = median_Fpreintfeces + median_Furine,
         mean_cumFpreintfeces = mean_Fpreintfeces + mean_Furine,
         q25_cumFpreintfeces = q25_Fpreintfeces + q25_Furine,
         q75_cumFpreintfeces = q75_Fpreintfeces + q75_Furine,
         median_cumFbile = median_Fpreintfeces + median_Furine + median_FbilEl,
         mean_cumFbile = mean_Fpreintfeces + mean_Furine + mean_FbilEl,
         q25_cumFbile = q25_Fpreintfeces + q25_Furine + q25_FbilEl,
         q75_cumFbile = q75_Fpreintfeces + q75_Furine + q75_FbilEl,
  ) |>  
  pivot_longer(median_Fa:q75_cumFbile, names_to = c("What", "Parameter"), values_to = "Value", names_sep = "_") |> 
  pivot_wider(names_from = "What", values_from = "Value") |> 
  filter(Parameter != "Ffeces") |> 
  mutate(Direction = ifelse(Parameter == "Fa", "In",
                            ifelse(Parameter == "Furine", "Urine",
                                   "Feces"))) 

FelFa_predicted_df <- full_join(FelFa_predicted_df, nC_MW_df)

Clean_FelFa <- FelFa_predicted_df |> 
  filter(Parameter %in% c("FbilEl", "Fpreintfeces", "Furine"))
Cum_FelFa <- FelFa_predicted_df |> 
  filter(Parameter %in% c("cumFbile", "cumFpreintfeces", "Furine")) |> 
  mutate(Parameter = str_replace_all(Parameter, "cum", ""),
         Parameter = str_replace_all(Parameter, "Fbile", "FbilEl"))

Fa_predicted_df <- TK_df |> 
  filter(Parameter %in% c("Fa")) 

Fa_predicted_df <- full_join(Fa_predicted_df, nC_MW_df)

Fa_plot <- ggplot()+
  geom_boxplot(data = Fa_predicted_df,
               aes(x = reorder(PFAS, MW), y = median,
                   ymin = q25, lower = q25, middle = median, upper = q75, ymax = q75,
                   fill = Category_OECD, colour = Category_OECD),
               stat = "identity",
               linewidth = 0.3,
               width = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.7) +
  scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
  theme_quick



Fe_plot <- ggplot() +
  geom_bar(data = Clean_FelFa,
           aes(x = reorder(PFAS, MW), y = median, 
               fill = Parameter),
           stat = "identity",
           position = "stack",
           width = 0.5) +
  scale_fill_manual(name = "Parameter",
                    values = c("FbilEl" = "#B198B9",
                               "Fpreintfeces" = "#8C6897",
                               "Furine" = "#89D499"),
                    labels = c("FbilEl" = "Feces: bile",
                               "Fpreintfeces" = "Feces: pre-systemic",
                               "Furine" = "Urine")) +
  geom_errorbar(data = Cum_FelFa,
                aes(x = reorder(PFAS, MW),
                    ymin = q25, ymax = q75,
                    colour = Parameter),
                position = "dodge",
                width = .3,
                linewidth = .4
  ) +
  scale_color_manual(name = "Parameter",
                     values = c("FbilEl" = "#683975",
                                "Fpreintfeces" = "#440154",
                                "Furine" = "#2A7741"),
                     labels = c("FbilEl" = "Feces: bile",
                                "Fpreintfeces" = "Feces: pre-systemic",
                                "Furine" = "Urine")) +
  labs(x = "", 
       y = "Cumulative fraction excreted") +
  theme_quick +
  theme(legend.position = "right")


Fe_plot <- Fe_plot + Patchwork_theme +
  theme(legend.position = "top",
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.8, "lines"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        axis.text.x = element_text(size = 35, face = "bold", margin = margin(0, 0, 0, 0)),
        axis.text.y = element_text(size = 35, margin = margin(0, 0, 0, 0)),   
        axis.title = element_text(size = 45) #, face = "bold"
  )


ggsave(plot = Fe_plot,
       file = here(OUTPUT,"Fe_plot.png"),
       dpi = 600,
       width = 16,
       height = 12,
       units = "cm")

## Urinary and fecal clearances ----

CLurCLfe_predicted_df <- TK_df |> 
  filter(Parameter %in% c("CLfe", 
                          "CLur")) 

CLurCLfe_predicted_df <- full_join(CLurCLfe_predicted_df, nC_MW_df)

RurRfe_predicted_df <- TK_df |> 
  filter(Parameter %in% c("RfeEl", 
                          "RurEl"))
RurRfe_predicted_df <- full_join(RurRfe_predicted_df, nC_MW_df)

CL_plot <- ggplot() +
  geom_boxplot(data = CLurCLfe_predicted_df,
               aes(x = reorder(PFAS, MW), y = median,
                   ymin = q25, lower = q25, middle = median, upper = q75, ymax = q75,
                   fill = Category_OECD, colour = Category_OECD),
               stat = "identity",
               linewidth = 0.3,
               width = 0.5) +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, alpha = 0.7) +
  scale_colour_viridis(option = "D", discrete = TRUE, end = 0.7) +
  facet_wrap(~ Parameter, nrow = 2, 
             labeller = as_labeller(c("CLfe" = "Biliary/ fecal",
                                      "CLur" = "Renal")),
             scales = "free_y") +
  theme_quick

# CL_plot <- CL_plot + Patchwork_theme +
#    theme(axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0)),
#          legend.position = "top"
#    )
#  
#  
#  ggsave(plot = CL_plot,
#         file = here(OUTPUT,"CL_plot.png"),
#         dpi = 600,
#         width = 12,
#         height = 11,
#         units = "cm")

# Final plot Uncertainty Analysis ----
UA_one <- CP_plot + 
  labs(title = "A.",
       x = "Time (years)",
       y = "Plasma concentration (μg/L)") +
  guides(fill = FALSE) + guides(colour = FALSE) + 
  Patchwork_theme + 
  theme(legend.position = "none")

UA_two <- Fa_plot + 
  labs(title = "B." ,
       x = "",
       y = "Fraction absorbed") +
  guides(colour = FALSE) + guides(linetype = FALSE) + 
  Patchwork_theme + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90))

UA_three <- CL_plot + 
  labs(title = "C." ,
       x = "",
       y = "Clearance (L/day)") +
  guides(fill = FALSE) + guides(colour = FALSE) + 
  Patchwork_theme + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90))

UA_four <- HL_predicted_plot + 
  labs(title = "D.",
       x = "",
       y = "Half-life log(years)") +
  Patchwork_theme + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))

UA_plot <- UA_one / (UA_two + UA_three + UA_four) +
  plot_layout(heights = c(0.7, 1),
              guides = "collect")  &
  theme(plot.margin = margin(0.1, 0.1, b = 0, 0.1, "cm"),
        legend.box.margin = margin(t = -0.5, 0, b = 0.1, 0, "cm"),
        legend.position='bottom')

ggsave(plot = UA_plot,
       file = here(OUTPUT,"UA_plot_n.png"),
       dpi = 600,
       width = 24,
       height = 16,
       units = "cm")

## Against HL measured ----


Measured_df <- read_excel(here("Input", "Evaluation_data.xlsx")) |>
  rename(PFAS = Acronym) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))

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


HL_measured_df <- full_join(HL_measured, nC_MW_df)

HL_MesPred_plot <- ggplot() +
  geom_jitter(data = HL_measured_df,
              aes(reorder(PFAS, MW),
                  log10(Value),
                  shape = Category_OECD),
              color = "grey60",
              width = 0.2,
              alpha = 0.9,
              size = .3) +
  geom_errorbar(data = HL_predicted_df,
                aes(x = reorder(PFAS, MW),
                    y = log10(median),
                    colour = Category_OECD,
                    ymin = log10(q25),
                    ymax = log10(q75)),
                width = .01,
                linewidth = .4) +
  geom_point(data = HL_predicted_df,
             aes(x = reorder(PFAS, MW),
                 y = log10(median),
                 shape = Category_OECD,
                 colour = Category_OECD),
             size = 1.5) +
  scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme_quick

HL_MesPred_plot

HL_MesPred_plot <- HL_MesPred_plot + Patchwork_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 30, face = "bold", margin = margin(0, 0, 0, 0))
        )

ggsave(plot = HL_MesPred_plot,
       file = here(OUTPUT,"HL_MesPred_plot.png"),
       dpi = 600,
       width = 16,
       height = 11,
       units = "cm")

### Plasma against predicted ----
# 
# Measured_Abraham_df <- Measured_df |>
#   filter(Reference == "Abraham_2024")
# 
# C_plasma_measured_df <- Measured_df |>
#   filter(Parameter == "C_plasma") |>
#   mutate(Time = as.numeric(Time_days)) |>
#   select(PFAS, Time, Value)
# C_plasma_measured_df <- full_join(C_plasma_measured_df, nC_MW_df) 
# 
# C_plasma_predicted_df <- C_df |>
#   filter(Tissue == "Plasma",
#          Time < 450)
# 
# C_plasma_predicted_df <- full_join(C_plasma_predicted_df, nC_MW_df) 
# 
# 
# ObsVsPred_PFCAs_plasma <- ggplot() +
#   geom_ribbon(data = filter(C_plasma_predicted_df, Category_OECD == "PFCAs"),
#               aes(x = Time/365,
#                   ymin = q25,
#                   ymax = q75,
#                   fill = as.factor(MW)),
#               alpha = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   geom_path(data = filter(C_plasma_predicted_df, Category_OECD == "PFCAs"),
#             aes(Time/365, median,
#                 colour = as.factor(MW))) +
#   geom_point(data = filter(C_plasma_measured_df, Category_OECD == "PFCAs"),
#              aes(Time/365, Value),
#              colour = "black",
#              size = 0.1) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.7) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA")
#   ),
#   nrow = 1) +
#   labs(x = "Time (years)",
#        y = "Plasma concentration (μg/L)") +
#   theme_quick
# 
# ObsVsPred_PFSAs_plasma <- ggplot() +
#   geom_ribbon(data = filter(C_plasma_predicted_df, Category_OECD == "PFSAs"),
#               aes(x = Time/365,
#                   ymin = q25,
#                   ymax = q75,
#                   fill = as.factor(MW)),
#               alpha = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.3) +
#   geom_path(data = filter(C_plasma_predicted_df, Category_OECD == "PFSAs"),
#             aes(Time/365, median,
#                 colour = as.factor(MW))) +
#   scale_color_viridis(option = "D", discrete = TRUE, end = 0.3) +
#   geom_point(data = filter(C_plasma_measured_df, Category_OECD == "PFSAs"),
#              aes(Time/365, Value),
#              colour = "black",
#              size = 0.1) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA"
#   )),
#   nrow = 1) +
#   labs(x = "Time (years)",
#        y = "Plasma concentration (μg/L)") +
#   theme_quick
# 
# ObsVsPred_PFECAs_plasma <- ggplot() +
#   geom_ribbon(data = filter(C_plasma_predicted_df, Category_OECD == "PFECAs"),
#               aes(x = Time/365,
#                   ymin = q25,
#                   ymax = q75,
#                   fill = as.factor(MW)),
#               alpha = 0.5) +
#   scale_fill_viridis(option = "D", discrete = TRUE, begin = 0.2, end = 0.3) +
#   geom_path(data = filter(C_plasma_predicted_df, Category_OECD == "PFECAs"),
#             aes(Time/365, median,
#                 colour = as.factor(MW))) +
#   scale_color_viridis(option = "D", discrete = TRUE, begin = 0.2, end = 0.3) +
#   geom_point(data = filter(C_plasma_measured_df, Category_OECD == "PFECAs"),
#              aes(Time/365, Value),
#              colour = "black",
#              size = 0.1) +
#   facet_wrap(~ factor(PFAS, c( #factor(n_C),
#     "PFBA",
#     "PFPeA",
#     "PFBS",
#     "PFHxA",
#     "HFPO-DA",
#     "PFPeS",
#     "PFHpA",
#     "ADONA",
#     "PFHxS",
#     "PFOA",
#     "PFHpS",
#     "PFNA",
#     "PFOS",
#     "PFDA",
#     "PFUnDA",
#     "PFDoDA",
#     "PFTrDA"
#   )),
#   nrow = 1) +
#   labs(x = "Time (years)",
#        y = "Plasma concentration (μg/L)") +
#   theme_quick
# 
# 
# 
# design <- "AAAAAAAAAA
#            BBBBB##          ###
#            ##CC###          ###"
# 
# ObsVsPred_one_plasma_UA <- ObsVsPred_PFCAs_plasma + labs(title = "PFCAs") +
#   theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5, "cm")
#   )
# 
# ObsVsPred_two_plasma_UA <- ObsVsPred_PFSAs_plasma + labs(title = "PFSAs") +
#   theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0.5, "cm")
#   )
# 
# ObsVsPred_three_plasma_UA <- ObsVsPred_PFECAs_plasma + labs(title = "PFECAs") +
#   theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0.5, "cm")
#   )
# 
# ObsVsPred_Plasma_plot <- (ObsVsPred_one_plasma_UA+Patchwork_theme) / 
#   (ObsVsPred_two_plasma_UA+Patchwork_theme) / 
#   (ObsVsPred_three_plasma_UA+Patchwork_theme) + 
#   plot_layout(design = design,
#               # nrow = 3,
#               guides = "collect",
#               axes = "collect",
#               axis_titles = "collect") &
#   Patchwork_theme +
#   theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, "cm"))
# 
# ggsave(plot = ObsVsPred_Plasma_plot,
#        file = here(OUTPUT, "ObsVsPred_Plasma_plot_UA.png"),
#        dpi = 600,
#        width = 24,
#        height = 12,
#        units = "cm")
# 

