# To find fup and GFR !


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


 
# Produce results ----


# functions ----

source(here("Script", "functions_calculate_results.R"))


# produce result files and plots for all PFAS ----
INPUT <- here("Output", "CleanInitialParams_UncertaintyAnalysis")
OUTPUT <- here("Output","Re_analyse_2ndUA")
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

Dose_df <- tibble(
  PFAS = PFAS_list,
  Dose = 0.1) 
# |>
#   mutate(Dose = ifelse(PFAS == "PFBS", 18.4, 
#                        ifelse(PFAS == "HFPO-DA" , 19.6,
#                               Dose)))

# PFAS_list <- Dose_df$PFAS

nPFAS <- length(PFAS_list)

for (i in 1:nPFAS) {
  
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}-n{i}-----------------------------------"))
  
  # Inputs ----
  
  # load results (called savey)
  load(here(INPUT, paste0(PFAS_i, "_LHS_MC_y.RData" )
  ))
  # load oral dose (different per PFAS if comparing results to Abraham)
  OD <- filter(Dose_df, PFAS == PFAS_i) |>
    select(Dose) |> pull()
  
  # TK results ----
  # calculate AUC, half-life, elimination rates, fraction absorbed etc
  # To see if there are any errors
  
  TK_results <- sapply(savey, check_error, OD)
  
  if(is.list(TK_results)){
    errored_iterations <- which(sapply(TK_results, function(x) length(x) < 3)) #which(is.na(TK_results[2,]))
    message(glue("Removing iterations n{errored_iterations} that give errors. \\"))
    savey <- savey[-errored_iterations] 
  }
  
  TK_results <- sapply(savey, calc_TK_results, OD) #savey
  
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
  
  errored_TK_result_iterations <- errored_TK_results$Iteration
  
  if(is.numeric(errored_TK_results$Iteration)){
    message(glue("Removing iterations n{errored_TK_result_iterations} that give unrealistic clearance values (CL > 1e5)"))
    write.csv(errored_TK_results, file = here(OUTPUT, paste0(PFAS_i, "errored_TK_results.csv")))
  }
  
  TK_summary_results <- calc_TK_summary_results(clean_TK_results) |>
    as.data.frame() |> mutate(PFAS = PFAS_i)
  write.csv(TK_summary_results, file = here(OUTPUT, paste0(PFAS_i, "_TK_sumstat.csv")))
  
  Results_df <- do.call(rbind, savey) 
  
  # Actual parameters ----
  
  Parameters_df <- Results_df |>
    select(GFR:CF_K_reabsorption) |>
    select(!starts_with("CF_")) |> 
    distinct() |> 
    mutate(fup_GFR = fup*GFR) |> 
    rename(Clear_reabsorption = Clear_net_reabsorption_unscaled) 
  
  write.csv(Parameters_df, file = here(OUTPUT, paste0(PFAS_i, "_Final_model_parameters.csv")))
  
  P_long_df <- Parameters_df |>
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
                          "tco",
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
                          "GFR",
                          "fup_GFR",
                          "Clear_OAT3",              
                          "Clear_OAT1",                      
                          "Clear_OAT4",              
                          "Clear_reabsorption",
                          "QT"
                        )), scales = "free") +
    labs(x = "Value",
         y = "Frequency") +
    coord_flip()+
    theme_quick
  
  Parameters_plot = Parameters_plot + Patchwork_theme +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(plot = Parameters_plot,
         file = here(OUTPUT, paste0(PFAS_i, "_final_model_parameters_plot.png")),
         dpi = 600,
         width = 24,
         height = 18,
         units = "cm")
  
}