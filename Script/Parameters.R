# --------------------------------------------------------------------------- #
# Script for deriving the model parameters
# CP, 30-06-2026
# --------------------------------------------------------------------------- #

# Pre-run ----
# Clear out the global environment
rm(list=ls()) 

# Load packages
library(readxl)
# library(xlsx)
library(here)
library(tidyverse)
library(showtext)
library(ggplot2)
library(purrr)
library(knitr)
library(patchwork)
library(dendroTools)
library(viridis)

# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'), #element_blank(), #
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), 
    plot.title = element_text(hjust = 0, size = 35, face = "bold"),
    plot.subtitle = element_text(hjust = 0, size = 30),
    plot.caption = element_text(hjust = 0, size = 30),
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
OUTPUT <- here("Output", "QSPR_results") 
dir.create(OUTPUT, recursive = TRUE)

# Initialise ----

# PFAS included in the study
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
               "PFUnDA")
PFAS_list <- sort(PFAS_list)

# PFAS OECD categories included in the study
OECD_cat_list <- c("PFCAs", "PFSAs", "PFECAs")

# Load in vitro and in silico input
All_input_df <- read_excel(here("Input", "PFAS_database.xlsx"), 
                           col_types = c("text", "text", "text", 
                                         "text", "text", "text", "text", "text", 
                                         "text", "text", "numeric", "numeric", 
                                         "text", "numeric", "text", "text", 
                                         "text", "text", "text", "text", "numeric", 
                                         "text", "text", "text", "text")
                                ) |> 
  mutate(Value = signif(Value, 3),
         SD = signif(SD, 3))

# Keep the input parameters for only the PFAS included in the study
Datab_df <- All_input_df %>% filter(Acronym %in% PFAS_list)

# Confirm that the database contains only the selected PFAS
Check <- unique(Datab_df$Acronym)
Check %in% PFAS_list 


# Toxicokinetic parameters ----

# Parameters that we are interested in for the study
Selected_params <- c("logD_AlbuminWater",
                     "logD_StructuralProteinWater",
                     "logD_MembraneLipidWater",
                     "logD_StorageLipidWater",
                     "logD_GlobulinWater",
                     "Papp",
                     "Vmax_OATP1B1","Vmax_OATP1B3","Vmax_OATP2B1","Vmax_NTCP","Vmax_OAT1","Vmax_OAT3","Vmax_OAT4",   
                     "Km_OAT1","Km_OAT3","Km_OAT4","Km_OATP1B1","Km_OATP1B3","Km_OATP2B1","Km_NTCP",   
                     "Ka_ASBT",
                     "Kd_FABP", "Kd_ASBT",
                     "UR_OAT1","UR_OAT2","UR_OAT3","UR_OAT4",
                     "Dbile_liver",
                     "freab_ASBT_NTCP",
                     "BPR")

# Molecular weight dataframe
MW_df <- Datab_df |> filter(Parameter == "MW") |> 
  select(Acronym, Category_OECD, Value) |> 
  rename(MW = Value, 
         PFAS = Acronym) |> 
  arrange(MW)


# Keep only the parameters we are interested in
Datab <- Datab_df %>% 
  filter(Parameter %in% Selected_params) |> 
  select(Acronym, Category_OECD, 
         PBK_mechanism, Parameter, Value, SD,  
         Parameter_specific, Parameter_source, Parameter_source_system, 
         Parameter_explanation, Parameter_comment,
         Reference) |> 
  rename(PFAS = Acronym)

# Add the MW column for the QSPR
Datab <- full_join(Datab, MW_df)

# Freabsorbed ----
CleanfReab <- Datab %>% 
  filter(Parameter == "freab_ASBT_NTCP")

CleanfReab |> ggplot(aes(x = MW, Value,
                        colour = Category_OECD)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = (Value - SD),
                    ymax = (Value + SD)),
                width = .2) +
  labs(y = "fraction reabsorbed",
       x = "MW (g/mol)") +
  theme_quick


QSPR_fReab <- CleanfReab 

y = QSPR_fReab$Value
x = QSPR_fReab$MW
data <- tibble(x = x, y = y)

fit_fReab <- lm(y ~ x, data = data)

Results_fReab <- summary(fit_fReab)
print(Results_fReab)

Formula <- "Intercept * X"
RSE <- signif(Results_fReab$sigma, 3)
R.squared <- signif(Results_fReab$r.squared, 3)
Intercept <- signif(Results_fReab$coefficients[1], 3)


Residuals_fReab <- tibble(
  x = x,
  y = residuals(fit_fReab)
)

Residuals_fReab_plot <- ggplot(Residuals_fReab) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_fReab_plot

predicted <- tibble(
  MW = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10),
  predicted_y = predict(fit_fReab, list(x = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

observed_fReab <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_fReab_plot <- QSPR_fReab |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = (Value - SD),
                    ymax = (Value + SD)),
                width = .2) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", RSE: ", RSE,
                        ", R.squared: ", R.squared,
                        ", Intercept: ", Intercept),
       y = "fReab", 
       x = "MW (g/mol)") +
  theme_quick
QSPR_fReab_plot

D <- QSPR_fReab_plot + CP_theme
R <- Residuals_fReab_plot + CP_theme

QSPR_fReab_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_fReab_plots,
       file = here(OUTPUT, "QSPR_fReab_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")

predicted_fReab <- tibble(
  MW = PFAS_MW_list,
  predicted_value = predict(fit_fReab, list(x = PFAS_MW_list))
) %>% arrange(MW)

collect_fReab <- predicted_fReab |> 
  rename(Value = predicted_value) |> 
  mutate(PBK_mechanism = "Transport",
         Parameter = "fractionReabsorbed", 
         SD = RSE, 
         SD_comment = NA,
         Specific = "predicted_linear",
         Distribution = "uniform", 
         Distribution_comment = NA
  )

write.csv(collect_fReab, here(OUTPUT,"PFAS_fReab.csv"))

# add PFAS 
predicted_OATP2B1 <- full_join(collect_OATP2B1, MW_df)

# Absorption ----

# Clean Papp values 

CleanPapp <- Datab %>% 
  filter(Parameter == "Papp") |> 
  mutate(Parameter_specific = str_replace_all(Parameter_specific, "NA", "in_vitro"),
         Parameter_source_system = str_replace_all(Parameter_source_system, " ", "-")) |> 
  unite(Specific, 
        Parameter_specific, Parameter_source_system,  
        sep = "__") |> 
  # removing the Papp__neutral as given the physiological pH and pKa of PFAS the predominant species in the human body is the ionised
  filter(Specific != "neutral__COSMOperm") |> 
  select(!c(Parameter_source, Parameter_explanation, Parameter_comment)) |> 
  mutate(Value = Value * 60*60*24) |>  # unit conversion, cm/s to cm/dau
  mutate(SD = SD * 60*60*24)   # unit conversion, cm/s to cm/dau

GlimpsePapp <- CleanPapp |> ggplot(aes(x = MW, Value,
           colour = Category_OECD,
           shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = (Value - SD),
                    ymax = (Value + SD)),
                width = .2) +
  labs(y = "Papp (cm/day)",
       x = "MW (g/mol)") +
  theme_quick
GlimpsePapp
# we can see in the graph that there is quite some variance between all these values
# two values seem to be outliers for PFOS, cosmo perm and equilibrium dialysis

GlimpsePapp_log <- CleanPapp |> ggplot(aes(log(MW), log(Value),
                                       colour = Category_OECD,
                                       shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 10, 8)) +
  geom_point() +
  geom_errorbar(aes(ymin = log(Value - SD),
                    ymax = log(Value + SD)),
                width = .002) +
  labs(y = "Papp ln(cm/day)",
       x = "MW ln(g/mol)") +
  theme_quick
GlimpsePapp_log
GlimpsePapp_log <- GlimpsePapp_log + CP_theme

ggsave(plot = GlimpsePapp_log,
       file = here(OUTPUT, "Papp_withCOSMOperm.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")
# the Papp values derived using COSMOperm seem to not follow the same trend as those derived in vitro so I will remove those

# remove the two COSMOperm values
CleanPapp <- CleanPapp |> filter(Specific != "ionic__COSMOperm")
# 
# GlimpsePapp <- CleanPapp |> ggplot(aes(x = MW, Value, 
#                                        colour = Category_OECD, 
#                                        shape = Specific)) +
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
#   geom_point() +
#   geom_errorbar(aes(ymin = (Value - SD), 
#                     ymax = (Value + SD)),
#                 width = .2) + 
#   labs(y = "Papp (cm/day)", 
#        x = "MW (g/mol)") +
#   theme_quick
# GlimpsePapp
# GlimpsePapp <- GlimpsePapp + CP_theme
# ggsave(plot = GlimpsePapp,
#        file = here(OUTPUT, "Papp.png"),
#        dpi = 600,
#        width = 12,
#        height = 9,
#        units = "cm")
# 
# GlimpsePapp_log <- CleanPapp |> ggplot(aes(log(MW), log(Value), 
#                                            colour = Category_OECD, 
#                                            shape = Specific)) +
#   scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
#   geom_point() +
#   geom_errorbar(aes(ymin = log(Value - SD), 
#                     ymax = log(Value + SD)),
#                 width = .002) + 
#   labs(y = "Papp ln(cm/day)", 
#        x = "MW ln(g/mol)") +
#   theme_quick
# GlimpsePapp_log
# These are the final values to be included for the Papp QSPR
# QSPR to be done in log-log scale
# GlimpsePapp_log <- GlimpsePapp_log + CP_theme
# ggsave(plot = GlimpsePapp_log, 
#        file = here(OUTPUT, "Papp_log.png"),
#        dpi = 600,
#        width = 12,
#        height = 9,
#        units = "cm")

# QSPR
glimpse(CleanPapp)

QSPR_CleanPapp <- CleanPapp |> 
  mutate(Value = log(Value), 
         MW = log(MW))

PFAS_MW_list <- MW_df$MW
PFAS_MW_list_log <- log(MW_df$MW)

y = QSPR_CleanPapp$Value
x = QSPR_CleanPapp$MW
data <- tibble(x = x, y = y)


fit_Papp <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
                 data = data,
                 start = list(Const = min(QSPR_CleanPapp$Value),
                              Asym = -2.8772,
                              xmid = 6.0187 ,
                              scal = -0.1451),
                 # control = nls.control(maxiter = 1000000),
                 algorithm='port')

Results_Papp <- summary(fit_Papp)
print(Results_Papp)

Formula <- Results_Papp$formula[3]
RSE <- signif(Results_Papp$sigma, 3)
Const <- signif(Results_Papp$coefficients[1], 3)
Asym <- signif(Results_Papp$coefficients[2], 3)
xmid <- signif(Results_Papp$coefficients[3], 3)
scal <- signif(Results_Papp$coefficients[4], 3)

Residuals_Papp <- tibble(
  x = x,
  y = residuals(fit_Papp)
)

Residuals_Papp_plot <- ggplot(Residuals_Papp) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick

predicted <- tibble(
  MW = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001),
  predicted_y = predict(fit_Papp, list(x = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
  ) %>% arrange(MW)

predicted_Papp <- tibble(
  MW = PFAS_MW_list_log,
  predicted_value = predict(fit_Papp, list(x = PFAS_MW_list_log))
  ) %>% arrange(MW)

observed_Papp <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_CleanPapp2 <- QSPR_CleanPapp |> 
  mutate(nValue = exp(Value))

QSPR_Papp_plot <- QSPR_CleanPapp2 |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = log(nValue - SD),
                    ymax = log(nValue + SD)),
                width = .002) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                         ", const: ", Const,
                         ", asym: ", Asym,
                         ", xmid: ", xmid,
                         ", scal: ", scal),
       y = "Papp ln(cm/day)", 
       x = "MW ln(g/mol)") +
  theme_quick

D <- QSPR_Papp_plot + CP_theme
R <- Residuals_Papp_plot + CP_theme

QSPR_Papp_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_Papp_plots,
       file = here(OUTPUT, "QSPR_Papp_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")

# collect results
convert_Papp <- predicted_Papp |> 
  rename(MWlog = MW,
         meanlog = predicted_value) |> 
  mutate(MW = signif(exp(MWlog),3),
         PBK_mechanism = "Absorption",
         Parameter = "Papp_SI", 
         Value = signif(meanlog,3), # keep in in log scale!
         SD = RSE, 
         SD_comment = "is sdlog",
         Specific = "predicted_sigmoid",
         Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and sdlog is SD") |> 
  select(! c(MWlog, meanlog))

# add PFAS 

convert_Papp <- full_join(convert_Papp, MW_df)

# merge dataframes
CleanPapp <- CleanPapp |> 
  mutate(Value = log(Value),
         SD = log(SD),
         Parameter = "Papp_SI",
         Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and SD of log Values is SD", 
         SD_comment = "NA")

# Final dataframe
FinalPapp <- full_join(CleanPapp, convert_Papp) 



# Transport ----

CleanTransport <- Datab %>% 
  filter(PBK_mechanism == "Transport") |> 
  mutate(Parameter_source_system = str_replace_all(Parameter_source_system, "Calculated based on binding affinity of TCA, assuming equal contribution of ASBT and NTCP", "MolecularDocking"),
         Parameter_source_system = str_replace_all(Parameter_source_system, "Calculated based on binding affinity of TCA_ assuming equal contribution of ASBT and NTCP", "MolecularDocking"),
         Parameter_source_system = str_replace_all(Parameter_source_system, "Molecular docking Random Forest", "MolecularDocking")) |> 
  unite(Specific, 
      Parameter_specific, Parameter_source_system,  
      sep = "__") |>
  select(!c(Parameter_explanation, Parameter_comment))


# Keep only in vitro values
CleanTransport <- CleanTransport |> filter(Parameter_source == "in vitro") |> 
  select(!Parameter_source)

# Keep only Vmax and Km
CleanTransport <- CleanTransport |> 
  separate(Parameter, 
           into = c("Parameter", "Transporter")) |> 
  filter(Parameter %in% c("Vmax", "Km")) 

# Calculate clearances
CleanTransport <- CleanTransport |> 
  rename(Mean = Value) 

MeanTransport <- CleanTransport |> select(!SD) |> 
  pivot_wider(values_from = Mean, names_from = Parameter)
SDTransport <- CleanTransport |> select(!Mean) |> 
  pivot_wider(values_from = SD, names_from = Parameter) |> 
  rename(SD_Km = Km, 
         SD_Vmax = Vmax)

# Vmax < LOQ for:
# OAT4: PFBS
# OAT3: PFBS, PFOS
# OAT1: PFBS, PFHxS, PFOS, PFNA
# To exclude from QSPR

CleanTransport <- full_join(MeanTransport, SDTransport) |> 
  separate(Specific, into = c("Delete", "Specific")) |> 
  #unit correction Vmax: nmol/mg protein/min -> umol/mg protein/day
  mutate(Vmax = Vmax*1e-3*60*24, 
         SD_Vmax = SD_Vmax*1e-3*60*24) |> 
  # CL = umol/mg protein/day / umol/L => L/day/mg protein
  mutate(CL = Vmax/Km) |> 
  # propagation of standard deviation (ΔCL = CL * sqrt((ΔVmax/Vmax)^2 + (ΔKm/Km)^2))
  mutate(SD_CL = Vmax/Km * sqrt((SD_Vmax/Vmax)**2 + (SD_Km/Km)**2)) |> 
  # in the database, values reported as < LOQ were assumed to be very small 1e-10, 
  # then converted 1e-13 due to the unit conversion above, here putting them to zero
  mutate(CL = ifelse(Vmax <= 1e-1, 0, CL), 
         SD_CL = ifelse(Vmax<=1e-1, 0, SD_CL)) |> 
  select(!Delete) |> 
  mutate(Transporter = factor(Transporter, levels = c("OATP2B1", "OATP1B1", "OATP1B3", "NTCP", "OAT1", "OAT3", "OAT4")))

GlimpseTransport <- CleanTransport |> ggplot(aes(x = MW, y = CL, 
                                       colour = Category_OECD, 
                                       shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = (CL - SD_CL), 
                    ymax = (CL + SD_CL)),
                width = .2) + 
  labs(y = "Clearance (L/day/mg protein)", 
       x = "MW (g/mol)") +
  facet_wrap(~Transporter, scales = "free") +
  theme_quick
GlimpseTransport

GlimpseTransport <- GlimpseTransport + CP_theme
ggsave(plot = GlimpseTransport,
       file = here(OUTPUT, "Transport.png"),
       dpi = 600,
       width = 18,
       height = 14,
       units = "cm")

# We can see that in general where we have values, PFSAs and PFCAs have similar values
# PFOA could be an outlier in regard to NTCP
# OAT1 seem to only be involved for PFCAs, but there is a large uncertainty
# There is also a large uncertainty for the OAT3 clearance as well. 

# Removing the values below LOQ
CleanTransport <- CleanTransport |> 
  filter(CL > 0) |>
  # Checking if the SD is higher than CL and putting it to equal to CL
  mutate(SD_CL = ifelse(SD_CL >= CL, CL-0.0000001, SD_CL)) 

# Calculate ymin and ymax
CleanTransport <- CleanTransport |> 
  mutate(ymin = CL-SD_CL,
         ymax = CL+SD_CL)

GlimpseTransport <- CleanTransport |> ggplot(aes(x = MW, y = CL, 
                                                 colour = Category_OECD, 
                                                 shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = ymin, 
                    ymax = ymax),
                width = .2) + 
  labs(y = "Clearance (L/day/mg protein)", 
       x = "MW (g/mol)") +
  facet_wrap(~Transporter, scales = "free") +
  theme_quick
GlimpseTransport

GlimpseTransport_log <- CleanTransport |> ggplot(aes(log(MW), log(CL), 
                                           colour = Category_OECD, 
                                           shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = log(ymin), 
                    ymax = log(ymax)),
                width = .002) + 
  labs(y = "Clearance ln(L/day/mg protein)", 
       x = "MW ln(g/mol)") +
  facet_wrap(~Transporter, scales = "free") +
  theme_quick
GlimpseTransport_log
# Log transformation doesn't really change the trends

GlimpseTransport_log <- GlimpseTransport_log + CP_theme
ggsave(plot = GlimpseTransport_log,
       file = here(OUTPUT, "Transport_log.png"),
       dpi = 600,
       width = 18,
       height = 14,
       units = "cm")

# Saving the clean dataframe with the clearances
CleanTransport <- CleanTransport |> 
  select(!c(ymin, ymax)) 

Mean2Transport <- CleanTransport |> 
  select(!c(SD_CL, SD_Vmax, SD_Km)) |> 
  pivot_longer(Km:CL, names_to = "Parameter", values_to = "Value")

SD2Transport <- CleanTransport |> 
  select(!c(CL, Vmax, Km)) |> 
  pivot_longer(SD_Km:SD_CL, names_to = "Parameter", values_to = "SD") |> 
  mutate(Parameter = str_replace_all(Parameter, "SD_", ""))

Clean2Transport <- full_join(Mean2Transport, SD2Transport)


## CL_OATP2B1 ----
QSPR_OATP2B1 <- Clean2Transport |>
  filter(Transporter == "OATP2B1", 
         Parameter == "CL")

rm(x)
rm(y)
rm(data)
y = QSPR_OATP2B1$Value
x = QSPR_OATP2B1$MW
data <- tibble(x = x, y = y)

fit_OATP2B1 <- lm(y ~ x, data = data)

Results_OATP2B1 <- summary(fit_OATP2B1)
print(Results_OATP2B1)

rm(Formula)
rm(RSE)
rm(Const)
rm(Asym)
rm(xmid)
rm(scal)

Formula <- "Intercept * X"
RSE <- signif(Results_OATP2B1$sigma, 3)
R.squared <- signif(Results_OATP2B1$r.squared, 3)
Intercept <- signif(Results_OATP2B1$coefficients[1], 3)


Residuals_OATP2B1 <- tibble(
  x = x,
  y = residuals(fit_OATP2B1)
)

Residuals_OATP2B1_plot <- ggplot(Residuals_OATP2B1) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_OATP2B1_plot
rm(predicted)
predicted <- tibble(
  MW = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10),
  predicted_y = predict(fit_OATP2B1, list(x = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

observed_OATP2B1 <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_OATP2B1_plot <- QSPR_OATP2B1 |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = (Value - SD),
                    ymax = (Value + SD)),
                width = .2) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", RSE: ", RSE,
                        ", R.squared: ", R.squared,
                        ", Intercept: ", Intercept),
       y = "OATP2B1 clearance (L/day/mg protein)", 
       x = "MW (g/mol)") +
  theme_quick
QSPR_OATP2B1_plot
D <- QSPR_OATP2B1_plot + CP_theme
R <- Residuals_OATP2B1_plot + CP_theme

QSPR_OATP2B1_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_OATP2B1_plots,
       file = here(OUTPUT, "QSPR_OATP2B1_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")


# collect results
predicted_OATP2B1 <- tibble(
  MW = PFAS_MW_list,
  predicted_value = predict(fit_OATP2B1, list(x = PFAS_MW_list))
) %>% arrange(MW) |> 
  mutate(predicted_value = ifelse(predicted_value <0, 0, predicted_value))

collect_OATP2B1 <- predicted_OATP2B1 |> 
  rename(Value = predicted_value) |> 
  mutate(PBK_mechanism = "Transport",
         Parameter = "CL_OATP2B1c", 
         SD = RSE, 
         SD_comment = NA,
         Specific = "predicted_linear",
         Distribution = "uniform", 
         Distribution_comment = NA
         )

# add PFAS 
predicted_OATP2B1 <- full_join(collect_OATP2B1, MW_df)

# merge dataframes
QSPR_OATP2B1 <- QSPR_OATP2B1 |> 
  mutate(Parameter = "CL_OATP2B1c",
         Distribution = "uniform", 
         Distribution_comment = NA, 
         SD_comment = NA) |> 
  select(!Transporter)

# Final dataframe
FinalOATP2B1 <- full_join(QSPR_OATP2B1, predicted_OATP2B1) 


## CL_OATP1B1 ----

QSPR_OATP1B1 <- Clean2Transport |>
  filter(Transporter == "OATP1B1", 
         Parameter == "CL")

rm(x)
rm(y)
rm(data)
y = QSPR_OATP1B1$Value
x = QSPR_OATP1B1$MW
data <- tibble(x = x, y = y)

fit_OATP1B1 <- lm(y ~ x, data = data)

Results_OATP1B1 <- summary(fit_OATP1B1)
print(Results_OATP1B1)

rm(Formula)
rm(RSE)
rm(Const)
rm(Asym)
rm(xmid)
rm(scal)

Formula <- "Intercept * X"
RSE <- signif(Results_OATP1B1$sigma, 3)
R.squared <- signif(Results_OATP1B1$r.squared, 3)
Intercept <- signif(Results_OATP1B1$coefficients[1], 3)


Residuals_OATP1B1 <- tibble(
  x = x,
  y = residuals(fit_OATP1B1)
)

Residuals_OATP1B1_plot <- ggplot(Residuals_OATP1B1) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_OATP1B1_plot
rm(predicted)
predicted <- tibble(
  MW = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10),
  predicted_y = predict(fit_OATP1B1, list(x = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

observed_OATP1B1 <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_OATP1B1_plot <- QSPR_OATP1B1 |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = (Value - SD),
                    ymax = (Value + SD)),
                width = .2) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", RSE: ", RSE,
                        ", R.squared: ", R.squared,
                        ", Intercept: ", Intercept),
       y = "OATP1B1 clearance (L/day/mg protein)", 
       x = "MW (g/mol)") +
  theme_quick
QSPR_OATP1B1_plot
D <- QSPR_OATP1B1_plot + CP_theme
R <- Residuals_OATP1B1_plot + CP_theme

QSPR_OATP1B1_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_OATP1B1_plots,
       file = here(OUTPUT, "QSPR_OATP1B1_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")


# collect results
predicted_OATP1B1 <- tibble(
  MW = PFAS_MW_list,
  predicted_value = predict(fit_OATP1B1, list(x = PFAS_MW_list))
) %>% arrange(MW) |> 
  mutate(predicted_value = ifelse(predicted_value <0, 0, predicted_value))

collect_OATP1B1 <- predicted_OATP1B1 |> 
  rename(Value = predicted_value) |> 
  mutate(PBK_mechanism = "Transport",
         Parameter = "CL_OATP1B1c", 
         SD = RSE, 
         SD_comment = NA,
         Specific = "predicted_linear",
         Distribution = "uniform", 
         Distribution_comment = NA)

# add PFAS 
predicted_OATP1B1 <- full_join(collect_OATP1B1, MW_df)

# merge dataframes
QSPR_OATP1B1 <- QSPR_OATP1B1 |> 
  mutate(Parameter = "CL_OATP1B1c",
         Distribution = "uniform", 
         Distribution_comment = NA, 
         SD_comment = NA) |> 
  select(!Transporter)

# Final dataframe
FinalOATP1B1 <- full_join(QSPR_OATP1B1, predicted_OATP1B1)




## CL_OATP1B3 ----

QSPR_OATP1B3 <- Clean2Transport |>
  filter(Transporter == "OATP1B3", 
         Parameter == "CL")

rm(x)
rm(y)
rm(data)
y = QSPR_OATP1B3$Value
x = QSPR_OATP1B3$MW
data <- tibble(x = x, y = y)

fit_OATP1B3 <- lm(y ~ x, data = data)

Results_OATP1B3 <- summary(fit_OATP1B3)
print(Results_OATP1B3)

rm(Formula)
rm(RSE)
rm(Const)
rm(Asym)
rm(xmid)
rm(scal)

Formula <- "Intercept * X"
RSE <- signif(Results_OATP1B3$sigma, 3)
R.squared <- signif(Results_OATP1B3$r.squared, 3)
Intercept <- signif(Results_OATP1B3$coefficients[1], 3)


Residuals_OATP1B3 <- tibble(
  x = x,
  y = residuals(fit_OATP1B3)
)

Residuals_OATP1B3_plot <- ggplot(Residuals_OATP1B3) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_OATP1B3_plot
rm(predicted)
predicted <- tibble(
  MW = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10),
  predicted_y = predict(fit_OATP1B3, list(x = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

observed_OATP1B3 <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_OATP1B3_plot <- QSPR_OATP1B3 |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = (Value - SD),
                    ymax = (Value + SD)),
                width = .2) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", RSE: ", RSE,
                        ", R.squared: ", R.squared,
                        ", Intercept: ", Intercept),
       y = "OATP1B3 clearance (L/day/mg protein)", 
       x = "MW (g/mol)") +
  theme_quick
QSPR_OATP1B3_plot
D <- QSPR_OATP1B3_plot + CP_theme
R <- Residuals_OATP1B3_plot + CP_theme

QSPR_OATP1B3_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_OATP1B3_plots,
       file = here(OUTPUT, "QSPR_OATP1B3_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")



#nls
fit_OATP1B3_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
                       data = data,
                       start = list(Const = 0,
                                    Asym = 0.04591,
                                    xmid = 349.85558  ,
                                    scal = 25.03187 ),
                       # control = nls.control(maxiter = 1000000),
                       algorithm='port')

Results_OATP1B3_nls <- summary(fit_OATP1B3_nls)
print(Results_OATP1B3_nls)

Formula_nls <- "Intercept * X"
RSE_nls <- signif(Results_OATP1B3_nls$sigma, 3)
Const_nls <- signif(Results_OATP1B3_nls$coefficients[1], 3)
Asym_nls <- signif(Results_OATP1B3_nls$coefficients[2], 3)
xmid_nls <- signif(Results_OATP1B3_nls$coefficients[3], 3)
scal_nls <- signif(Results_OATP1B3_nls$coefficients[4], 3)

Residuals_OATP1B3_nls <- tibble(
  x = x,
  y = residuals(fit_OATP1B3_nls)
)

Residuals_OATP1B3_nls_plot <- ggplot(Residuals_OATP1B3_nls) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE_nls)) +
  theme_quick
Residuals_OATP1B3_nls_plot
rm(predicted)
predicted_nls <- tibble(
  MW = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10),
  predicted_y = predict(fit_OATP1B3_nls, list(x = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10))),
  min = predicted_y - RSE_nls, 
  max = predicted_y + RSE_nls
) %>% arrange(MW)


observed_OATP1B3 <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_OATP1B3_nls_plot <- QSPR_OATP1B3 |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  # geom_errorbar(aes(MW, Value, 
  #                   colour = Category_OECD,
  #                   ymin = (Value - SD),
  #                   ymax = (Value + SD)),
  #               width = .2) +
  geom_path(data = predicted_nls, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted_nls, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula_nls, 
                        ", const: ", Const_nls,
                        ", asym: ", Asym_nls,
                        ", xmid: ", xmid_nls,
                        ", scal: ", scal_nls),,
       y = "OATP1B3 clearance (L/day/mg protein)", 
       x = "MW (g/mol)") +
  theme_quick
QSPR_OATP1B3_nls_plot
D <- QSPR_OATP1B3_nls_plot + CP_theme
R <- Residuals_OATP1B3_nls_plot + CP_theme

QSPR_OATP1B3_nls_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_OATP1B3_nls_plots,
       file = here(OUTPUT, "QSPR_OATP1B3_nls_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")


# collect results

predicted_OATP1B3 <- tibble(
  MW = PFAS_MW_list,
  predicted_value = predict(fit_OATP1B3, list(x = PFAS_MW_list))
) %>% arrange(MW) |> 
  mutate(predicted_value = ifelse(predicted_value <0, 0, predicted_value))


collect_OATP1B3 <- predicted_OATP1B3 |>  
  rename(Value = predicted_value) |> 
  mutate(PBK_mechanism = "Transport",
         Parameter = "CL_OATP1B3c", 
         Value = ifelse(Value <0, 0, Value),
         SD = RSE,  
         SD_comment = NA,
         Specific = "predicted",
         Model = "linear",
         Distribution = "uniform", 
         Distribution_comment = NA)


# add PFAS 
predicted_OATP1B3 <- full_join(collect_OATP1B3, MW_df)


predicted_OATP1B3_nls <- tibble(
  MW = PFAS_MW_list,
  predicted_value = predict(fit_OATP1B3_nls, list(x = PFAS_MW_list))
) %>% arrange(MW) |> 
  mutate(predicted_value = ifelse(predicted_value <0, 0, predicted_value))


collect_OATP1B3_nls <- predicted_OATP1B3_nls |> 
  rename(Value = predicted_value) |> 
  mutate(PBK_mechanism = "Transport",
         Parameter = "CL_OATP1B3c", 
         Value = ifelse(Value <0, 0, Value),
         SD = 0.04, # based on the experimental data, not the prediction , 
         SD_comment = NA,
         Specific = "predicted",
         Model = "sigmoid",
         Distribution = "uniform", 
         Distribution_comment = NA)

# add PFAS 
predicted_OATP1B3_nls <- full_join(collect_OATP1B3_nls, MW_df)

# merge predicted dfs
predicted_OATP1B3 <- rbind(predicted_OATP1B3, predicted_OATP1B3_nls) |> 
  unite(Specific, 
        Specific, Model, sep = "_")

# merge dataframes
QSPR_OATP1B3 <- QSPR_OATP1B3 |> 
  mutate(Parameter = "CL_OATP1B3c",
         Distribution = "uniform", 
         Distribution_comment = NA, 
         SD_comment = NA) |> 
  select(!Transporter)


# Final dataframe
FinalOATP1B3 <- full_join(QSPR_OATP1B3, predicted_OATP1B3)


## CL_NTCP ----

QSPR_NTCP <- Clean2Transport |>
  filter(Transporter == "NTCP", 
         Parameter == "CL")

rm(x)
rm(y)
rm(data)
y = QSPR_NTCP$Value
x = QSPR_NTCP$MW
data <- tibble(x = x, y = y)

fit_NTCP <- lm(y ~ x, data = data)

Results_NTCP <- summary(fit_NTCP)
print(Results_NTCP)

rm(Formula)
rm(RSE)
rm(Const)
rm(Asym)
rm(xmid)
rm(scal)

Formula <- "Intercept * X"
RSE <- signif(Results_NTCP$sigma, 3)
R.squared <- signif(Results_NTCP$r.squared, 3)
Intercept <- signif(Results_NTCP$coefficients[1], 3)


Residuals_NTCP <- tibble(
  x = x,
  y = residuals(fit_NTCP)
)

Residuals_NTCP_plot <- ggplot(Residuals_NTCP) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_NTCP_plot
rm(predicted)
predicted <- tibble(
  MW = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10),
  predicted_y = predict(fit_NTCP, list(x = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

observed_NTCP <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_NTCP_plot <- QSPR_NTCP |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = (Value - SD),
                    ymax = (Value + SD)),
                  width = .2) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", RSE: ", RSE,
                        ", R.squared: ", R.squared,
                        ", Intercept: ", Intercept),
       y = "NTCP clearance (L/day/mg protein)", 
       x = "MW (g/mol)") +
  theme_quick
QSPR_NTCP_plot
D <- QSPR_NTCP_plot + CP_theme
R <- Residuals_NTCP_plot + CP_theme

QSPR_NTCP_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_NTCP_plots,
       file = here(OUTPUT, "QSPR_NTCP_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")


# collect results
predicted_NTCP <- tibble(
  MW = PFAS_MW_list,
  predicted_value = predict(fit_NTCP, list(x = PFAS_MW_list))
) %>% arrange(MW) |> 
  mutate(predicted_value = ifelse(predicted_value <0, 0, predicted_value))

collect_NTCP <- predicted_NTCP |> 
  rename(Value = predicted_value) |> 
  mutate(PBK_mechanism = "Transport",
         Parameter = "CL_NTCPc", 
         SD = RSE, 
         SD_comment = NA,
         Specific = "predicted_linear",
         Distribution = "uniform", 
         Distribution_comment = NA)

# add PFAS 
predicted_NTCP <- full_join(collect_NTCP, MW_df)

# merge dataframes
QSPR_NTCP <- QSPR_NTCP |> 
  mutate(Parameter = "CL_NTCPc",
         Distribution = "uniform", 
         Distribution_comment = NA, 
         SD_comment = NA) |> 
  select(!Transporter)

# Final dataframe
FinalNTCP <- full_join(QSPR_NTCP, predicted_NTCP)



## CL_OAT1 ----
OAT1_df <- Clean2Transport |>
  filter(Transporter == "OAT1",
         Parameter == "CL") |> 
  mutate(Specific = "measured", 
         Parameter = "CL_OAT1c",
         Distribution = "uniform", 
         Distribution_comment = NA, 
         SD_comment = NA) |> 
  select(!Transporter)

Pred_OAT1 <- MW_df |>
  mutate(Specific = "predicted_assumed", 
         PBK_mechanism = "Transport",
         Parameter = "CL_OAT1c",
         Value = 0.024,
         SD = 0.024, 
         Distribution = "uniform",
         Distribution_comment = NA, 
         SD_comment = NA)

FinalOAT1 <- full_join(OAT1_df, Pred_OAT1)

## CL_OAT3 ----
OAT3_df <- Clean2Transport |>
  filter(Transporter == "OAT3",
         Parameter == "CL") |> 
  mutate(Specific = "measured", 
         Parameter = "CL_OAT3c",
         Distribution = "uniform", 
         Distribution_comment = NA, 
         SD_comment = NA) |> 
  select(!Transporter)

Pred_OAT3 <- MW_df |>
mutate(Specific = "predicted_assumed", 
       PBK_mechanism = "Transport",
       Parameter = "CL_OAT3c",
       Value = 0.019,
       SD = 0.019, 
       Distribution = "uniform",
       Distribution_comment = NA, 
       SD_comment = NA)


FinalOAT3 <- full_join(OAT3_df, Pred_OAT3)


## CL_OAT4 ----

QSPR_OAT4 <- Clean2Transport |>
  filter(Transporter == "OAT4", 
         Parameter == "CL")

rm(x)
rm(y)
rm(data)
y = QSPR_OAT4$Value
x = QSPR_OAT4$MW
data <- tibble(x = x, y = y)

fit_OAT4 <- lm(y ~ x, data = data)

Results_OAT4 <- summary(fit_OAT4)
print(Results_OAT4)

rm(Formula)
rm(RSE)
rm(Const)
rm(Asym)
rm(xmid)
rm(scal)

Formula <- "Intercept * X"
RSE <- signif(Results_OAT4$sigma, 3)
R.squared <- signif(Results_OAT4$r.squared, 3)
Intercept <- signif(Results_OAT4$coefficients[1], 3)


Residuals_OAT4 <- tibble(
  x = x,
  y = residuals(fit_OAT4)
)

Residuals_OAT4_plot <- ggplot(Residuals_OAT4) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_OAT4_plot
rm(predicted)
predicted <- tibble(
  MW = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10),
  predicted_y = predict(fit_OAT4, list(x = seq(min(PFAS_MW_list), max(PFAS_MW_list), by = 10))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

observed_OAT4 <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_OAT4_plot <- QSPR_OAT4 |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = (Value - SD),
                    ymax = (Value + SD)),
                width = .2) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", RSE: ", RSE,
                        ", R.squared: ", R.squared,
                        ", Intercept: ", Intercept),
       y = "OAT4 clearance (L/day/mg protein)", 
       x = "MW (g/mol)") +
  theme_quick
QSPR_OAT4_plot
D <- QSPR_OAT4_plot + CP_theme
R <- Residuals_OAT4_plot + CP_theme

QSPR_OAT4_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_OAT4_plots,
       file = here(OUTPUT, "QSPR_OAT4_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")


# collect results
predicted_OAT4 <- tibble(
  MW = PFAS_MW_list,
  predicted_value = predict(fit_OAT4, list(x = PFAS_MW_list))
) %>% arrange(MW) |> 
  mutate(predicted_value = ifelse(predicted_value <0, 0, predicted_value))

collect_OAT4 <- predicted_OAT4 |> 
  rename(Value = predicted_value) |> 
  mutate(PBK_mechanism = "Transport",
         Parameter = "CL_OAT4c", 
         SD = RSE, 
         SD_comment = NA,
         Specific = "predicted_sigmoid",
         Distribution = "uniform", 
         Distribution_comment = NA)

# add PFAS 
predicted_OAT4 <- full_join(collect_OAT4, MW_df)

# merge dataframes
QSPR_OAT4 <- QSPR_OAT4 |> 
  mutate(Parameter = "CL_OAT4c",
         Distribution = "uniform", 
         Distribution_comment = NA, 
         SD_comment = NA) |> 
  select(!Transporter)

# Final dataframe
FinalOAT4 <- full_join(QSPR_OAT4, predicted_OAT4)

## CL_BSEP ----
Vmax_BSEP <- mean(c(5.8*60*24, 7*60*24, 8.4*60*24))                             # umol/min->day/mg BSEP, Average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011 
SD_Vmax_BSEP <- sd(c(5.8*60*24, 7*60*24, 8.4*60*24))

Km_BSEP <- mean(c(4.3, 25.9, 19))                                               # umol/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
SD_Km_BSEP <- sd(c(4.3, 25.9, 19))

CL_BSEP <- Vmax_BSEP/Km_BSEP                                                    # L/day
SD_CL_BSEP <- Vmax_BSEP/Km_BSEP * sqrt((SD_Vmax_BSEP/Vmax_BSEP)**2 + (SD_Km_BSEP/Km_BSEP)**2)

FinalBSEP <- tibble(
  PFAS = MW_df$PFAS,
  Category_OECD = MW_df$Category_OECD,
  PBK_mechanism = "Transport",
  Specific = "predicted_assumed",
  Reference = "De Bruijn et al. 2023",
  MW = MW_df$MW,
  Parameter = "CL_BSEPc",
  Value = CL_BSEP, 
  SD = SD_CL_BSEP,
  Distribution = "uniform", 
  Distribution_comment = NA,
  SD_comment = NA
)

## Scaling factors ----

# Input transporter protein abundance in different tissues
ProteinAb <-  read_excel(here("Input", "ProteinAbundance.xlsx"))
ProteinAb_df <- ProteinAb %>% mutate(Reference = factor(ProteinAb$Reference),
                                     Transporter = factor(ProteinAb$Transporter),
                                     Localisation = factor(ProteinAb$Localisation),
                                     min = Protein_abundance-SD,
                                     max = Protein_abundance+SD) %>% 
  mutate(min = pmax(min, 0.000)) 

# Summarise protein abundance
Meta_df <- ProteinAb_df %>%
  group_by(Transporter, Localisation) %>%
  summarise(.groups = "drop",
            minimum = min(min, na.rm = TRUE), 
            maximum = max(max, na.rm = TRUE),
            mean = mean(Protein_abundance, na.rm = TRUE), 
            sd = sd(Protein_abundance, na.rm = TRUE)) |> 
  filter(Localisation %in% c("intestine", "liver", "kidney"), 
         Transporter != "ASBT") 

# Calculate scaling factors
Calc_SF_df <- Meta_df %>% 
  mutate(minimum = pmax(minimum, 0.001)) %>% 
  mutate(vitro = c(1001.82, 1.202, 2.596, 0.711, 0.063, 0.120, 0.719, 0.812, 0.812),
         sd_vitro = vitro, # standard deviation not reported therefore assuming that sd equals the mean
         mgProtein.mgT = c(0.18, 0.18, 0.17, 0.17, 0.17, 0.18, 0.18, 0.13, 0.18),
         sd_mgProtein.mgT = c(0.06, 0.06, 0.03, 0.03, 0.03, 0.06, 0.06, 0.02, 0.06)) |>  # only minimum and maximum available, therefore taking the maximum deviation 
  mutate(sd = ifelse(is.na(sd), mean, sd)) %>%
  rename(sd_vivo = sd,
         vivo = mean) |>
  select(!c(minimum, maximum)) |> 
  mutate(REF = vivo/vitro,
         SF = mgProtein.mgT*REF*1e6) |> 
  mutate(sd_SF = SF * sqrt((sd_vivo/vivo)**2 + (sd_vitro/vitro)**2 + (sd_mgProtein.mgT/mgProtein.mgT)**2)) |> 
  mutate(Transporter = factor(Transporter, levels = c("OATP2B1", "OATP1B1", "OATP1B3", "NTCP", "BSEP", "OAT1", "OAT3", "OAT4")))


SF_plot <-  Calc_SF_df |> 
  mutate(SF = SF*1e-07, 
         sd_SF = sd_SF*1e-07) |> 
  ggplot(aes(x = Transporter, y = SF)) +
  geom_col(stat = "identity",
           fill = "grey50", 
           linewidth = 1,
           width = 0.5) +
  geom_errorbar(aes(ymin = SF-sd_SF,
                    ymax = SF+sd_SF),
                  # size = 0.3,
                  width = 0.0000000001
                  ) +
  scale_color_viridis(option = "H", discrete = TRUE) +
  facet_wrap(~ Localisation, scales = "free"
             ) +
  labs(x = "",
       y = "Scaling Factor (1e+07)") +
  theme_quick


BSEP_SF_plot <- Calc_SF_df |> 
  mutate(SF = SF*1e-07, 
         sd_SF = sd_SF*1e-07) |> 
  filter(Transporter == "BSEP") |> 
  ggplot(aes(x = Transporter, y = SF)) +
  geom_col(stat = "identity",
           fill = "grey50", 
           linewidth = 1,
           width = 0.5) +
  geom_errorbar(aes(ymin = SF-sd_SF,
                    ymax = SF+sd_SF),
                # size = 0.3,
                width = 0.0000000001
  ) +
  scale_color_viridis(option = "H", discrete = TRUE) +
  facet_wrap(~ Localisation, scales = "free"
  ) +
  labs(x = "",
       y = "Scaling Factor (1e+07)") +
  theme_quick


save_SF_plot <- SF_plot + CP_theme +
  theme(axis.text.x = element_text(size = 35, angle = 30)
 )
save_BSEP_SF_plot <- BSEP_SF_plot + CP_theme +
  theme(axis.text.x = element_text(size = 35, angle = 30)
  )

ggsave(plot = save_SF_plot,
       file = here(OUTPUT,"ScalingFactors.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

ggsave(plot = save_BSEP_SF_plot,
       file = here(OUTPUT,"BSEPScalingFactors.png"),
       dpi = 600,
       width = 4,
       height = 9,
       units = "cm")


SaveSF <- Calc_SF_df |> 
  select(Transporter, Localisation, SF, sd_SF) |> 
  rename(Value = SF, 
         SD = sd_SF, 
         Parameter = Localisation) |> 
  mutate(Parameter = paste0("SF_", Parameter)) |> 
  unite(Parameter,
        Parameter, Transporter, 
        sep = "_") |> 
  mutate(Distribution = "uniform", 
         Specific = "calculated",
         PFAS = list(PFAS_list)) |> 
  unnest(PFAS) |> 
  mutate(PBK_mechanism = "Distribution",
         Distribution_comment = NA,
         SD_comment = NA,
         Reference = "NA") |> 
  mutate(Parameter = str_replace_all(Parameter, "SF_liver_BSEP", "SF_BSEP"),
         Parameter = str_replace_all(Parameter, "SF_liver_NTCP", "SF_NTCP"),
         Parameter = str_replace_all(Parameter, "SF_kidney_OAT1", "SF_OAT1"),
         Parameter = str_replace_all(Parameter, "SF_kidney_OAT3", "SF_OAT3"),
         Parameter = str_replace_all(Parameter, "SF_kidney_OAT4", "SF_OAT4"),
         Parameter = str_replace_all(Parameter, "SF_liver_OATP1B1", "SF_OATP1B1"),
         Parameter = str_replace_all(Parameter, "SF_liver_OATP1B3", "SF_OATP1B3"),
         Parameter = str_replace_all(Parameter, "SF_intestine_OATP2B1", "SF_OATP2B1_I"),
         Parameter = str_replace_all(Parameter, "SF_liver_OATP2B1", "SF_OATP2B1"))

FinalSF <- full_join(SaveSF, MW_df)



# Distribution ----

# Create df
Distribution <- Datab %>%
  filter(PBK_mechanism == "Distribution") %>%
  select(PFAS, Category_OECD, 
         PBK_mechanism, Parameter, 
         Value, SD, MW, Parameter_source_system, Reference) |> 
  mutate(SD = as.numeric(SD), 
         # MW = ln(MW),# to log as values are also in log scale
         Parameter_source_system = str_replace_all(Parameter_source_system, " ", "-")) |> 
  filter(!Parameter %in% c("BPR", "Kd_FABP")) |> 
  rename(Specific = Parameter_source_system) |> 
  mutate(Value = log(10^Value),
         SD = log(10^SD))

Distribution <- Distribution %>%
  mutate(Title = case_when(
    Parameter == "logD_AlbuminWater" ~ "logD (Albumin:Water)",
    Parameter == "logD_GlobulinWater" ~ "logD (Globulin:Water)",
    Parameter == "logD_StructuralProteinWater" ~ "logD (Structural Protein:Water)",
    Parameter == "logD_StorageLipidWater" ~ "logD (Storage Lipid:Water)",
    Parameter == "logD_MembraneLipidWater" ~ "logD (Membrane Lipid:Water)"
    ))

GlimpseDistribution <- Distribution |> ggplot(aes(x = log(MW), y = Value, 
                                                 colour = Category_OECD, 
                                                 shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = (Value - SD), 
                    ymax = (Value + SD)),
                width = .0002) + 
  labs(y = "ln(D)", 
       x = "ln(MW)") +
  facet_wrap(~Title, scales = "free") +
  theme_quick
GlimpseDistribution

GlimpseDistribution <- GlimpseDistribution + CP_theme
ggsave(plot = GlimpseDistribution,
       file = here(OUTPUT, "Distribution.png"),
       dpi = 600,
       width = 18,
       height = 14,
       units = "cm")

# We can see that in general where we have values, PFSAs and PFCAs have similar values


## Albumin Water ----

CleanAlbumin <- Distribution %>%
  filter(Parameter == "logD_AlbuminWater") 

GlimpseAlbumin <- CleanAlbumin |> ggplot(aes(x = log(MW), # to log as values are also in log scale
                                             y = Value, 
                                             colour = Category_OECD, 
                                             shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = Value - SD, 
                ymax = Value + SD)) +
  labs(y = paste0(unique(CleanAlbumin$Title)), 
       x = "MW ln(g/mol)") +
  theme_quick
  
GlimpseAlbumin <- GlimpseAlbumin + CP_theme
ggsave(plot = GlimpseAlbumin,
       file = here(OUTPUT, "Albumin.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

# QSPR
glimpse(CleanAlbumin)

QSPR_CleanAlbumin <- CleanAlbumin |> 
  mutate(MW = log(MW)) # to log as values are also in log scale

y = QSPR_CleanAlbumin$Value
x = QSPR_CleanAlbumin$MW
data <- tibble(x = x, y = y)

fit_Albumin <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
                data = data,
                start = list(Const = min(QSPR_CleanAlbumin$Value),
                             Asym = 12.1127,
                             xmid = 5.5517  ,
                             scal = 0.3041 ),
                # control = nls.control(maxiter = 1000000),
                algorithm='port')

Results_Albumin <- summary(fit_Albumin)
print(Results_Albumin)

Formula <- Results_Albumin$formula[3]
RSE <- signif(Results_Albumin$sigma, 3)
Const <- signif(Results_Albumin$coefficients[1], 3)
Asym <- signif(Results_Albumin$coefficients[2], 3)
xmid <- signif(Results_Albumin$coefficients[3], 3)
scal <- signif(Results_Albumin$coefficients[4], 3)

Residuals_Albumin <- tibble(
  x = x,
  y = residuals(fit_Albumin)
)

Residuals_Albumin_plot <- ggplot(Residuals_Albumin) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_Albumin_plot

predicted <- tibble(
  MW = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001),
  predicted_y = predict(fit_Albumin, list(x = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

predicted_Albumin <- tibble(
  MW = PFAS_MW_list_log,
  predicted_value = predict(fit_Albumin, list(x = PFAS_MW_list_log))
) %>% arrange(MW)

observed_Albumin <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_CleanAlbumin2 <- QSPR_CleanAlbumin |> 
  mutate(nValue = exp(Value), 
         nSD = exp(SD))

QSPR_Albumin_plot <- QSPR_CleanAlbumin |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = Value - SD,
                    ymax = Value + SD),
                width = .002) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", const: ", Const,
                        ", asym: ", Asym,
                        ", xmid: ", xmid,
                        ", scal: ", scal),
       y = "ln(D_albumin/water)", 
       x = "MW ln(g/mol)") +
  theme_quick
QSPR_Albumin_plot

D <- QSPR_Albumin_plot + CP_theme
R <- Residuals_Albumin_plot + CP_theme

QSPR_Albumin_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_Albumin_plots,
       file = here(OUTPUT, "QSPR_Albumin_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")


# collect results
convert_Albumin <- predicted_Albumin |> 
  rename(MWlog = MW,
         meanlog = predicted_value) |> 
  mutate(MW = signif(exp(MWlog),3),
         PBK_mechanism = "Distribution",
         Parameter = unique(CleanAlbumin$Parameter), 
         Value = signif(meanlog,3),
         SD = RSE, 
         SD_comment = "is sdlog",
         Specific = "predicted_sigmoid",
         Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and sdlog is SD") |> 
  select(! c(MWlog, meanlog))

  
# add PFAS 
convert_Albumin <- full_join(convert_Albumin, MW_df)

# merge dataframes
CleanAlbumin <- CleanAlbumin |>  
  mutate(Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and SD of log Values is SD", 
         SD_comment = "NA") 

# Final dataframe
FinalAlbumin <- full_join(CleanAlbumin, convert_Albumin) |> 
  select(!Title)

## Globulin Water ----

CleanGlobulin <- Distribution %>%
  filter(Parameter == "logD_GlobulinWater") 

GlimpseGlobulin <- CleanGlobulin |> ggplot(aes(x = log(MW), # to log as values are also in log scale
                                             y = Value, 
                                             colour = Category_OECD, 
                                             shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = Value - SD, 
                    ymax = Value + SD)) +
  labs(y = paste0(unique(CleanGlobulin$Title)), 
       x = "MW ln(g/mol)") +
  theme_quick
GlimpseGlobulin
GlimpseGlobulin <- GlimpseGlobulin + CP_theme
ggsave(plot = GlimpseGlobulin,
       file = here(OUTPUT, "Globulin.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

# QSPR
glimpse(CleanGlobulin)

QSPR_CleanGlobulin <- CleanGlobulin |> 
  mutate(MW = log(MW)) # to log as values are also in log scale

y = QSPR_CleanGlobulin$Value
x = QSPR_CleanGlobulin$MW
data <- tibble(x = x, y = y)

fit_Globulin <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
                   data = data,
                   start = list(Const = min(QSPR_CleanGlobulin$Value),
                                Asym = 33.7336,
                                xmid = 6.8628,
                                scal = 0.5285),
                   # control = nls.control(maxiter = 1000000),
                   algorithm='port')

Results_Globulin <- summary(fit_Globulin)
print(Results_Globulin)

Formula <- Results_Globulin$formula[3]
RSE <- signif(Results_Globulin$sigma, 3)
Const <- signif(Results_Globulin$coefficients[1], 3)
Asym <- signif(Results_Globulin$coefficients[2], 3)
xmid <- signif(Results_Globulin$coefficients[3], 3)
scal <- signif(Results_Globulin$coefficients[4], 3)

Residuals_Globulin <- tibble(
  x = x,
  y = residuals(fit_Globulin)
)

Residuals_Globulin_plot <- ggplot(Residuals_Globulin) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_Globulin_plot

predicted <- tibble(
  MW = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001),
  predicted_y = predict(fit_Globulin, list(x = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

predicted_Globulin <- tibble(
  MW = PFAS_MW_list_log,
  predicted_value = predict(fit_Globulin, list(x = PFAS_MW_list_log))
) %>% arrange(MW)

observed_Globulin <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_CleanGlobulin2 <- QSPR_CleanGlobulin |> 
  mutate(nValue = exp(Value), 
         nSD = exp(SD))

QSPR_Globulin_plot <- QSPR_CleanGlobulin |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = Value - SD,
                    ymax = Value + SD),
                width = .002) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(title = "QSPR Globulin",
       caption = paste0("Formula: ", Formula, 
                        ", const: ", Const,
                        ", asym: ", Asym,
                        ", xmid: ", xmid,
                        ", scal: ", scal),
       y = "ln(D_globulin/water)", 
       x = "MW ln(g/mol)") +
  theme_quick
QSPR_Globulin_plot

D <- QSPR_Globulin_plot + CP_theme
R <- Residuals_Globulin_plot + CP_theme

QSPR_Globulin_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_Globulin_plots,
       file = here(OUTPUT, "QSPR_Globulin_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")

# collect results
convert_Globulin <- predicted_Globulin |> 
  rename(MWlog = MW,
         meanlog = predicted_value) |> 
  mutate(MW = signif(exp(MWlog),3),
         PBK_mechanism = "Distribution",
         Parameter = unique(CleanGlobulin$Parameter), 
         Value = signif(meanlog,3),
         SD = RSE, 
         SD_comment = "is sdlog",
         Specific = "predicted_sigmoid",
         Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and sdlog is SD") |> 
  select(! c(MWlog, meanlog))

# add PFAS 
convert_Globulin <- full_join(convert_Globulin, MW_df)

# merge dataframes
CleanGlobulin <- CleanGlobulin |> 
  mutate(Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and SD of log Values is SD", 
         SD_comment = "NA") |> 
  select(!Title)


# Final dataframe
FinalGlobulin <- full_join(CleanGlobulin, convert_Globulin) 


## StructuralProtein Water ----

CleanStructuralProtein <- Distribution %>%
  filter(Parameter == "logD_StructuralProteinWater") 

GlimpseStructuralProtein <- CleanStructuralProtein |> ggplot(aes(x = log(MW), # to log as values are also in log scale
                                             y = Value, 
                                             colour = Category_OECD, 
                                             shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = Value - SD, 
                    ymax = Value + SD)) +
  labs(y = paste0(unique(CleanStructuralProtein$Title)), 
       x = "MW ln(g/mol)") +
  theme_quick
GlimpseStructuralProtein
GlimpseStructuralProtein <- GlimpseStructuralProtein + CP_theme
ggsave(plot = GlimpseStructuralProtein,
       file = here(OUTPUT, "StructuralProtein.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

# QSPR
glimpse(CleanStructuralProtein)

QSPR_CleanStructuralProtein <- CleanStructuralProtein |> 
  mutate(MW = log(MW)) # to log as values are also in log scale

y = QSPR_CleanStructuralProtein$Value
x = QSPR_CleanStructuralProtein$MW
data <- tibble(x = x, y = y)

fit_StructuralProtein <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
                   data = data,
                   start = list(Const = min(QSPR_CleanStructuralProtein$Value),
                                Asym = 12.7597,
                                xmid = 6.2135,
                                scal = 0.2395),
                   # control = nls.control(maxiter = 1000000),
                   algorithm='port')

Results_StructuralProtein <- summary(fit_StructuralProtein)
print(Results_StructuralProtein)

Formula <- Results_StructuralProtein$formula[3]
RSE <- signif(Results_StructuralProtein$sigma, 3)
Const <- signif(Results_StructuralProtein$coefficients[1], 3)
Asym <- signif(Results_StructuralProtein$coefficients[2], 3)
xmid <- signif(Results_StructuralProtein$coefficients[3], 3)
scal <- signif(Results_StructuralProtein$coefficients[4], 3)

Residuals_StructuralProtein <- tibble(
  x = x,
  y = residuals(fit_StructuralProtein)
)

Residuals_StructuralProtein_plot <- ggplot(Residuals_StructuralProtein) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_StructuralProtein_plot

predicted <- tibble(
  MW = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001),
  predicted_y = predict(fit_StructuralProtein, list(x = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

predicted_StructuralProtein <- tibble(
  MW = PFAS_MW_list_log,
  predicted_value = predict(fit_StructuralProtein, list(x = PFAS_MW_list_log))
) %>% arrange(MW)

observed_StructuralProtein <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_CleanStructuralProtein2 <- QSPR_CleanStructuralProtein |> 
  mutate(nValue = exp(Value), 
         nSD = exp(SD))

QSPR_StructuralProtein_plot <- QSPR_CleanStructuralProtein |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = Value - SD,
                    ymax = Value + SD),
                width = .002) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", const: ", Const,
                        ", asym: ", Asym,
                        ", xmid: ", xmid,
                        ", scal: ", scal),
       y = "ln(D_structural protein/water)", 
       x = "MW ln(g/mol)") +
  theme_quick
QSPR_StructuralProtein_plot

D <- QSPR_StructuralProtein_plot + CP_theme
R <- Residuals_StructuralProtein_plot + CP_theme

QSPR_StructuralProtein_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_StructuralProtein_plots,
       file = here(OUTPUT, "QSPR_StructuralProtein_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")

# collect results
convert_StructuralProtein <- predicted_StructuralProtein |> 
  rename(MWlog = MW,
         meanlog = predicted_value) |> 
  mutate(MW = signif(exp(MWlog),3),
         PBK_mechanism = "Distribution",
         Parameter = unique(CleanStructuralProtein$Parameter), 
         Value = signif(meanlog,3),
         SD = RSE, 
         SD_comment = "is sdlog",
         Specific = "predicted_sigmoid",
         Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and sdlog is SD") |> 
  select(! c(MWlog, meanlog))

# add PFAS 
convert_StructuralProtein <- full_join(convert_StructuralProtein, MW_df)

# merge dataframes
CleanStructuralProtein <- CleanStructuralProtein |> 
  mutate(Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and SD of log Values is SD", 
         SD_comment = "NA") |> 
  select(!Title)


# Final dataframe
FinalStructuralProtein <- full_join(CleanStructuralProtein, convert_StructuralProtein) 


## StorageLipid Water ----

CleanStorageLipid <- Distribution %>%
  filter(Parameter == "logD_StorageLipidWater") 

GlimpseStorageLipid <- CleanStorageLipid |> ggplot(aes(x = log(MW), # to log as values are also in log scale
                                             y = Value, 
                                             colour = Category_OECD, 
                                             shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = Value - SD, 
                    ymax = Value + SD)) +
  labs(y = paste0(unique(CleanStorageLipid$Title)), 
       x = "MW ln(g/mol)") +
  theme_quick
GlimpseStorageLipid
GlimpseStorageLipid <- GlimpseStorageLipid + CP_theme
ggsave(plot = GlimpseStorageLipid,
       file = here(OUTPUT, "StorageLipid.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

# QSPR
glimpse(CleanStorageLipid)

QSPR_CleanStorageLipid <- CleanStorageLipid |> 
  mutate(MW = log(MW)) # to log as values are also in log scale

y = QSPR_CleanStorageLipid$Value
x = QSPR_CleanStorageLipid$MW
data <- tibble(x = x, y = y)

fit_StorageLipid <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
                   data = data,
                   start = list(Const = min(QSPR_CleanStorageLipid$Value),
                                Asym = -3.2391,
                                xmid = 6.1866,
                                scal = -0.1057 ),
                   # control = nls.control(maxiter = 1000000),
                   algorithm='port')

Results_StorageLipid <- summary(fit_StorageLipid)
print(Results_StorageLipid)

Formula <- Results_StorageLipid$formula[3]
RSE <- signif(Results_StorageLipid$sigma, 3)
Const <- signif(Results_StorageLipid$coefficients[1], 3)
Asym <- signif(Results_StorageLipid$coefficients[2], 3)
xmid <- signif(Results_StorageLipid$coefficients[3], 3)
scal <- signif(Results_StorageLipid$coefficients[4], 3)

Residuals_StorageLipid <- tibble(
  x = x,
  y = residuals(fit_StorageLipid)
)

Residuals_StorageLipid_plot <- ggplot(Residuals_StorageLipid) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_StorageLipid_plot

predicted <- tibble(
  MW = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001),
  predicted_y = predict(fit_StorageLipid, list(x = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

predicted_StorageLipid <- tibble(
  MW = PFAS_MW_list_log,
  predicted_value = predict(fit_StorageLipid, list(x = PFAS_MW_list_log))
) %>% arrange(MW)

observed_StorageLipid <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_CleanStorageLipid2 <- QSPR_CleanStorageLipid |> 
  mutate(nValue = exp(Value), 
         nSD = exp(SD))

QSPR_StorageLipid_plot <- QSPR_CleanStorageLipid |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = Value - SD,
                    ymax = Value + SD),
                width = .002) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", const: ", Const,
                        ", asym: ", Asym,
                        ", xmid: ", xmid,
                        ", scal: ", scal),
       y = "ln(D_storage lipid/water)", 
       x = "MW ln(g/mol)") +
  theme_quick
QSPR_StorageLipid_plot

D <- QSPR_StorageLipid_plot + CP_theme
R <- Residuals_StorageLipid_plot + CP_theme

QSPR_StorageLipid_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_StorageLipid_plots,
       file = here(OUTPUT, "QSPR_StorageLipid_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")

# collect results
convert_StorageLipid <- predicted_StorageLipid |> 
  rename(MWlog = MW,
         meanlog = predicted_value) |> 
  mutate(MW = signif(exp(MWlog),3),
         PBK_mechanism = "Distribution",
         Parameter = unique(CleanStorageLipid$Parameter), 
         Value = signif(meanlog,3),
         SD = RSE, 
         SD_comment = "is sdlog",
         Specific = "predicted_sigmoid",
         Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and sdlog is SD") |> 
  select(! c(MWlog, meanlog))

# add PFAS 
convert_StorageLipid <- full_join(convert_StorageLipid, MW_df)

# merge dataframes
CleanStorageLipid <- CleanStorageLipid |> 
  mutate(Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and SD of log Values is SD", 
         SD_comment = "NA") |> 
  select(!Title)

# Final dataframe
FinalStorageLipid <- full_join(CleanStorageLipid, convert_StorageLipid) 


## MembraneLipid Water ----

### PFCAs and PFECAs ---- 
CleanMembraneLipid <- Distribution %>%
  filter(Parameter == "logD_MembraneLipidWater", 
         Category_OECD != "PFSAs") 

GlimpseMembraneLipid <- CleanMembraneLipid |> ggplot(aes(x = log(MW), # to log as values are also in log scale
                                             y = Value, 
                                             colour = Category_OECD, 
                                             shape = Specific)) +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.5) +
  geom_point() +
  geom_errorbar(aes(ymin = Value - SD, 
                    ymax = Value + SD)) +
  labs(y = paste0(unique(CleanMembraneLipid$Title)), 
       x = "MW ln(g/mol)") +
  theme_quick
GlimpseMembraneLipid
GlimpseMembraneLipid <- GlimpseMembraneLipid + CP_theme
ggsave(plot = GlimpseMembraneLipid,
       file = here(OUTPUT, "MembraneLipid.png"),
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

# QSPR
glimpse(CleanMembraneLipid)

QSPR_CleanMembraneLipid <- CleanMembraneLipid |> 
  mutate(MW = log(MW)) # to log as values are also in log scale

y = QSPR_CleanMembraneLipid$Value
x = QSPR_CleanMembraneLipid$MW
data <- tibble(x = x, y = y)

fit <- nls(y ~ SSlogis(x, Asym, xmid, scal),
           data = data,
           control = nls.control(maxiter = 1000000, minFactor = 0.0000000000002))

fit_MembraneLipid <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
                   data = data,
                   start = list(Const = min(QSPR_CleanMembraneLipid$Value),
                                Asym = 29.4370 ,
                                xmid = 6.6275 ,
                                scal = 0.6137 ),
                   # control = nls.control(maxiter = 1000000),
                   algorithm='port')

Results_MembraneLipid <- summary(fit_MembraneLipid)
print(Results_MembraneLipid)

Formula <- Results_MembraneLipid$formula[3]
RSE <- signif(Results_MembraneLipid$sigma, 3)
Const <- signif(Results_MembraneLipid$coefficients[1], 3)
Asym <- signif(Results_MembraneLipid$coefficients[2], 3)
xmid <- signif(Results_MembraneLipid$coefficients[3], 3)
scal <- signif(Results_MembraneLipid$coefficients[4], 3)

Residuals_MembraneLipid <- tibble(
  x = x,
  y = residuals(fit_MembraneLipid)
)

Residuals_MembraneLipid_plot <- ggplot(Residuals_MembraneLipid) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", RSE)) +
  theme_quick
Residuals_MembraneLipid_plot

predicted <- tibble(
  MW = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001),
  predicted_y = predict(fit_MembraneLipid, list(x = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001))),
  min = predicted_y - RSE, 
  max = predicted_y + RSE
) %>% arrange(MW)

predicted_MembraneLipid <- tibble(
  MW = PFAS_MW_list_log,
  predicted_value = predict(fit_MembraneLipid, list(x = PFAS_MW_list_log))
) %>% arrange(MW)

observed_MembraneLipid <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_MembraneLipid_plot <- QSPR_CleanMembraneLipid |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.5) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = Value - SD,
                    ymax = Value + SD),
                width = .002) +
  geom_path(data = predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", Formula, 
                        ", const: ", Const,
                        ", asym: ", Asym,
                        ", xmid: ", xmid,
                        ", scal: ", scal),
       y = "ln(D_membrane lipid/water)", 
       x = "MW ln(g/mol)") +
  theme_quick
QSPR_MembraneLipid_plot

D <- QSPR_MembraneLipid_plot + CP_theme
R <- Residuals_MembraneLipid_plot + CP_theme

QSPR_MembraneLipid_plots <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_MembraneLipid_plots,
       file = here(OUTPUT, "QSPR_MembraneLipid_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")

# collect results
convert_MembraneLipid <- predicted_MembraneLipid |> 
  rename(MWlog = MW,
         meanlog = predicted_value) |> 
  mutate(MW = signif(exp(MWlog),3),
         PBK_mechanism = "Distribution",
         Parameter = unique(CleanMembraneLipid$Parameter), 
         Value = signif(meanlog,3),
         SD = RSE, 
         SD_comment = "is sdlog",
         Specific = "predicted_sigmoid",
         Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and sdlog is SD") |> 
  select(! c(MWlog, meanlog))

# add PFAS 
convert_MembraneLipid <- full_join(convert_MembraneLipid, MW_df)


CleanMembraneLipid <- Distribution %>%
  filter(Parameter == "logD_MembraneLipidWater", 
         Category_OECD != "PFSAs")

# merge dataframes
CleanMembraneLipid <- CleanMembraneLipid |> 
  mutate(Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and SD of log Values is SD", 
         SD_comment = "NA") |> 
  select(!Title)

# Final dataframe
InitialMembraneLipid <- full_join(CleanMembraneLipid, convert_MembraneLipid) |>
  mutate(method = "PFCAs")
CleanInitialMembraneLipid <- InitialMembraneLipid |> filter(Category_OECD != "PFSAs")

### PFSAs ----

CleanMembraneLipid <- Distribution %>%
  filter(Parameter == "logD_MembraneLipidWater", 
         Category_OECD == "PFSAs")
QSPR_CleanMembraneLipid <- CleanMembraneLipid |> mutate(MW = log(MW))

y = QSPR_CleanMembraneLipid$Value 
x = QSPR_CleanMembraneLipid$MW
data <- tibble(x = x, y = y)

PFSAsfit_MembraneLipid <- nls(y ~ 3.95 + 7.46/(1+exp((xmid-x)/0.135)),
                              data = data,
                              start = list(xmid = 5.98),
                              control = nls.control(maxiter = 1000000),
                              algorithm='port')

Results_PFSAsMembraneLipid <- summary(PFSAsfit_MembraneLipid)
print(Results_PFSAsMembraneLipid)

PFSA_Formula <- Results_PFSAsMembraneLipid$formula[3]
PFSA_RSE <- signif(Results_PFSAsMembraneLipid$sigma, 3)
PFSA_xmid <- signif(Results_PFSAsMembraneLipid$coefficients[1], 3)

ResidualsPFSA__MembraneLipid <- tibble(
  x = x,
  y = residuals(PFSAsfit_MembraneLipid)
)

ResidualsPFSA_MembraneLipid_plot <- ggplot(ResidualsPFSA__MembraneLipid) +
  geom_point(aes(x, y)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks=seq(1,13))+
  labs(y = "Residuals", x = "MW ln(g/mol)",
       title =  paste0("RSE: ", PFSA_RSE)) +
  theme_quick
ResidualsPFSA_MembraneLipid_plot

PFSAs_predicted <- tibble(
  MW = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001),
  predicted_y = predict(PFSAsfit_MembraneLipid, list(x = seq(min(PFAS_MW_list_log), max(PFAS_MW_list_log), by = 0.000001))),
  min = predicted_y - PFSA_RSE, 
  max = predicted_y + PFSA_RSE
) %>% arrange(MW)

PFSAs_predicted_MembraneLipid <- tibble(
  MW = PFAS_MW_list_log,
  predicted_value = predict(PFSAsfit_MembraneLipid, list(x = PFAS_MW_list_log))
) %>% arrange(MW)

observed_MembraneLipid <- data %>% 
  mutate(MW = x, Value = y) |> 
  arrange(MW)

QSPR_PFSAs_MembraneLipid_plot <- QSPR_CleanMembraneLipid |> ggplot() +
  scale_color_viridis(option = "H", discrete = TRUE, end = 0.8) +
  geom_point(aes(MW, Value, 
                 colour = Category_OECD)) +
  geom_errorbar(aes(MW, Value, 
                    colour = Category_OECD,
                    ymin = Value - SD,
                    ymax = Value + SD),
                width = .002) +
  geom_path(data = PFSAs_predicted, aes(MW, predicted_y)) +
  geom_ribbon(data = PFSAs_predicted, aes(x = MW, ymin = min, ymax = max), alpha = .15) +
  labs(caption = paste0("Formula: ", PFSA_Formula, 
                        ", xmid: ", PFSA_xmid),
       y = "ln(D_membrane lipid/water)", 
       x = "MW ln(g/mol)") +
  theme_quick
QSPR_PFSAs_MembraneLipid_plot

D <- QSPR_PFSAs_MembraneLipid_plot + CP_theme
R <- ResidualsPFSA_MembraneLipid_plot + CP_theme

QSPR_PFSAs_MembraneLipid_plot <- D / R + 
  plot_layout(heights = c(5,1))

ggsave(plot = QSPR_PFSAs_MembraneLipid_plot,
       file = here(OUTPUT, "QSPR_PFSAs_MembraneLipid_plots.png"),
       dpi = 600,
       width = 12,
       height = 12,
       units = "cm")

# collect results
PFSAs_convert_MembraneLipid <- PFSAs_predicted_MembraneLipid |> 
  rename(MWlog = MW,
         meanlog = predicted_value) |> 
  mutate(MW = signif(exp(MWlog),3),
         PBK_mechanism = "Distribution",
         Parameter = unique(CleanMembraneLipid$Parameter), 
         Value = signif(meanlog,3),
         SD = RSE, 
         SD_comment = "is sdlog",
         Specific = "predicted_sigmoid",
         Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and sdlog is SD") |> 
  select(! c(MWlog, meanlog))

# add PFAS 
PFSAs_convert_MembraneLipid <- full_join(PFSAs_convert_MembraneLipid, MW_df)

# merge dataframes
CleanMembraneLipid <- Distribution %>%
filter(Parameter == "logD_MembraneLipidWater", 
       Category_OECD == "PFSAs")

CleanMembraneLipid <- CleanMembraneLipid |> 
  mutate(Distribution = "lognormal", 
         Distribution_comment = "meanlog is ln(Value) and SD of log Values is SD", 
         SD_comment = "NA") |> 
  select(!Title)

PFSAs_FinalMembraneLipid <- full_join(CleanMembraneLipid, PFSAs_convert_MembraneLipid) |> 
  mutate(method = "PFSAs")
CleanPFSAs_FinalMembraneLipid <- PFSAs_FinalMembraneLipid |> 
  filter(Category_OECD == "PFSAs")

FinalMembraneLipid <- rbind(CleanInitialMembraneLipid, 
                            CleanPFSAs_FinalMembraneLipid) |> 
  select(!method)

# Collect all dfs ----

FinalParameters <- rbind(FinalPapp |>  mutate(Unit = "ln(cm/day)"),
                         FinalAlbumin |> mutate(Unit = "ln()"),
                         FinalGlobulin |> mutate(Unit = "ln()"),
                         FinalStructuralProtein |> mutate(Unit = "ln()"),
                         FinalStorageLipid |> mutate(Unit = "ln()"),
                         FinalMembraneLipid |> mutate(Unit = "ln()"),
                         FinalBSEP |> mutate(Unit = "L/day/mg protein"),
                         FinalNTCP |> mutate(Unit = "L/day/mg protein"),
                         FinalOATP2B1 |> mutate(Unit = "L/day/mg protein"),
                         FinalOATP1B1 |> mutate(Unit = "L/day/mg protein"),
                         FinalOATP1B3 |> mutate(Unit = "L/day/mg protein"),
                         FinalOAT1 |> mutate(Unit = "L/day/mg protein"),
                         FinalOAT3 |> mutate(Unit = "L/day/mg protein"),
                         FinalOAT4 |> mutate(Unit = "L/day/mg protein"),
                         FinalSF |> mutate(Unit = "mg protein/kg tissue")) |> 
  mutate(
  Predicted = ifelse(str_detect(Specific, "predicted_"), 
                     "yes", "no"))

write.xlsx(FinalParameters, here(OUTPUT,"FinalParameters.xlsx"))
FinalParameters <- read_xlsx(here(OUTPUT, "FinalParameters.xlsx"))

glimpseParameters <- FinalParameters |> 
  # filter(!str_detect(Parameter, 'SF_')) |> 
  ggplot(aes(PFAS, Value, colour = Predicted,
             shape = Predicted)) +
  scale_color_viridis(option = "B", discrete = TRUE, 
                      begin = 0.2, end = 0.8, alpha = 0.7,
                      direction = -1) +
  # geom_pointrange(aes(ymin = Value-SD,
  #                     ymax = Value+SD),
  #                 position=position_jitter(width=0.3),
  #                 # linetype='dotted',
  #                 size = 0.3,
  #                 linewidth = 0.5) +
  geom_jitter() +
  # scale_color_manual(values = c("yes" = "purple",
  #                               "no" = "grey50")) +
  facet_wrap(~ Parameter, scales = "free") +
  coord_flip() +
  labs(title = "PFAS-specific parameters",
       x = "", 
       y = "") +
  theme_quick

glimpseParameters <- glimpseParameters + CP_theme + 
  theme(legend.title = element_text(size = 30),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 25, face = "bold", angle = 15)
  )

ggsave(plot = glimpseParameters,
       file = here(OUTPUT, "AllParameters.png"),
       dpi = 600,
       width = 20,
       height = 24,
       units = "cm")
