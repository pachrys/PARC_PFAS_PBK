# --------------------------------------------------------------------------- #
# Script for running the MCMC using eFAST
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
library(sensitivity)
library(ggrepel)
library(scales)


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
    
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 30, margin = margin(0, 0, 0, 0)),
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
    
    strip.text.x = element_text(size = 30, face = "bold", margin = margin(0,0,0,0))
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

Patchwork_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 50, face = "bold", margin = margin(0, 0, b = 0.5, 0)),
    plot.subtitle = element_text(hjust = 0, size = 35, margin = margin(0, 0, b = 0.5, 0)),
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


# Set output storage directory
OUTPUT <- here("Output", "test_eFAST") 
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

Dose_df <- tibble(
  PFAS = PFAS_list,
  Dose = 0.1) 

nPFAS <- length(PFAS_list)


## Tissue constituents for partition coefficient calculation
Physio.data <- read_csv(here("Input", "TissueComposition.csv")) %>% 
  filter(Species == "Human") %>%
  filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
  select(- c(Comment, Species)) %>% 
  mutate(across(-Tissue, as.numeric))


## All PBK parameters
SA_Physio_params <- read_xlsx(here("Input", "Physio_Parameters.xlsx")) 
SA_TK_params <- read_xlsx(here("Input", "TK_Parameters.xlsx")) |>
  filter(Predicted == "yes") 

# load functions
source(here("Script", "functions_NG_PBK_model.R"))


# Function to run the eFAST design
run_eFAST <- function(parm.GSA, parm.c){
  
  
  # Final parameters per test
  names(parm.GSA) <- eFAST_parms
  parm.final <- parm.c
  parm.final[names(parm.GSA)] <- parm.GSA
  
  # print(parm.final[names(parm.GSA)])

  Dose_ug_day = 0.1
  
  A_init <- c(AIL = Dose_ug_day,
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
              InMB = Dose_ug_day)
  
  
  pars <- GSA_Calc_final_params(parms = parm.final)
  
  # save parameters to check when analysing results
  run_pars <- as.data.frame(pars)
  
  PBK.out <- lsoda(y = A_init,
                   times = seq(0, 365*3, by = 1), 
                   func = PBK.model,
                   parms = pars,
                   rtol = 1e-4,        
                   atol = 1e-6,
                   hmax = 0.5,         
                   maxsteps = 100000)
  
  PBK_OUTPUT <- as.data.frame(PBK.out) 
  PBK_OUTPUT[,-1] <- pmax(PBK_OUTPUT[,-1], 0)
  
  data <- PBK_OUTPUT |> drop_na() 
  
  # AUC ----
  AUC <- trapz(data$time, data$CP) # ng * day/ml
  
  # Cmax, Tmax, Tlast ----
  Cmax <- max(data$CP) # ng/ml
  Tmax <- filter(data, CP == max(data$CP)) %>% pull(time) # days
  Tlast <- filter(data, CP > 1e-9) |> select(time) |> last() |> pull() # 1e-9 ng/ml is already a very small number
  
  # calculate half life ----
  elimination_phase <- data |>
    filter(time >= Tmax & time <= Tlast) |>
    mutate(logCP = log(CP))
  
  # linear regression to find the slope
  elimination_model <- lm(logCP ~ time, data = elimination_phase)
  
  # calculate elimination constant (i.e. slope of the line)
  kel <- -elimination_model$coefficients[[2]]
  
  # calculate the half-life
  half_life <- (log(2) / kel) / 365 # unit correction to years
  
  return( c(Cmax=Cmax, 
            half_life = half_life, 
            AUC = AUC))
}

for (i in 1:nPFAS) {
  
  
  
  PFAS_i <- PFAS_list[i]
  message(glue("{PFAS_i}------------------------------------"))
  
  # physiological parameters; 
  Physio_params <- tibble(
    Parameter = SA_Physio_params$Abbreviation,
    Value = SA_Physio_params$Value,
    binf = SA_Physio_params$Binf,
    bsup = SA_Physio_params$Bsup
  )
  
  Dose_ug_day <- filter(Dose_df, PFAS == PFAS_i) |>
    select(Dose) |> pull()
  
  PhysChemTK_params <- SA_TK_params |>
    filter(PFAS == PFAS_i) |>
    mutate(Value = ifelse(PBK_mechanism == "Distribution", exp(Value),
                          ifelse(PBK_mechanism == "Absorption", exp(Value),
                                 Value)),
    binf = min,
    bsup = max
    ) |>
    select(Parameter, Value, binf, bsup)
  
  Model_params <- rbind(Physio_params, PhysChemTK_params)
  
  eFAST_parms <- c(
    "logD_AlbuminWater",
    "SF_OAT4",
    "BW",
    "logD_MembraneLipidWater",
    "SF_OATP2B1_I",
    "SF_NTCP",
    "SF_BSEP",
    "logD_GlobulinWater",
    "SF_OAT3",
    "Calb_P",
    "logD_StructuralProteinWater",
    "GFR",
    "VAc",
    "Qmenses",
    "SF_OATP1B1",
    "Clear_NTCPc",
    "Papp_SI"
    )
  
  params_qlnorm <- c(
    "logD_AlbuminWater",
    "logD_MembraneLipidWater",
    "logD_GlobulinWater",
    "logD_StructuralProteinWater",
    "Clear_NTCPc",
    "Papp_SI"
  )
  
  eFAST_TK_params <- SA_TK_params |> 
    filter(PFAS == PFAS_i,
           Parameter %in% eFAST_parms) |> 
    select(Parameter, Value, RSE, Distribution) |> 
    mutate(Value = ifelse(Parameter == "Clear_NTCPc", log(Value), Value)) |> 
    transmute(Parameter = Parameter, 
              Value = Value,
              meanlog = ifelse(Parameter %in% params_qlnorm, Value, NA),
              sdlog = ifelse(Parameter %in% params_qlnorm, RSE, NA),
              min = ifelse(Parameter %in% params_qlnorm, NA, 1),
              max = ifelse(Parameter %in% params_qlnorm, NA, Value),
              Distribution = ifelse(Distribution == "lognormal", "qlnorm", "qunif"))  
     
  eFAST_phys_params <- Physio_params |> 
    filter(Parameter %in% eFAST_parms) |> 
    transmute(Parameter = Parameter, 
              Value = Value, 
              meanlog = ifelse(Parameter %in% params_qlnorm, Value, NA),
              sdlog = ifelse(Parameter %in% params_qlnorm, RSE, NA),
              min = ifelse(Parameter %in% params_qlnorm, NA, binf),
              max = ifelse(Parameter %in% params_qlnorm, NA, bsup),
              Distribution = "qunif")  
  eFAST_Ps <- full_join(eFAST_TK_params, eFAST_phys_params)
  eFAST_Ps <- full_join(eFAST_Ps, eFAST_TK_params)
  
  P <- eFAST_parms
  parm.c <- as.list(setNames(Model_params$Value, Model_params$Parameter))
  
  ## eFAST experiment design ----

  eFAST_parms <- eFAST_Ps$Parameter
  eFAST_factors <- length(eFAST_parms)
  
  set.seed(10010) # set seeds for reproducibility
  
  # Define q: q needs to be a list of character strings, giving the names of the quantile functions
  q <- eFAST_Ps$Distribution
  
  # Define q.arg: q.arg needs to be a list of lists
  
  q.arg <- lapply(1:eFAST_factors, function(i) { # apply the function to eachone of the rows of P
    if (eFAST_Ps$Distribution[i] == "qunif") { #if the distribution in uniform
      list(min = eFAST_Ps$min[i], max = eFAST_Ps$max[i]) #create a list containing min (Bing of row i) and max (Bsup of row i) per row
    } else {
      list(meanlog = eFAST_Ps$meanlog[i], sdlog = eFAST_Ps$sdlog[i]) #if not then create a list containing mean and sdlog
    }
  })
  
  names(q.arg) <- eFAST_Ps$Parameter
  length(q) == length(q.arg)
  
  set.seed(10010) # set seeds for reproducibility
  
  # eFAST test function
  Experience <- fast99(
    model = NULL, 
    factors = eFAST_parms, # These are the parameters that will be varying in the model
    n = 1000,# Integer giving the sample size, i.e. the length of the discretization of the s-space
    q = q,
    q.arg = q.arg
  )
  
  # create the final eFAST design and remove infinite values if any
  
  dim(Experience$X) # dataframe of: n factors(number of parameters) * n n(1000) number of observations, of n factors(number of parameters) number of variables
  eFAST_design <- Experience$X |>
    as.data.frame() |>
    mutate(across(everything(), ~na_if(., Inf))) |>
    drop_na()
  colnames(eFAST_design) <- eFAST_parms

  # visualise and save design
  Design_plot <- eFAST_design |>
    pivot_longer(everything(),
                 names_to = "Parameter",
                 values_to = "Value") |>
    ggplot(aes(x=Value))+
    geom_histogram()+
    facet_wrap(~Parameter, scales = "free_y")+
    coord_flip()+
    labs(title = paste0(PFAS_i," eFAST test parameter distributions"),
         y = "Frequency",
         x = "Value (different unit per parameter)") +
    theme_quick
  Design_plot
  
  Design_plot <- Design_plot + Patchwork_theme +
    theme(axis.text.x = element_text(angle = 0))
  
  ggsave(plot = Design_plot,
         file = here(OUTPUT, paste0(PFAS_i,"_eFAST_ExpDesign.png")),
         dpi = 600,
         width = 34,
         height = 24,
         units = "cm")
  
  write.csv(eFAST_design,
            file = here(OUTPUT, paste0(PFAS_i,"_eFAST_ExpDesign.csv")),
            row.names = FALSE)
  
  results <- apply(eFAST_design, 1, function(row) {
    run_eFAST(row, parm.c = parm.c)  
  })
  
  y <- results
  y <- as.data.frame(y) 
  savey = y
  save(savey, file = here(OUTPUT, paste0(PFAS_i,"_eFAST_y.RData")))
  
  # load(here(OUTPUT, paste0(PFAS_i,"_eFAST_y.RData")))
  # y <- savey
  # y <- t(y) %>% na.omit() %>% t()
  sim.results.eFAST <- as.matrix(y)
  
  # Tell eFAST y for each model output
  for (i in (1:length(sim.results.eFAST[,1])))
  {
    tell(Experience, sim.results.eFAST[i,])
  }
  
  plot(Experience)
  print(PFAS_i)
  print(Experience)
  
  Variance <- Experience$V # Total variance
  Done <- Experience$D1
  Dt <- Experience$Dt
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  total_order <- 1 - (Dt / Variance)  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance

  names(Variance) <- names(Done) <- names(Dt) <- names(first_order) <- names(total_order) <- eFAST_parms

  lowry_data <- data.frame(
    Parameter = eFAST_parms,
    Main.Effect = first_order,
    Interaction = total_order)
  
  write.csv(lowry_data, here(OUTPUT, paste0(PFAS_i, "_eFAST_lowry_data.csv")), row.names = FALSE)
  
  ordered_data <- lowry_data %>%
    arrange(desc(Main.Effect)) %>% 
    mutate(
      Parameter = factor(Parameter, levels = Parameter),
      Total.effect = Main.Effect + Interaction,
      Cumulative.main = cumsum(Main.Effect),
      Cumulative.total = cumsum(Total.effect)) %>% 
    mutate(
      Main.Effect = Main.Effect/max(Cumulative.total),
      Interaction = Interaction/max(Cumulative.total),
      Total.effect = Main.Effect + Interaction,
      Cumulative.total = cumsum(Total.effect),
      ymin = cumsum(Main.Effect), 
      ymax = pmin(Cumulative.total, 1)
    )
  
  long_data <- ordered_data %>% 
    pivot_longer(c(Main.Effect, Interaction))
  long_data$name <- factor(long_data$name, levels = c("Main.Effect", "Interaction"))
  
  lowry_plot <- 
    ggplot(long_data) +
    geom_col(aes(x = Parameter, y = value, fill = name),
             position = position_stack(reverse = TRUE),
             width = 0.8,
             alpha = 0.8) +
    geom_ribbon(
      data = ordered_data,
      aes(x = as.numeric(Parameter),
          ymin = ymin, 
          ymax = ymax),
      fill = "grey50", alpha = 0.3, color = "grey50", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Main.Effect" = "grey5", "Interaction" = "grey50"),
      labels = c("Main Effect", "Interaction")) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      name = "Sensitivity Index") +
    labs(
      x = "Parameter",
      title = paste0(PFAS_i, " eFAST Sensitivity Analysis")) +
    CP_theme +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(size = 30, angle = 45, hjust = 1,),
      axis.line = element_line(size = 0.1)
    )
  
  # lowry_plot
  ggsave(
    file = here(OUTPUT, paste0(PFAS_i, "_eFAST_plot.png")),
    plot = lowry_plot,
    dpi = 600, 
    width = 16, 
    height = 9, 
    units = "cm")
  
}

