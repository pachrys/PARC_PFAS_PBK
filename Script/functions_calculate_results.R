# Post run results functions

calc_TK_results <- function(PBK_OUTPUT, OD) {
  
  data <- PBK_OUTPUT |> drop_na()
  
  # AUC ----
  AUC <- trapz(data$time, data$CP) # ng * day/ml
  
  # Cmax, Tmax, Tlast ----
  Cmax <- max(data$CP) # ng/ml
  Tmax <- filter(data, CP == max(data$CP)) %>% pull(time) # days
  Tlast <- filter(data, CP > 0) |> select(time) |> last() |> pull()
  
  # calculate half life ----
  # check that there is enough data points to run linear regression
  # if (is_empty(Tlast)){
  #   
  #   print("Tlast empty")
  #   elimination_phase <- NA
  #   kel <- NA
  #   half_life <- NA
  #   
  # } else if(Tlast <= Tmax){
  #   
  #   print("not enough data-points to calculate half life")
  #   elimination_phase <- NA
  #   kel <- NA
  #   half_life <- NA
  #   
  # } else {
    
    # elimination phase, i.e. after the Tmax
    
    elimination_phase <- data |>
      filter(time >= Tmax & time <= Tlast) |>
      mutate(logCP = log(CP))
    
    # linear regression to find the slope
    elimination_model <- lm(logCP ~ time, data = elimination_phase)
    
    # calculate elimination constant (i.e. slope of the line)
    kel <- -elimination_model$coefficients[[2]]
    
    # calculate the half-life
    half_life <- (log(2) / kel) / 365 # unit correction to years
    
    
  # }
  
 
  
  # fractions absorbed and eliminated ----

  # all amounts in tissues at Tmax
  Amounts <- data |> select(time, starts_with("A"), -AIn) # AIn is the input amount (used to calculate MB)
  
  ## Fa ----
  A_Tmax <- Amounts %>% filter(time == Tmax) |>
    select(starts_with("A"), -c(AIL, AFe))
  ATotal_Tmax <- rowSums(A_Tmax) # total amount absorbed at Tmax
  Fa <- ATotal_Tmax/OD *100 # includes urine
  
  ## Ffeces and Furine ----
  # Cumulative amounts 
  A_feces <- Amounts |> select(time, AFe) 
  A_urine <- Amounts |> select(time, AUr) 
  
  # Fractions excreted
  feEl <- last(A_feces$AFe)
  Ffeces <- feEl/OD * 100
  
  urEl <- last(A_urine$AUr) 
  Furine <- urEl/OD *100
  
  FtotalEl = Ffeces + Furine
  RfeEl = Ffeces/FtotalEl
  RurEl = Furine/FtotalEl
  
  # split Ffe as total Ffe = Ffe pre-systemic + Ffe left in intestine during absorption phase + Ffe bile 

  # Cumulative amount in feces before it reached the systemic circulation (before Tmax so during absorption phase)
  feEl_preTmax <- filter(A_feces, time <= Tmax) |> 
    select(!time) |> last() |> pull()
  Ffeces_presystemic <- feEl_preTmax/OD *100   # Fraction eliminated in feces before reaching he systemic circ
  FinIL_presystemic <- 100 - Fa - Ffeces_presystemic # amount remaining in the intestinal lumen at Tmax
  
  feEl_duringelimination <- feEl - feEl_preTmax # the total amount in feces minus what was removed before reaching the systemic
  Ffeces_duringelimination <- feEl_duringelimination/OD *100
  FbilEl <- Ffeces_duringelimination - FinIL_presystemic # is the total amount found in feces minous what was already removed before reaching the systemic circulation minus what was remaining in the intestinal lumen at that moment
  
  ## fecal and urinary clearances ----
  times_days = seq(round(Tmax,0)+1, round(Tlast,0)-1, by = 1) # it should be after Tmax as I'm using the pre-systemic first

  # daily fecal elimination during elimination phase
  dailyFe_elimination <- A_feces |>
    filter(time %in% times_days) |>
    mutate(daily = AFe - lag(AFe, n = 1, default = feEl_preTmax)) |>
    filter(daily > 0) |>
    select(daily) |> pull() # ug

  urEl_preTmax <- filter(A_urine, time <= Tmax) |>
    select(!time) |> last() |> pull()

  # daily urinary elimination during elimination phase
  dailyUr_elimination <- A_urine |>
    filter(time %in% times_days) |>
    mutate(daily = AUr - lag(AUr, n = 1, default = urEl_preTmax)) |>
    filter(daily > 0) |>
    select(daily) |> pull() # ug

  # daily plasma concentration
  dailyCP <- data |> select(time, CP) |>
    filter(time %in% times_days) |>
    select(CP) |> pull()

  CLfe = dailyFe_elimination/dailyCP # ug/day / ug/L => L/day
  median_CLfe <- median(CLfe) # ug/day / ug/L => L/day
  CLur = dailyUr_elimination/dailyCP # ug/day / ug/L => L/day
  median_CLur <- median(CLur) # ug/day / ug/L => L/day

  return( c(AUC = AUC,
           Cmax = Cmax,
           # Tmax = Tmax/365,
           # Tlast = Tlast/365,
           # kel = kel,
           HL = half_life,
           
           Fa = Fa/100,
           FelTotal = FtotalEl/100,
           
           Ffeces = Ffeces/100,
           FbilEl = FbilEl/100, 
           # Ffeces_presystemic = Ffeces_presystemic,
           # Ffeces_duringelimination = Ffeces_duringelimination,
           
           Furine = Furine/100,
           
           RfeEl = RfeEl,
           RurEl = RurEl,
           
           CLfe = median_CLfe,
           CLur = median_CLur
           )
  )
}

# quick function to find errors...
check_error <- function(PBK_OUTPUT, OD) {
  
  data <- PBK_OUTPUT |> drop_na()
  Tmax <- filter(data, CP == max(data$CP)) %>% pull(time) # days
  Tlast <- filter(data, CP > 0) |> select(time) |> last() |> pull()
  Iteration <- unique(data$Iteration)
  
  return(c(Tmax = Tmax, 
           Tlast = Tlast, 
           Iteration = Iteration))}



calc_TK_summary_results <- function(raw_TK_results){
  
  Results_all <- raw_TK_results |>  #do.call(rbind, raw_TK_results) |>
    na.omit() |> t() |> 
    as.data.frame() |> 
    summarise(across(where(is.numeric), .fns =
                       list(median = ~median(.,na.rm = TRUE),
                            mean = ~mean(.,na.rm = TRUE),
                            q25 = ~quantile(., 0.25, na.rm = TRUE),
                            q75 = ~quantile(.,0.75, na.rm = TRUE)))) |>
    pivot_longer(everything(), names_to = "Variable",
                 values_to = "Value")  |>
    separate(Variable,
             into = c("Parameter", "Variable"))|>
    pivot_wider(names_from = Variable, values_from = Value)
  
  return(Results_all)
}

calc_C_summary <- function(Results_df, Tissue_j){
  
  C_df <- Results_df |> filter(Tissue == Tissue_j)
  C_wide_df <- C_df |> select(!PFAS) |>
    mutate(Iteration = paste0("n_", Iteration)) |>
    pivot_wider(names_from = Iteration, values_from = Concentration)
  C_summary <- C_wide_df |>
    select(!c(time, Tissue)) |> t() |>
    as.data.frame() |>
    summarise(across(where(is.numeric), .fns =
                       list(median = ~median(.,na.rm = TRUE),
                            mean = ~mean(.,na.rm = TRUE),
                            q25 = ~quantile(., 0.25, na.rm = TRUE),
                            q75 = ~quantile(.,0.75, na.rm = TRUE))))  |>
    pivot_longer(everything(), names_to = "Variable",
                 values_to = "Value") |>
    separate(Variable,
             into = c("Time", "Variable")) |>
    mutate(Time = as.numeric(str_replace(Time, "V", ""))) |>
    # mutate(Time = (Time*10)/365) |>  #unit correction, to years
    pivot_wider(names_from = Variable, values_from = Value)
  
  return(C_summary)
}

plot_per_tissue <- function(C_sum){
  
  Tissue_i <- unique(C_sum$Tissue)
  
  C_plot <- ggplot(C_sum, aes(Time/365, median)) +
    geom_ribbon(aes(ymin = q25,
                    ymax = q75),
                fill = "grey70") +
    geom_path()+
    labs(title = paste0(PFAS_i),
         x = "Time (years)",
         y = paste0(Tissue_i, " Concentration (ng/ml)")) +
    theme_quick
  
  return(C_plot)
  
}

plot_all_tissues <- function(C_ts){
  
  C_ts_plot <- ggplot(C_ts, aes(Time, median)) +
    geom_ribbon(aes(ymin = q25,
                    ymax = q75),
                fill = "grey70") +
    geom_path() +
    facet_wrap(~ Tissue, scales = "free", nrow = 2) +
    labs(title = paste0(PFAS_i),
         x = "Time (years)",
         y = "Concentration (ng/ml)") +
    theme_quick
  
  return(C_ts_plot)
  
}