# find parameters that gave errors ----


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


nC_MW_df <- read.csv(here("Input", "nC_MW_df.csv")) |>
  mutate(PFAS = str_replace_all(PFAS, "_", "-"))


## load errored files ----
LOC_ER_FILES <- here("Output",
                     "UA_with_updated_Ps"
)

errors_df <-
  list.files(path = here(LOC_ER_FILES),
             pattern = "\\_errored_iterations.csv$") |>
  map_df(~ read_csv(file = here(LOC_ER_FILES, .),
                    id = "name")) 

errored_iterations_df <- errors_df |>  
  # select(!"...2") |> 
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/Admin/Desktop/PhD_Chrysa/CP_L_R/PARC_PFAS_PBK/Output/UA_with_updated_Ps/",
                                "")) |> 
  mutate(PFAS = str_replace_all(PFAS,
                                "_errored_iterations.csv","")) |> 
  select(!name) |> 
  rename(Iteration = x) |>
  unite("PFAS_iteration", PFAS:Iteration, sep = "_") |> 
  pull()

# load final parameters ----
LOC_FINALPARAMETER_FILES <- LOC_ER_FILES

actual_simulated_parameters_df <-
  list.files(path = here(LOC_FINALPARAMETER_FILES),
             pattern = "\\_parameters.csv$") |>
  # imap_df(~ read_csv(file = here(location, .x)) |>
  #           mutate(file = .y, row_num = row_number()))
  # 
  map_df(~ read_csv(file = here(LOC_FINALPARAMETER_FILES, .),
                    id = "name"))|> 
  rename(Iteration = "...2") |> 
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/Admin/Desktop/PhD_Chrysa/CP_L_R/PARC_PFAS_PBK/Output/UA_with_updated_Ps/",
                                "")) |> 
  mutate(PFAS = str_replace_all(PFAS,
                                "_final_model_parameters.csv",""))|>
  select(!name) |> 
  unite("PFAS_iteration", c(PFAS,Iteration), sep = "_") 

errored_simulated_parameters <- actual_simulated_parameters_df |> 
  filter(PFAS_iteration %in% errored_iterations_df$PFAS_iteration)

errored_simulated_parameters_long <- errored_simulated_parameters |> 
  pivot_longer(!PFAS_iteration, names_to = "Parameter", values_to = "Value") |> 
  separate(PFAS_iteration, into = c("PFAS", "Iteration"), sep = "_")

errored_simulated_parameters_df <- full_join(errored_simulated_parameters_long, 
                                             nC_MW_df) 

good_simulated_parameters <- actual_simulated_parameters_df |> 
  filter_out(PFAS_iteration %in% errored_iterations_df$PFAS_iteration)

good_simulated_parameters_long <- good_simulated_parameters |> 
  pivot_longer(!PFAS_iteration, names_to = "Parameter", values_to = "Value") |> 
  separate(PFAS_iteration, into = c("PFAS", "Iteration"), sep = "_")

good_simulated_parameters_df <- full_join(good_simulated_parameters_long, 
                                          nC_MW_df) 
# load experimental designs ----
LOC_DESIGN_FILES <- LOC_ER_FILES

experimental_design_df <-
  list.files(path = here(LOC_DESIGN_FILES),
             pattern = "\\MC_ExpDesign.csv$") |>
  map_df(~ read_csv(file = here(LOC_DESIGN_FILES, .),
                    id = "name")) |> 
  mutate(PFAS = str_replace_all(name,
                                "C:/Users/Admin/Desktop/PhD_Chrysa/CP_L_R/PARC_PFAS_PBK/Output/UA_with_updated_Ps/",
                                "")) |> 
  mutate(PFAS = str_replace_all(PFAS,
                                "_LHS_MC_ExpDesign.csv",""),
         Iteration = rep(seq(1,1000, by = 1), 17))|>
  select(!name) |> 
  unite("PFAS_iteration", c(PFAS,Iteration), sep = "_") 

experimental_design_long <- experimental_design_df |> 
  pivot_longer(!PFAS_iteration, names_to = "Parameter", values_to = "Value") |> 
  separate(PFAS_iteration, into = c("PFAS", "Iteration"), sep = "_") #|> 
  # mutate(Value = ifelse(str_detect(Parameter, "logD_"), log(Value), 
  #                       ifelse(Parameter == "Papp_SI", log(Value), 
  #                              Value))
  # ) 
experimental_design_long <- full_join(experimental_design_long, 
                                      nC_MW_df) 


errored_experimental_design_df <- experimental_design_df |> 
  filter(PFAS_iteration %in% errored_iterations_df$PFAS_iteration)
errored_designs_long <- errored_experimental_design_df |> 
  pivot_longer(!PFAS_iteration, names_to = "Parameter", values_to = "Value") |> 
  separate(PFAS_iteration, into = c("PFAS", "Iteration"), sep = "_")

errored_designs_long_df <- full_join(errored_designs_long, 
                                     nC_MW_df) 

good_experimental_design_df <- experimental_design_df |> 
  filter_out(PFAS_iteration %in% errored_iterations_df$PFAS_iteration)
good_experimental_design_long <- good_experimental_design_df |> 
  pivot_longer(!PFAS_iteration, names_to = "Parameter", values_to = "Value") |> 
  separate(PFAS_iteration, into = c("PFAS", "Iteration"), sep = "_") #|> 
  # mutate(Value = ifelse(str_detect(Parameter, "logD_"), log(Value), 
  #                       ifelse(Parameter == "Papp_SI", log(Value), 
  #                              Value))
  # ) 
good_experimental_design_df <- full_join(good_experimental_design_long, 
                                         nC_MW_df) 




### ADONA ----

# Design
ADONA_error_params <- errored_designs_long_df |> 
  filter(PFAS == "ADONA") |> mutate(Result = "error")
ADONA_good_params <- good_experimental_design_df |> 
  filter(PFAS == "ADONA") |> mutate(Result = "good")
ADONA_design <- rbind(ADONA_error_params, ADONA_good_params)

ADONA_ErrorGood_D_plot <- ggplot(data = ADONA_design, 
       aes(x = Value, fill = factor(Result, c("good", "error")))
      ) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "ADONA Design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "right")
ADONA_ErrorGood_D_plot

ADONA_Error_D_plot <- ggplot(data = ADONA_error_params, 
       aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "ADONA Error design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
ADONA_Error_D_plot

ADONA_Good_D_plot <- ggplot(data = ADONA_good_params, 
       aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "ADONA Good design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
ADONA_Good_D_plot

# Simulated parameters
ADONA_errored_simulated_parameters <- errored_simulated_parameters_df |> 
  filter(PFAS == "ADONA") |> mutate(Result = "error")
ADONA_good_simulated_parameters <- good_simulated_parameters_df |> 
  filter(PFAS == "ADONA") |> mutate(Result = "good")

ADONA_parameters <- rbind(ADONA_errored_simulated_parameters, ADONA_good_simulated_parameters)

ADONA_ErrorGood_P_plot <- ggplot(data = ADONA_parameters, 
       aes(x = Value, fill = factor(Result, c("good", "error")))) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "ADONA Final parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "right")
ADONA_ErrorGood_P_plot

ADONA_Error_P_plot <- ggplot(data = ADONA_errored_simulated_parameters, 
       aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "ADONA Error parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
ADONA_Error_P_plot

ADONA_Good_P_plot <- ggplot(data = ADONA_good_simulated_parameters, 
       aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "ADONA Good parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
ADONA_Good_P_plot

P1 <- ADONA_ErrorGood_D_plot + Patchwork_theme
P2 <- ADONA_Error_D_plot
P3 <- ADONA_Good_D_plot

P4 <- ADONA_ErrorGood_P_plot + Patchwork_theme
P5 <- ADONA_Error_P_plot
P6 <- ADONA_Good_P_plot

ggsave(plot = P1,
       file = here(OUTPUT, "ADONA_ErrorGood_Design_plot.png"),
       dpi = 600,
       width = 24,
       height = 24,
       units = "cm")

Good_Error_D_plot <- (P2 + P3) &
  Patchwork_theme
ggsave(plot = Good_Error_D_plot,
       file = here(OUTPUT, "ADONA_GoodandError_Design_plot.png"),
       dpi = 600,
       width = 48,
       height = 12,
       units = "cm")

ggsave(plot = P4,
       file = here(OUTPUT, "ADONA_ErrorGood_FinalP_plot.png"),
       dpi = 600,
       width = 24,
       height = 24,
       units = "cm")

Good_Error_P_plot <- (P5 + P6) &
  Patchwork_theme
ggsave(plot = Good_Error_P_plot,
       file = here(OUTPUT, "ADONA_GoodandError_FinalP_plot.png"),
       dpi = 600,
       width = 48,
       height = 12,
       units = "cm")


### HFPO_DA ----

# Design
HFPO_DA_error_params <- errored_designs_long_df |> 
  filter(PFAS == "HFPO-DA") |> mutate(Result = "error")
HFPO_DA_good_params <- good_experimental_design_df |> 
  filter(PFAS == "HFPO-DA") |> mutate(Result = "good")
HFPO_DA_design <- rbind(HFPO_DA_error_params, HFPO_DA_good_params)

HFPO_DA_ErrorGood_D_plot <- ggplot(data = HFPO_DA_design, 
                                 aes(x = Value, fill = factor(Result, c("good", "error")))
) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "HFPO_DA Design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "right")
HFPO_DA_ErrorGood_D_plot

HFPO_DA_Error_D_plot <- ggplot(data = HFPO_DA_error_params, 
                             aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "HFPO_DA Error design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
HFPO_DA_Error_D_plot

HFPO_DA_Good_D_plot <- ggplot(data = HFPO_DA_good_params, 
                            aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "HFPO_DA Good design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
HFPO_DA_Good_D_plot

# Simulated parameters
HFPO_DA_errored_simulated_parameters <- errored_simulated_parameters_df |> 
  filter(PFAS == "HFPO-DA") |> mutate(Result = "error")
HFPO_DA_good_simulated_parameters <- good_simulated_parameters_df |> 
  filter(PFAS == "HFPO-DA") |> mutate(Result = "good")

HFPO_DA_parameters <- rbind(HFPO_DA_errored_simulated_parameters, HFPO_DA_good_simulated_parameters)

HFPO_DA_ErrorGood_P_plot <- ggplot(data = HFPO_DA_parameters, 
                                 aes(x = Value, fill = factor(Result, c("good", "error")))) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "HFPO_DA Final parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "right")
HFPO_DA_ErrorGood_P_plot

HFPO_DA_Error_P_plot <- ggplot(data = HFPO_DA_errored_simulated_parameters, 
                             aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "HFPO_DA Error parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
HFPO_DA_Error_P_plot

HFPO_DA_Good_P_plot <- ggplot(data = HFPO_DA_good_simulated_parameters, 
                            aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "HFPO_DA Good parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
HFPO_DA_Good_P_plot

P1 <- HFPO_DA_ErrorGood_D_plot + Patchwork_theme
P2 <- HFPO_DA_Error_D_plot
P3 <- HFPO_DA_Good_D_plot

P4 <- HFPO_DA_ErrorGood_P_plot + Patchwork_theme
P5 <- HFPO_DA_Error_P_plot
P6 <- HFPO_DA_Good_P_plot

ggsave(plot = P1,
       file = here(OUTPUT, "HFPO_DA_ErrorGood_Design_plot.png"),
       dpi = 600,
       width = 24,
       height = 24,
       units = "cm")

Good_Error_D_plot <- (P2 + P3) &
  Patchwork_theme
ggsave(plot = Good_Error_D_plot,
       file = here(OUTPUT, "HFPO_DA_GoodandError_Design_plot.png"),
       dpi = 600,
       width = 48,
       height = 12,
       units = "cm")

ggsave(plot = P4,
       file = here(OUTPUT, "HFPO_DA_ErrorGood_FinalP_plot.png"),
       dpi = 600,
       width = 24,
       height = 24,
       units = "cm")

Good_Error_P_plot <- (P5 + P6) &
  Patchwork_theme
ggsave(plot = Good_Error_P_plot,
       file = here(OUTPUT, "HFPO_DA_GoodandError_FinalP_plot.png"),
       dpi = 600,
       width = 48,
       height = 12,
       units = "cm")



### PFBS ----

# Design
PFBS_error_params <- errored_designs_long_df |> 
  filter(PFAS == "PFBS") |> mutate(Result = "error")
PFBS_good_params <- good_experimental_design_df |> 
  filter(PFAS == "PFBS") |> mutate(Result = "good")
PFBS_design <- rbind(PFBS_error_params, PFBS_good_params)

PFBS_ErrorGood_D_plot <- ggplot(data = PFBS_design, 
                                 aes(x = Value, fill = factor(Result, c("good", "error")))
) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBS Design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "right")
PFBS_ErrorGood_D_plot

PFBS_Error_D_plot <- ggplot(data = PFBS_error_params, 
                             aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBS Error design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
PFBS_Error_D_plot

PFBS_Good_D_plot <- ggplot(data = PFBS_good_params, 
                            aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBS Good design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
PFBS_Good_D_plot

# Simulated parameters
PFBS_errored_simulated_parameters <- errored_simulated_parameters_df |> 
  filter(PFAS == "PFBS") |> mutate(Result = "error")
PFBS_good_simulated_parameters <- good_simulated_parameters_df |> 
  filter(PFAS == "PFBS") |> mutate(Result = "good")

PFBS_parameters <- rbind(PFBS_errored_simulated_parameters, PFBS_good_simulated_parameters)

PFBS_ErrorGood_P_plot <- ggplot(data = PFBS_parameters, 
                                 aes(x = Value, fill = factor(Result, c("good", "error")))) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBS Final parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "right")
PFBS_ErrorGood_P_plot

PFBS_Error_P_plot <- ggplot(data = PFBS_errored_simulated_parameters, 
                             aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBS Error parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
PFBS_Error_P_plot

PFBS_Good_P_plot <- ggplot(data = PFBS_good_simulated_parameters, 
                            aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBS Good parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
PFBS_Good_P_plot

P1 <- PFBS_ErrorGood_D_plot + Patchwork_theme
P2 <- PFBS_Error_D_plot
P3 <- PFBS_Good_D_plot

P4 <- PFBS_ErrorGood_P_plot + Patchwork_theme
P5 <- PFBS_Error_P_plot
P6 <- PFBS_Good_P_plot

ggsave(plot = P1,
       file = here(OUTPUT, "PFBS_ErrorGood_Design_plot.png"),
       dpi = 600,
       width = 24,
       height = 24,
       units = "cm")

Good_Error_D_plot <- (P2 + P3) &
  Patchwork_theme
ggsave(plot = Good_Error_D_plot,
       file = here(OUTPUT, "PFBS_GoodandError_Design_plot.png"),
       dpi = 600,
       width = 48,
       height = 12,
       units = "cm")

ggsave(plot = P4,
       file = here(OUTPUT, "PFBS_ErrorGood_FinalP_plot.png"),
       dpi = 600,
       width = 24,
       height = 24,
       units = "cm")

Good_Error_P_plot <- (P5 + P6) &
  Patchwork_theme
ggsave(plot = Good_Error_P_plot,
       file = here(OUTPUT, "PFBS_GoodandError_FinalP_plot.png"),
       dpi = 600,
       width = 48,
       height = 12,
       units = "cm")



### PFBA ----

# Design
PFBA_error_params <- errored_designs_long_df |> 
  filter(PFAS == "PFBA") |> mutate(Result = "error")
PFBA_good_params <- good_experimental_design_df |> 
  filter(PFAS == "PFBA") |> mutate(Result = "good")
PFBA_design <- rbind(PFBA_error_params, PFBA_good_params)

PFBA_ErrorGood_D_plot <- ggplot(data = PFBA_design, 
                                 aes(x = Value, fill = factor(Result, c("good", "error")))
) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBA Design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "right")
PFBA_ErrorGood_D_plot

PFBA_Error_D_plot <- ggplot(data = PFBA_error_params, 
                             aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBA Error design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
PFBA_Error_D_plot

PFBA_Good_D_plot <- ggplot(data = PFBA_good_params, 
                            aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBA Good design",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
PFBA_Good_D_plot

# Simulated parameters
PFBA_errored_simulated_parameters <- errored_simulated_parameters_df |> 
  filter(PFAS == "PFBA") |> mutate(Result = "error")
PFBA_good_simulated_parameters <- good_simulated_parameters_df |> 
  filter(PFAS == "PFBA") |> mutate(Result = "good")

PFBA_parameters <- rbind(PFBA_errored_simulated_parameters, PFBA_good_simulated_parameters)

PFBA_ErrorGood_P_plot <- ggplot(data = PFBA_parameters, 
                                 aes(x = Value, fill = factor(Result, c("good", "error")))) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBA Final parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "right")
PFBA_ErrorGood_P_plot

PFBA_Error_P_plot <- ggplot(data = PFBA_errored_simulated_parameters, 
                             aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBA Error parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
PFBA_Error_P_plot

PFBA_Good_P_plot <- ggplot(data = PFBA_good_simulated_parameters, 
                            aes(x = Value)) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_y", nrow = 5) +
  coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "Result") +
  labs(title = "PFBA Good parameters",
       x = "Parameter value",
       y = "Frequency") +
  theme_quick +
  theme(legend.position = "none")
PFBA_Good_P_plot

P1 <- PFBA_ErrorGood_D_plot + Patchwork_theme
P2 <- PFBA_Error_D_plot
P3 <- PFBA_Good_D_plot

P4 <- PFBA_ErrorGood_P_plot + Patchwork_theme
P5 <- PFBA_Error_P_plot
P6 <- PFBA_Good_P_plot

ggsave(plot = P1,
       file = here(OUTPUT, "PFBA_ErrorGood_Design_plot.png"),
       dpi = 600,
       width = 24,
       height = 24,
       units = "cm")

Good_Error_D_plot <- (P2 + P3) &
  Patchwork_theme
ggsave(plot = Good_Error_D_plot,
       file = here(OUTPUT, "PFBA_GoodandError_Design_plot.png"),
       dpi = 600,
       width = 48,
       height = 12,
       units = "cm")

ggsave(plot = P4,
       file = here(OUTPUT, "PFBA_ErrorGood_FinalP_plot.png"),
       dpi = 600,
       width = 24,
       height = 24,
       units = "cm")

Good_Error_P_plot <- (P5 + P6) &
  Patchwork_theme
ggsave(plot = Good_Error_P_plot,
       file = here(OUTPUT, "PFBA_GoodandError_FinalP_plot.png"),
       dpi = 600,
       width = 48,
       height = 12,
       units = "cm")















# rest

ggplot(data = experimental_design_long, 
       aes(x = Value, fill = as.factor(n_C))) + 
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free_x") +
  # coord_flip() +
  scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
  labs(title = "Uncertainty Analysis experimental design",
       x = "Parameter value",
       y = "Frequency of occurence in errored simulations") +
  theme_quick +
  theme(legend.position = "right")

# 
# cleanexperimental_design_long <- full_join(cleanexperimental_design_long, 
#                                            nC_MW_df) 
# 
# UA_experimental_design_plot <- ggplot() + 
#   geom_histogram(data = experimental_design_long, 
#                  aes(x = Value), 
#                  fill = "grey70") +
#   facet_wrap(~ Parameter, scales = "free_x") +
#   geom_histogram(data = cleanexperimental_design_long, 
#                  aes(x = Value, fill = as.factor(n_C))) +
#   facet_wrap(~ Parameter, scales = "free_x") +
#   # coord_flip() +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
#   labs(title = "Uncertainty analysis experimental design",
#        x = "Parameter value",
#        y = "Frequency") +
#   theme_quick +
#   theme(legend.position = "right")
# 
# UA_experimental_design_plot <- UA_experimental_design_plot + Patchwork_theme +
#   theme(legend.position = "right",
#         legend.title = element_text(size = 35, face = "bold"),
#         axis.text = element_text(size = 35))
# 
# ggsave(plot = UA_experimental_design_plot,
#        file = here(OUTPUT, "UA_experimental_design_plot.png"),
#        dpi = 600,
#        width = 24,
#        height = 18,
#        units = "cm")
# 
# 
# 
# errored_simulated_parameters_long <- full_join(errored_simulated_parameters_long, 
#                                                nC_MW_df) 
# errored_designs_long <- full_join(errored_designs_long, 
#                                   nC_MW_df) 
# 
# loged_errored_designs_long <- errored_designs_long |> 
#   mutate(Value = ifelse(str_detect(Parameter, "logD_"), log(Value), 
#                         ifelse(Parameter == "Papp_SI", log(Value), 
#                                Value))
#   ) 
# 
# Errored_designs_plot <- ggplot(data = loged_errored_designs_long, 
#                                aes(x = Value, fill = as.factor(n_C))) + 
#   geom_histogram() +
#   facet_wrap(~ Parameter, scales = "free_x") +
#   # coord_flip() +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
#   labs(title = "Uncertainty analysis experimental design that caused errors",
#        x = "Parameter value",
#        y = "Frequency") +
#   theme_quick +
#   theme(legend.position = "right")
# 
# Errored_designs_plot <- Errored_designs_plot + Patchwork_theme +
#   theme(legend.position = "right",
#         legend.title = element_text(size = 35, face = "bold"),
#         axis.text = element_text(size = 35))
# 
# ggsave(plot = Errored_designs_plot,
#        file = here(OUTPUT, "Errored_designs_plot.png"),
#        dpi = 600,
#        width = 24,
#        height = 18,
#        units = "cm")
# 
# 
# 
# errored_final_simulated_parameters_long <- errored_simulated_parameters_long |> 
#   filter(Parameter %in% c("Clear_uptake_I",
#                           "Clear_uptake_L",
#                           "Clear_excretion_L",
#                           "Clear_net_reabsorption",
#                           "fup"))
# 
# Errored_final_simulated_parameters_plot <- ggplot(data = errored_final_simulated_parameters_long, 
#                                                   aes(x = Value, fill = as.factor(n_C))) + 
#   geom_histogram() +
#   facet_wrap(~ factor(Parameter,
#                       c("Clear_uptake_I",
#                         "Clear_uptake_L",
#                         "Clear_excretion_L",
#                         "Clear_net_reabsorption",
#                         "fup")), 
#              scales = "free_x", nrow = 1) +
#   # coord_flip() +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
#   labs(title = "Calculated model parameter distribution in that caused errors",
#        x = "Parameter value",
#        y = "Frequency") +
#   theme_quick +
#   theme(legend.position = "right"
#   ) 
# 
# # Errored_final_simulated_parameters_plot <- Errored_final_simulated_parameters_plot + Patchwork_theme +
# #   theme(legend.position = "right",
# #         legend.title = element_text(size = 35, face = "bold"),
# #         axis.text = element_text(size = 35))
# # 
# # ggsave(plot = Errored_final_simulated_parameters_plot,
# #        file = here(OUTPUT, "Errored_final_simulated_parameters_plot.png"),
# #        dpi = 600,
# #        width = 24,
# #        height = 6,
# #        units = "cm")
# 
# 
# 
# CL_out_L_design_long <- loged_errored_designs_long |> 
#   filter(Parameter %in% c("SF_BSEP", "Clear_BSEPc"
#   ))
# 
# CL_out_L_design_plot <- ggplot(data = CL_out_L_design_long, 
#                                aes(x = Value, fill = as.factor(n_C))) + 
#   geom_histogram() +
#   facet_wrap(~ Parameter, 
#              scales = "free_x") +
#   # coord_flip() +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
#   labs(title = "Parameters related to: Clear_excretion_L",
#        x = "Parameter value",
#        y = "Frequency") +
#   theme_quick +
#   theme(legend.position = "right")
# 
# 
# # CL_out_L_design_plot <- CL_out_L_design_plot + Patchwork_theme +
# #   theme(legend.position = "none",
# #         legend.title = element_text(size = 35, face = "bold"),
# #         axis.text = element_text(size = 35))
# # 
# # 
# # ggsave(plot = CL_out_L_design_plot,
# #        file = here(OUTPUT, "CL_out_L_design_plot.png"),
# #        dpi = 600,
# #        width = 13,
# #        height = 11,
# #        units = "cm")
# 
# 
# 
# CL_reab_errored_design_long <- loged_errored_designs_long |> 
#   filter(Parameter %in% c("logD_AlbuminWater",
#                           "logD_GlobulinWater",
#                           "SF_OAT1", "Clear_OAT1c",
#                           "SF_OAT3", "Clear_OAT3c",
#                           "SF_OAT4", "Clear_OAT4c"
#   ))
# 
# CL_reab_errored_design_plot <- ggplot(data = CL_reab_errored_design_long, 
#                                       aes(x = Value, fill = as.factor(n_C))) + 
#   geom_histogram() +
#   facet_wrap(~ Parameter, 
#              scales = "free_x") +
#   # coord_flip() +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
#   labs(title = "Parameters related to: Clear_net_reabsorption",
#        x = "Parameter value",
#        y = "Frequency") +
#   theme_quick +
#   theme(legend.position = "right")
# 
# # CL_reab_errored_design_plot <- CL_reab_errored_design_plot + CP_theme +
# #     theme(legend.position = "none")
# # 
# # 
# # ggsave(plot = CL_reab_errored_design_plot,
# #        file = here(OUTPUT, "CL_reab_errored_design_plot.png"),
# #        dpi = 600,
# #        width = 16,
# #        height = 11,
# #        units = "cm")
# 
# 
# 
# 
# CL_up_I_errored_design_long <- loged_errored_designs_long |> 
#   filter(Parameter %in% c("Papp_SI",
#                           "SF_OATP2B1_I", "Clear_OATP2B1c"
#   ))
# 
# CL_up_I_errored_design_plot <- ggplot(data = CL_up_I_errored_design_long, 
#                                       aes(x = Value, fill = as.factor(n_C))) + 
#   geom_histogram() +
#   facet_wrap(~ Parameter, 
#              scales = "free_x") +
#   # coord_flip() +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
#   labs(title = "Parameters related to: Clear_uptake_I",
#        x = "Parameter value",
#        y = "Frequency") +
#   theme_quick +
#   theme(legend.position = "right")
# 
# # CL_up_I_errored_design_plot <- CL_up_I_errored_design_plot + CP_theme +
# #   theme(legend.position = "none")
# # 
# # ggsave(plot = CL_up_I_errored_design_plot,
# #        file = here(OUTPUT, "CL_up_I_errored_design_plot.png"),
# #        dpi = 600,
# #        width = 11,
# #        height = 11,
# #        units = "cm")
# 
# 
# CL_up_L_errored_design_long <- loged_errored_designs_long |> 
#   filter(Parameter %in% c("logD_AlbuminWater",
#                           "logD_GlobulinWater",
#                           "SF_OATP1B1", "Clear_OATP1B1c",
#                           "SF_OATP1B3", "Clear_OATP1B3c",
#                           "SF_OATP2B1", "Clear_OATP2B1c",
#                           "SF_NTCP", "Clear_NTCPc"
#   )
#   )
# 
# CL_up_L_errored_design_plot <- ggplot(data = CL_up_L_errored_design_long, 
#                                       aes(x = Value, fill = as.factor(n_C))) + 
#   geom_histogram() +
#   facet_wrap(~ Parameter, 
#              scales = "free_x") +
#   # coord_flip() +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
#   labs(title = "Parameters related to: Clear_uptake_L",
#        x = "Parameter value",
#        y = "Frequency") +
#   theme_quick +
#   theme(legend.position = "right")
# 
# 
# # CL_up_L_errored_design_plot <- CL_up_L_errored_design_plot + CP_theme +
# #   theme(legend.position = "none")
# # 
# # ggsave(plot = CL_up_L_errored_design_plot,
# #        file = here(OUTPUT, "CL_up_L_errored_design_plot.png"),
# #        dpi = 600,
# #        width = 16,
# #        height = 11,
# #        units = "cm")
# 
# 
# 
# fup_L_errored_design_long <- loged_errored_designs_long |> 
#   filter(Parameter %in% c("logD_AlbuminWater",
#                           "logD_GlobulinWater")
#   )
# 
# fup_L_errored_design_plot <- ggplot(data = fup_L_errored_design_long, 
#                                     aes(x = Value, fill = as.factor(n_C))) + 
#   geom_histogram() +
#   facet_wrap(~ Parameter, 
#              scales = "free_x") +
#   # coord_flip() +
#   scale_fill_viridis(option = "D", discrete = TRUE, end = 0.7, name = "nC") +
#   labs(title = "Parameters related to: fup",
#        x = "Parameter value",
#        y = "Frequency") +
#   theme_quick +
#   theme(legend.position = "right")
# 
# # fup_L_errored_design_plot <- fup_L_errored_design_plot + CP_theme +
# #   theme(legend.position = "none")
# # 
# # ggsave(plot = fup_L_errored_design_plot,
# #        file = here(OUTPUT, "fup_L_errored_design_plot.png"),
# #        dpi = 600,
# #        width = 11,
# #        height = 11,
# #        units = "cm")
# 
# #   facet_wrap(~ factor(Parameter, c("SF_OATP2B1_I",
# #                                    "SF_OATP2B1", "SF_OATP1B3", "SF_OATP1B1", "SF_NTCP", "SF_BSEP",
# #                                    "SF_OAT1", "SF_OAT3", "SF_OAT4",
# #                                    "fup", "fu_L_ec", "fu_PTL", "fu_PT_ec",
# #                                    "Clear_apparent_uptake_I", "Clear_OATP2B1_uptake_I", 
# #                                    "Clear_OATP1B1_uptake_L", "Clear_OATP1B3_uptake_L", "Clear_OATP2B1_uptake_L", "Clear_NTCP_uptake_L", 
# #                                    "Clear_BSEP_elimination_L", 
# #                                    "Clear_OAT3_elimination_K", "Clear_OAT1_elimination_K", "Clear_total_elimination_K", "Clear_OAT4_reabsorption_K", 
# #                                    "Clear_total_uptake_I", "Clear_total_uptake_L", "Clear_total_elimination_L", "Clear_net_reabsorption_K"
# #                                    )), scales = "free") +
# #   labs(x = "Value (Clear: L/day, SF: unitless, fu: unitless)",
# #        y = "Frequency") +
# #   coord_flip()+
# #   theme_quick
# # 
# 
# 
# 
