# Packages
library(readxl)
library(xlsx)
library(here)
library(tidyverse)
library(showtext)
library(ggplot2)
library(purrr)
library(knitr)
library(patchwork)
library(dendroTools)

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

# Transporter protein abudance ----

ProteinAb <-  read_excel(here("Input", "ProteinAbundance.xlsx"))
ProteinAb_df <- ProteinAb %>% mutate(Reference = factor(ProteinAb$Reference),
                                     Transporter = factor(ProteinAb$Transporter),
                                     Localisation = factor(ProteinAb$Localisation),
                                     min = Protein_abundance-SD,
                                     max = Protein_abundance+SD)

Meta_df <- ProteinAb_df %>%
  group_by(Transporter, Localisation) %>%
  summarise(.groups = "drop",
            minimum = min(min, na.rm = TRUE), 
            maximum = max(max, na.rm = TRUE),
            mean = mean(Protein_abundance, na.rm = TRUE), 
            sd = sd(Protein_abundance, na.rm = TRUE)) 

Calc_SF_df <- Meta_df %>% 
  filter(! Localisation %in% c("ileum", "duodenum", "jejunum", "colon"), 
         Transporter != "ASBT") %>% 
  mutate(minimum = pmax(minimum, 0.001)) %>% 
  mutate(vitro = c(1001.82, 1.202, 2.596, 0.711, 0.063, 0.120, 0.719, 0.812, 0.812),
         mgProtein.mgT = c(0.35, 0.35, 0.21, 0.21, 0.21, 0.35, 0.35, 0.13, 0.35)) %>% 
  pivot_longer(minimum:mean, names_to = "id", values_to = "vivo") %>% 
  select(!sd) %>% 
  mutate(REF = vivo/vitro,
         SF = mgProtein.mgT*REF*1e6)

Calc_SF_wide <- Calc_SF_df %>% 
  select(Transporter, Localisation, id, SF) %>% 
  pivot_wider(names_from = id, values_from = SF) %>% 
  mutate("Mult100" = mean*100,
         "Div100" = mean/100) %>% 
  select(!mean) %>% 
  pivot_longer(minimum:Div100, names_to = "id", values_to = "SF")

Tests <- Calc_SF_wide %>% 
  select(Transporter, Localisation, id, SF) %>% 
  mutate(What = "SF_",
         Localisation = ifelse(Localisation == "intestine", "_I", "")) %>% 
  unite(P, What,Transporter, Localisation, sep="") %>% 
  mutate(Test = paste0(id, "_", P)) %>% 
  select(!id)

write.csv(Tests, "LSA_parameters.csv", row.names = FALSE)


Liver_SF <- ggplot(filter(Calc_SF_wide, Localisation == "liver"), 
       aes(x = Transporter, y = log10(mean))) +
  # geom_jitter(aes(Transporter, Protein_abundance, colour = Reference))+
  geom_pointrange(aes(ymin = log10(minimum), 
                      ymax = log10(maximum)),
                  position=position_jitter(width=0.3), 
                  # linetype='dotted',
                  size = 0.3,
                  linewidth = 0.5) +
  # scale_y_continuous(breaks=seq(0,30, by = 5)) +
  labs(title = "Liver", 
       x = "",
       y = "Scaling Factor (log10)") +
  CP_theme

ggsave(plot = Liver_SF,
       file = "Liver_SF.png",
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")


Kidney_SF <- ggplot(filter(Calc_SF_wide, Localisation == "kidney"), 
                   aes(x = Transporter, y = log10(mean))) +
  # geom_jitter(aes(Transporter, Protein_abundance, colour = Reference))+
  geom_pointrange(aes(ymin = log10(minimum), 
                      ymax = log10(maximum)),
                  position=position_jitter(width=0.3), 
                  # linetype='dotted',
                  size = 0.3,
                  linewidth = 0.5) +
  # scale_y_continuous(breaks=seq(0,30, by = 5)) +
  labs(title = "Kidney", 
       x = "",
       y = "Scaling Factor (log10)") +
  CP_theme

ggsave(plot = Kidney_SF,
       file = "Kidney_SF.png",
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

Intestine_SF <- ggplot(filter(Calc_SF_wide, Localisation == "intestine"), 
                    aes(x = Transporter, y = log10(mean))) +
  # geom_jitter(aes(Transporter, Protein_abundance, colour = Reference))+
  geom_pointrange(aes(ymin = log10(minimum), 
                      ymax = log10(maximum)),
                  position=position_jitter(width=0.3), 
                  # linetype='dotted',
                  size = 0.3,
                  linewidth = 0.5) +
  # scale_y_continuous(breaks=seq(0,30, by = 5)) +
  labs(title = "Intestine", 
       x = "",
       y = "Scaling Factor (log10)") +
  CP_theme

ggsave(plot = Intestine_SF,
       file = "Intestine_SF.png",
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")


Liver_p <- ggplot(filter(ProteinAb_df, Localisation == "liver"), 
       aes(x = Transporter, y = Protein_abundance, colour = Reference)) +
  # geom_jitter(aes(Transporter, Protein_abundance, colour = Reference))+
  geom_pointrange(aes(ymin = Protein_abundance-SD, 
                      ymax = Protein_abundance+SD),
                  position=position_jitter(width=0.3), 
                  # linetype='dotted',
                  size = 0.3,
                  linewidth = 0.5) +
  scale_y_continuous(breaks=seq(0,30, by = 5)) +
  labs(title = "Liver", 
       x = "",
       y = "Protein abundance (pmol/mg membrane protein)") +
  CP_theme
  
ggsave(plot = Liver_p,
       file = "Liver_p.png",
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

Intestine_p <- ggplot(filter(ProteinAb_df, Localisation != "liver", Localisation != "kidney"), 
                  aes(x = Transporter, y = Protein_abundance, colour = Reference, shape = Localisation)) +
  # geom_jitter(aes(Transporter, Protein_abundance, colour = Reference))+
  geom_pointrange(aes(ymin = Protein_abundance-SD, 
                      ymax = Protein_abundance+SD),
                  position=position_jitter(width=0.3), 
                  # linetype='dotted',
                  size = 0.3,
                  linewidth = 0.5) +
  scale_y_continuous(breaks=seq(0,60, by = 5)) +
  labs(title = "Intestine", 
       x = "",
       y = "Protein abundance (pmol/mg membrane protein)") +
  CP_theme
Intestine_p
ggsave(plot = Intestine_p,
       file = "Intestine_p.png",
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

Intestine_p_log <- ggplot(filter(ProteinAb_df, Localisation != "liver", Localisation != "kidney"), 
                      aes(x = Transporter, y = log10(Protein_abundance), colour = Reference, shape = Localisation)) +
  # geom_jitter(aes(Transporter, Protein_abundance, colour = Reference))+
  geom_pointrange(aes(ymin = log10(Protein_abundance-SD), 
                      ymax = log10(Protein_abundance+SD)),
                  position=position_jitter(width=0.3), 
                  # linetype='dotted',
                  size = 0.3,
                  linewidth = 0.5) +
  # scale_y_continuous(breaks=seq(0,60, by = 5)) +
  labs(title = "Intestine", 
       x = "",
       y = "Protein abundance log10(pmol/mg membrane protein)") +
  CP_theme
Intestine_p_log
ggsave(plot = Intestine_p_log,
       file = "Intestine_p_log.png",
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

Kidney_p <-   ggplot(filter(ProteinAb_df, Localisation == "kidney"), 
                  aes(x = Transporter, y = Protein_abundance, colour = Reference)) +
  # geom_jitter(aes(Transporter, Protein_abundance, colour = Reference))+
  geom_pointrange(aes(ymin = Protein_abundance-SD, 
                      ymax = Protein_abundance+SD),
                  position=position_jitter(width=0.3), 
                  # linetype='dotted',
                  size = 0.3,
                  linewidth = 0.5) +
  scale_y_continuous(breaks=seq(0,25, by = 5)) +
  labs(title = "Kidney", 
       x = "",
       y = "Protein abundance (pmol/mg membrane protein)") +
  CP_theme

ggsave(plot = Kidney_p,
       file = "Kidney_p.png",
       dpi = 600,
       width = 12,
       height = 9,
       units = "cm")

# Load data ----
DataValidation <- read_excel(here("DeriveParameters", "Data_availability.xlsx"), sheet = "ValidationData")

Datab <- read_excel(here("DeriveParameters", "PFAS_database.xlsx"), 
     col_types = c("text", "text", "text", 
          "text", "text", "text", "text", 
          "text", "numeric", "numeric", "text", 
          "text", "text", "text", "numeric", 
          "text", "text", "text", "text", "text", 
          "text", "numeric", "text", "text",  
          "text")) 

DataAv <- read_excel(here("DeriveParameters", "Data_availability.xlsx"), sheet = "Selected_NAM_Data_Summarized")

Include <- read_excel(here("DeriveParameters", "Data_availability.xlsx"), sheet = "CombinedDataAvailability")
IncludedPFAS <- filter(Include, ToInclude == "yes")$Acronym
IncludedPFAS <- unique(IncludedPFAS)
IncludedOECDCat <- filter(Include, ToInclude == "yes")$Category_OECD
IncludedOECDCat <- unique(IncludedOECDCat)

Datab <- Datab %>% filter(Acronym %in% IncludedPFAS)
Check <- unique(Datab$Acronym)
Check %in% IncludedPFAS 
IncludedPFAS <- sort(IncludedPFAS)

DataAv <- DataAv %>% pivot_longer(!c(PBK_mechanism, Parameter), 
          names_to = "Category_OECD",
          values_to = "Occurence")%>% 
  filter(Category_OECD %in% IncludedOECDCat)
Check <- unique(DataAv$Category_OECD)
Check %in% IncludedOECDCat 
DataAv <- DataAv %>% pivot_wider(names_from = Category_OECD, 
         values_from = Occurence)

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

MW <- Datab %>% filter(Parameter == "MW")
SelectedPFAS_Parms_df <- Datab %>% 
  filter(Parameter %in% Selected_params) 
write.csv(MW, "MW.csv")

RmPappTransport <- SelectedPFAS_Parms_df %>% 
  filter(Parameter != "Papp",
   PBK_mechanism != "Transport") %>% 
  select(!Parameter_specific)

CleanPapp <- SelectedPFAS_Parms_df %>% 
  filter(Parameter == "Papp") %>% 
  filter(Parameter_source != "in silico") %>% 
  mutate(Parameter_specific = str_replace_all(Parameter_specific,
             "NA", "in_vitro")) %>% 
  unite(Parameter, 
  Parameter,  Parameter_specific,
  sep = "__")

CleanTransport <- SelectedPFAS_Parms_df %>% 
  filter(PBK_mechanism == "Transport") %>% 
  unite(Parameter, 
  Parameter,  Parameter_specific,
  sep = "__")

Final_df <- rbind(RmPappTransport, CleanPapp, CleanTransport) 

nC_PFCAs <- unique(filter(Final_df, 
        Category_OECD =="PFCAs")$n_C) %>% sort()
nC_PFSAs <- unique(filter(Final_df, 
        Category_OECD =="PFSAs")$n_C) %>% sort()
nC_PFECAs <- unique(filter(Final_df, 
         Category_OECD =="PFECAs")$n_C) %>% sort()

df <- Final_df %>% select("Category_OECD", "Parameter", "n_C", "Value") %>% 
  unite(Var, 
  Category_OECD, Parameter, sep = "___") %>% 
  filter(!Var %in% c("PFCAs___Km__NTCP",
      "PFCAs___Km__OATP1B1",
      "PFCAs___Km__OATP1B3",
      "PFCAs___Km__OATP2B1",
      "PFCAs___Vmax__NTCP",
      "PFCAs___Vmax__OATP1B1",
      "PFCAs___Vmax__OATP1B3",
      "PFCAs___Vmax__OATP2B1",
      "PFSAs___Km__OAT1",
      "PFSAs___Vmax__OAT1"))
# write.csv(Final_df, here("NotOnGit","For_Regression_df.csv"), row.names = FALSE)

# Transport parameters ----

## Create dfs ----
Transport <- Final_df %>%
  filter(PBK_mechanism == "Transport") %>%
  select(Acronym, Category_OECD, Parameter, 
   Value, SD, n_C, Parameter_source, Parameter_source_system, Reference) %>%
  separate(Parameter, sep = "__", into = c("Parameter", "Delete")) %>% 
  select(!Delete)

Transport$SD <- as.numeric(Transport$SD)

Transport <- Transport %>%
  mutate(Title = case_when(
 Parameter == "ka_ASBT" ~ "Binding Energy (kcal/mol) ASBT iV",
 Parameter == "Ka_ASBT" ~ "Binding Energy (kcal/mol) ASBT iV",
 Parameter == "Kd_ASBT" ~ "ASBT RF score",
 
 Parameter == "freab_ASBT_NTCP" ~ "fraction Reabsorbed ASBT/NTCP",
 
 Parameter == "Vmax_OATP2B1" ~ "Vmax OATP2B1",
 Parameter == "Vmax_OATP1B1" ~ "Vmax OATP1B1",
 Parameter == "Vmax_OATP1B3" ~ "Vmax OATP1B3",
 Parameter == "Vmax_NTCP" ~ "Vmax NTCP",
 Parameter == "Km_OATP1B1" ~ "Km OATP1B1",
 Parameter == "Km_OATP1B3" ~ "Km OATP1B3",
 Parameter == "Km_OATP2B1" ~ "Km OATP2B1",
 Parameter == "Km_NTCP" ~ "Km NTCP",
 
 Parameter == "Vmax_OAT1" ~ "Vmax OAT1",
 Parameter == "Vmax_OAT3" ~ "Vmax OAT3",
 Parameter == "Vmax_OAT4" ~ "Vmax OAT4",
 Parameter == "Km_OAT1" ~ "Km OAT1",
 Parameter == "Km_OAT3" ~ "Km OAT3",
 Parameter == "Km_OAT4" ~ "Km OAT4",
 Parameter == "UR_OAT1" ~ "Uptake Rate OAT1",
 Parameter == "UR_OAT2" ~ "Uptake Rate OAT2",
 Parameter == "UR_OAT3" ~ "Uptake Rate OAT3",
 Parameter == "UR_OAT4" ~ "Uptake Rate OAT4",
 Parameter == "ka_OAT4" ~ "Binding energy OAT4")) 


Transport$Title <- factor(Transport$Title,
        levels = c("Uptake Rate OAT1", "Uptake Rate OAT2", "Uptake Rate OAT3", 
          "Uptake Rate OAT4", "Binding energy OAT4",
          "Vmax OAT1", "Vmax OAT3", "Vmax OAT4",
          "Km OAT1", "Km OAT3", "Km OAT4",
          
          "Vmax NTCP", "Vmax OATP2B1", "Vmax OATP1B1", "Vmax OATP1B3",
          "Km NTCP", "Km OATP2B1", "Km OATP1B1", "Km OATP1B3",
          
          "Binding Energy (kcal/mol) ASBT iV",
          "ASBT RF score",
          "fraction Reabsorbed ASBT/NTCP"
        ))
Transport$Category_OECD <- factor(Transport$Category_OECD,
          levels = c("PFCAs", "PFSAs", "PFECAs"))


# Transport_Parameters <- unique(Transport$Title) %>% print()
# # 
# Parameter_i <- NULL 
# Plot_i <- NULL
# Data_i <- NULL
# Plot_list <- list()
# 
# for (i in Transport_Parameters) {
#   
#   Parameter_i <- i
#   Data_i <- filter(Transport, Title == Parameter_i) 
#   
#   Plot_i <- ggplot(Data_i, 
#     aes(y = Value,
#      x = n_C, 
#      colour = Category_OECD,
#      shape = Category_OECD)
#   ) +
#  geom_errorbar(data = Data_i, 
#       aes(ymin = Value - SD, ymax = Value + SD,
#        colour = Category_OECD),
#       width = 0.25, alpha = 0.8) + 
#  scale_x_continuous(breaks=seq(1,13))+
#  labs(x = "nC",
#    y = "Parameter value") 
#   
#   Plot_i <- Plot_i + 
#  geom_point(size = 1.5, stroke = 1, alpha = 0.9) + 
#  scale_shape_manual(values = c("PFCAs" = 19,
#           "PFECAs"  = 15,
#           "PFSAs" = 18)) +
#  scale_colour_manual(values = c("PFCAs" = "#4ea9cf",
#         "PFECAs"  = "#55b483",
#         "PFSAs" = "#b45c73")) +
#  labs(title = paste0(i),
#    colour = "C. OECD",   
#    shape = "C. OECD") +
#  CP_theme +
#  theme(
#    legend.position = "none"
#  ) 
#   
#   Plot_list[[i]] <- Plot_i
#   
# }
# 
# for(i in 1:Transport_Parameters){
#   ggsave(plot = Plot_list[[i]], file = paste("Transport",i,".png",sep=""),
#    dpi = 600,
#    width = 9, 
#    height = 10, 
#    units = "cm")
# }
 Vmax <- Transport %>% filter(Parameter %in% c("Vmax_NTCP", "Vmax_OATP2B1","Vmax_OATP1B1", "Vmax_OATP1B3", "Vmax_OAT1", "Vmax_OAT3", "Vmax_OAT4"))
 Km <- Transport %>% filter(Parameter %in% c("Km_NTCP", "Km_OATP2B1", "Km_OATP1B1", "Km_OATP1B3", "Km_OAT1", "Km_OAT3", "Km_OAT4"))
 
 
 Vmax_All <- ggplot(Vmax,
       aes(y = Value,
        x = n_C,
        colour = Category_OECD, 
        shape = Parameter_source_system)) +
   geom_point(size = 1.5, stroke = 1, alpha = 0.9) +
   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD,
      colour = Category_OECD),
     width = 0.25, alpha = 0.8) + 
   scale_colour_manual(values = c("PFCAs" = "#4ea9cf",
          "PFECAs"  = "#55b483",
          "PFSAs" = "#b45c73")) +
   scale_x_continuous(breaks=seq(1,13)) +
   labs(x = "nC",
  y = "Vmax (nmol/min/mg protein)") +
   facet_grid(. ~ Title, scales = "free") +
   CP_theme +
   theme(legend.title = element_blank(),
   legend.position = "right") # c(0.9,0.9))
 
 ggsave(plot = Vmax_All,
  file = "Vmax_All.png",
  dpi = 600,
  width = 20,
  height = 5,
  units = "cm")
 
 Km <- Km %>% mutate(Value = ifelse(Value == 1.0000e-10, 1.0000e+10, Value))
 Km_All <- ggplot(Km,
      aes(y = 1/Value,
        x = n_C,
        colour = Category_OECD, 
        shape = Parameter_source_system)) +
   geom_point(size = 1.5, stroke = 1, alpha = 0.9) +
   geom_errorbar(aes(ymin = (1/Value) - (SD* (-1/Value**2)), 
                     ymax = (1/Value) + (SD* (-1/Value**2)), #ymin = Value - SD, ymax = Value + SD, #
      colour = Category_OECD),
     width = 0.25, alpha = 0.8) + 
   scale_colour_manual(values = c("PFCAs" = "#4ea9cf",
          "PFECAs"  = "#55b483",
          "PFSAs" = "#b45c73")) +
   scale_x_continuous(breaks=seq(1,13)) +
   labs(x = "nC",
  y = "Km-1 (μM-1)") +
   facet_grid(. ~ Title, scales = "free") +
   CP_theme +
   theme(legend.title = element_blank(),
   legend.position = "right") # c(0.9,0.9))
 
 ggsave(plot = Km_All,
  file = "Km_All.png",
  dpi = 600,
  width = 20,
  height = 5,
  units = "cm")
 
 
 TTransport <- Transport %>%
   filter(Parameter_source == "in vitro") %>% 
   select(Acronym, Category_OECD, Parameter, Parameter_source_system, Value, Reference, n_C, SD) %>% 
   unite(Parameter, 
   Parameter, Parameter_source_system,
   sep = "__") %>% 
   unite(Parameter, 
   Parameter, Reference,
   sep = "___") %>% 
   filter(! Parameter %in% c(
  "UR_OAT1__HEK293cells___S.Ryu et al. 2024",
  "UR_OAT2__HEK293cells___S.Ryu et al. 2024",
  "UR_OAT3__HEK293cells___S.Ryu et al. 2024",
  "UR_OAT4__HEK293cells___S.Ryu et al. 2024",
  "Ka_ASBT__ASBTinPBS___Y.Liang et al. 2025" 
   )) %>% 
   mutate(Parameter = str_replace_all(Parameter, " et al. 2024", ""),
    Parameter = str_replace_all(Parameter, " et al. 2023", ""),
    Parameter = str_replace_all(Parameter, " et al. 2017", ""),
    Parameter = str_replace_all(Parameter, " et al. 2015", ""),
    Parameter = str_replace_all(Parameter, " et al. 2007", ""))
 
 Transport_wide <- TTransport %>% select(!SD) %>% 
   pivot_wider(names_from = Parameter, values_from = Value)
 
 TransportSD_wide <- TTransport %>% select(!Value) %>% 
   mutate(Parameter = paste0("SD_", Parameter)) %>% 
   pivot_wider(names_from = Parameter, values_from = SD)
 
 Transport_SD_clearances <- cbind(Transport_wide, select(TransportSD_wide, !c(Acronym, Category_OECD, n_C)))
 
 Transport_clearances <- Transport_wide %>% 
   transmute(Acronym = Acronym, 
    Category_OECD = Category_OECD, 
    n_C = n_C, 
    CL_OAT1__HEK293cells___Louisse = (Vmax_OAT1__HEK293cells___Louisse*1e-3)/Km_OAT1__HEK293cells___Louisse,
    CL_OAT3__HEK293cells___Louisse = (Vmax_OAT3__HEK293cells___Louisse*1e-3)/Km_OAT3__HEK293cells___Louisse,
    CL_OAT4__HEK293cells___Louisse = (Vmax_OAT4__HEK293cells___Louisse*1e-3)/Km_OAT4__HEK293cells___Louisse,
    CL_OAT1__HEK293cells___Hatsuki = (Vmax_OAT1__HEK293cells___Hatsuki*1e-3)/Km_OAT1__HEK293cells___Hatsuki,
    CL_OAT3__HEK293cells___Hatsuki = (Vmax_OAT3__HEK293cells___Hatsuki*1e-3)/Km_OAT3__HEK293cells___Hatsuki,
    CL_OATP1B3__CHOcells___Zhao = (Vmax_OATP1B3__CHOcells___Zhao*1e-3)/Km_OATP1B3__CHOcells___Zhao, 
    CL_OATP1B1__CHOcells___Zhao = (Vmax_OATP1B1__CHOcells___Zhao*1e-3)/Km_OATP1B1__CHOcells___Zhao, 
    CL_OATP2B1__CHOcells___Zhao = (Vmax_OATP2B1__CHOcells___Zhao*1e-3)/Km_OATP2B1__CHOcells___Zhao,
    CL_NTCP__CHOcells___Zhao =  (Vmax_NTCP__CHOcells___Zhao*1e-3)/Km_NTCP__CHOcells___Zhao,
    CL_OATP1B3__CHOcells___Lin = (Vmax_OATP1B3__CHOcells___Lin*1e-3)/Km_OATP1B3__CHOcells___Lin,
    CL_OATP1B1__CHOcells___Lin = (Vmax_OATP1B1__CHOcells___Lin*1e-3)/Km_OATP1B1__CHOcells___Lin,
    CL_OATP2B1__CHOcells___Lin = (Vmax_OATP2B1__CHOcells___Lin*1e-3)/Km_OATP2B1__CHOcells___Lin,
    CL_NTCP__HEK293cells___Lin =  (Vmax_NTCP__HEK293cells___Lin*1e-3)/Km_NTCP__HEK293cells___Lin
   ) %>% 
   pivot_longer(CL_OAT1__HEK293cells___Louisse:CL_NTCP__HEK293cells___Lin,
    names_to = "Parameter", values_to = "Value") %>% 
   separate(Parameter, into = c("Parameter", "Reference"), sep = "___")  %>% 
   separate(Parameter, into = c("Parameter", "Cells"), sep = "__") %>% 
   separate(Parameter, into = c("Delete", "Parameter"), sep = "_") %>% 
   select(!Delete) %>% 
   mutate(Parameter = factor(Parameter, levels = c("OAT1", "OAT3", "OAT4", "NTCP", "OATP1B1", "OATP1B3", "OATP2B1")))
 
 Transport_clearances_SD <- Transport_SD_clearances %>% 
   transmute(Acronym = Acronym, 
    Category_OECD = Category_OECD, 
    n_C = n_C,
    SD_CL_OAT1__HEK293cells___Louisse = sqrt((((SD_Vmax_OAT1__HEK293cells___Louisse*1e-3)**2)*((1/Km_OAT1__HEK293cells___Louisse)**2)) + 
              ((SD_Km_OAT1__HEK293cells___Louisse**2)*(-(Vmax_OAT1__HEK293cells___Louisse*1e-3)/(Km_OAT1__HEK293cells___Louisse**2))**2)),
    SD_CL_OAT3__HEK293cells___Louisse = sqrt((((SD_Vmax_OAT3__HEK293cells___Louisse*1e-3)**2)*((1/Km_OAT3__HEK293cells___Louisse)**2)) + 
              ((SD_Km_OAT3__HEK293cells___Louisse**2)*(-(Vmax_OAT3__HEK293cells___Louisse*1e-3)/(Km_OAT3__HEK293cells___Louisse**2))**2)),
    SD_CL_OAT4__HEK293cells___Louisse = sqrt((((SD_Vmax_OAT4__HEK293cells___Louisse*1e-3)**2)*((1/Km_OAT4__HEK293cells___Louisse)**2)) + 
              ((SD_Km_OAT4__HEK293cells___Louisse**2)*(-(Vmax_OAT4__HEK293cells___Louisse*1e-3)/(Km_OAT4__HEK293cells___Louisse**2))**2)),
    SD_CL_OAT1__HEK293cells___Hatsuki = sqrt((((SD_Vmax_OAT1__HEK293cells___Hatsuki*1e-3)**2)*((1/Km_OAT1__HEK293cells___Hatsuki)**2)) + 
              ((SD_Km_OAT1__HEK293cells___Hatsuki**2)*(-(Vmax_OAT1__HEK293cells___Hatsuki*1e-3)/(Km_OAT1__HEK293cells___Hatsuki**2))**2)),
    SD_CL_OAT3__HEK293cells___Hatsuki = sqrt((((SD_Vmax_OAT3__HEK293cells___Hatsuki*1e-3)**2)*((1/Km_OAT3__HEK293cells___Hatsuki)**2)) + 
              ((SD_Km_OAT3__HEK293cells___Hatsuki**2)*(-(Vmax_OAT3__HEK293cells___Hatsuki*1e-3)/(Km_OAT3__HEK293cells___Hatsuki**2))**2)),
    SD_CL_OATP1B3__CHOcells___Zhao = sqrt((((SD_Vmax_OATP1B3__CHOcells___Zhao*1e-3)**2)*((1/Km_OATP1B3__CHOcells___Zhao)**2)) +
              ((SD_Km_OATP1B3__CHOcells___Zhao**2)*(-(Vmax_OATP1B3__CHOcells___Zhao*1e-3)/(Km_OATP1B3__CHOcells___Zhao**2))**2)), 
    SD_CL_OATP1B1__CHOcells___Zhao =  sqrt((((SD_Vmax_OATP1B1__CHOcells___Zhao*1e-3)**2)*((1/Km_OATP1B1__CHOcells___Zhao)**2)) + 
               ((SD_Km_OATP1B1__CHOcells___Zhao**2)*(-(Vmax_OATP1B1__CHOcells___Zhao*1e-3)/(Km_OATP1B1__CHOcells___Zhao**2))**2)), 
    SD_CL_OATP2B1__CHOcells___Zhao = sqrt((((SD_Vmax_OATP2B1__CHOcells___Zhao*1e-3)**2)*((1/Km_OATP2B1__CHOcells___Zhao)**2)) + 
              ((SD_Km_OATP2B1__CHOcells___Zhao**2)*(-(Vmax_OATP2B1__CHOcells___Zhao*1e-3)/(Km_OATP2B1__CHOcells___Zhao**2))**2)),
    SD_CL_NTCP__CHOcells___Zhao =  sqrt((((SD_Vmax_NTCP__CHOcells___Zhao*1e-3)**2)*((1/Km_NTCP__CHOcells___Zhao)**2)) + 
               ((SD_Km_NTCP__CHOcells___Zhao**2)*(-(Vmax_NTCP__CHOcells___Zhao*1e-3)/(Km_NTCP__CHOcells___Zhao**2))**2)),
    SD_CL_OATP1B3__CHOcells___Lin = sqrt((((SD_Vmax_OATP1B3__CHOcells___Lin*1e-3)**2)*((1/Km_OATP1B3__CHOcells___Lin)**2)) + 
             ((SD_Km_OATP1B3__CHOcells___Lin**2)*(-(Vmax_OATP1B3__CHOcells___Lin*1e-3)/(Km_OATP1B3__CHOcells___Lin**2))**2)),
    SD_CL_OATP1B1__CHOcells___Lin =  sqrt((((SD_Vmax_OATP1B1__CHOcells___Lin*1e-3)**2)*((1/Km_OATP1B1__CHOcells___Lin)**2)) + 
              ((SD_Km_OATP1B1__CHOcells___Lin**2)*(-(Vmax_OATP1B1__CHOcells___Lin*1e-3)/(Km_OATP1B1__CHOcells___Lin**2))**2)),
    SD_CL_OATP2B1__CHOcells___Lin = sqrt((((SD_Vmax_OATP2B1__CHOcells___Lin*1e-3)**2)*((1/Km_OATP2B1__CHOcells___Lin)**2)) + 
             ((SD_Km_OATP2B1__CHOcells___Lin**2)*(-(Vmax_OATP2B1__CHOcells___Lin*1e-3)/(Km_OATP2B1__CHOcells___Lin**2))**2)),
    SD_CL_NTCP__HEK293cells___Lin =  sqrt((((SD_Vmax_NTCP__HEK293cells___Lin*1e-3)**2)*((1/Km_NTCP__HEK293cells___Lin)**2)) + 
              ((SD_Km_NTCP__HEK293cells___Lin**2)*(-(Vmax_NTCP__HEK293cells___Lin*1e-3)/(Km_NTCP__HEK293cells___Lin**2))**2))
   ) %>% 
   pivot_longer(SD_CL_OAT1__HEK293cells___Louisse:SD_CL_NTCP__HEK293cells___Lin,
    names_to = "Parameter", values_to = "SD") %>% 
   separate(Parameter, into = c("Parameter", "Reference"), sep = "___")  %>% 
   separate(Parameter, into = c("Parameter", "Cells"), sep = "__") %>% 
   separate(Parameter, into = c("Delete", "Delete2","Parameter"), sep = "_") %>% 
   select(!c(Delete, Delete2)) %>% 
   mutate(Parameter = factor(Parameter, levels = c("NTCP", "OATP1B1", "OATP1B3", "OATP2B1","OAT1", "OAT3", "OAT4"))) 
 
 
 Transport2_clearances <- left_join(Transport_clearances, Transport_clearances_SD)
 Transport3_clearances <- Transport2_clearances %>% 
   mutate(Value = ifelse(Value == 1.000000e-03, 1.000000e-10, Value))
 
 Clearances_SD <-   ggplot(Transport3_clearances,
         aes(y = Value*1e3,
          x = n_C,
          colour = Category_OECD, 
          shape = Cells)) +
   geom_point(size = 1.5, stroke = 1, alpha = 0.9) +
   geom_errorbar(aes(ymin = Value*1e3 - SD*1e3, ymax = Value*1e3 + SD*1e3,
      colour = Category_OECD),
     width = 0.25, alpha = 0.8) + 
   scale_colour_manual(values = c("PFCAs" = "#4ea9cf",
          "PFECAs"  = "#55b483",
          "PFSAs" = "#b45c73")) +
   scale_y_continuous(breaks=seq(0,0.3, by = 0.05)) +
   scale_x_continuous(breaks=seq(1,13)) +
   labs(x = "nC",
  y = "Clearance (ml/min/mg protein)") +
   facet_grid(. ~ Parameter, scales = "free") +
   CP_theme +
   theme(legend.title = element_blank(),
   legend.position = "right") # c(0.9,0.9))
 
 ggsave(plot = Clearances_SD,
  file = "Clearances_SD.png",
  dpi = 600,
  width = 20,
  height = 9,
  units = "cm")
 
 ## Derive values ----
 
 ### OAT1 PFCAs ----
 
 Clearance <- Transport3_clearances %>% 
   filter(Category_OECD == "PFCAs") %>%
   filter(Parameter == "OAT1") %>% 
   mutate(Value = Value*1e3, 
          SD = SD*1e3) # unit change from L to ml
 
 OAT1_PFCAs <- Clearance %>%
   select(!SD) %>% na.omit()
 
 OAT1_df <- OAT1_PFCAs

 y = OAT1_df$Value
 x = OAT1_df$n_C
 a = 0
 data <- tibble(x = x, y = y, a=a)


 fit_nls <- lm(y ~ x, data = data)

 Results <- summary(fit_nls)

 Residuals_OAT1 <- tibble(
   x = x,
   y = residuals(fit_nls)
   )
 
  R <- ggplot(Residuals_OAT1) +
    geom_point(aes(x, y)) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks=seq(1,13))+
    labs(y = "Residuals", x = "n_C") +
    CP_theme

  predicted <- tibble(
    n_C = seq(4, 13, by = 0.01),
    predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_y = ifelse(predicted_y < 0, 0, predicted_y))
  

  predicted_OAT1_PFCAs <- tibble(
    n_C = seq(4, 13, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  

  observed_OAT1_PFCAs <- OAT1_df %>% arrange(n_C)

  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)

  D <- ggplot() +
    geom_point(data = filter(Clearance, Category_OECD == "PFCAs"), 
               aes(y = Value,  x = n_C),
    size = 1.5, stroke = 1, alpha = 0.9, colour =  "#4ea9cf") +
    geom_errorbar(data = filter(Clearance, Category_OECD == "PFCAs"),
      aes(x = n_C, ymin = Value - SD, ymax = Value + SD),
      width = 0.25, alpha = 0.8, colour =  "#4ea9cf") +
    geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
    geom_point(data = predicted_OAT1_PFCAs, aes(n_C, predicted_value), shape = 19) +
    scale_x_continuous(breaks=seq(1,13)) +
    scale_y_continuous(breaks=seq(0,0.06, by = 0.5e-02)) +
    annotate("text", x=3.5, y=0.06, label=paste0("y = ", Formula), size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.0585, label=paste0("Intercept = ", Intercept),
             size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.0570, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
             size = 10, hjust = 0) +
    labs(y = "Clearance OAT1 ml/min/mg protein")   +
    CP_theme +
    theme(legend.title = element_blank(),
    legend.position = "right")

  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("OAT1_PFCAs", ".png",sep=""),
   dpi = 600,
   width = 9,
   height = 12,
   units = "cm")

 ### OAT1 PFSAs ----
  observed_OAT1_PFSAs <- Transport3_clearances %>% 
    filter(Category_OECD == "PFSAs") %>%
    filter(Parameter == "OAT1") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) %>%# unit change from L to ml
    select(!SD) %>% 
    na.omit() 
  predicted_OAT1_PFSAs <- predicted_OAT1_PFCAs %>% filter(n_C <= 8)
  

 ### OAT3 PFCAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFCAs") %>%
    filter(Parameter == "OAT3") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OAT3_PFCAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  OAT3_df <- OAT3_PFCAs
  
  y = OAT3_df$Value
  x = OAT3_df$n_C
  data <- tibble(x = x, y = y)
  
  
  fit_nls <- lm(y ~ x, data = data)
  
  Results <- summary(fit_nls)
  
  Residuals_OAT3 <- tibble(
    x = x,
    y = residuals(fit_nls)
  )
  
  R <- ggplot(Residuals_OAT3) +
    geom_point(aes(x, y)) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks=seq(1,13))+
    labs(y = "Residuals", x = "n_C") +
    CP_theme
  
  predicted <- tibble(
    n_C = seq(4, 13, by = 0.01),
    predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_y = ifelse(predicted_y < 0, 0, predicted_y))
  
  predicted_OAT3_PFCAs <- tibble(
    n_C = seq(4, 13, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  observed_OAT3_PFCAs <- OAT3_df %>% arrange(n_C)
  
  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)
  
  D <- ggplot() +
    geom_point(data = filter(Clearance, Category_OECD == "PFCAs"), 
               aes(y = Value,  x = n_C),
               size = 1.5, stroke = 1, alpha = 0.9, colour =  "#4ea9cf") +
    geom_errorbar(data = filter(Clearance, Category_OECD == "PFCAs"),
                  aes(x = n_C, ymin = Value - SD, ymax = Value + SD),
                  width = 0.25, alpha = 0.8, colour =  "#4ea9cf") +
    geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
    geom_point(data = predicted_OAT3_PFCAs, aes(n_C, predicted_value), shape = 19) +
    scale_x_continuous(breaks=seq(1,13)) +
    scale_y_continuous(breaks=seq(0,0.05, by = 0.005)) +
    annotate("text", x=3.5, y=0.05, label=paste0("y = ", Formula), size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.0485, label=paste0("Intercept = ", Intercept),
             size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.0470, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
             size = 10, hjust = 0) +
    labs(y = "Clearance OAT3 ml/min/mg protein") +
    CP_theme +
    theme(legend.title = element_blank(),
          legend.position = "right")
  
  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("OAT3_PFCAs", ".png",sep=""),
         dpi = 600,
         width = 9,
         height = 12,
         units = "cm")
  
  

 ### OAT3 PFSAs ----
 
  observed_OAT3_PFSAs <- Transport3_clearances %>% 
    filter(Category_OECD == "PFSAs") %>%
    filter(Parameter == "OAT3") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) %>%# unit change from L to ml
    select(!SD) %>% 
    na.omit() 
  predicted_OAT3_PFSAs <- predicted_OAT3_PFCAs %>% filter(n_C <= 8)
  
  

  ### OAT4 PFCAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFCAs") %>%
    filter(Parameter == "OAT4") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OAT4_PFCAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  OAT4_df <- OAT4_PFCAs
  
  y = OAT4_df$Value
  x = OAT4_df$n_C
  data <- tibble(x = x, y = y)
  
  
  fit_nls <- lm(y ~ x, data = data)
  
  Results <- summary(fit_nls)
  
  Residuals_OAT4 <- tibble(
    x = x,
    y = residuals(fit_nls)
  )
  
  R <- ggplot(Residuals_OAT4) +
    geom_point(aes(x, y)) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks=seq(1,13))+
    labs(y = "Residuals", x = "n_C") +
    CP_theme
  
  predicted <- tibble(
    n_C = seq(4, 13, by = 0.01),
    predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_y = ifelse(predicted_y < 0, 0, predicted_y))
  
  predicted_OAT4_PFCAs <- tibble(
    n_C = seq(4, 13, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  observed_OAT4_PFCAs <- OAT4_df %>% arrange(n_C)
  
  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)
  
  D <- ggplot() +
    geom_point(data = filter(Clearance, Category_OECD == "PFCAs"), 
               aes(y = Value,  x = n_C),
               size = 1.5, stroke = 1, alpha = 0.9, colour =  "#4ea9cf") +
    geom_errorbar(data = filter(Clearance, Category_OECD == "PFCAs"),
                  aes(x = n_C, ymin = Value - SD, ymax = Value + SD),
                  width = 0.25, alpha = 0.8, colour =  "#4ea9cf") +
    geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
    geom_point(data = predicted_OAT4_PFCAs, aes(n_C, predicted_value), shape = 19) +
    scale_x_continuous(breaks=seq(1,13)) +
    scale_y_continuous(breaks=seq(0,0.25, by = 0.05)) +
    annotate("text", x=3.5, y=0.25, label=paste0("y = ", Formula), size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.24, label=paste0("Intercept = ", Intercept),
             size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.23, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
             size = 10, hjust = 0) +
    labs(y = "Clearance OAT4 ml/min/mg protein") +
    CP_theme +
    theme(legend.title = element_blank(),
          legend.position = "right")
  
  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("OAT4_PFCAs", ".png",sep=""),
         dpi = 600,
         width = 9,
         height = 12,
         units = "cm")
  
  
  
  ### OAT4 PFSAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFSAs") %>%
    filter(Parameter == "OAT4") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OAT4_PFSAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  OAT4_df <- OAT4_PFSAs
  
  y = OAT4_df$Value
  x = OAT4_df$n_C
  data <- tibble(x = x, y = y)
  
  
  fit_nls <- lm(y ~ x, data = data)
  
  Results <- summary(fit_nls)
  
  Residuals_OAT4 <- tibble(
    x = x,
    y = residuals(fit_nls)
  )
  
  R <- ggplot(Residuals_OAT4) +
    geom_point(aes(x, y)) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks=seq(1,8))+
    labs(y = "Residuals", x = "n_C") +
    CP_theme
  
  predicted <- tibble(
    n_C = seq(4, 8, by = 0.01),
    predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_y = ifelse(predicted_y < 0, 0, predicted_y))
  
  predicted_OAT4_PFSAs <- tibble(
    n_C = seq(4, 8, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  observed_OAT4_PFSAs <- OAT4_df %>% arrange(n_C)
  
  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)
  
  D <-  ggplot() +
    geom_point(data = filter(Clearance, Category_OECD == "PFSAs"), 
               aes(y = Value,  x = n_C),
               size = 1.5, stroke = 1, alpha = 0.9, colour =  "#b45c73") +
    geom_errorbar(data = filter(Clearance, Category_OECD == "PFSAs"),
                  aes(x = n_C, ymin = Value - SD, ymax = Value + SD),
                  width = 0.25, alpha = 0.8, colour =  "#b45c73") +
    geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
    geom_point(data = predicted_OAT4_PFSAs, aes(n_C, predicted_value), shape = 19) +
    scale_x_continuous(breaks=seq(1,8)) +
    scale_y_continuous(breaks=seq(0,0.11, by = 0.01)) +
    annotate("text", x=3.5, y=0.11, label=paste0("y = ", Formula), size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.104, label=paste0("Intercept = ", Intercept),
             size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.098, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
             size = 10, hjust = 0) +
    labs(y = "Clearance OAT4 ml/min/mg protein")  +
    CP_theme +
    theme(legend.title = element_blank(),
          legend.position = "right")
  
  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("OAT4_PFSAs", ".png",sep=""),
         dpi = 600,
         width = 9,
         height = 12,
         units = "cm")
  
  predicted_OAT4_PFSAs_asPFCAs <- predicted_OAT4_PFCAs %>% filter(n_C <= 8)
  
  
  
  ### NTCP PFSAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFSAs") %>%
    filter(Parameter == "NTCP") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  NTCP_PFSAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  NTCP_df <- NTCP_PFSAs
  
  y = NTCP_df$Value
  x = NTCP_df$n_C
  data <- tibble(x = x, y = y)
  
  
  fit_nls <- lm(y ~ x, data = data)
  
  Results <- summary(fit_nls)
  
  Residuals_NTCP <- tibble(
    x = x,
    y = residuals(fit_nls)
  )
  
  R <- ggplot(Residuals_NTCP) +
    geom_point(aes(x, y)) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks=seq(1,8))+
    labs(y = "Residuals", x = "n_C") +
    CP_theme
  
  predicted <- tibble(
    n_C = seq(4, 8, by = 0.01),
    predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_y = ifelse(predicted_y < 0, 0, predicted_y))
  
  predicted_NTCP_PFSAs <- tibble(
    n_C = seq(4, 8, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  observed_NTCP_PFSAs <- NTCP_df %>% arrange(n_C)
  
  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)
  
  D <- ggplot() +
    geom_point(data = filter(Clearance, Category_OECD == "PFSAs"), 
               aes(y = Value,  x = n_C),
               size = 1.5, stroke = 1, alpha = 0.9, colour =  "#b45c73") +
    geom_errorbar(data = filter(Clearance, Category_OECD == "PFSAs"),
                  aes(x = n_C, ymin = Value - SD, ymax = Value + SD),
                  width = 0.25, alpha = 0.8, colour =  "#b45c73") +
    geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
    geom_point(data = predicted_NTCP_PFSAs, aes(n_C, predicted_value), shape = 19) +
    scale_x_continuous(breaks=seq(1,8)) +
    scale_y_continuous(breaks=seq(0,0.3, by = 0.05)) +
    annotate("text", x=3.5, y=0.3, label=paste0("y = ", Formula), size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.29, label=paste0("Intercept = ", Intercept),
             size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.28, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
             size = 10, hjust = 0) +
    labs(y = "Clearance NTCP ml/min/mg protein")  +
    CP_theme +
    theme(legend.title = element_blank(),
          legend.position = "right")
  
  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("NTCP_PFSAs", ".png",sep=""),
         dpi = 600,
         width = 9,
         height = 12,
         units = "cm")
  
  ### NTCP PFCAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFCAs") %>%
    filter(Parameter == "NTCP") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  NTCP_PFCAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  observed_NTCP_PFCAs <- NTCP_PFCAs %>% arrange(n_C)
  
  predicted_NTCP_PFCAs <- tibble(
    n_C = seq(4, 13, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  
  ### OATP1B1 PFSAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFSAs") %>%
    filter(Parameter == "OATP1B1") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OATP1B1_PFSAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  OATP1B1_df <- OATP1B1_PFSAs
  
  y = OATP1B1_df$Value
  x = OATP1B1_df$n_C
  data <- tibble(x = x, y = y)
  
  
  fit_nls <- lm(y ~ x, data = data)
  
  Results <- summary(fit_nls)
  
  Residuals_OATP1B1 <- tibble(
    x = x,
    y = residuals(fit_nls)
  )
  
  R <- ggplot(Residuals_OATP1B1) +
    geom_point(aes(x, y)) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks=seq(1,8))+
    labs(y = "Residuals", x = "n_C") +
    CP_theme
  
  predicted <- tibble(
    n_C = seq(4, 8, by = 0.01),
    predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_y = ifelse(predicted_y < 0, 0, predicted_y))
  
  predicted_OATP1B1_PFSAs <- tibble(
    n_C = seq(4, 8, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  observed_OATP1B1_PFSAs <- OATP1B1_df %>% arrange(n_C)
  
  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)
  
  D <- ggplot() +
    geom_point(data = filter(Clearance, Category_OECD == "PFSAs"), 
               aes(y = Value,  x = n_C),
               size = 1.5, stroke = 1, alpha = 0.9, colour =  "#b45c73") +
    geom_errorbar(data = filter(Clearance, Category_OECD == "PFSAs"),
                  aes(x = n_C, ymin = Value - SD, ymax = Value + SD),
                  width = 0.25, alpha = 0.8, colour =  "#b45c73") +
    geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
    geom_point(data = predicted_OATP1B1_PFSAs, aes(n_C, predicted_value), shape = 19) +
    scale_x_continuous(breaks=seq(1,8)) +
    scale_y_continuous(breaks=seq(0,0.045, by = 0.005)) +
    annotate("text", x=3.5, y=0.045, label=paste0("y = ", Formula), size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.043, label=paste0("Intercept = ", Intercept),
             size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.041, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
             size = 10, hjust = 0) +
    labs(y = "Clearance OATP1B1 ml/min/mg protein")  +
    CP_theme +
    theme(legend.title = element_blank(),
          legend.position = "right")
  
  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("OATP1B1_PFSAs", ".png",sep=""),
         dpi = 600,
         width = 9,
         height = 12,
         units = "cm")
  
  ### OATP1B1 PFCAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFCAs") %>%
    filter(Parameter == "OATP1B1") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OATP1B1_PFCAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  observed_OATP1B1_PFCAs <- OATP1B1_PFCAs %>% arrange(n_C)
  
  predicted_OATP1B1_PFCAs <- tibble(
    n_C = seq(4, 13, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  
  
  
  
  ### OATP1B3 PFSAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFSAs") %>%
    filter(Parameter == "OATP1B3") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OATP1B3_PFSAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  OATP1B3_df <- OATP1B3_PFSAs
  
  y = OATP1B3_df$Value
  x = OATP1B3_df$n_C
  data <- tibble(x = x, y = y)
  
  
  fit_nls <- lm(y ~ x, data = data)
  
  Results <- summary(fit_nls)
  
  Residuals_OATP1B3 <- tibble(
    x = x,
    y = residuals(fit_nls)
  )
  
  R <- ggplot(Residuals_OATP1B3) +
    geom_point(aes(x, y)) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks=seq(1,8))+
    labs(y = "Residuals", x = "n_C") +
    CP_theme
  
  predicted <- tibble(
    n_C = seq(4, 8, by = 0.01),
    predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_y = ifelse(predicted_y < 0, 0, predicted_y))
  
  predicted_OATP1B3_PFSAs <- tibble(
    n_C = seq(4, 8, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  observed_OATP1B3_PFSAs <- OATP1B3_df %>% arrange(n_C)
  
  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)
  
  D <- ggplot() +
    geom_point(data = filter(Clearance, Category_OECD == "PFSAs"), 
               aes(y = Value,  x = n_C),
               size = 1.5, stroke = 1, alpha = 0.9, colour =  "#b45c73") +
    geom_errorbar(data = filter(Clearance, Category_OECD == "PFSAs"),
                  aes(x = n_C, ymin = Value - SD, ymax = Value + SD),
                  width = 0.25, alpha = 0.8, colour =  "#b45c73") +
    geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
    geom_point(data = predicted_OATP1B3_PFSAs, aes(n_C, predicted_value), shape = 19) + scale_y_continuous(breaks=seq(0,0.045, by = 0.005)) +
    annotate("text", x=3.5, y=0.055, label=paste0("y = ", Formula), size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.053, label=paste0("Intercept = ", Intercept),
             size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.051, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
             size = 10, hjust = 0) +
    labs(y = "Clearance OATP1B3 ml/min/mg protein")  +
    CP_theme +
    theme(legend.title = element_blank(),
          legend.position = "right")
  
  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("OATP1B3_PFSAs", ".png",sep=""),
         dpi = 600,
         width = 9,
         height = 12,
         units = "cm")
  
  ### OATP1B3 PFCAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFCAs") %>%
    filter(Parameter == "OATP1B3") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OATP1B3_PFCAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  observed_OATP1B3_PFCAs <- OATP1B3_PFCAs %>% arrange(n_C)
  
  predicted_OATP1B3_PFCAs <- tibble(
    n_C = seq(4, 13, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  
  
  
  
  
  ### OATP2B1 PFSAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFSAs") %>%
    filter(Parameter == "OATP2B1") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OATP2B1_PFSAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  OATP2B1_df <- OATP2B1_PFSAs
  
  y = OATP2B1_df$Value
  x = OATP2B1_df$n_C
  data <- tibble(x = x, y = y)
  
  
  fit_nls <- lm(y ~ x, data = data)
  
  Results <- summary(fit_nls)
  
  Residuals_OATP2B1 <- tibble(
    x = x,
    y = residuals(fit_nls)
  )
  
  R <- ggplot(Residuals_OATP2B1) +
    geom_point(aes(x, y)) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks=seq(1,8))+
    labs(y = "Residuals", x = "n_C") +
    CP_theme
  
  predicted <- tibble(
    n_C = seq(4, 8, by = 0.01),
    predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
  ) %>% arrange(n_C) %>% 
    mutate(predicted_y = ifelse(predicted_y < 0, 0, predicted_y))
  
  predicted_OATP2B1_PFSAs <- tibble(
    n_C = seq(4, 8, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  observed_OATP2B1_PFSAs <- OATP2B1_df %>% arrange(n_C)
  
  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)
  
  D <- ggplot() +
    geom_point(data = filter(Clearance, Category_OECD == "PFSAs"), 
               aes(y = Value,  x = n_C),
               size = 1.5, stroke = 1, alpha = 0.9, colour =  "#b45c73") +
    geom_errorbar(data = filter(Clearance, Category_OECD == "PFSAs"),
                  aes(x = n_C, ymin = Value - SD, ymax = Value + SD),
                  width = 0.25, alpha = 0.8, colour =  "#b45c73") +
    geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
    geom_point(data = predicted_OATP2B1_PFSAs, aes(n_C, predicted_value), shape = 19) +
    scale_x_continuous(breaks=seq(1,8)) +
    scale_y_continuous(breaks=seq(0,0.15, by = 0.05)) +
    annotate("text", x=3.5, y=0.15, label=paste0("y = ", Formula), size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.145, label=paste0("Intercept = ", Intercept),
             size = 10, hjust = 0) +
    annotate("text", x=3.5, y=0.140, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
             size = 10, hjust = 0) +
    labs(y = "Clearance OATP2B1 ml/min/mg protein")  +
    CP_theme +
    theme(legend.title = element_blank(),
          legend.position = "right")
  
  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("OATP2B1_PFSAs", ".png",sep=""),
         dpi = 600,
         width = 9,
         height = 12,
         units = "cm")
  
  ### OATP2B1 PFCAs ----
  
  Clearance <- Transport3_clearances %>% 
    filter(Category_OECD == "PFCAs") %>%
    filter(Parameter == "OATP2B1") %>% 
    mutate(Value = Value*1e3, 
           SD = SD*1e3) # unit change from L to ml
  
  OATP2B1_PFCAs <- Clearance %>%
    select(!SD) %>% na.omit()
  
  observed_OATP2B1_PFCAs <- OATP2B1_PFCAs %>% arrange(n_C)
  
  predicted_OATP2B1_PFCAs <- tibble(
    n_C = seq(4, 13, by = 1),
    predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
  ) %>% arrange(n_C)  %>% 
    mutate(predicted_value = ifelse(predicted_value < 0, 0, predicted_value))
  
  
  ## Save Results ----

  Observed <- rbind(observed_OAT1_PFSAs,
                    observed_OAT1_PFCAs,
                    observed_OAT3_PFSAs,
                    observed_OAT3_PFCAs,
                    observed_OAT4_PFSAs,
                    observed_OAT4_PFCAs,
                    observed_NTCP_PFSAs,
                    observed_NTCP_PFCAs,
                    observed_OATP1B1_PFSAs
        ) %>% mutate(Source = "Observed") 

  predicted_OAT1_PFSAs <- predicted_OAT1_PFSAs %>%
  mutate(Parameter = "OAT1", Category_OECD = "PFSAs")
  predicted_OAT1_PFCAs <- predicted_OAT1_PFCAs %>%
    mutate(Parameter = "OAT1", Category_OECD = "PFCAs")
  predicted_OAT3_PFSAs <- predicted_OAT3_PFSAs %>%
    mutate(Parameter = "OAT3", Category_OECD = "PFSAs")
  predicted_OAT3_PFCAs <- predicted_OAT3_PFCAs %>%
    mutate(Parameter = "OAT3", Category_OECD = "PFCAs")
  predicted_OAT4_PFSAs <- predicted_OAT4_PFSAs %>%
    mutate(Parameter = "OAT4", Category_OECD = "PFSAs")
  predicted_OAT4_PFCAs <- predicted_OAT4_PFCAs %>%
    mutate(Parameter = "OAT4", Category_OECD = "PFCAs")
  predicted_NTCP_PFSAs <- predicted_NTCP_PFSAs %>%
    mutate(Parameter = "NTCP", Category_OECD = "PFSAs")
  predicted_NTCP_PFCAs <- predicted_NTCP_PFCAs %>%
    mutate(Parameter = "NTCP", Category_OECD = "PFCAs")
  predicted_OATP1B1_PFSAs <- predicted_OATP1B1_PFSAs %>%
    mutate(Parameter = "OATP1B1", Category_OECD = "PFSAs")
  predicted_OATP1B1_PFCAs <- predicted_OATP1B1_PFCAs %>%
    mutate(Parameter = "OATP1B1", Category_OECD = "PFCAs")
  predicted_OATP1B3_PFSAs <- predicted_OATP1B3_PFSAs %>%
    mutate(Parameter = "OATP1B3", Category_OECD = "PFSAs")
  predicted_OATP1B3_PFCAs <- predicted_OATP1B3_PFCAs %>%
    mutate(Parameter = "OATP1B3", Category_OECD = "PFCAs")
  predicted_OATP2B1_PFSAs <- predicted_OATP2B1_PFSAs %>%
    mutate(Parameter = "OATP2B1", Category_OECD = "PFSAs")
  predicted_OATP2B1_PFCAs <- predicted_OATP2B1_PFCAs %>%
    mutate(Parameter = "OATP2B1", Category_OECD = "PFCAs")
  

  Predicted <- rbind(  predicted_OAT1_PFSAs,
                       predicted_OAT1_PFCAs,
                       predicted_OAT3_PFSAs,
                       predicted_OAT3_PFCAs,
                       predicted_OAT4_PFSAs,
                       predicted_OAT4_PFCAs,
                       predicted_NTCP_PFSAs,
                       predicted_NTCP_PFCAs,
                       predicted_OATP1B1_PFSAs,
                       predicted_OATP1B1_PFCAs,
                       predicted_OATP1B3_PFSAs,
                       predicted_OATP1B3_PFCAs,
                       predicted_OATP2B1_PFSAs,
                       predicted_OATP2B1_PFCAs) %>%
    rename("Value" = "predicted_value") %>%
    mutate(Source = "Predicted",
     Acronym = case_when(
    n_C == 4  & Category_OECD == "PFCAs" ~ "PFBA" ,
    n_C == 5  & Category_OECD == "PFCAs" ~ "PFPeA",
    n_C == 6  & Category_OECD == "PFCAs" ~ "PFHxA",
    n_C == 7  & Category_OECD == "PFCAs" ~ "PFHpA",
    n_C == 8  & Category_OECD == "PFCAs" ~ "PFOA",
    n_C == 9  & Category_OECD == "PFCAs" ~ "PFNA",
    n_C == 10  & Category_OECD == "PFCAs" ~ "PFDA",
    n_C == 11  & Category_OECD == "PFCAs" ~ "PFUnDA",
    n_C == 12  & Category_OECD == "PFCAs" ~ "PFDoDA",
    n_C == 13  & Category_OECD == "PFCAs" ~ "PFTrDA",
    n_C == 4  & Category_OECD == "PFSAs" ~ "PFBS" ,
    n_C == 5  & Category_OECD == "PFSAs" ~ "PFPeS",
    n_C == 6  & Category_OECD == "PFSAs" ~ "PFHxS",
    n_C == 7  & Category_OECD == "PFSAs" ~ "PFHpS",
    n_C == 8  & Category_OECD == "PFSAs" ~ "PFOS"),
    Cells = NA,
    Reference = NA)

  ClearanceData <- rbind(Observed, Predicted)
  write.csv(ClearanceData, "ClearanceData.csv", row.names = FALSE)

  
# Absorption parameters ----

## Create dfs ----
Absorption <- Final_df %>%
  filter(PBK_mechanism == "Absorption",
   Parameter != c("Papp__ionic")) %>%
  select(Acronym, Category_OECD, Parameter, Parameter_source_system,
   Value, SD, n_C) %>%
  unite(Parameter,
  Parameter,  Parameter_source_system,
  sep = "__")
Absorption$SD <- as.numeric(Absorption$SD)

Absorption <- Absorption %>%
  mutate(Title = case_when(
 Parameter == "Papp__6.5__MDCKII-LEcells" ~ "pH 6.5, MDCKII-LE cells",
 Parameter == "Papp__7.4__MDCKII-LEcells" ~ "pH 7.4, MDCKII-LE cells",
 Parameter == "Papp__in_vitro__hiPSC-derivedIEC" ~ "pH 7.4, hiPSC-derivedIEC",
 Parameter == "Papp__in_vitro__EpiIntestinal" ~ "pH 7.4, EpiIntestinal",
 Parameter == "Papp__in_vitro__Caco-2" ~ "pH 7.4, Caco-2"
  ))

Absorption$Title <- factor(Absorption$Title,
         levels = c("pH 6.5, MDCKII-LE cells",
           "pH 7.4, MDCKII-LE cells",
           "pH 7.4, hiPSC-derivedIEC",
           "pH 7.4, Caco-2",
           "pH 7.4, EpiIntestinal"
           ))
Absorption$Category_OECD <- factor(Absorption$Category_OECD,
        levels = c("PFCAs", "PFSAs", "PFECAs"))


  Papp <- ggplot(Absorption,
     aes(y = Value,
      x = n_C,
      colour = Title)) +
 geom_point(size = 1.5, stroke = 1, alpha = 0.9) +
 geom_errorbar(
   aes(ymin = Value - SD, ymax = Value + SD),
    width = 0.25, alpha = 0.8) +
 scale_x_continuous(breaks=seq(1,13)) +
 labs(x = "nC",
   y = "Papp cm/s") +
 facet_grid(. ~ Category_OECD, scales = "free") +
 CP_theme +
 theme(legend.title = element_blank(),
    legend.position = "right") # c(0.9,0.9))

  ggsave(plot = Papp,
   file = "Papp.png",
   dpi = 600,
   width = 12,
   height = 10,
   units = "cm")


## Derive values ----
### Papp ----
#### PFCAs ----
  Papp__7.4__MDCKII_nls <- Absorption %>%
 filter(Category_OECD == "PFCAs")
  # ,
  # Parameter %in% c("Papp__7.4__MDCKII-LEcells", "Papp__6.5__MDCKII-LEcells"))

  Papp_df <- Papp__7.4__MDCKII_nls[-20, ]  #-14

  y = Papp_df$Value*10e5
  x = Papp_df$n_C
  data <- tibble(x = x, y = y)

  fit_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
     data = data,
     start = list(Const = 1.03, Asym = 14, xmid = 9.5, scal = 0.42))

  Results <- summary(fit_nls)

  Residuals_Papp <- tibble(
 x = x,
 y = residuals(fit_nls)
  )

  R <- ggplot(Residuals_Papp) +
 geom_point(aes(x, y)) +
 geom_hline(yintercept = 0) +
 scale_x_continuous(breaks=seq(1,13))+
 labs(y = "Residuals", x = "n_C") +
 CP_theme

  predicted <- tibble(
 n_C = seq(4, 13, by = 0.01),
 predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
  ) %>% arrange(n_C)

  predicted_Papp_PFCAs <- tibble(
 n_C = seq(4, 13, by = 1),
 predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
  ) %>% arrange(n_C)

  observed_Papp_PFCAs <- Papp_df %>% arrange(n_C)
  # write.csv(observed_Papp_PFCAs, "Papp_PFCAs.csv_Observed.csv", row.names = FALSE)
  # write.csv(predicted_Papp_PFCAs, "Papp_PFCAs.csv", row.names = FALSE)


  Formula <- Results$formula[3]
  RSE <- round(Results$sigma, 3)
  Const <- round(Results$coefficients[1], 3)
  Asym <- round(Results$coefficients[2], 3)
  xmid <- round(Results$coefficients[3], 3)
  scal <- round(Results$coefficients[4], 3)

  D <- ggplot() +
 geom_point(data = filter(Absorption, Category_OECD == "PFCAs"), aes(y = Value*10e5,  x = n_C),
 size = 1.5, stroke = 1, alpha = 0.9, colour =  "#4ea9cf") +
 geom_errorbar(data = filter(Absorption, Category_OECD == "PFCAs"),
      aes(x = n_C, ymin = Value*10e5 - SD*10e5, ymax = Value*10e5 + SD*10e5),
      width = 0.25, alpha = 0.8, colour =  "#4ea9cf") +
 geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
 geom_point(data = predicted_Papp_PFCAs, aes(n_C, predicted_value), shape = 19) +
 scale_x_continuous(breaks=seq(1,13)) +
 annotate("text", x=3.5, y=16, label=paste0("y = ", Formula), size = 10, hjust = 0) +
 annotate("text", x=3.5, y=15.7, label=paste0("Const = ", Const, ", Asym = ", Asym,
            ", xmid = ", xmid, ", scal = ", scal),
    size = 10, hjust = 0) +
 annotate("text", x=3.5, y=15.4, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
 labs(y = "Papp 10e-5 cm/s") +
 CP_theme +
 theme(legend.title = element_blank(),
    legend.position = "right")

  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("Papp_nn_PFCAs", ".png",sep=""),
   dpi = 600,
   width = 9,
   height = 12,
   units = "cm")


#### PFSAs ----
  Papp__7.4__MDCKII_nls <- Absorption %>%
 filter(Category_OECD == "PFSAs") %>%
 filter(Parameter != "Papp__6.5__MDCKII-LEcells",
  Parameter != "Papp__in_vitro__Caco-2",
  Parameter != "Papp__in_vitro__hiPSC-derivedIEC")

  Papp_df <- Papp__7.4__MDCKII_nls

  y = Papp_df$Value*10e5
  x = Papp_df$n_C
  data <- tibble(x = x, y = y)

  fit_nls <- lm(y ~ x, data = data)

  Results <- summary(fit_nls)

  Residuals_Papp <- tibble(
 x = x,
 y = residuals(fit_nls)
  )

  R <- ggplot(Residuals_Papp) +
 geom_point(aes(x, y)) +
 geom_hline(yintercept = 0) +
 scale_x_continuous(breaks=seq(1,13))+
 labs(y = "Residuals", x = "n_C") +
 CP_theme

  predicted <- tibble(
 n_C = seq(4, 8, by = 0.01),
 predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
  ) %>% arrange(n_C)

  predicted_Papp_PFSAs <- tibble(
 n_C = seq(4, 8, by = 1),
 predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
  ) %>% arrange(n_C)

  observed_Papp_PFSAs <- Papp_df %>% arrange(n_C)

  Formula <- "Intercept * X"
  RSE <- round(Results$sigma, 3)
  R.squared <- round(Results$r.squared, 3)
  Intercept <- round(Results$coefficients[1], 3)

  D <- ggplot() +
 geom_point(data = filter(Absorption, Category_OECD == "PFSAs"), aes(y = Value*10e5,  x = n_C),
 size = 1.5, stroke = 1, alpha = 0.9, colour =  "#b45c73") +
 geom_errorbar(data = filter(Absorption, Category_OECD == "PFSAs"),
      aes(x = n_C, ymin = Value*10e5 - SD*10e5, ymax = Value*10e5 + SD*10e5),
      width = 0.25, alpha = 0.8, colour =  "#b45c73") +
 geom_path(data = predicted, aes(x = n_C, y = predicted_y)) +
 geom_point(data = predicted_Papp_PFSAs, aes(n_C, predicted_value), shape = 19) +
 scale_x_continuous(breaks=seq(1,8)) +
   annotate("text", x=3.5, y=16, label=paste0("y = ", Formula), size = 10, hjust = 0) +
   annotate("text", x=3.5, y=15.7, label=paste0("Intercept = ", Intercept),
   size = 10, hjust = 0) +
   annotate("text", x=3.5, y=15.4, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared),
   size = 10, hjust = 0) +
 labs(y = "Papp 10e-5 cm/s") +
 CP_theme +
 theme(legend.title = element_blank(),
    legend.position = "right")

  SavePlots <- D / R + plot_layout(heights = c(5,1))
  ggsave(plot = SavePlots, file = paste("Papp_n_PFSAs", ".png",sep=""),
   dpi = 600,
   width = 9,
   height = 12,
   units = "cm")


### Save results ----
  Observed <- rbind(observed_Papp_PFCAs,
     observed_Papp_PFSAs) %>%
 mutate(Source = "Observed")

  predicted_Papp_PFCAs <- predicted_Papp_PFCAs %>%
  mutate(Title = "Papp", Category_OECD = "PFCAs")
  predicted_Papp_PFSAs <- predicted_Papp_PFCAs %>%
 mutate(Title = "Papp", Category_OECD = "PFSAs")

  Predicted <- rbind(predicted_Papp_PFCAs,
      predicted_Papp_PFSAs) %>%
 rename("Value" = "predicted_value") %>%
 mutate(Source = "Predicted",
  Acronym = case_when(
    n_C == 4  & Category_OECD == "PFCAs" ~ "PFBA" ,
    n_C == 5  & Category_OECD == "PFCAs" ~ "PFPeA",
    n_C == 6  & Category_OECD == "PFCAs" ~ "PFHxA",
    n_C == 7  & Category_OECD == "PFCAs" ~ "PFHpA",
    n_C == 8  & Category_OECD == "PFCAs" ~ "PFOA",
    n_C == 9  & Category_OECD == "PFCAs" ~ "PFNA",
    n_C == 10  & Category_OECD == "PFCAs" ~ "PFDA",
    n_C == 11  & Category_OECD == "PFCAs" ~ "PFUnDA",
    n_C == 12  & Category_OECD == "PFCAs" ~ "PFDoDA",
    n_C == 13  & Category_OECD == "PFCAs" ~ "PFTrDA",
    n_C == 4  & Category_OECD == "PFSAs" ~ "PFBS" ,
    n_C == 5  & Category_OECD == "PFSAs" ~ "PFPeS",
    n_C == 6  & Category_OECD == "PFSAs" ~ "PFHxS",
    n_C == 7  & Category_OECD == "PFSAs" ~ "PFHpS",
    n_C == 8  & Category_OECD == "PFSAs" ~ "PFOS"),
  Parameter = case_when(
    Title == "Papp" ~ "Papp"),
  SD = NA)

  AbsorptionData <- rbind(Observed, Predicted)
  write.csv(AbsorptionData, "AbsorptionData.csv", row.names = FALSE)



# # Distribution parameters ----
# 
# ## Create dfs ----
# Distribution <- Final_df %>% 
#   filter(PBK_mechanism == "Distribution") %>% 
#   select(Acronym, Category_OECD, Parameter, Value, SD, n_C)
# Distribution$SD <- as.numeric(Distribution$SD)
# 
# Distribution <- Distribution %>%  
#   mutate(Title = case_when(
#  Parameter == "logD_AlbuminWater" ~ "logD (Albumin:Water)",
#  Parameter == "logD_GlobulinWater" ~ "logD (Globulin:Water)",
#  Parameter == "BPR" ~ "BPR",
#  Parameter == "Kd_FABP" ~ "Kd FABP",
#  Parameter == "logD_StructuralProteinWater" ~ "logD (Structural Protein:Water)",
#  Parameter == "logD_StorageLipidWater" ~ "logD (Storage Lipid:Water)",
#  Parameter == "logD_MembraneLipidWater" ~ "logD (Membrane Lipid:Water)"
#   )) %>%   
#   mutate(Value = case_when(
#  Parameter == "Kd_FABP" ~ log10(Value),
#  Parameter != "Kd_FABP" ~ Value
#   ))
# 
# Distribution_Parameters <- unique(Distribution$Title) %>% print()
# 
# DistributionAv <- DataAv %>% 
#   filter(PBK_mechanism == "Distribution")
# 
# DistributionAv <- DistributionAv %>%  
#   mutate(Title = case_when(
#  Parameter == "logD_AlbuminWater" ~ "logD (Albumin:Water)",
#  Parameter == "logD_GlobulinWater" ~ "logD (Globulin:Water)",
#  Parameter == "BPR" ~ "BPR",
#  Parameter == "Kd_FABP" ~ "Kd FABP",
#  Parameter == "logD_StructuralProteinWater" ~ "logD (Structural Protein:Water)",
#  Parameter == "logD_StorageLipidWater" ~ "logD (Storage Lipid:Water)",
#  Parameter == "logD_MembraneLipidWater" ~ "logD (Membrane Lipid:Water)"
#   ))
# 
# ## Derive values ----
# ### Albumin Water ----
# #### PFCAs ----
# logD_AlbuminWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFCAs", 
#    Parameter == "logD_AlbuminWater") 
# 
# y = logD_AlbuminWater_nls$Value
# x = logD_AlbuminWater_nls$n_C
# data <- tibble(x = x, y = y) 
# 
# fit_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)), 
#    data = data,
#    start = list(Const = 1.8, Asym = 5, xmid = 9.5, scal = 1.8))
# 
# Results <- summary(fit_nls)
# 
# Residuals_logD_AlbuminWater_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_logD_AlbuminWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,13))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 13, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_AlbuminWater_PFCAs <- tibble(
#   n_C = seq(4, 13, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_AlbuminWater_PFCAs <- logD_AlbuminWater_nls %>% arrange(n_C)
# # write.csv(observed_AlbuminWater_PFCAs, "logD_AlbuminWater_nls_PFCAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_AlbuminWater_PFCAs, "logD_AlbuminWater_nls_PFCAs.csv", row.names = FALSE)
# 
# 
# Formula <- Results$formula[3]
# RSE <- round(Results$sigma, 3)
# Const <- round(Results$coefficients[1], 3)
# Asym <- round(Results$coefficients[2], 3)
# xmid <- round(Results$coefficients[3], 3)
# scal <- round(Results$coefficients[4], 3)
# 
# D <- ggplot(data = observed_AlbuminWater_PFCAs, aes(n_C, Value)) +
#   geom_point(shape = 19, color = "#4ea9cf") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#4ea9cf", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_AlbuminWater_PFCAs, aes(n_C, predicted_value), shape = 19) +
#   scale_x_continuous(breaks=seq(1,13)) +
#   annotate("text", x=3.5, y=5.4, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.3, label=paste0("Const = ", Const, ", Asym = ", Asym, 
#              ", xmid = ", xmid, ", scal = ", scal),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.2, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
#   labs(y = "logD (Albumin:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_AlbuminWater_PFCAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# #### PFSAs ----
# logD_AlbuminWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFSAs", 
#    Parameter == "logD_AlbuminWater")
# 
# y = logD_AlbuminWater_nls$Value
# x = logD_AlbuminWater_nls$n_C
# data <- tibble(x = x, y = y)
# 
# fit_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)), 
#    data = data,
#    start = list(Const = 2.0738, Asym = 2.6764, xmid = 4.2436, scal = 0.8869))
# 
# Results <- summary(fit_nls)
# 
# Residuals_logD_AlbuminWater_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_logD_AlbuminWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,8))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 8, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_PFSAs <- tibble(
#   n_C = seq(4, 8, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
# ) %>% arrange(n_C)
# 
# logD_AlbuminWater_nls <- logD_AlbuminWater_nls %>% arrange(n_C)
# write.csv(logD_AlbuminWater_nls, "logD_AlbuminWater_nls_PFSAs_Observed.csv", row.names = FALSE)
# write.csv(predicted_PFSAs, "logD_AlbuminWater_nls_PFSAs.csv", row.names = FALSE)
# 
# 
# predicted_AlbuminWater_PFSAs <- tibble(
#   n_C = seq(4, 8, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_AlbuminWater_PFSAs <- logD_AlbuminWater_nls %>% arrange(n_C)
# # write.csv(observed_AlbuminWater_PFSAs, "logD_AlbuminWater_nls_PFSAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_AlbuminWater_PFSAs, "logD_AlbuminWater_nls_PFSAs.csv", row.names = FALSE)
# 
# Formula <- Results$formula[3]
# RSE <- round(Results$sigma, 3)
# Const <- round(Results$coefficients[1], 3)
# Asym <- round(Results$coefficients[2], 3)
# xmid <- round(Results$coefficients[3], 3)
# scal <- round(Results$coefficients[4], 3)
# 
# D <- ggplot(data = observed_AlbuminWater_PFSAs, aes(n_C, Value)) +
#   geom_point(shape = 18, color = "#b45c73") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#b45c73", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_AlbuminWater_PFSAs, aes(n_C, predicted_value), shape = 18) +
#   scale_x_continuous(breaks=seq(1,8)) +
#   annotate("text", x=3.5, y=5.8, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.7, label=paste0("Const = ", Const, ", Asym = ", Asym, 
#              ", xmid = ", xmid, ", scal = ", scal),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.6, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
#   labs(y = "logD (Albumin:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_AlbuminWater_PFSAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12, 
#  units = "cm")
# 
# 
# 
# 
# 
# ### Globulin Water ----
# #### PFCAs ----
# logD_GlobulinWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFCAs", 
#    Parameter == "logD_GlobulinWater")
# 
# y = logD_GlobulinWater_nls$Value
# x = logD_GlobulinWater_nls$n_C
# data <- tibble(x = x, y = y)
# 
# fit_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)), 
#    data = data,
#    start = list(Const = 1.8, Asym = 5, xmid = 9.5, scal = 1.8))
# 
# Results <- summary(fit_nls)
# 
# Residuals_logD_GlobulinWater_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_logD_GlobulinWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,13))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 13, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_GlobulinWater_PFCAs <- tibble(
#   n_C = seq(4, 13, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_GlobulinWater_PFCAs <- logD_GlobulinWater_nls %>% arrange(n_C)
# # write.csv(observed_GlobulinWater_PFCAs, "logD_GlobulinWater_nls_PFCAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_GlobulinWater_PFCAs, "logD_GlobulinWater_nls_PFCAs.csv", row.names = FALSE)
# 
# 
# Formula <- Results$formula[3]
# RSE <- round(Results$sigma, 3)
# Const <- round(Results$coefficients[1], 3)
# Asym <- round(Results$coefficients[2], 3)
# xmid <- round(Results$coefficients[3], 3)
# scal <- round(Results$coefficients[4], 3)
# 
# D <- ggplot(data = observed_GlobulinWater_PFCAs, aes(n_C, Value)) +
#   geom_point(shape = 19, color = "#4ea9cf") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#4ea9cf", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_GlobulinWater_PFCAs, aes(n_C, predicted_value), shape = 19) +
#   scale_x_continuous(breaks=seq(1,13)) +
#   annotate("text", x=3.5, y=5.4, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.3, label=paste0("Const = ", Const, ", Asym = ", Asym, 
#              ", xmid = ", xmid, ", scal = ", scal),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.2, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
#   labs(y = "logD (Globulin:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_GlobulinWater_PFCAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# 
# 
# 
# #### PFSAs ----
# logD_GlobulinWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFSAs", 
#    Parameter == "logD_GlobulinWater")
# 
# y = logD_GlobulinWater_nls$Value
# x = logD_GlobulinWater_nls$n_C
# data <- tibble(x = x, y = y)
# 
# fit_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)), 
#    data = data,
#    start = list(Const = 1.8, Asym = 5, xmid = 9.5, scal = 1.8))
# 
# Results <- summary(fit_nls)
# 
# Residuals_logD_GlobulinWater_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_logD_GlobulinWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,8))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 8, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_GlobulinWater_PFSAs <- tibble(
#   n_C = seq(4, 8, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_GlobulinWater_PFSAs <- logD_GlobulinWater_nls %>% arrange(n_C)
# # write.csv(observed_GlobulinWater_PFSAs, "logD_GlobulinWater_nls_PFSAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_GlobulinWater_PFSAs, "logD_GlobulinWater_nls_PFSAs.csv", row.names = FALSE)
# 
# Formula <- Results$formula[3]
# RSE <- round(Results$sigma, 3)
# Const <- round(Results$coefficients[1], 3)
# Asym <- round(Results$coefficients[2], 3)
# xmid <- round(Results$coefficients[3], 3)
# scal <- round(Results$coefficients[4], 3)
# 
# D <- ggplot(data = observed_GlobulinWater_PFSAs, aes(n_C, Value)) +
#   geom_point(shape = 18, color = "#b45c73") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#b45c73", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_GlobulinWater_PFSAs, aes(n_C, predicted_value), shape = 18) +
#   scale_x_continuous(breaks=seq(1,8)) +
#   annotate("text", x=3.5, y=5.4, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.3, label=paste0("Const = ", Const, ", Asym = ", Asym, 
#              ", xmid = ", xmid, ", scal = ", scal),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.2, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
#   labs(y = "logD (Globulin:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_GlobulinWater_PFSAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# ### Structural Protein ----
# #### PFCAs ----
# logD_StructuralProteinWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFCAs", 
#    Parameter == "logD_StructuralProteinWater") 
# 
# y = logD_StructuralProteinWater_nls$Value
# x = logD_StructuralProteinWater_nls$n_C
# data <- tibble(x = x, y = y) 
# 
# fit_nls <- lm(y ~ x, data = data)
# 
# Results <- summary(fit_nls)
# 
# Residuals_logD_StructuralProteinWater_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_logD_StructuralProteinWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,13))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 13, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_StructuralProteinWater_PFCAs <- tibble(
#   n_C = seq(4, 13, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_StructuralProteinWater_PFCAs <- logD_StructuralProteinWater_nls %>% arrange(n_C)
# # write.csv(observed_StructuralProteinWater_PFCAs, "logD_StructuralProteinWater_nls_PFCAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_StructuralProteinWater_PFCAs, "logD_StructuralProteinWater_nls_PFCAs.csv", row.names = FALSE)
# 
# Formula <- "Intercept * X"
# RSE <- round(Results$sigma, 3)
# R.squared <- round(Results$r.squared, 3)
# Intercept <- round(Results$coefficients[1], 3)
# 
# D <- ggplot(data = observed_StructuralProteinWater_PFCAs, aes(n_C, Value)) +
#   geom_point(shape = 19, color = "#4ea9cf") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#4ea9cf", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_StructuralProteinWater_PFCAs, aes(n_C, predicted_value), shape = 19) +
#   scale_x_continuous(breaks=seq(1,13)) +
#   annotate("text", x=3.5, y=5.4, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.2, label=paste0("Intercept = ", Intercept),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.0, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared), 
#   size = 10, hjust = 0) +
#   labs(y = "logD (Structural Protein:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_StructuralProteinWater_PFCAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# #### PFSAs ----
# logD_StructuralProteinWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFSAs", 
#    Parameter == "logD_StructuralProteinWater") 
# 
# y = logD_StructuralProteinWater_nls$Value
# x = logD_StructuralProteinWater_nls$n_C
# data <- tibble(x = x, y = y) 
# 
# fit_nls <- lm(y ~ x, data = data)
# 
# Results <- summary(fit_nls)
# 
# Residuals_logD_StructuralProteinWater_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_logD_StructuralProteinWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,8))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 8, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_StructuralProteinWater_PFSAs <- tibble(
#   n_C = seq(4, 8, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_StructuralProteinWater_PFSAs <- logD_StructuralProteinWater_nls %>% arrange(n_C)
# # write.csv(observed_StructuralProteinWater_PFSAs, "logD_StructuralProteinWater_nls_PFSAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_StructuralProteinWater_PFSAs, "logD_StructuralProteinWater_nls_PFSAs.csv", row.names = FALSE)
# 
# Formula <- "Intercept * X"
# RSE <- round(Results$sigma, 3)
# R.squared <- round(Results$r.squared, 3)
# Intercept <- round(Results$coefficients[1], 3)
# 
# D <- ggplot(data = observed_StructuralProteinWater_PFSAs, aes(n_C, Value)) +
#   geom_point(shape = 18, color = "#b45c73") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#b45c73", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_StructuralProteinWater_PFSAs, aes(n_C, predicted_value), shape = 18) +
#   scale_x_continuous(breaks=seq(1,8)) +
#   annotate("text", x=3.5, y=5.6, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.4, label=paste0("Intercept = ", Intercept),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.2, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared), 
#   size = 10, hjust = 0) +
#   labs(y = "logD (Structural Protein:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_StructuralProteinWater_PFSAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# 
# 
# ### Storage Lipid ----
# #### PFCAs ----
# logD_StorageLipidWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFCAs", 
#    Parameter == "logD_StorageLipidWater") 
# 
# y = c(-1.55, -0.55, -1.50, -1.50, -0.73, -1.37, -0.45)
# x = c(4,  10,  7,  6,  9,  8, 11)
# data <- tibble(x = x, y = y) 
# 
# fit_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
#    data = data,
#    start = list(Const = -1.5, Asym = 0.86767, xmid = 8.5, scal = 0.26))
# 
# Results <- summary(fit_nls)
# 
# # Manually calculate residuals as the model was not fit to the exact observed data
# Residuals_logD_StorageLipidWater_nls = tibble(
#   x = c(4, 6, 7, 8, 9, 10, 11), 
#   observed = c(-1.40, -1.50, -1.10, -1.37, -0.73, -0.55, -0.32),
#   predicted = c(-1.52, -1.52, -1.51, -1.36, -0.739, -0.506, -0.488)
# ) %>% mutate(y = observed-predicted)
# 
# #calculate the number of model parameters - 1
# k=length(fit_nls$coefficients)-1
# 
# #calculate sum of squared residuals
# SSE=sum(Residuals_logD_StorageLipidWater_nls$y**2)
# 
# #calculate total observations in dataset
# n=length(Residuals_logD_StorageLipidWater_nls$y)
# 
# #calculate residual standard error
# RSE = sqrt(SSE/(n-(1+k)))
# 
# R <- ggplot(Residuals_logD_StorageLipidWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,13))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 13, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_StorageLipidWater_PFCAs <- tibble(
#   n_C = seq(4, 13, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_StorageLipidWater_PFCAs <- logD_StorageLipidWater_nls %>% arrange(n_C)
# # 
# # write.csv(observed_StorageLipidWater_PFCAs, "logD_StorageLipidWater_nls_PFCAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_StorageLipidWater_PFCAs, "logD_StorageLipidWater_nls_PFCAs.csv", row.names = FALSE)
# #  
# Formula <- Results$formula[3]
# RSE <- round(RSE, 3)
# Const <- round(Results$coefficients[1],3)
# Asym <- round(Results$coefficients[2], 3)
# xmid <- round(Results$coefficients[3], 3)
# scal <- round(Results$coefficients[4], 3)
# 
# D <- ggplot(data = observed_StorageLipidWater_PFCAs, aes(n_C, Value)) +
#   geom_point(shape = 19, color = "#4ea9cf") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#4ea9cf", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_StorageLipidWater_PFCAs, aes(n_C, predicted_value), shape = 19) +
#   scale_x_continuous(breaks=seq(1,13)) +
#   annotate("text", x=3.5, y=-0.33, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=-0.36, label=paste0("Const = ", Const, ", Asym = ", Asym, 
#              ", xmid = ", xmid, ", scal = ", scal),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=-0.39, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
#   labs(y = "logD (Storage Lipid:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_StorageLipidWater_PFCAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# 
# #### PFSAs ----
# logD_StorageLipidWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFSAs", 
#    Parameter == "logD_StorageLipidWater") 
# 
# y = c(-1.63, -0.63, -1.57, -1.57, -0.8, -1.44, -0.52)
# x = c(2.5,  8.5,  5.5,  4.5,  7.5,  6.5, 9.5)
# data <- tibble(x = x, y = y) 
# 
# fit_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)),
#    data = data,
#    start = list(Const = -1.6, Asym = 1.036, xmid = 6.4, scal = 0.353))
# 
# Results <- summary(fit_nls)
# 
# # Manually calculate residuals as the model was not fit to the exact observed data
# Residuals_logD_StorageLipidWater_nls = tibble(
#   x = c(6, 8), 
#   observed = c(-0.56, -1.58),
#   predicted = c(-0.636, -1.55)
# ) %>% mutate(y = observed-predicted)
# 
# #calculate the number of model parameters - 1
# k=length(fit_nls$coefficients)-1
# 
# #calculate sum of squared residuals
# SSE=sum(Residuals_logD_StorageLipidWater_nls$y**2)
# 
# #calculate total observations in dataset
# n=length(Residuals_logD_StorageLipidWater_nls$y)
# 
# #calculate residual standard error
# RSE = sqrt(SSE/(n-(1+k)))
# 
# R <- ggplot(Residuals_logD_StorageLipidWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,8))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 8, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_StorageLipidWater_PFSAs <- tibble(
#   n_C = seq(4, 8, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_StorageLipidWater_PFSAs <- logD_StorageLipidWater_nls %>% arrange(n_C)
# 
# # write.csv(observed_StorageLipidWater_PFSAs, "logD_StorageLipidWater_nls_PFSAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_StorageLipidWater_PFSAs, "logD_StorageLipidWater_nls_PFSAs.csv", row.names = FALSE)
# 
# Formula <- Results$formula[3]
# RSE <- round(RSE, 3)
# Const <- round(Results$coefficients[1],3)
# Asym <- round(Results$coefficients[2], 3)
# xmid <- round(Results$coefficients[3], 3)
# scal <- round(Results$coefficients[4], 3)
# 
# D <- ggplot(data = observed_StorageLipidWater_PFSAs, aes(n_C, Value)) +
#   geom_point(shape = 18, color = "#b45c73") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#b45c73", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_StorageLipidWater_PFSAs, aes(n_C, predicted_value), shape = 18) +
#   scale_x_continuous(breaks=seq(1,8)) +
#   annotate("text", x=3.5, y=-0.5, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=-0.53, label=paste0("Const = ", Const, ", Asym = ", Asym,
#             ", xmid = ", xmid, ", scal = ", scal),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=-0.56, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
#   labs(y = "logD (Storage Lipid:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_StorageLipidWater_PFSAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# 
# 
# ### Membrane Lipid Water ----
# #### PFCAs ----
# logD_MembraneLipidWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFCAs", 
#    Parameter == "logD_MembraneLipidWater")
# 
# y = logD_MembraneLipidWater_nls$Value
# x = logD_MembraneLipidWater_nls$n_C
# data <- tibble(x = x, y = y)
# 
# fit_nls <- nls(y ~ Const + Asym/(1+exp((xmid-x)/scal)), 
#    data = data,
#    start = list(Const = 1.8, Asym = 5, xmid = 9.5, scal = 1.8))
# 
# Results <- summary(fit_nls)
# 
# Residuals_logD_MembraneLipidWater_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_logD_MembraneLipidWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,13))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 13, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_MembraneLipid_PFCAs <- tibble(
#   n_C = seq(4, 13, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_MembraneLipid_PFCAs <- logD_MembraneLipidWater_nls %>% arrange(n_C)
# # write.csv(observed_MembraneLipid_PFCAs, "logD_MembraneLipidWater_nls_PFCAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_MembraneLipid_PFCAs, "logD_MembraneLipidWater_nls_PFCAs.csv", row.names = FALSE)
# 
# 
# Formula <- Results$formula[3]
# RSE <- round(Results$sigma, 3)
# Const <- round(Results$coefficients[1], 3)
# Asym <- round(Results$coefficients[2], 3)
# xmid <- round(Results$coefficients[3], 3)
# scal <- round(Results$coefficients[4], 3)
# 
# D <- ggplot(data = observed_MembraneLipid_PFCAs, aes(n_C, Value)) +
#   geom_point(shape = 19, color = "#4ea9cf") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#4ea9cf", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_MembraneLipid_PFCAs, aes(n_C, predicted_value), shape = 19) +
#   scale_x_continuous(breaks=seq(1,13)) +
#   annotate("text", x=3.5, y=5.4, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.3, label=paste0("Const = ", Const, ", Asym = ", Asym, 
#              ", xmid = ", xmid, ", scal = ", scal),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.2, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
#   labs(y = "logD (Membrane Lipid:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_MembraneLipidWater_PFCAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# 
# 
# 
# #### PFSAs ----
# logD_MembraneLipidWater_nls <- Distribution %>% 
#   filter(Category_OECD == "PFSAs", 
#    Parameter == "logD_MembraneLipidWater")
# 
# y = logD_MembraneLipidWater_nls$Value
# x = logD_MembraneLipidWater_nls$n_C
# 
# data <- tibble(x = x, y = y)
# 
# fit_nls <- lm(y ~ x)
# 
# Results <- summary(fit_nls)
# 
# Residuals_logD_MembraneLipidWater_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_logD_MembraneLipidWater_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,8))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 8, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_MembraneLipid_PFSAs <- tibble(
#   n_C = seq(4, 8, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_MembraneLipid_PFSAs <- logD_MembraneLipidWater_nls %>% arrange(n_C)
# # write.csv(observed_MembraneLipid_PFSAs, "logD_MembraneLipidWater_nls_PFSAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_MembraneLipid_PFSAs, "logD_MembraneLipidWater_nls_PFSAs.csv", row.names = FALSE)
# 
# Formula <- "Intercept * X"
# RSE <- round(Results$sigma, 3)
# R.squared <- round(Results$r.squared, 3)
# Intercept <- round(Results$coefficients[1], 3)
# 
# D <- ggplot(data = observed_MembraneLipid_PFSAs, aes(n_C, Value)) +
#   geom_point(shape = 18, color = "#b45c73") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#b45c73", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_MembraneLipid_PFSAs, aes(n_C, predicted_value), shape = 18) +
#   scale_x_continuous(breaks=seq(1,8)) +
#   annotate("text", x=3.5, y=5.6, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.4, label=paste0("Intercept = ", Intercept),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=5.2, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared), 
#   size = 10, hjust = 0) +
#   labs(y = "logD (Membrane Lipid:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("logD_MembraneLipidWater_PFSAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# ### FABP Distribution coefficient ----
# #### PFCAs ----
# Kd_FABP_nls <- Distribution %>% 
#   filter(Category_OECD == "PFCAs", 
#    Parameter == "Kd_FABP") 
# 
# y = Kd_FABP_nls$Value
# x = Kd_FABP_nls$n_C
# data <- tibble(x = x, y = y) 
# 
# fit_nls <- nls(y ~ SSlogis(x, Asym, xmid, scal),
#    data = data,
#    control = nls.control(maxiter = 1000000, minFactor = 1e-20))
# 
# 
# Results <- summary(fit_nls)
# 
# Residuals_Kd_FABP_nls <- tibble(
#   x = te,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_Kd_FABP_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,13))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 13, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 13, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_KdFABP_PFCAs <- tibble(
#   n_C = seq(4, 13, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 13, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_KdFABP_PFCAs <- Kd_FABP_nls %>% arrange(n_C)
# # write.csv(observed_KdFABP_PFCAs, "Kd_FABP_nls_PFCAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_KdFABP_PFCAs, "Kd_FABP_nls_PFCAs.csv", row.names = FALSE)
# 
# 
# Formula <- Results$formula[3]
# RSE <- round(Results$sigma, 3)
# Asym <- round(Results$coefficients[2], 3)
# xmid <- round(Results$coefficients[3], 3)
# scal <- round(Results$coefficients[4], 3)
# 
# D <- ggplot(data = observed_KdFABP_PFCAs, aes(n_C, Value)) +
#   geom_point(shape = 19, color = "#4ea9cf") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#4ea9cf", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_KdFABP_PFCAs, aes(n_C, predicted_value), shape = 19) +
#   scale_x_continuous(breaks=seq(1,13)) +
#   annotate("text", x=3.5, y=-0.9, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=-0.6, label=paste0( "Asym = ", Asym, 
#              ", xmid = ", xmid, ", scal = ", scal),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=-0.3, label=paste0("RSE = ", RSE), size = 10, hjust = 0) +
#   labs(y = "logD (Albumin:Water)") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("Kd_FABP_PFCAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12,
#  units = "cm")
# 
# #### PFSAs ----
# Kd_FABP_nls <- Distribution %>% 
#   filter(Category_OECD == "PFSAs", 
#    Parameter == "Kd_FABP")
# 
# y = Kd_FABP_nls$Value
# x = Kd_FABP_nls$n_C
# data <- tibble(x = x, y = y)
# 
# fit_nls <- lm(y ~ x)
# 
# Results <- summary(fit_nls)
# 
# Residuals_Kd_FABP_nls <- tibble(
#   x = x,
#   y = residuals(fit_nls)
# )
# 
# R <- ggplot(Residuals_Kd_FABP_nls) +
#   geom_point(aes(x, y)) +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(breaks=seq(1,8))+
#   labs(y = "Residuals", x = "n_C") +
#   CP_theme 
# 
# predicted <- tibble(
#   n_C = seq(4, 8, by = 0.01),
#   predicted_y = predict(fit_nls, list(x = seq(4, 8, by = 0.01)))
# ) %>% arrange(n_C)
# 
# predicted_PFSAs <- tibble(
#   n_C = seq(4, 8, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
# ) %>% arrange(n_C)
# 
# Kd_FABP_nls <- Kd_FABP_nls %>% arrange(n_C)
# write.csv(Kd_FABP_nls, "Kd_FABP_nls_PFSAs_Observed.csv", row.names = FALSE)
# write.csv(predicted_PFSAs, "Kd_FABP_nls_PFSAs.csv", row.names = FALSE)
# 
# 
# predicted_KdFABP_PFSAs <- tibble(
#   n_C = seq(4, 8, by = 1),
#   predicted_value = predict(fit_nls, list(x = seq(4, 8, by = 1)))
# ) %>% arrange(n_C)
# 
# observed_KdFABP_PFSAs <- Kd_FABP_nls %>% arrange(n_C)
# # write.csv(observed_KdFABP_PFSAs, "Kd_FABP_nls_PFSAs.csv_Observed.csv", row.names = FALSE)
# # write.csv(predicted_KdFABP_PFSAs, "Kd_FABP_nls_PFSAs.csv", row.names = FALSE)
# 
# Formula <- "Intercept * X"
# RSE <- round(Results$sigma, 3)
# R.squared <- round(Results$r.squared, 3)
# Intercept <- round(Results$coefficients[1], 3)
# 
# 
# D <- ggplot(data = observed_KdFABP_PFSAs, aes(n_C, Value)) +
#   geom_point(shape = 18, color = "#b45c73") +
#   geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
#     colour =  "#b45c73", width = 0.25, alpha = 0.8) +
#   geom_path(data = predicted, aes(n_C, predicted_y)) +
#   geom_point(data = predicted_KdFABP_PFSAs, aes(n_C, predicted_value), shape = 18) +
#   scale_x_continuous(breaks=seq(1,8)) +
#   annotate("text", x=3.5, y=-0.9, label=paste0("y = ", Formula), size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=-0.6, label=paste0("Intercept = ", Intercept),
#   size = 10, hjust = 0) +
#   annotate("text", x=3.5, y=-0.3, label=paste0("RSE = ", RSE, ", R.squared = ", R.squared), 
#   size = 10, hjust = 0) +
#   labs(y = "Kd FABP") +
#   CP_theme 
# 
# SavePlots <- D / R + plot_layout(heights = c(5,1))
# ggsave(plot = SavePlots, file = paste("Kd_FABP_PFSAs", ".png",sep=""),
#  dpi = 600,
#  width = 9, 
#  height = 12, 
#  units = "cm")
# 
# 
# # Save results ----
# observed_AlbuminWater_PFECAs <- Distribution %>% 
#   filter(Category_OECD == "PFECAs", 
#    Parameter == "logD_AlbuminWater") %>% arrange(n_C)
# observed_StorageLipidWater_PFECAs <- Distribution %>% 
#   filter(Category_OECD == "PFECAs", 
#    Parameter == "logD_StorageLipidWater") %>% arrange(n_C)
# observed_StructuralProteinWater_PFECAs <- Distribution %>% 
#   filter(Category_OECD == "PFECAs", 
#    Parameter == "logD_StructuralProteinWater") %>% arrange(n_C)
# observed_KdFABP_PFECAs <- Distribution %>% 
#   filter(Category_OECD == "PFECAs", 
#    Parameter == "Kd_FABP") %>% arrange(n_C)
# observed_BPR_PFECAs <- Distribution %>% 
#   filter(Category_OECD == "PFECAs", 
#    Parameter == "BPR") %>% arrange(n_C)
# observed_BPR_PFCAs <- Distribution %>% 
#   filter(Category_OECD == "PFCAs", 
#    Parameter == "BPR") %>% arrange(n_C)
# observed_BPR_PFSAs <- Distribution %>% 
#   filter(Category_OECD == "PFSAs", 
#    Parameter == "BPR") %>% arrange(n_C)
# observed_MembraneLipidWater_PFECAs <- Distribution %>% 
#   filter(Category_OECD == "PFECAs", 
#    Parameter == "logD_MembraneLipidWater") %>% arrange(n_C)
# observed_GlobulinWater_PFECAs <- Distribution %>% 
#   filter(Category_OECD == "PFECAs", 
#    Parameter == "logD_GlobulinWater") %>% arrange(n_C)
# 
# Observed <- rbind(observed_AlbuminWater_PFCAs, 
#       observed_AlbuminWater_PFSAs, 
#       observed_AlbuminWater_PFECAs,
#       observed_GlobulinWater_PFCAs, 
#       observed_GlobulinWater_PFSAs,
#       observed_GlobulinWater_PFECAs,
#       observed_StructuralProteinWater_PFCAs, 
#       observed_StructuralProteinWater_PFSAs,
#       observed_StructuralProteinWater_PFECAs,
#       observed_StorageLipidWater_PFCAs,
#       observed_StorageLipidWater_PFSAs,
#       observed_StorageLipidWater_PFECAs,
#       observed_MembraneLipid_PFCAs, 
#       observed_MembraneLipid_PFSAs, 
#       observed_MembraneLipidWater_PFECAs,
#       observed_KdFABP_PFCAs, 
#       observed_KdFABP_PFSAs,
#       observed_KdFABP_PFECAs,
#       observed_BPR_PFCAs,
#       observed_BPR_PFSAs,
#       observed_BPR_PFECAs
#       ) %>% mutate(Source = "Observed")
# 
# predicted_AlbuminWater_PFCAs <- predicted_AlbuminWater_PFCAs %>% 
# mutate(Title = "logD (Albumin:Water)", Category_OECD = "PFCAs")
# predicted_AlbuminWater_PFSAs <- predicted_AlbuminWater_PFSAs %>% 
#   mutate(Title = "logD (Albumin:Water)", Category_OECD = "PFSAs")
# predicted_GlobulinWater_PFCAs <- predicted_GlobulinWater_PFCAs %>% 
#   mutate(Title = "logD (Globulin:Water)", Category_OECD = "PFCAs")
# predicted_GlobulinWater_PFSAs <- predicted_GlobulinWater_PFSAs %>% 
#   mutate(Title = "logD (Globulin:Water)", Category_OECD = "PFSAs")
# predicted_StructuralProteinWater_PFCAs <- predicted_StructuralProteinWater_PFCAs %>% 
#   mutate(Title = "logD (Structural Protein:Water)", Category_OECD = "PFCAs") 
# predicted_StructuralProteinWater_PFSAs <- predicted_StructuralProteinWater_PFSAs %>% 
#   mutate(Title = "logD (Structural Protein:Water)", Category_OECD = "PFSAs") 
# predicted_StorageLipidWater_PFCAs <- predicted_StorageLipidWater_PFCAs %>% 
#   mutate(Title = "logD (Storage Lipid:Water)", Category_OECD = "PFCAs") 
# predicted_StorageLipidWater_PFSAs <- predicted_StorageLipidWater_PFSAs %>% 
#   mutate(Title = "logD (Storage Lipid:Water)", Category_OECD = "PFSAs") 
# predicted_MembraneLipid_PFCAs <- predicted_MembraneLipid_PFCAs %>% 
#   mutate(Title = "logD (Membrane Lipid:Water)", Category_OECD = "PFCAs") 
# predicted_MembraneLipid_PFSAs <- predicted_MembraneLipid_PFSAs %>% 
#   mutate(Title = "logD (Membrane Lipid:Water)", Category_OECD = "PFSAs") 
# predicted_KdFABP_PFCAs <- predicted_KdFABP_PFCAs %>% 
#   mutate(Title = "Kd FABP", Category_OECD = "PFCAs") 
# predicted_KdFABP_PFSAs <- predicted_KdFABP_PFSAs %>% 
#   mutate(Title = "Kd FABP", Category_OECD = "PFSAs") 
# 
# Predicted <- rbind(predicted_AlbuminWater_PFCAs, 
#     predicted_AlbuminWater_PFSAs,
#     predicted_GlobulinWater_PFCAs,
#     predicted_GlobulinWater_PFSAs,
#     predicted_StructuralProteinWater_PFCAs, 
#     predicted_StructuralProteinWater_PFSAs,
#     predicted_StorageLipidWater_PFCAs,
#     predicted_StorageLipidWater_PFSAs,
#     predicted_MembraneLipid_PFCAs,
#     predicted_MembraneLipid_PFSAs,
#     predicted_KdFABP_PFCAs,
#     predicted_KdFABP_PFSAs) %>% 
#   rename("Value" = "predicted_value") %>% 
#   mutate(Source = "Predicted", 
#    Acronym = case_when(
#   n_C == 4  & Category_OECD == "PFCAs" ~ "PFBA" ,
#   n_C == 5  & Category_OECD == "PFCAs" ~ "PFPeA",  
#   n_C == 6  & Category_OECD == "PFCAs" ~ "PFHxA",  
#   n_C == 7  & Category_OECD == "PFCAs" ~ "PFHpA",  
#   n_C == 8  & Category_OECD == "PFCAs" ~ "PFOA",  
#   n_C == 9  & Category_OECD == "PFCAs" ~ "PFNA",  
#   n_C == 10  & Category_OECD == "PFCAs" ~ "PFDA",  
#   n_C == 11  & Category_OECD == "PFCAs" ~ "PFUnDA",  
#   n_C == 12  & Category_OECD == "PFCAs" ~ "PFDoDA",  
#   n_C == 13  & Category_OECD == "PFCAs" ~ "PFTrDA",  
#   n_C == 4  & Category_OECD == "PFSAs" ~ "PFBS" ,
#   n_C == 5  & Category_OECD == "PFSAs" ~ "PFPeS",  
#   n_C == 6  & Category_OECD == "PFSAs" ~ "PFHxS",  
#   n_C == 7  & Category_OECD == "PFSAs" ~ "PFHpS",  
#   n_C == 8  & Category_OECD == "PFSAs" ~ "PFOS"), 
#    Parameter = case_when(
#   Title == "logD (Albumin:Water)" ~ "logD_AlbuminWater",  
#   Title == "logD (Globulin:Water)" ~ "logD_GlobulinWater",  
#   Title == "logD (Structural Protein:Water)" ~ "logD_StructuralProteinWater",  
#   Title == "logD (Storage Lipid:Water)" ~ "logD_StorageLipidWater",  
#   Title == "logD (Membrane Lipid:Water)" ~ "logD_MembraneLipidWater",  
#   Title == "Kd FABP" ~ "Kd_FABP"), 
#    SD = NA)
# 
# DistributionData <- rbind(Observed, Predicted)
# write.csv(DistributionData, "DistributionData.csv", row.names = FALSE)
# 
# 
# 
# # fit_nls <- lm(y ~ x, data = data)
# # fit_nls <- nls(y ~ SSlogis(x, Asym, xmid, scal),
# #    data = data,
# #    control = nls.control(maxiter = 1000000, minFactor = 1e-20))
# 
# 
# 

# Initial distribution plots ####
# 
# Parameter_i <- NULL 
# Plot_i <- NULL
# PlotAv_i <- NULL
# Data_i <- NULL
# DataAv_i <- NULL
# Plot_list <- list()
# 
# for (i in Distribution_Parameters) {
#   
#   Parameter_i <- i
#   Data_i <- filter(Distribution, Title == Parameter_i) 
#   DataAv_i <- filter(DistributionAv, Title == Parameter_i) %>% 
#  pivot_longer("PFCAs":"PFSAs", names_to = "Category_OECD", values_to = "Percent")
#   
#   Plot_i <- ggplot(Data_i, 
#     aes(y = Value,
#      x = n_C, 
#      colour = Category_OECD,
#      shape = Category_OECD)
#   ) +
#  geom_errorbar(data = Data_i, 
#       aes(ymin = Value - SD, ymax = Value + SD,
#        colour = Category_OECD),
#       width = 0.25, alpha = 0.8) + 
#  scale_x_continuous(breaks=seq(1,13))+
#  labs(x = "nC",
#    y = "Parameter value") 
#   
#   PlotAv_i <- ggplot(DataAv_i, 
#       aes(y = Percent, 
#        x = Category_OECD,
#        fill = Category_OECD)) +
#  geom_col() +
#  scale_fill_manual(values = c("PFCAs" = "#4ea9cf",
#          "PFECAs"  = "#55b483",
#          "PFSAs" = "#b45c73")
#  )+
#  coord_flip() +
#  labs(fill = "",
#    x = "",
#    y = "") +
#  CP_theme + #bw() +
#  theme(#legend.text = element_text(size = 22), #
#    #legend.title = element_text(size = 20),
#    legend.position = "none",
#    axis.text = element_text(size = 30)) 
#   
#   Plot_i <- Plot_i + 
#  geom_point(size = 1.5, stroke = 1, alpha = 0.9) + 
#  scale_shape_manual(values = c("PFCAs" = 19,
#           "PFECAs"  = 15,
#           "PFSAs" = 18)) +
#  scale_colour_manual(values = c("PFCAs" = "#4ea9cf",
#         "PFECAs"  = "#55b483",
#         "PFSAs" = "#b45c73")) +
#  labs(title = paste0(i),
#    colour = "C. OECD",   
#    shape = "C. OECD") +
#  CP_theme +# theme_bw() +
#  theme(#plot.title = element_text(size = 28, face = "bold"),
#    #plot.title.position = "plot", 
#    #legend.text = element_text(size = 22), #
#    legend.position = "none"#,
#    #axis.title = element_text(size = 22, face = "bold"),
#    #axis.text = element_text(size = 20, face = "bold"),
#    #legend.key.spacing = unit(0.2, "cm"),
#    #legend.key.size = unit(1, "lines"),
#    #legend.margin = margin(0,0,0,0)
#  ) 
#   
#   Plot_I <- Plot_i + PlotAv_i + plot_layout(widths = c(5,1))
#   
#   Plot_list[[i]] <- Plot_I
#   
# }
# 
# for(i in 1:7){
#   ggsave(plot = Plot_list[[i]], file = paste("file",i,".png",sep=""),
#    dpi = 600,
#    width = 9, 
#    height = 10, 
#    units = "cm")
# }