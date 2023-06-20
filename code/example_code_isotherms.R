# Library ----
library(data.table)
library(readxl)
library(latex2exp)
library(ggtext)
library(scales)
library(writexl)
library(tidyverse)
library(ggtext)
library(glue)
library(plotrix)
library(moderndive)
library(broom)
library(knitr)


# Data manipulation sorption isotherms ----
Sorption_BC <- read_excel("data_raw/160222_sorption_rawdata.xlsx") %>% 
  na.omit() %>% 
  filter(Conc_point != 1) %>% 
  mutate(SoilLogic = as.logical(Soil_binary),
         mixLogic = as.logical(mix_binary),
         isothermLogic = as.logical(isotherm_binary),
         Kd = Cs / Cw) %>% 
  filter(mixLogic == FALSE) %>% 
  select(Conc_point, Compound, Biochar, type, SoilLogic, mixLogic, isothermLogic, 
         Cs, Cw, Kd, Ci) 

Sorption_soil <- read_excel("data_raw/010322_sorption_rawdata_soil.xlsx") %>% 
  drop_na(log_Cs) %>% 
  mutate(BClogic = if_else(Biochar == "no",
                           TRUE,
                           FALSE),
         SoilLogic = as.logical(Soil_binary),
         mixLogic = as.logical(mix_binary),
         isothermLogic = as.logical(isotherm_binary),
         Kd = Cs / Cw) %>% 
  filter(Biochar != "no",
         #Compound %in% c("PFOA", "PFNA", "PFDA")
  ) %>% 
  select(Conc_point, Compound, Biochar, type, SoilLogic, mixLogic, isothermLogic, 
         BClogic, Cs, Cw, K_ds, Kd)

Sorption_BC_soil <- full_join(Sorption_BC, Sorption_soil)

Sorption_BC_soil %>% filter(isothermLogic == F) %>% write_xlsx("R/data_manipulated/090822_rawdata_joined.xlsx")

regression_glance <- Sorption_BC_soil %>% 
  group_by(Compound, Biochar, type) %>%
  do(fit_isotherms = glance(lm(log10(Cs) ~ log10(Cw), data = .))) %>% 
  unnest(fit_isotherms)

regression_tidy <- Sorption_BC_soil %>% 
  group_by(Compound, Biochar, type) %>%
  do(fit_isotherms = tidy(lm(log10(Cs) ~ log10(Cw), data = .))) %>% 
  unnest(fit_isotherms) %>% 
  pivot_wider(names_from = term,
              values_from = c(estimate, std.error, statistic, p.value))

regression_statistics_BC_soil <- full_join(regression_glance, 
                                           regression_tidy) %>% 
  rename(K_F = "estimate_(Intercept)",
         n_F = "estimate_log10(Cw)",
         K_F_se = "std.error_(Intercept)",
         n_F_se = "std.error_log10(Cw)"
  ) %>% 
  select(Compound, Biochar, type, r.squared, p.value, K_F, n_F, K_F_se, n_F_se, nobs)

write_xlsx(regression_statistics_BC_soil, "R/data_manipulated/150622_Freundlich_coefficients.xlsx")

# Sorption isotherm plots ----
Sorption_isotherms_nonsigremoved <- Sorption_BC_soil %>% 
  filter(SoilLogic == FALSE,
         Compound != "PFPeA" | Biochar != "DSL",
         Compound != "PFPeA" | Biochar != "CWC",
         Compound != "PFHxA" | Biochar != "CWC") %>% 
  mutate(Compound = factor(Compound, 
                           levels = c("PFDA", "PFNA", "PFOA", 
                                      "PFHpA", "PFHxA","PFPeA")
  )) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  geom_point(mapping = aes(x = log10(Cw), 
                           y = log10(Cs), 
                           color = factor(Biochar)
  ), 
  #color = "grey",
  size = 3,
  alpha = 0.6) + 
  geom_smooth(mapping = aes(x = log10(Cw), y = log10(Cs), color = factor(Biochar)), 
              formula = y ~ x, 
              method=lm, 
              se=F, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~Compound,
             #scales = "free_x",
             ncol = 2) +
  theme_bw() +
  theme(text = element_text(size = 30),
        #legend.text=element_text(size=20),
        panel.spacing = unit(0.2, "cm"),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,0,-10),
        axis.title.x = element_text(margin = margin(10,0,5,0))) +
  scale_color_manual(labels = c("WCBC", "SSBC1", "SSBC2"),
                     breaks = c("CWC", "ULS", "DSL"),
                     values = c("#FFB547FF","#4E9C81","#40E0CF")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
Sorption_isotherms_nonsigremoved
ggsave(filename="R/figs/article/Sorption_isotherms_single_nonsigremoved.pdf")