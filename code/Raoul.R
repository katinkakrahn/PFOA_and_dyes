# Library ----
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
library(smplot2)
library(ggpubr)
library(RColorBrewer)

# raw data dyes ----
dyes_raw_data <- read_excel("data_raw/190123_raw_data_dyes.xlsx") 

dyes <- dyes_raw_data %>% filter(comment == "good",
                                 sample_type == "sample") %>% 
  mutate(Cw_mgL = ifelse(dye == "RB", absorbance * 25.953-1.9764, absorbance * 5.4014-0.0348),
         Cs_mgkg = (50-Cw_mgL)*0.05/mass_BC_g*1000,
         Cw_mgL = ifelse(Cw_mgL > 0, Cw_mgL, 0.01),
         Kd_Lkg = Cs_mgkg/Cw_mgL,
         log_Kd = log10(Kd_Lkg))

dye_summary <- dyes %>% 
  filter(comment == "good") %>% 
  group_by(sample, sample_ID, sample_name, dye, sample_type, date_meas) %>% 
  summarise(Cw_mean = mean(Cw_mgL),
            Cs_mean = mean(Cs_mgkg),
            Kd_mean = mean(Cs_mgkg / Cw_mgL),
            mean_log_Kd = mean(log10(Kd_mean)),
            sd_log_Kd = sd(log_Kd),
            n = n())

write_xlsx(dye_summary, "data_manipulated/200123_dye_summary.xlsx")

# Raw data PFOA ----
PFOA_raw_data <- read_excel("data_raw/18112022_raw_data_PFOA.xlsx") %>% 
  filter(sample_type == "sample")

PFOA_samplesonly <- PFOA_raw_data %>% 
  drop_na(log_Kd) %>%
  filter(!(sample %in% c("SB", "NB")))


# PFOA Freundlich isotherm statistics ----
regression_glance <- PFOA_samplesonly %>% 
  group_by(sample, sample_ID, sample_name, replicate, biochar, temperature, residence_time, source) %>%
  do(fit_isotherms = glance(lm(log10(Cs_ugkg) ~ log10(Cw_ugL), data = .))) %>% 
  unnest(fit_isotherms)

regression_tidy <- PFOA_samplesonly %>% 
  group_by(sample, sample_ID, sample_name, replicate, biochar, temperature, residence_time, source) %>%
  do(fit_isotherms = tidy(lm(log10(Cs_ugkg) ~ log10(Cw_ugL), data = .))) %>% 
  unnest(fit_isotherms) %>% 
  pivot_wider(names_from = term,
              values_from = c(estimate, std.error, statistic, p.value))

PFOA_regression_statistics <- full_join(regression_glance, 
                                        regression_tidy) %>% 
  rename(log_KF = "estimate_(Intercept)",
         n_F = "estimate_log10(Cw_ugL)",
         log_KF_se = "std.error_(Intercept)",
         n_F_se = "std.error_log10(Cw_ugL)"
  ) %>% 
  select(sample, sample_ID, sample_name, replicate, 
         biochar, temperature, residence_time, source,
         r.squared, p.value, 
         log_KF, n_F, log_KF_se, n_F_se, nobs)

write_xlsx(PFOA_regression_statistics, "data_manipulated/181122_Freundlich_regression_statistics.xlsx")

# join dyes and isotherms ----
PFOA_dyes <- full_join(PFOA_regression_statistics, dye_summary)

PFOA_dyes %>% 
  filter(sample_type == "sample") %>% 
  select(sample, log_KF, log_KF_se, log_Kd_1mgL, mean_log_Kd, sd_log_Kd, dye) %>%
  write_xlsx("data_manipulated/111222_PFOA_dyes.xlsx")

# correlation plots ----
correlation_MB <- PFOA_dyes %>% 
  na.omit() %>% 
  filter(dye == "MB") %>% 
  filter(!(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>% 
  ggplot(mapping = aes(x = mean_log_Kd, 
                       y = log_KF)) +
  geom_point(color = "blue",
             size = 3) + 
  labs(x = TeX(r'($K_{d}~(L/kg)~dyes$)'), 
       y = TeX(r'($K_{F}~(L/kg)~PFOA~isotherms$)')) +
  sm_corr_theme() +
  sm_statCorr(corr_method = 'spearman',
              fit.params = list(color = 'blue',
                                linetype = 'dashed'))
correlation_MB
ggsave(filename = "figs/correlation_MB.jpeg")

correlation_RB <- PFOA_dyes %>% 
  na.omit() %>% 
  filter(dye == "RB") %>% 
  filter(!(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>% 
  ggplot(mapping = aes(x = mean_log_Kd, 
                       y = log_KF)) +
  geom_point(fill = "red",
             color = "red",
             size = 3) + 
  labs(x = TeX(r'($K_{d}~(L/kg)~dyes$)'), 
       y = TeX(r'($K_{F}~(L/kg)~PFOA~isotherms$)')) +
  sm_corr_theme() +
  sm_statCorr(corr_method = 'spearman',
              fit.params = list(color = 'red',
                                linetype = 'dashed'))
correlation_RB
ggsave(filename = "figs/correlation_RB.jpeg")

# linear model ----
lm_MB <- PFOA_dyes %>% 
  na.omit() %>% 
  filter(dye == "MB",
         !(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>%  
  lm(log_KF~mean_log_Kd,.)
summary(lm_MB)

lm_RB <- PFOA_dyes %>% 
  na.omit() %>% 
  filter(dye == "RB",
         !(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>%  
  lm(log_KF~mean_log_Kd,.)
summary(lm_RB)
# får ikke samme p-verdi her som i ggplot