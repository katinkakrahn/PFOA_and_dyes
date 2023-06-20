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
library(patchwork)

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

write_xlsx(PFOA_regression_statistics, "data/dye_PFOA/181122_Freundlich_regression_statistics.xlsx")

# join dyes and isotherms ----
PFOA_dyes <- full_join(PFOA_regression_statistics, dye_summary)

PFOA_dyes %>% 
  filter(sample_type == "sample") %>% 
  # select(sample, log_KF, log_KF_se, log_Kd_1mgL, mean_log_Kd, sd_log_Kd, dye) %>%
  select(sample, log_KF, log_KF_se, mean_log_Kd, sd_log_Kd, dye) %>%
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

################################################################################
################################################################################
################################################################################

# ---- RWo ----
# careful with the use of na.omit() too early in your pipelines!
# https://lindeloev.github.io/tests-as-linear/
# Spearman correlation != "normal" linear regression
#    rank(y) ~ rank(x) != y ~ x

## ---- MB
lm_MB_data_RWo <- PFOA_dyes %>% 
  filter(dye == "MB",
         !(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>%  
  select(log_KF, mean_log_Kd) |> 
  na.omit()

lm_MB_RWo <- lm_MB_data_RWo %>% 
  lm(log_KF~mean_log_Kd,.)

summary(lm_MB_RWo)

lm_MB_pred_RWo <- predict(lm_MB_RWo, 
                          newdata = data.frame("mean_log_Kd" = seq(min(lm_MB_data_RWo$mean_log_Kd), 
                                                                   max(lm_MB_data_RWo$mean_log_Kd), 
                                                                   length.out = 100)),
                          interval = "prediction",
                          level = 0.95) |> 
  as.data.frame() |> 
  cbind("mean_log_Kd" = seq(min(lm_MB_data_RWo$mean_log_Kd), 
                            max(lm_MB_data_RWo$mean_log_Kd), 
                            length.out = 100))

rsq <- summary(lm_MB_RWo)$r.squared
pval <- summary(lm_MB_RWo)$coefficients[2, 4]

lm_MB_data_RWo %>% 
  ggplot() +
  geom_point(aes(x = mean_log_Kd, y = log_KF), alpha = 0.5, size = 2) +
  geom_ribbon(aes(x = mean_log_Kd, ymin = lwr, ymax = upr), alpha = 0.3, data = lm_MB_pred_RWo) +
  geom_line(aes(x = mean_log_Kd, y = fit), data = lm_MB_pred_RWo, linewidth = 2, color = "blue") +
  labs(x = expression(paste("log K"[d]," dye")), y = expression(paste("log K"[F]," PFOA")), 
       title = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  annotate("text", x = min(lm_MB_data_RWo$mean_log_Kd) + 0, y = max(lm_MB_data_RWo$log_KF) + 0.4,
           label = paste("R-squared = ", round(rsq, 3)), size = 6, hjust = 0) +
  annotate("text", x = min(lm_MB_data_RWo$mean_log_Kd) + 0, y = max(lm_MB_data_RWo$log_KF) - 0,
           label = paste("p-value = ", round(pval, 5)), size = 6, hjust = 0) +
  annotate("text", x = max(lm_MB_data_RWo$mean_log_Kd), y = min(lm_MB_data_RWo$log_KF) - 0.8,
           label = expression(bold("a")), size = 10, hjust = 1, vjust = 1)
ggsave("figs/lm_MB.jpeg")

# ---- RB
lm_RB_data_RWo <- PFOA_dyes %>% 
  filter(dye == "RB",
         !(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>%  
  select(log_KF, mean_log_Kd) |> 
  na.omit()

lm_RB_RWo <- lm_RB_data_RWo %>% 
  lm(log_KF~mean_log_Kd,.)

summary(lm_RB_RWo)

lm_RB_pred_RWo <- predict(lm_RB_RWo, 
                          newdata = data.frame("mean_log_Kd" = seq(min(lm_RB_data_RWo$mean_log_Kd), 
                                                                   max(lm_RB_data_RWo$mean_log_Kd), 
                                                                   length.out = 100)),
                          interval = "prediction",
                          level = 0.95) |> 
  as.data.frame() |> 
  cbind("mean_log_Kd" = seq(min(lm_RB_data_RWo$mean_log_Kd), 
                            max(lm_RB_data_RWo$mean_log_Kd), 
                            length.out = 100)) %>% 
  mutate(delta_upr = upr-fit,
         delta_lwr =fit-lwr)

rsq_RB <- summary(lm_RB_RWo)$r.squared
pval_RB <- summary(lm_RB_RWo)$coefficients[2, 4]

lm_RB_data_RWo %>% 
  ggplot() +
  geom_point(aes(x = mean_log_Kd, y = log_KF), alpha = 0.5, size = 2) +
  geom_ribbon(aes(x = mean_log_Kd, ymin = lwr, ymax = upr), alpha = 0.3, data = lm_RB_pred_RWo) +
  geom_line(aes(x = mean_log_Kd, y = fit), data = lm_RB_pred_RWo, linewidth = 2, color = "#FB1894") +
  labs(x = expression(paste("log K"[d]," dye")), y = expression(paste("log K"[F]," PFOA")), 
       title = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  annotate("text", x = min(lm_RB_data_RWo$mean_log_Kd) + 0, y = max(lm_RB_data_RWo$log_KF) + 1.4,
           label = paste("R-squared = ", round(rsq_RB, 3)), size = 6, hjust = 0) +
  annotate("text", x = min(lm_RB_data_RWo$mean_log_Kd) + 0, y = max(lm_RB_data_RWo$log_KF) + 0.9,
           label = paste("p-value = ", round(pval_RB, 5)), size = 6, hjust = 0) +
  annotate("text", x = max(lm_MB_data_RWo$mean_log_Kd), y = min(lm_MB_data_RWo$log_KF) - 0.2,
           label = expression(bold("b")), size = 10, hjust = 1, vjust = 1)
ggsave("figs/lm_RB.jpeg")

####

# Reconstruction of Spearman correlations (based on wrong data subset)
# https://lindeloev.github.io/tests-as-linear/

# wrong_data_MB <- PFOA_dyes %>% 
#   na.omit() %>% 
#   filter(dye == "MB") %>% 
#   filter(!(sample %in% c("CWC-BC-600", "WT-BC-600)")))

cor.test(x = lm_MB_data_RWo$mean_log_Kd, 
         y = lm_MB_data_RWo$log_KF, 
         method = "spearman")

lm_MB_data_RWo |> 
  lm(rank(log_KF) ~ rank(mean_log_Kd), data = _) |> 
  summary()
          

# wrong_data_RB <- PFOA_dyes %>% 
#   na.omit() %>% 
#   filter(dye == "RB") %>% 
#   filter(!(sample %in% c("CWC-BC-600", "WT-BC-600)")))

cor.test(x = lm_RB_data_RWo$mean_log_Kd, 
         y = lm_RB_data_RWo$log_KF, 
         method = "spearman")

lm_RB_data_RWo |> 
  lm(rank(log_KF) ~ 1 + rank(mean_log_Kd), data = _) |> 
  summary()

# plotting ----
plot1 <- lm_MB_data_RWo %>% 
  ggplot() +
  geom_point(aes(x = mean_log_Kd, y = log_KF), alpha = 0.5, size = 2) +
  geom_ribbon(aes(x = mean_log_Kd, ymin = lwr, ymax = upr), alpha = 0.3, data = lm_MB_pred_RWo) +
  geom_line(aes(x = mean_log_Kd, y = fit), data = lm_MB_pred_RWo, linewidth = 2, color = "blue") +
  labs(x = expression(paste("log K"[d]," dye")), y = expression(paste("log K"[F]," PFOA")), 
       title = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  annotate("text", x = min(lm_MB_data_RWo$mean_log_Kd) + 0, y = max(lm_MB_data_RWo$log_KF) + 0.4,
           label = paste("R-squared = ", round(rsq, 3)), size = 6, hjust = 0) +
  annotate("text", x = min(lm_MB_data_RWo$mean_log_Kd) + 0, y = max(lm_MB_data_RWo$log_KF) - 0,
           label = paste("p-value = ", round(pval, 5)), size = 6, hjust = 0) +
  annotate("text", x = max(lm_MB_data_RWo$mean_log_Kd), y = min(lm_MB_data_RWo$log_KF) - 0.8,
           label = expression(bold("a")), size = 10, hjust = 1, vjust = 1)

plot2 <- lm_RB_data_RWo %>% 
  ggplot() +
  geom_point(aes(x = mean_log_Kd, y = log_KF), alpha = 0.5, size = 2) +
  geom_ribbon(aes(x = mean_log_Kd, ymin = lwr, ymax = upr), alpha = 0.3, data = lm_RB_pred_RWo) +
  geom_line(aes(x = mean_log_Kd, y = fit), data = lm_RB_pred_RWo, linewidth = 2, color = "#FB1894") +
  labs(x = expression(paste("log K"[d]," dye")), y = expression(paste("log K"[F]," PFOA")), 
       title = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  annotate("text", x = min(lm_RB_data_RWo$mean_log_Kd) + 0, y = max(lm_RB_data_RWo$log_KF) + 1.4,
           label = paste("R-squared = ", round(rsq_RB, 3)), size = 6, hjust = 0) +
  annotate("text", x = min(lm_RB_data_RWo$mean_log_Kd) + 0, y = max(lm_RB_data_RWo$log_KF) + 0.9,
           label = paste("p-value = ", round(pval_RB, 5)), size = 6, hjust = 0) +
  annotate("text", x = max(lm_MB_data_RWo$mean_log_Kd), y = min(lm_MB_data_RWo$log_KF) - 0.2,
           label = expression(bold("b")), size = 10, hjust = 1, vjust = 1)

grid.arrange(plot1, plot2, ncol=2)
