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


# Raw data dyes ----
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


# dye summary statistics ----
dye_summary <- dyes %>% 
  filter(comment == "good") %>% 
  group_by(sample, sample_ID, sample_name, dye, sample_type, date_meas) %>% 
  summarise(Cw_mean = mean(Cw_mgL),
            Cs_mean = mean(Cs_mgkg),
            Kd_mean = mean(Cs_mgkg / Cw_mgL),
            mean_log_Kd = mean(log_Kd),
            sd_log_Kd = sd(log_Kd),
            n = n())

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
         log_KF, n_F, log_KF_se, n_F_se, nobs) %>% 
  mutate(log_Kd_1mgL = log10(10^log_KF*1^n_F/1)-3*(1-n_F),
         sd_logKd_1mgL = sqrt(((log_KF_se))/log_KF)^2+n_F_se^2,
         log_Kd_1ugL = log10(10^log_KF*1^n_F/1)-3*(1-n_F),
         sd_logKd_1ugL = sqrt(((log_KF_se))/log_KF)^2+n_F_se^2)

write_xlsx(PFOA_regression_statistics, "data_manipulated/181122_Freundlich_regression_statistics.xlsx")
PFOA_regression_statistics_sludge <- PFOA_regression_statistics %>% 
  filter(source == "sludge")

PFOA_regression_statistics %>% 
  filter(source == "sludge") %>% 
  write_xlsx("data_manipulated/Freundlich_regression_statistics_sludge.xlsx")

# join dyes and isotherms ----
PFOA_dyes <- full_join(PFOA_regression_statistics, dye_summary)

PFOA_dyes %>% 
  filter(sample_type == "sample") %>% 
  select(sample, log_KF, log_KF_se, log_Kd_1mgL, mean_log_Kd, sd_log_Kd, dye) %>%
  write_xlsx("data_manipulated/111222_PFOA_dyes.xlsx")

# Sorption isotherms non-linear ----
Sorption_isotherms_nonlinear <- PFOA_samplesonly %>% 
  group_by(sample, replicate) %>% 
  filter(!(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>% 
  ggplot() +
  geom_point(mapping = aes(x = Cw_ugL, 
                           y = Cs_ugkg, 
                           color = factor(sample)
  ), 
  size = 3,
  alpha = 0.6) +
  geom_smooth(mapping = aes(x = Cw_ugL, 
                            y = Cs_ugkg, 
                            color = factor(sample)), 
              formula = y ~ x, 
              method=lm, 
              se=F, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~temperature,
             #scales = "free_x",
             ncol = 2) +
  theme_bw()
Sorption_isotherms_nonlinear
ggsave(filename="figs/Sorption_isotherms_nonlinear.jpeg")

Sorption_isotherms_nonlinear_sludge <- PFOA_samplesonly %>% 
  filter(source == "sludge") %>% 
  group_by(sample, replicate) %>% 
  ggplot() +
  geom_point(mapping = aes(x = Cw_ugL, 
                           y = Cs_ugkg, 
                           color = factor(sample)
  ), 
  size = 3,
  alpha = 0.6) +
  geom_smooth(mapping = aes(x = Cw_ugL, 
                            y = Cs_ugkg, 
                            color = factor(sample)), 
              formula = y ~ x, 
              method=lm, 
              se=F, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~temperature,
             scales = "free_x",
             ncol = 2) +
  theme_bw()
Sorption_isotherms_nonlinear_sludge
ggsave(filename="figs/Sorption_isotherms_nonlinear_sludge.jpeg")

# Sorption isotherms Freundlich ----
temp.labs <- c("500 \u00B0C", "600 \u00B0C", "700 \u00B0C", "750 \u00B0C", "800 \u00B0C")
names(temp.labs) <- c("500", "600", "700", "750", "800")

Sorption_isotherms_Freundlich <- PFOA_samplesonly %>% 
  group_by(sample_ID, replicate) %>% 
  filter(!(sample %in% c("CWC-BC-600", "WT-BC-600", "GW-MAP"))) %>% 
  mutate(biochar = factor(biochar, 
                          levels = c("DWSS", "DSS-1", "DSS-2", "LSS", 
                                     "FWR",
                                     "WT", "GW", "CWC")
  )) %>% 
  ggplot(aes(x = log10(Cw_ugL), 
             y = log10(Cs_ugkg), 
             color = factor(biochar)
  )) +
  geom_point( 
    size = 3,
    alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  geom_smooth(formula = y ~ x, 
              method=lm, 
              se=F, 
              fullrange = F) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~temperature,
             #scales = "free_x",
             ncol = 2,
             labeller = labeller(temperature = temp.labs)) +
  theme_bw() +
  guides(color=guide_legend(ncol=2)) +
  theme(panel.grid = element_blank(),
        legend.position=c(0.75,0.18)) +
  scale_color_brewer(palette = "Dark2")
Sorption_isotherms_Freundlich
ggsave(filename="figs/Sorption_isotherms_Freundlich.jpeg")

Sorption_isotherms_Freundlich_sludge <- PFOA_samplesonly %>% 
  filter(source == "sludge") %>% 
  group_by(sample, replicate) %>% 
  ggplot() +
  geom_point(mapping = aes(x = log10(Cw_ugL), 
                           y = log10(Cs_ugkg), 
                           color = factor(biochar)
  ), 
  size = 3,
  alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  geom_smooth(mapping = aes(x = log10(Cw_ugL), 
                            y = log10(Cs_ugkg), 
                            color = factor(biochar)), 
              formula = y ~ x, 
              method=lm, 
              se=F, 
              fullrange = FALSE) + 
  labs(x = TeX(r'($log~C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($log~C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  facet_wrap(~temperature,
             #scales = "free_x",
             ncol = 2) +
  sm_corr_theme()
Sorption_isotherms_Freundlich_sludge
ggsave(filename="figs/Sorption_isotherms_Freundlich_sludge.jpeg")

# correlation plot Kd dyes and KF isotherms ----
correlation <- PFOA_dyes %>% 
  na.omit() %>% 
  group_by(sample, replicate, dye) %>% 
  filter(!(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>% 
  ggplot(mapping = aes(x = mean_log_Kd, 
                       y = log_KF, 
                       color = factor(dye),
                       shape = factor(temperature))) +
  geom_point(size = 3,
             alpha = 0.8) + 
  labs(x = TeX(r'($K_{d}~(L/kg)~dyes$)'), 
       y = TeX(r'($K_{F}~(L/kg)~PFOA~isotherms$)'), 
       color = "dye",
       shape = "pyrolysis temperature") +
  sm_corr_theme()
correlation
ggsave(filename = "figs/correlation.jpeg")  

correlation2 <- PFOA_dyes %>% 
  na.omit() %>% 
  group_by(sample, replicate, dye, biochar) %>% 
  filter(!(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>% 
  ggplot(mapping = aes(x = mean_log_Kd, 
                       y = log_KF, 
                       color = factor(dye),
                       size = factor(biochar),
                       alpha = factor(biochar))) +
  geom_point() + 
  labs(x = TeX(r'($K_{d}~(L/kg)~dyes$)'), 
       y = TeX(r'($K_{F}~(L/kg)~PFOA~isotherms$)'), 
       color = "dye",
       size = "biochar",
       alpha = "biochar") +
  sm_corr_theme()
correlation2
ggsave(filename = "figs/correlation2.jpeg")  

correlation3 <- PFOA_dyes %>% 
  na.omit() %>% 
  group_by(sample, replicate, dye, biochar) %>% 
  filter(!(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>% 
  ggplot(mapping = aes(x = mean_log_Kd, 
                       y = log_KF,
                       color = factor(dye))) +
  geom_point() + 
  geom_smooth() +
  labs(x = TeX(r'($K_{d}~(L/kg)~dyes$)'), 
       y = TeX(r'($K_{F}~(L/kg)~PFOA~isotherms$)'), 
       color = "dye",
       size = "biochar",
       alpha = "biochar") +
  sm_corr_theme()
correlation3
ggsave(filename = "figs/correlation3.jpeg")  

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

filter_PFOA_KF <- PFOA_regression_statistics %>% 
  drop_na(p.value) %>% 
  select(sample, sample_ID, sample_name, biochar, temperature, residence_time, log_KF, log_KF_se, log_Kd_1mgL, sd_logKd_1mgL)

# PFOA, dyes line chart ----
PFOA_dyes_plot <- merge(dye_summary, filter_PFOA_KF, by = c("sample","sample_ID","sample_name"), all = TRUE)

PFOA_dyes_plot2 <- PFOA_dyes_plot %>% 
  ggplot(aes(x = reorder(sample_name, log_KF, max))) +
  geom_point(aes(y=mean_log_Kd, color = dye),
             size = 4,
             alpha = 0.5) +
  geom_line(aes(y=mean_log_Kd, group = dye, color = dye)) +
  geom_errorbar(aes(ymin = mean_log_Kd-sd_log_Kd,
                    ymax = mean_log_Kd+sd_log_Kd,
                    color = dye),
                width = .3,
                show.legend = F) +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        plot.margin = margin(t = 0.5,
                             r = 0.5,
                             b = 0.5,
                             l = 0.5,
                             unit = "cm"),
        legend.position = c(0.05, 0.85),
        panel.grid = element_blank()) +
  guides(fill=guide_legend(title=""),
         color=guide_legend(title=""),
         shape=guide_legend(title="")) +
  geom_errorbar(aes(ymin = log_KF-log_KF_se,
                    ymax = log_KF+log_KF_se),
                color = "black",
                width = .3,
                show.legend = F) +
  geom_line(aes(y=log_KF, group = sample_type), color = "black", size = 2) +
  geom_point(aes(y=log_KF,
                 color = "PFOA"),
             size = 4,
             alpha = 0.5) + 
  scale_color_manual(values = c("blue", "#FB1894", "black"),
                     breaks = c("MB", "RB", "PFOA")) +
  scale_y_continuous(
    name = expression("log K"[F]~"PFOA"),
    sec.axis = sec_axis(~., name = expression("log K"[d]~"dyes")))
PFOA_dyes_plot2
ggsave(filename="figs/PFOA_dyes_plot2.jpeg")  

PFOA_dyes_plot3 <- PFOA_dyes_plot %>% 
  ggplot(aes(x = reorder(sample_name, log_Kd_1mgL, max))) +
  geom_point(aes(y=mean_log_Kd, color = dye),
             size = 2,
             alpha = 0.5) +
  geom_point(aes(y=log_Kd_1mgL), 
             color = "black",
             size = 2,
             alpha = 0.5) + 
  geom_line(aes(y=mean_log_Kd, group = dye, color = dye)) + 
  geom_errorbar(aes(ymin = mean_log_Kd-sd_log_Kd,
                    ymax = mean_log_Kd+sd_log_Kd,
                    color = dye),
                width = .3,
                show.legend = F) +
  geom_errorbar(aes(ymin = log_Kd_1mgL-sd_logKd_1mgL,
                    ymax = log_Kd_1mgL+sd_logKd_1mgL),
                color = "grey45",
                width = .3,
                show.legend = F) +
  scale_color_manual(values = c("blue", "#FB1894")) +
  geom_line(aes(y=log_Kd_1mgL, group = sample_type), color = "black", linewidth = 1) +
  scale_y_continuous(name = expression("log K"[d])) +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.margin = margin(t = 0.5,
                             r = 0.5,
                             b = 0,
                             l = 0.7,
                             unit = "cm"),
        legend.position = c(0.05, 0.88)) #+
  #guides(color=guide_legend(title=""))
PFOA_dyes_plot3
ggsave(filename="figs/PFOA_dyes_plot3.pdf")
ggsave(filename="figs/PFOA_dyes_plot3.jpeg")   
