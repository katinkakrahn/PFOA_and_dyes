library(corrr)

# Multiple linear regression ----
PFOA_dyes_wider <- PFOA_dyes %>% 
  select(sample, K_F, dye, mean_log_Kd) %>% 
  pivot_wider(names_from = dye, values_from = mean_log_Kd, names_prefix = "log_Kd_") %>% 
  na.omit()

res <- lm(K_F ~ log_Kd_MB + log_Kd_RB,
          data = PFOA_dyes_wider)
summary(res)

res <- lm(K_F ~ log_Kd_MB * log_Kd_RB,
          data = PFOA_dyes_wider)

library(car)

avPlots(res)

# - The x-axis displays a single predictor variable and the y-axis displays the response variable.
# - The blue line shows the association between the predictor variable and the response variable, 
#   while holding the value of all other predictor variables constant.
# - The points that are labelled in each plot represent the 2 observations with the largest residuals 
#   and the 2 observations with the largest partial leverage.

# Correlation Spearman ----
PFOA_MB <- PFOA_dyes %>% 
  drop_na(c(log_KF, mean_log_Kd)) %>% 
  filter(dye == "MB")

PFOA_RB <- PFOA_dyes %>% 
  drop_na(c(log_KF, mean_log_Kd)) %>% 
  filter(dye == "RB")

cor.test(PFOA_MB$log_KF, PFOA_MB$mean_log_Kd, 
         method = "spearman", 
         alternative = "greater",
         exact = FALSE)
cor.test(PFOA_RB$log_KF, PFOA_RB$mean_log_Kd, 
         method = "spearman", 
         alternative = "greater",
         exact = FALSE)
cor.test(PFOA_dyes$K_F, PFOA_dyes$mean_log_Kd, 
         method = "spearman", 
         alternative = "greater",
         exact = FALSE)

# Raoul's sorption model ----
Sorption_isotherms_nonlinear <- PFOA_samplesonly %>% 
  filter(source == "sludge") %>% 
  group_by(sample, replicate) %>% 
  ggplot(aes(x = Cw_ugL, 
             y = Cs_ugkg)) +
  geom_point(mapping = aes(color = factor(sample)), 
             size = 3,
             alpha = 0.6) +
  facet_wrap(~sample,
             ncol = 2,
             scales = "free") +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  labs(x = TeX(r'($C_{w}~(\mu g/L)$)'), 
       y = TeX(r'($C_{s}~(\mu g/kg)$)'), 
       color = "") + 
  theme_bw()
Sorption_isotherms_nonlinear
ggsave(filename="figs/Sorption_isotherms_nonlinear_sludge.jpeg")

regression_glance <- PFOA_samplesonly %>% 
  filter(!(sample %in% c("CWC-BC-600", 
                         "WT-BC-600)",
                         "GW-BC-500",
                         "GW-BC-600",
                         "MS-BC-500",
                         "ULS-BC-800-40",
                         "VS-BC-600",
                         "CWC-BC-600",
                         "WT-BC-600"))) %>% 
  mutate(log_Cw = log10(Cw_ugL),
         log_Cs = log10(Cs_ugkg)) %>% 
  group_by(sample) %>%
  do(fit_isotherms = glance(lm(log_Cs ~ poly(log_Cw, 2), data = .))) %>% 
  unnest(fit_isotherms)

regression_tidy <- PFOA_samplesonly %>% 
  filter(!(sample %in% c("CWC-BC-600", 
                         "WT-BC-600)",
                         "GW-BC-500",
                         "GW-BC-600",
                         "MS-BC-500",
                         "ULS-BC-800-40",
                         "VS-BC-600",
                         "CWC-BC-600",
                         "WT-BC-600"))) %>% 
  mutate(log_Cw = log10(Cw_ugL),
         log_Cs = log10(Cs_ugkg)) %>% 
  group_by(sample) %>%
  do(fit_isotherms = tidy(lm(log_Cs ~ poly(log_Cw, 2), data = .))) %>% 
  unnest(fit_isotherms) %>% 
  pivot_wider(names_from = term,
              values_from = c(estimate, std.error, statistic, p.value))

regression_statistics_Raoul <- full_join(regression_glance, 
                                         regression_tidy)
write_xlsx(regression_statistics_Raoul, "data_manipulated/091222_regression_stats_Raoul.xlsx")
