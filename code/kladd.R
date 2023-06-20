lm_RB_RWo_ranked <- lm_RB_data_RWo %>% 
  lm(rank(log_KF) ~ rank(mean_log_Kd),.)

summary(lm_RB_RWo_ranked)

lm_RB_data_RWo_ranked <- PFOA_dyes %>%
  filter(dye == "RB", !(sample %in% c("CWC-BC-600", "WT-BC-600"))) %>%
  select(log_KF, mean_log_Kd) %>%
  na.omit() %>%
  mutate(rank_log_KF = rank(log_KF),
         rank_mean_log_Kd = rank(mean_log_Kd))

lm_RB_pred_RWo_ranked <- predict(lm_RB_data_RWo_ranked, 
                          newdata = data.frame("rank_mean_log_Kd" = seq(min(lm_RB_data_RWo_ranked$rank_mean_log_Kd), 
                                                                   max(lm_RB_data_RWo_ranked$rank_mean_log_Kd), 
                                                                   length.out = 100)),
                          interval = "prediction",
                          level = 0.95) |> 
  as.data.frame() |> 
  cbind("rank_mean_log_Kd" = seq(min(lm_RB_data_RWo_ranked$rank_mean_log_Kd), 
                            max(lm_RB_data_RWo_ranked$rank_mean_log_Kd), 
                            length.out = 100))


lm_RB_spearman <- cor.test(x = lm_RB_data_RWo$mean_log_Kd, 
                           y = lm_RB_data_RWo$log_KF, 
                           method = "spearman")

rho_RB <- lm_RB_spearman$estimate
pval_spearman_RB <- lm_RB_spearman$p.value

lm_RB_data_RWo_ranked %>% 
  ggplot() +
  geom_point(aes(x = rank_mean_log_Kd, y = rank_log_KF), alpha = 0.5, size = 2) +
  geom_ribbon(aes(x = rank_mean_log_Kd, ymin = lwr, ymax = upr), alpha = 0.3, data = lm_RB_pred_RWo_ranked) +
  geom_line(aes(x = rank_mean_log_Kd, y = fit), data = lm_RB_pred_RWo_ranked, linewidth = 2, color = "#FB1894") +
  labs(x = expression(paste("log K"[d]," dye")), y = expression(paste("ranked log K"[F]," PFOA")), 
       title = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14)) +
  annotate("text", x = min(lm_RB_data_RWo$mean_log_Kd) + 0, y = max(lm_RB_data_RWo$log_KF) + 1.2,
           label = paste("R-squared = ", round(rho_RB, 3)), size = 4, hjust = 0) +
  annotate("text", x = min(lm_RB_data_RWo$mean_log_Kd) + 0, y = max(lm_RB_data_RWo$log_KF) + 0.9,
           label = paste("p-value = ", round(pval_spearman_RB, 5)), size = 4, hjust = 0)































## ---- MB
lm_MB_data_ranked <- PFOA_dyes %>% 
  filter(dye == "MB",
         !(sample %in% c("CWC-BC-600", "WT-BC-600)"))) %>%  
  select(log_KF, mean_log_Kd) |> 
  na.omit() %>% 
  mutate(rank_log_KF = rank(log_KF),
         rank_mean_log_Kd = rank(mean_log_Kd))

lm_MB_ranked <- lm_MB_data_ranked %>% 
  lm(rank(log_KF)~rank(mean_log_Kd),.)

summary(lm_MB_ranked)

lm_MB_pred_ranked <- predict(lm_MB_ranked, 
                          newdata = data.frame("mean_log_Kd" = seq(min(lm_MB_data_ranked$rank_mean_log_Kd), 
                                                                   max(lm_MB_data_ranked$rank_mean_log_Kd), 
                                                                   length.out = 100)),
                          interval = "prediction",
                          level = 0.95) |> 
  as.data.frame() |> 
  cbind("mean_log_Kd" = seq(min(lm_MB_data_ranked$rank_mean_log_Kd), 
                            max(lm_MB_data_ranked$rank_mean_log_Kd), 
                            length.out = 100))

lm_MB_spearman <- cor.test(x = lm_MB_data_RWo$mean_log_Kd, 
                           y = lm_MB_data_RWo$log_KF, 
                           method = "spearman")

rho_MB <- lm_RB_spearman$estimate
pval_spearman_MB <- lm_RB_spearman$p.value

lm_MB_data_ranked %>% 
  ggplot() +
  geom_point(aes(x = rank_mean_log_Kd, y = rank_log_KF), alpha = 0.5, size = 2) +
  geom_ribbon(aes(x = mean_log_Kd, ymin = lwr, ymax = upr), alpha = 0.3, data = lm_MB_pred_ranked) +
  geom_line(aes(x = mean_log_Kd, y = fit), data = lm_MB_pred_ranked, linewidth = 2, color = "blue") +
  # scale_x_continuous(breaks = seq(0, max(lm_MB_data_ranked$rank_mean_log_Kd), by = 1),
  #                    limits = c(0, max(lm_MB_data_ranked$rank_mean_log_Kd))) +
  labs(x = expression(paste("log K"[d]," dye")), y = expression(paste("log K"[F]," PFOA")), 
       title = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14)) +
  annotate("text", x = min(lm_MB_data_ranked$rank_mean_log_Kd) + 0, y = max(lm_MB_data_ranked$rank_log_KF) + 0.3,
           label = paste("R-squared = ", round(rho_MB, 3)), size = 4, hjust = 0) +
  annotate("text", x = min(lm_MB_data_ranked$rank_mean_log_Kd) + 0, y = max(lm_MB_data_ranked$rank_log_KF) - 0,
           label = paste("p-value = ", round(pval_spearman_MB, 5)), size = 4, hjust = 0)
ggsave("figs/lm_MB.jpeg")

cor.test(x = lm_MB_data_RWo$mean_log_Kd, 
         y = lm_MB_data_RWo$log_KF, 
         method = "spearman")

lm_MB_data_RWo |> 
  lm(rank(log_KF) ~ rank(mean_log_Kd), data = _) |> 
  summary()