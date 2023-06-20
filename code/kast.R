PFOA_dyes %>% 
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
        legend.position = "top_left",
        panel.grid = element_blank(),
        text = element_text(size = 14)
  ) +
  guides(fill=guide_legend(title=""),
         color=guide_legend(title=""),
         shape=guide_legend(title="")) +
  geom_errorbar(aes(ymin = log_KF-log_KF_se,
                    ymax = log_KF+log_KF_se),
                color = "black",
                width = .3,
                show.legend = F) +
  geom_line(aes(y=log_KF, group = sample_type), color = "black", linewidth = 2) +
  geom_point(aes(y=log_KF,
                 color = "PFOA"),
             size = 4,
             alpha = 0.5) + 
  scale_color_manual(values = c("blue", "#FB1894", "black"),
                     breaks = c("MB", "RB", "PFOA")) +
  scale_y_continuous(
    name = expression("log K"[F]~"PFOA"),
    sec.axis = sec_axis(~., name = expression("log K"[d]~"dyes"))
  )

PFOA_dyes %>% 
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
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
        legend.position = c(0.02, 0.99),
        legend.justification = c(0, 1),
        panel.grid = element_blank(),
        text = element_text(size = 14)
  ) +
  guides(fill=guide_legend(title=""),
         color=guide_legend(title=""),
         shape=guide_legend(title="")) +
  geom_errorbar(aes(ymin = log_KF-log_KF_se,
                    ymax = log_KF+log_KF_se),
                color = "black",
                width = .3,
                show.legend = T) +
  geom_line(aes(y=log_KF, group = sample_type), color = "black", linewidth = 2) +
  geom_point(aes(y=log_KF,
                 color = "PFOA"),
             size = 4,
             alpha = 0.5) + 
  scale_color_manual(values = c("blue", "#FB1894", "black"),
                     breaks = c("MB", "RB", "PFOA")) +
  scale_y_continuous(
    name = expression("log K"[F]~"PFOA"),
    sec.axis = sec_axis(~., name = expression("log K"[d]~"dyes"))
  )

library(ggplot2)

# Set theme for the plot
theme_set(theme_classic(base_size = 16))

# Plot the data points and regression line
ggplot(lm_MB_data_RWo, aes(x = mean_log_Kd, y = log_KF)) +
  geom_point(color = "blue") +
  geom_line(data = lm_MB_pred_RWo, aes(x = mean_log_Kd, y = fit), color = "red", linewidth = 1) +
  
  # Add the prediction interval ribbon
  geom_ribbon(data = lm_MB_pred_RWo, aes(x = mean_log_Kd, ymin = lwr, ymax = upr), fill = "grey80", alpha = 0.5) +
  
  # Add axis labels and title
  labs(x = expression(paste("log K"[d])), y = expression(paste("log K"[F])), 
       title = "") +
  
  # Customize axis and legend labels
  theme(axis.line = element_line(size = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  
  # Adjust the plot margins
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))