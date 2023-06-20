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
library(tidytext)

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

Kd_MB_RB <- dye_summary %>% 
  na.omit() %>% 
  ggplot(aes(x = reorder(sample, mean_log_Kd), y = mean_log_Kd, group = dye, fill = dye)) +
  geom_bar(stat = "identity", 
           color = "black",
           position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_log_Kd-sd_log_Kd, 
                    ymax = mean_log_Kd+sd_log_Kd),
                width = .4,
                position = position_dodge(.9)) +
  labs(x = "", 
       y = "log Kd") +
  theme_bw() +
  scale_fill_manual(values = c("blue", "#FB1894")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=1)) +
  theme(plot.margin = margin(t = 0.5,
                             r = 0,
                             b = 0,
                             l = 0,
                             unit = "cm"))
Kd_MB_RB
ggsave(filename="figs/Kd_MB_RB.pdf")
ggsave(filename="figs/Kd_MB_RB.jpeg")

Kd_MB_RB_point <- dye_summary %>% 
  na.omit() %>% 
  ggplot(aes(x = reorder(sample, mean_log_Kd), 
             y = mean_log_Kd,
             color = dye,
             fill = dye,
             shape = dye)) +
  geom_errorbar(aes(ymin = mean_log_Kd-sd_log_Kd,
                    ymax = mean_log_Kd+sd_log_Kd),
                width = .3,
                show.legend = F) +
  geom_point(size = 4,
             alpha = 0.5) +
  scale_fill_manual(values = c("blue", "#FB1894"),
                    labels = c("methylene blue",
                               "rose bengal")) +
  scale_color_manual(values = c("blue", "#FB1894"),
                     labels = c("methylene blue",
                                "rose bengal")) +
  scale_shape_manual(values = c(21, 22),
                     labels = c("methylene blue",
                                "rose bengal")) +
  labs(x = "", 
       y = TeX(r'($log~K_{d}~(L/kg)$)')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=1),
        plot.margin = margin(t = 0.5,
                             r = 0,
                             b = 0,
                             l = 0,
                             unit = "cm"),
        legend.position = c(0.1, 0.9))
Kd_MB_RB_point
ggsave(filename="figs/Kd_MB_RB_point.pdf")
ggsave(filename="figs/Kd_MB_RB_point.jpeg")

# Kd PFOA and dyes together ----
Kd_PFOA_dyes <- read_excel("data_manipulated/200123_PFOA_dyes_long.xlsx")

Kd_PFOA_dyes_point <- Kd_PFOA_dyes %>% 
  na.omit() %>% 
  ggplot(aes(x = reorder(sample_name, log_Kd, max), 
             y = log_Kd,
             color = type,
             fill = type,
             shape = type)) +
  # geom_errorbar(aes(ymin = mean_log_Kd-sd_log_Kd,
  #                   ymax = mean_log_Kd+sd_log_Kd),
  #               width = .3,
  #               show.legend = F) +
  geom_point(size = 4,
             alpha = 0.5) +
  scale_fill_manual(values = c("blue", "#449e48", "#FB1894"),
                    labels = c("methylene blue",
                               "PFOA",
                               "rose bengal")) +
  scale_color_manual(values = c("blue", "#449e48", "#FB1894"),
                     labels = c("methylene blue",
                                "PFOA",
                                "rose bengal")) +
  scale_shape_manual(values = c(21, 22, 23),
                     labels = c("methylene blue",
                                "PFOA",
                                "rose bengal")) +
  labs(x = "", 
       y = TeX(r'($log~K_{d}~(L/kg)$)')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        plot.margin = margin(t = 0.5,
                             r = 0.5,
                             b = 0.5,
                             l = 0.5,
                             unit = "cm"),
        legend.position = c(0.1, 0.86)) +
  guides(fill=guide_legend(title=""),
         color=guide_legend(title=""),
         shape=guide_legend(title=""))
Kd_PFOA_dyes_point
ggsave(filename="figs/Kd_MB_RB_point.pdf")
ggsave(filename="figs/Kd_MB_RB_point.jpeg")

Kd_PFOA_dyes_plot <- Kd_PFOA_dyes %>% 
  ggplot(aes(x = reorder(sample, log_Kd, max))) +
  geom_col(aes(y=log_Kd, 
               fill = type), 
           color = "black", 
           position = "dodge") + 
  scale_fill_manual(values = c("#006ee6", "#449e48", "#FB1894")) +
  geom_errorbar(aes(ymin = log_Kd-sd_log_Kd, 
                    ymax = log_Kd+sd_log_Kd,
                    group = type),
                width = .4,
                position = position_dodge(1)) +
  labs(x = "") +
  scale_y_continuous(
    name = expression("log K"[F]~"PFOA"),
    sec.axis = sec_axis(~., name = expression("log K"[d]~"dyes"))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=1),
        axis.title.y = element_text(color = "black", size=11),
        axis.title.y.right = element_text(color = "black", size=11),
        plot.margin = margin(t = 0.5,
                             r = 0,
                             b = 0,
                             l = 0,
                             unit = "cm"))
Kd_PFOA_dyes_plot
ggsave(filename="figs/Kd_PFOA_dyes_plot.pdf")
ggsave(filename="figs/Kd_PFOA_dyes_plot.jpeg")


# Individual dyes ----
Kd_MB <- raw_data_summary %>% 
  na.omit() %>% 
  filter(dye == "MB") %>% 
  ggplot(mapping = aes(x = reorder(sample, mean_log_Kd), y = mean_log_Kd, group = dye, fill = dye)) +
  geom_bar(stat = "identity", 
           color = "black",
           position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_log_Kd-sd_log_Kd, 
                    ymax = mean_log_Kd+sd_log_Kd),
                width = .4,
                position = position_dodge(.9)) +
  labs(x = "", 
       y = "log Kd") +
  theme_bw() +
  scale_fill_manual(values = "blue") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=1)) +
  theme(plot.margin = margin(t = 0.5,
                             r = 0,
                             b = 0,
                             l = 0,
                             unit = "cm"))
Kd_MB
ggsave(filename="figs/Kd_MB.pdf")
ggsave(filename="figs/Kd_MB.jpeg")

Kd_RB <- raw_data_summary %>% 
  na.omit() %>% 
  filter(dye == "RB") %>% 
  ggplot(mapping = aes(x = reorder(sample, mean_log_Kd), y = mean_log_Kd, group = dye, fill = dye)) +
  geom_bar(stat = "identity", 
           color = "black",
           position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_log_Kd-sd_log_Kd, 
                    ymax = mean_log_Kd+sd_log_Kd),
                width = .4,
                position = position_dodge(.9)) +
  labs(x = "", 
       y = "log Kd") +
  theme_bw() +
  scale_fill_manual(values = "#FB1894") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=1)) +
  theme(plot.margin = margin(t = 0.5,
                             r = 0,
                             b = 0,
                             l = 0,
                             unit = "cm"))
Kd_RB
ggsave(filename="figs/Kd_RB.pdf")
ggsave(filename="figs/Kd_RB.jpeg")

           