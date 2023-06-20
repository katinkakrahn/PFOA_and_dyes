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

raw_data <- read_excel("/VOW_fargetester_total_041122.xlsx")

raw_data_summary <- raw_data %>% 
  group_by(sample, dye) %>% 
  summarise(Cw_mean = mean(Cw_mgL),
            Cs_mean = mean(Cs_mgkg),
            Kd_mean = mean(Kd_Lkg),
            mean_log_Kd = mean(log_Kd))

raw_data_summary %>% 
  na.omit() %>% 
  ggplot() +
  geom_bar(mapping = aes(x = reorder(sample, mean_log_Kd), y = mean_log_Kd, group = dye, fill = dye), 
           stat = "identity", 
           position = "dodge") +
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
           