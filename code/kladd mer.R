filter_PFOA_KF <- PFOA_regression_statistics %>% 
  drop_na(p.value) %>% 
  select(sample, sample_ID, sample_name, biochar, temperature, residence_time, log_KF, log_KF_se, log_Kd_1mgL, sd_logKd_1mgL)

PFOA_dyes_plot <- merge(dye_summary, filter_PFOA_KF, by = c("sample","sample_ID","sample_name"), all = TRUE)


PFOA_dyes_plot2 <- PFOA_dyes_plot %>% 
  ggplot(aes(x = reorder(sample_name, log_KF, max))) +
  geom_point(aes(y=log_KF), 
           color = "red",
           size = 4,
           alpha = 0.5) + 
  geom_point(aes(y=mean_log_Kd, color = dye),
             size = 4,
             alpha = 0.5) +
  geom_line(aes(y=mean_log_Kd, group = dye, color = dye)) + 
  scale_color_manual(values = c("blue", "#FB1894")) +
  geom_line(aes(y=log_KF, group = sample_type), color = "red", size = 2) +
  scale_y_continuous(
    name = expression("log K"[F]~"PFOA"),
    sec.axis = sec_axis(~., name = expression("log K"[d]~"dyes"))) +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        plot.margin = margin(t = 0.5,
                             r = 0.5,
                             b = 0.5,
                             l = 0.5,
                             unit = "cm"),
        legend.position = c(0.05, 0.89)) +
  guides(fill=guide_legend(title=""),
         color=guide_legend(title=""),
         shape=guide_legend(title=""))
PFOA_dyes_plot2
ggsave(filename="figs/PFOA_dyes_plot2.pdf")
ggsave(filename="figs/PFOA_dyes_plot2.jpeg")  



breaks_ID <- c("CWC-F", "CWC","CWC-500" , "CWC-600" , "CWC-700" , "CWC-750" , 
               "WT-F", "WT","WT-500" , "WT-600" , "WT-700",  "WT-800", 
               "GW-F", "GW","GW-500", "GW-600",  "GW-800", 
               "DSS-1-F", "DSS-1","DSS-1-500","DSS-1-600" ,"DSS-1-700", "DSS-1-800", 
               "DSS-2-F", "DSS-2","DSS-2-500", "DSS-2-600" ,"DSS-2-700", "DSS-2-800" ,
               "LSS-F","LSS","LSS-600" , "LSS-750" ,
               "FWR-F", "FWR","FWR-600","FWR-800",
               "DWSS")

values_ID <- c('#8c510a','#8c510a','#8c510aB3','#8c510a99','#8c510a66','#8c510a4D',
               '#d8b365','#d8b365','#d8b365B3','#d8b36599','#d8b36566','#d8b3654D',
               '#a1d76a','#a1d76a','#a1d76aB3','#a1d76a99',            '#a1d76a4D',
               '#2166ac','#2166ac','#2166acB3','#2166ac99','#2166ac66','#2166ac4D',
               '#5ab4ac','#5ab4ac','#5ab4acB3','#5ab4ac99','#5ab4ac66','#5ab4ac4D',
               '#1b7837','#1b7837',            '#1b783799',            '#1b78374D',
               '#762a83','#762a83',            '#762a8399',            '#762a834D',
               '#01665e')


breaks_ID <- c("CWC", "GW","WT", 
               "FWR",
               "DSS-1","DSS-2","LSS", 
               "DWSS")

values_ID <- c('#8c510a', '#d8b365','#a1d76a', 
               '#762a83',
               '#5ab4ac','#01665e','#1b7837',
               '#2166ac')

temp.labs <- c("500 \u00B0C", "600 \u00B0C", "700 \u00B0C", "750 \u00B0C", "800 \u00B0C", "MAP")
names(temp.labs) <- c("500", "600", "700", "750", "800", "MAP")

PFOA_samplesonly %>% 
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
    size = 2,
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
             ncol = 2,
             labeller = labeller(temperature = temp.labs)) +
  theme_bw() +
  scale_color_manual(values = values_ID,
                    breaks = breaks_ID,
                    name = "") +
  #guides(color=guide_legend(ncol=2)) +
  theme(panel.grid = element_blank(),
        legend.position="bottom")
ggsave(filename="figs/Sorption_isotherms_Freundlich_newcolors.jpeg")


breaks_ID <- c("raw sludge", 
               "digested sludge", 
               "food waste reject", 
               "wood-based")

values_ID <- c('#1b7837',
               '#5ab4ac',
               '#2166ac',
               '#d8b365')


PFOA_samplesonly %>% 
  group_by(sample_ID, replicate) %>% 
  filter(!(sample %in% c("CWC-BC-600", "WT-BC-600", "GW-MAP"))) %>% 
  mutate(source2 = factor(source2, 
                          levels = c("raw sludge", "digested sludge", "food waste reject", "wood-based")
  )) %>% 
  ggplot(aes(x = log10(Cw_ugL), 
             y = log10(Cs_ugkg), 
             color = factor(source2),
             group = factor(biochar)
  )) +
  geom_point( 
    size = 2,
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
             ncol = 2,
             labeller = labeller(temperature = temp.labs)) +
  theme_bw() +
  scale_color_manual(values = values_ID,
                     breaks = breaks_ID,
                     name = "") +
  theme(panel.grid = element_blank(),
        legend.position="bottom")
ggsave(filename="figs/Sorption_isotherms_Freundlich_source.jpeg")
























