pHcond <- read_excel("R/data_raw/pH_conductivity_R.xlsx")
as.data.table(pHcond)
pHcond <- as.data.table(pHcond)

#Mean and standard deviation pH and cond
pHcondSummary <- pHcond[, .(mean_ph = mean(pH), 
                            mean_cond = mean(Conductivity),
                            se_cond = std.error(Conductivity),
                            se_ph = std.error(pH)
                            ),
                        keyby = .(Sample)]

pHcondSummary <- pHcondSummary[order(mean_ph),]
pHcondSummary <- pHcondSummary[,c(1,2,5,3,4)]
target <- c("ULS", "BRL", "CWC", "ULS+S", "BRL+S", "CWC+S", "S")

# Can't remember what this code did and now it's not working...
#pHcondSummary <- pHcondSummary[match(target, pHcondsummary$Sample),]

pHcond_table <- kable(pHcondSummary, "latex", digits = 2, booktabs = TRUE)
pHcond_table

# pH plot ----
pH <- ggplot(data = pHcondSummary, aes(x = reorder(Sample, mean_ph), y = mean_ph)) + 
  geom_point()+ 
  geom_errorbar(aes(ymin=mean_ph-se_ph, ymax=mean_ph+se_ph), width=.2,
                       position=position_dodge(.9))+ 
  labs(x = "", y = "pH") + 
  theme_bw()
pH
ggsave(filename="R/figs/pH.pdf")

# conductivity plot ----
cond <- ggplot(data = pHcondSummary, aes(x = reorder(Sample,mean_cond), y = mean_cond)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=mean_cond-se_cond, ymax=mean_cond+se_cond), width=.2,
                         position=position_dodge(.9)) + 
  labs(x = "", y = TeX(r'(Conductivity $(\mu S/cm)$)')) + 
  theme_bw()
cond
ggsave(filename = "R/figs/conductivity.pdf")


