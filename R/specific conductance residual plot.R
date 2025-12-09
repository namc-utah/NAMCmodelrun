Residual=setNames(as.data.frame(ranfor.mod$y-ranfor.mod$predicted),c("Residual"))
library(ggplot2)
p=ggplot(Residual, aes(x=Residual)) +
  geom_histogram(fill="white", color="black")+
  geom_vline(aes(xintercept=106.7, color="gray"),
             linetype="dashed",linewidth=1.5)+
  labs(title="",x=paste0("Observed - Predicted Specific Conductance μS/cm"), y = "Count")+
  theme_classic()+
  theme(legend.position = "none",text = element_text(size = 14))

windows()
p
