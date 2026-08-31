library(ggplot2)
library(ggpubr)

# -----------------------------------------------------------------------------
# Combine IFN and NFkB figures in a grid
# -----------------------------------------------------------------------------

p1 <- ifn_violin_plot
p2 <- seap_violin_plot

combined_plot <- ggarrange(p1, p2, ncol=1, common.legend = FALSE, labels = c("A", "B")) 

ggsave("figures/combined_ifn_nfkb_bioactivity_plots.png", combined_plot, width=30, height=15, units=c("cm"))
