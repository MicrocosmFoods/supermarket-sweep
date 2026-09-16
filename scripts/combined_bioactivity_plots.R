library(ggplot2)
library(ggpubr)

# -----------------------------------------------------------------------------
# combine viability unfiltered and filtered plots in a grid
# -----------------------------------------------------------------------------

p1 <- unfiltered_viability_violin_plot
p2 <- filtered_viability_violin_plot

combined_plot <- ggarrange(p1, p2, ncol=1, common.legend = FALSE, labels = c("A", "B")) 

ggsave("figures/combined_viability_plots.png", combined_plot, width=30, height=15, units=c("cm"))


# -----------------------------------------------------------------------------
# combine NFKB unfiltered and filtered plots in a grid
# -----------------------------------------------------------------------------

p1 <- unfiltered_seap_violin_plot
p2 <- filtered_seap_violin_plot

combined_plot <- ggarrange(p1, p2, ncol=1, common.legend = FALSE, labels = c("A", "B")) 

ggsave("figures/combined_seap_plots.png", combined_plot, width=30, height=15, units=c("cm"))
