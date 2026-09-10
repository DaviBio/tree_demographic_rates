# ==============================================================================
# 07 - Supplementary Figure S1: phylogenetic tree with caption
# ==============================================================================
# Combines the phylogenetic tree panel produced in 06_phylogenetic_signal.R
# (object `phylo_panel`) with a formatted caption and saves the final PDF.
# Run this script right after 06_phylogenetic_signal.R, in the same session.
#
# Output: supplementary/figure_S1_phylogenetic_tree_with_caption.pdf
# ==============================================================================

if (!require(cowplot)) install.packages("cowplot")

library(tidyverse)
library(patchwork)
library(cowplot)

dir.create("supplementary", showWarnings = FALSE)

# ==============================================================================
# 1. FORMATTED CAPTION TEXT
# ==============================================================================
raw_caption <- "Supporting Information 9: Phylogenetic tree for the 31 tree species sampled in the study, constructed using the V.PhyloMaker2 package (scenario S3) based on the vascular plant megaphylogeny. The species *Siphoneugena reitzii* was manually included in the phylogeny as a sister lineage to *Myrciaria floribunda*, in accordance with recent phylogenies proposed for Myrtaceae in the Neotropical region (The Neotropical Myrtaceae Working Group, 2024). Phylogenetic signal tests indicated that the categorical predictor variables do not exhibit strong phylogenetic structuring or structuring different from random expectation (Fritz-Purvis D-statistic for Biogeographic Origin: D = 0.75, p_random = 0.24; for Leaf Habit/Deciduousness: D = 0.59, p_random = 0.15). Additionally, the Phylogenetic Moran's I test applied to the mean residuals of the generalized linear model (ordbeta) confirmed the complete absence of phylogenetic autocorrelation in the model error (I_observed = -0.05; p = 0.57), validating the statistical independence of the species-level analysis."

# wrap the text at ~130 characters per line
wrapped_caption <- str_wrap(raw_caption, width = 130)

# turn the text into a ggplot drawing object
caption_plot <- ggdraw() +
  draw_label(
    wrapped_caption,
    fontface = "plain",
    size = 11,
    x = 0.05,          # left margin
    y = 0.95,          # aligned to the top of the caption block
    hjust = 0,
    vjust = 1,
    lineheight = 1.3,
    color = "#333333"
  )

# ==============================================================================
# 2. COMBINE THE FIGURE WITH THE CAPTION BELOW IT
# ==============================================================================

full_figure <- plot_grid(
  phylo_panel,
  caption_plot,
  ncol = 1,
  rel_heights = c(1, 0.3)
)

# save the final PDF
ggsave(
  filename = "supplementary/figure_S1_phylogenetic_tree_with_caption.pdf",
  plot = full_figure,
  width = 12,
  height = 11.5,
  units = "in",
  dpi = 300
)
