# ==============================================================================
# 00 - Packages and setup
# ==============================================================================
# Loads the packages used throughout the analysis pipeline (scripts 01-07).
# Run this script once before running the numbered analysis scripts, or simply
# make sure the packages below are installed.
# ==============================================================================

required_packages <- c(
  "dplyr", "tidyverse", "readr", "tidyr",
  "ggplot2", "ggthemes", "patchwork", "ggeffects",
  "data.table", "formattable",
  "ggpubr", "rstatix",
  "glmmTMB", "betareg", "performance", "DHARMa",
  "broom.mixed", "purrr",
  "colorspace", "visreg", "scales", "ggsignif",
  "ape", "caper", "cowplot",
  "remotes"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))

# V.PhyloMaker2 is only available on GitHub (used in script 06)
if (!requireNamespace("V.PhyloMaker2", quietly = TRUE)) {
  remotes::install_github("jinyizju/V.PhyloMaker2")
}
library(V.PhyloMaker2)

# record package versions for reproducibility
# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
