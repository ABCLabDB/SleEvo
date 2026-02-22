############################################################
# 00_setup.R
#
# Install and load R packages required for Sleep Evolution
# analyses. Run once from the repository root.
#
# Required working directory: Sleep_Evolution/
############################################################

required_packages <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "DescTools"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(DescTools)
})
