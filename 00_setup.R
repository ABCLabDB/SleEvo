## scripts/00_setup.R

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

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DescTools)
