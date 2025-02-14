rm(list = ls())

# Setup -------------------------------------------------------------------

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggfortify)
library(broom)
library(fmsb)

# All additional info can be found in the README file

# Mouse mitocarta proteomics ----------------------------------------------
## Figure 2A-B: Heatmap & PCA
source("Code/mouse.R")


# Human protein atlas -----------------------------------------------------
## Figure 2D-F: Heatmap & PCA
source("Code/hpa.R")


# mtPPS conceptual figure ---------------------------------------------
## Supplemental Figure 5: Example mtPPS calculations using Complex I in cortex and liver
source("Code/hpa_mtPPS.R")


# mtPPS Fibroblasts - controls only ---------------------------------------
## Figures 6C-E: PCA, Heatmap with passage correlations, scatterplots
source("Code/fibroblasts_mtPPS_ctrls.R")


# mtPPS Tissues and Fibroblasts - everything combined ---------------------
## Figures 6G+H, Supplemental Figures S8-S9: Heatmap, radar charts, dynamic ranges
## !! Takes longer and requires more memory 
source("Code/hpa_fibroblasts_combined.R") 

# Session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
