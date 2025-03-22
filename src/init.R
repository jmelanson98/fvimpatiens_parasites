## models
library(brms)
library(bayestestR)

## general
library(tidyverse)
library(dplyr)

## plotting
library(ggplot2)
library(tidybayes)
library(ggthemes)
library(ggtext)
library(stringr)
library(viridis)
library(gridExtra)
library(cowplot)
library(modelr)
library(grid)


## variograms + other spatial stuff
library(raster)
library(ggspatial)
library(terra)
library(sf)
library(lme4)
library(gstat)
library(sp)
library(glme)
library(patchwork)

## trees
library(phytools)
library(ggtree)
library(ape)

save.dir0 <- "saved/"
if(!dir.exists(save.dir0)) {
  dir.create(save.dir0, showWarnings = FALSE)
}

save.dir <- "saved/tables"
if(!dir.exists(save.dir)) {
  dir.create(save.dir, showWarnings = FALSE)
}

fig.dir <- "figures"
if(!dir.exists(fig.dir)) {
  dir.create(fig.dir, showWarnings = FALSE)
}


fig.dir2 <- "figures/diagnostics"
if(!dir.exists(fig.dir2)) {
  dir.create(fig.dir2, showWarnings = FALSE)
}

