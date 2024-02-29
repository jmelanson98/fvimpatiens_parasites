## models
library(brms)
library(bayesplot)
library(tidybayes)
library(ggplot2)


## plotting
library(tidyr)
library(dplyr)
library(viridis)
library(tidybayes)
library(gridExtra)
library(grid)
library(scales)
library(RColorBrewer)

library(rstantools)
library(performance)
library(bayestestR)
library(see)


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

