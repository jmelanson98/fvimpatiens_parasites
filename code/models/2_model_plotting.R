## plotting based on tutorial:
## https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#different-kinds-of-averCanopyBin-predictions-with-multilevel-models

setwd('/Users/jenna1/Documents/UBC/bombus_project/fvimpatiens_parasites')
rm(list=ls())
library(marginaleffects)
source("code/src/ggplotThemes.R")
source("code/src/init.R")
source("code/src/misc.R")
source("code/src/posterior_manip.R")

## load model results and data
load(file="saved/Base500m_32610.Rdata")
load(file="saved/NativePar500m_32610.Rdata")
load(file="saved/ImpatiensPar500m_32610.Rdata")


# compute & save and/or load conditional effects
base.cond.effects <- conditional_effects(fit.bombus.base)
native.cond.effects <- conditional_effects(fit.par.native)
impatiens.cond.effects <- conditional_effects(fit.par.impatiens)

save(base.cond.effects, native.cond.effects, impatiens.cond.effects,
     file="saved/conditional_effects_32610.Rdata")
load(file="saved/conditional_effects_32610.Rdata")

# load landscape metrics
landscapemetrics = read.csv("data/landscapemetrics_32610.csv")

## *********************************************************************************
## Prepping axis labels -- veg & native bee abundance models
## *********************************************************************************
new.net = fvimp_brmsdf[fvimp_brmsdf$Subset == TRUE, ]
new.orig = orig.spec[fvimp_brmsdf$Subset == TRUE, ]

#create axis values for standardized variables
labs.doy = (pretty(new.orig$julian_date, n=6))
axis.doy =  standardize.axis(labs.doy,
                             new.orig$julian_date)

labs.fdiv = (pretty(new.orig$floral_diversity, n=6))
axis.fdiv =  standardize.axis(labs.fdiv,
                              new.orig$floral_diversity)

labs.fabun = (pretty(new.orig$floral_abundance, n=6))
axis.fabun =  standardize.axis(labs.fabun,
                               new.orig$floral_abundance)

labs.blueberry = (pretty(new.orig$prop_blueberry_500, n=6))
axis.blueberry =  standardize.axis(labs.blueberry,
                                   new.orig$prop_blueberry_500)

labs.edge = (pretty(new.orig$prop_edge_500, n=6))
axis.edge =  standardize.axis(labs.edge,
                              new.orig$prop_edge_500)


## ***********************************************************************
## Native bombus abundance ~ floral abundance (Fig 3a)
## ***********************************************************************
babun = base.cond.effects[["nativebeeabundance.nativebeeabundance_floral_abundance"]]

#ggplot
babun.fabun =
  #plot raw data
  ggplot(new.net, aes(x = floral_abundance, y = native_bee_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x="", y="Native *Bombus* abundance") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.fabun,
    labels =  labs.fabun) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  #plot model prediction with credible interval
  geom_line(data = babun, aes(x = floral_abundance, y=estimate__)) +
  geom_ribbon(data = babun, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red")

babun.fabun

# We computed average marginal effects on the outcome scale by 
# marginalizing over the observed distribution of covariates with
# avg_comparisons() from the marginaleffects package (Arel-Dalmasse, 2024)

# AME for 1 unit increase in log-scaled floral abundance
floral_sd = sd(new.orig$floral_abundance, na.rm = TRUE)
new_data = insight::get_data(fit.bombus.base)
subset_data = new_data[new_data$Subset == TRUE,]
floral_babun_ame = avg_comparisons(fit.bombus.base, 
                                   resp = "nativebeeabundance",
                                   variables = list(floral_abundance = c(0, 1/floral_sd)),
                                   newdata = subset_data)

## ***********************************************************************
## native bombus abundance ~ floral diversity (Fig 3b)
## ***********************************************************************
babun = base.cond.effects[["nativebeeabundance.nativebeeabundance_floral_diversity"]]

babun.fdiv =
  ggplot(new.net, aes(x = floral_diversity, y = native_bee_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x="", y="") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.fdiv,
    labels =  labs.fdiv) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) 
babun.fdiv


## ***********************************************************************
## native bombus abundance ~ proportion blueberry (Fig 3c)
## ***********************************************************************
babun = base.cond.effects[["nativebeeabundance.nativebeeabundance_prop_blueberry_500"]]

babun.bberry =
  ggplot(new.net, aes(x = prop_blueberry_500, y = native_bee_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x="", y="") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.blueberry,
    labels =  labs.blueberry) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  #plot model prediction with credible interval
  geom_line(data = babun, aes(x = prop_blueberry_500, y=estimate__)) +
  geom_ribbon(data = babun, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red")
babun.bberry

# AME for 10% increase in blueberry area
bb_sd = sd(new.orig$prop_blueberry_500, na.rm = TRUE)
bb_babun_ame = avg_comparisons(fit.bombus.base, 
                               resp = "nativebeeabundance",
                               variables = list(prop_blueberry_500 = c(0, 0.1/bb_sd)),
                               newdata = subset_data)
bb_babun_ame

## ***********************************************************************
## native bombus abundance ~ edge density (Fig 3d)
## ***********************************************************************
babun = base.cond.effects[["nativebeeabundance.nativebeeabundance_prop_edge_500"]]

#ggplot
babun.edge =
  ggplot(new.net, aes(x = prop_edge_500, y = native_bee_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.edge,
    labels =  labs.edge) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

babun.edge

## ***********************************************************************
## native bombus abundance ~ doy (Fig S3a)
## ***********************************************************************
babun = base.cond.effects[["nativebeeabundance.nativebeeabundance_julian_date"]]

#ggplot
babun.doy =
  #plot raw data
  ggplot(new.net, aes(x = julian_date, y = native_bee_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.05, width = 0.05) +
  labs(x = "Julian date", y = "Native *Bombus* abundance") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.doy,
    labels =  labs.doy) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = babun, aes(x = julian_date, y=estimate__)) +
  geom_ribbon(data = babun, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red") 
babun.doy


## ***********************************************************************
## prepping axis labels -- Bombus richness and impatiens abundance models
## ***********************************************************************
new.net = fvimp_brmsdf[fvimp_brmsdf$impSubset == TRUE, ]
new.orig = orig.spec[fvimp_brmsdf$impSubset == TRUE, ]

#create axis values for standardized variables
labs.doy = (pretty(new.orig$julian_date, n=6))
axis.doy =  standardize.axis(labs.doy,
                             new.orig$julian_date)

labs.fdiv = (pretty(new.orig$floral_diversity, n=6))
axis.fdiv =  standardize.axis(labs.fdiv,
                              new.orig$floral_diversity)

labs.fabun = (pretty(new.orig$floral_abundance, n=6))
axis.fabun =  standardize.axis(labs.fabun,
                               new.orig$floral_abundance)

labs.blueberry = (pretty(new.orig$prop_blueberry_500, n=6))
axis.blueberry =  standardize.axis(labs.blueberry,
                                   new.orig$prop_blueberry_500)

labs.edge = (pretty(new.orig$prop_edge_500, n=6))
axis.edge =  standardize.axis(labs.edge,
                              new.orig$prop_edge_500)

## ***********************************************************************
## bombus richness ~ floral abundance (Fig 3i)
## ***********************************************************************
brich = base.cond.effects[["bombusrichness.bombusrichness_floral_abundance"]]

brich.fabund =
  #plot raw data
  ggplot(new.net, aes(x = floral_abundance, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "", y = "*Bombus* species richness") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.fabun,
    labels =  labs.fabun) +
  scale_y_continuous() +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = brich, aes(x = floral_abundance, y=estimate__)) +
  geom_ribbon(data = brich, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red") #+

brich.fabund

# AME for 1 unit increase in log-scaled floral abundance
floral_sd = sd(new.orig$floral_abundance, na.rm = TRUE)
new_data = insight::get_data(fit.bombus.base)
subset_data = new_data[new_data$impSubset == TRUE,]
floral_brich_ame = avg_comparisons(fit.bombus.base, 
                                   resp = "bombusrichness",
                                   variables = list(floral_abundance = c(0, 1/floral_sd)),
                                   newdata = subset_data)

## ***********************************************************************
## bombus richness ~ floral diversity (Fig 3j)
## ***********************************************************************
brich = base.cond.effects[["bombusrichness.bombusrichness_floral_diversity"]]

#ggplot
brich.fdiv =
  #plot raw data
  ggplot(new.net, aes(x = floral_diversity, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.fdiv,
    labels =  labs.fdiv) +
  scale_y_continuous() +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = brich, aes(x = floral_diversity, y=estimate__)) +
  geom_ribbon(data = brich, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red") #+

brich.fdiv

# AME for 1 SD increase in floral diversity
new_data = insight::get_data(fit.bombus.base)
subset_data = new_data[new_data$impSubset == TRUE,]
fdiv_brich_ame = avg_comparisons(fit.bombus.base, 
                                 resp = "bombusrichness",
                                 variables = list(floral_diversity = c(0, 1)),
                                 newdata = subset_data)

## ***********************************************************************
## bombus richness ~ proportion blueberry (Figure 3k)
## ***********************************************************************
brich = base.cond.effects[["bombusrichness.bombusrichness_prop_blueberry_500"]]

#ggplot
brich.bberry =
  #plot raw data
  ggplot(new.net, aes(x = prop_blueberry_500, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.blueberry,
    labels =  labs.blueberry) +
  scale_y_continuous() +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = brich, aes(x = prop_blueberry_500, y=estimate__)) +
  geom_ribbon(data = brich, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red") #+

brich.bberry

## ***********************************************************************
## bombus richness ~ edge density (figure 3l)
## ***********************************************************************
brich = base.cond.effects[["bombusrichness.bombusrichness_prop_edge_500"]]

#ggplot
brich.edge =
  #plot raw data
  ggplot(new.net, aes(x = prop_edge_500, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.edge,
    labels =  labs.edge) +
  scale_y_continuous() +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = brich, aes(x = prop_edge_500, y=estimate__)) +
  geom_ribbon(data = brich, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "blue") #+

brich.edge

## ***********************************************************************
## bombus richness ~ julian date (Figure S3c)
## ***********************************************************************
brich = base.cond.effects[["bombusrichness.bombusrichness_julian_date"]]

#ggplot
brich.doy =
  #plot raw data
  ggplot(new.net, aes(x = julian_date, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.2, width = 0.2) +
  labs(x = "Julian date", y = "*Bombus* species richness") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.doy,
    labels =  labs.doy) +
  scale_y_continuous() +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = brich, aes(x = julian_date, y=estimate__)) +
  geom_ribbon(data = brich, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red")

brich.doy

## ***********************************************************************
## impatiens abundance ~ floral abundance (figure 3e)
## ***********************************************************************
iabund = base.cond.effects[["impatiensabundance.impatiensabundance_floral_abundance"]]

#ggplot
iabund.fabund =
  #plot raw data
  ggplot(new.net, aes(x = floral_abundance, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Floral abundance (log)", y = "*B. impatiens* abundance") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.fabun,
    labels =  labs.fabun) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = iabund, aes(x = floral_abundance, y=estimate__)) +
  geom_ribbon(data = iabund, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.5), fill = "red") 

iabund.fabund

## ***********************************************************************
## impatiens abundance ~ floral diversity (figure 3f)
## ***********************************************************************
iabund = base.cond.effects[["impatiensabundance.impatiensabundance_floral_diversity"]]

#ggplot
iabund.fdiv <- 
  
  #plot raw data
  ggplot(new.net, aes(x = floral_diversity, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Floral diversity", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.fdiv,
    labels =  labs.fdiv) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) 

iabund.fdiv

## ***********************************************************************
## impatiens abundance ~ proportion blueberry (fig 3g)
## ***********************************************************************
iabund = base.cond.effects[["impatiensabundance.impatiensabundance_prop_blueberry_500"]]

#ggplot
iabund.bberry <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_blueberry_500, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Proportion blueberry", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.blueberry,
    labels =  labs.blueberry) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

iabund.bberry

## ***********************************************************************
## impatiens abundance ~ edge density (fig 3h)
## ***********************************************************************
iabund =base.cond.effects[["impatiensabundance.impatiensabundance_prop_edge_500"]]

#ggplot
iabund.edge <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_edge_500, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Edge density", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.edge,
    labels =  labs.edge) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) 

iabund.edge

## ***********************************************************************
## impatiens abundance ~ julian date (figure S3b)
## ***********************************************************************
iabund = base.cond.effects[["impatiensabundance.impatiensabundance_julian_date"]]

#ggplot
iabund.doy <- 
  
  #plot raw data
  ggplot(new.net, aes(x = julian_date, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.05, width = 0.05) +
  labs(x = "Julian date", y = "*B*. *impatiens* abundance") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.doy,
    labels =  labs.doy) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = iabund, aes(x = julian_date, y=estimate__)) +
  geom_ribbon(data = iabund, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.8), fill = "red")
iabund.doy



## ***********************************************************************
## Parasitism Models
## ***********************************************************************
data.par.native = fvimp_brmsdf[fvimp_brmsdf$subsetNativePar == TRUE, ]
data.par.impatiens = fvimp_brmsdf[fvimp_brmsdf$subsetImpatiensPar == TRUE, ]


# Set color scheme for plots
faded_pale = "#D2E4D4"
faded_light = "#B6D3B8"
faded_medium = "#609F65"
faded_strong = "#4E8353"
faded_green = "#355938"
faded_dark = "#1B2D1C"


light_gold = "#F7EAC0"
lm_gold = "#F2DC97"
medium_gold = "#ECCF6F"
gold = "#E2B41D"
dark_gold = "#B99318"
darker_gold = "#907313"

## ***********************************************************************
## has crithidia ~ impatiens abundance
## ***********************************************************************
hascrithnative = native.cond.effects[["hascrithidia.hascrithidia_impatiens_abundance"]]
hascrithimpatiens = impatiens.cond.effects[["hascrithidia.hascrithidia_impatiens_abundance"]]

#ggplot
hascrith.imp =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hascrithimpatiens, aes(x = impatiens_abundance, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hascrithimpatiens, aes(x = impatiens_abundance, y=estimate__), 
            color = dark_gold, linetype = "dashed") +
  geom_ribbon(data = hascrithnative, aes(x = impatiens_abundance, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "plum") +
  geom_line(data = hascrithnative, aes(x = impatiens_abundance, y=estimate__), color = "darkorchid4") +
  
  # labels
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(-0.5,11)) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hascrith.imp


## ***********************************************************************
## has crithidia ~ doy
## ***********************************************************************
hascrithnative = native.cond.effects[["hascrithidia.hascrithidia_julian_date"]]
hascrithimpatiens = impatiens.cond.effects[["hascrithidia.hascrithidia_julian_date"]]

#ggplot
hascrith.doy =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hascrithimpatiens, aes(x = julian_date, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hascrithimpatiens, aes(x = julian_date, y=estimate__), 
            color = dark_gold, linetype = "dashed") +
  geom_ribbon(data = hascrithnative, aes(x = julian_date, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "plum") +
  geom_line(data = hascrithnative, aes(x = julian_date, y=estimate__), 
            color = "darkorchid4") +
  
  # labels
  labs(x = "Julian date", y = "*Crithidia spp.* prevalence") +
  scale_x_continuous(breaks = axis.doy,
                     labels =  labs.doy) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hascrith.doy

## ***********************************************************************
## has crithidia ~ bombus richness
## ***********************************************************************
hascrithnative = native.cond.effects[["hascrithidia.hascrithidia_bombus_richness"]]
hascrithimpatiens = impatiens.cond.effects[["hascrithidia.hascrithidia_bombus_richness"]]

#ggplot
hascrith.brich =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hascrithimpatiens, aes(x = bombus_richness, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hascrithimpatiens, aes(x = bombus_richness, y=estimate__), 
            color = dark_gold, linetype = "dashed") +
  geom_ribbon(data = hascrithnative, aes(x = bombus_richness, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hascrithnative, aes(x = bombus_richness, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  
  # labels
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hascrith.brich

## ***********************************************************************
## has crithidia ~ prop blueberry
## ***********************************************************************
hascrithnative = native.cond.effects[["hascrithidia.hascrithidia_prop_blueberry_500"]]
hascrithimpatiens = impatiens.cond.effects[["hascrithidia.hascrithidia_prop_blueberry_500"]]

#ggplot
hascrith.bberry =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hascrithimpatiens, aes(x = prop_blueberry_500, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hascrithimpatiens, aes(x = prop_blueberry_500, y=estimate__), 
            color = dark_gold, linetype = "dashed") +
  geom_ribbon(data = hascrithnative, aes(x = prop_blueberry_500, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "plum") +
  geom_line(data = hascrithnative, aes(x = prop_blueberry_500, y=estimate__), 
            color = "darkorchid4") +
  
  
  # labels
  labs(x = "", y = "*Crithidia spp.* prevalence") +
  scale_x_continuous(breaks = axis.blueberry,
                     labels =  labs.blueberry) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hascrith.bberry


## ***********************************************************************
## has crithidia ~ floral diversity
## ***********************************************************************
hascrithnative = native.cond.effects[["hascrithidia.hascrithidia_floral_diversity"]]
hascrithimpatiens = impatiens.cond.effects[["hascrithidia.hascrithidia_floral_diversity"]]

#ggplot
hascrith.fdiv =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hascrithimpatiens, aes(x = floral_diversity, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hascrithimpatiens, aes(x = floral_diversity, y=estimate__), 
            color = dark_gold, linetype = "dashed") +
  geom_ribbon(data = hascrithnative, aes(x = floral_diversity, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hascrithnative, aes(x = floral_diversity, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  # labels
  labs(x = "", y = "") +
  scale_x_continuous(breaks = axis.fdiv,
                     labels =  labs.fdiv) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hascrith.fdiv



## ***********************************************************************
## has apicystis ~ impatiens abundance
## ***********************************************************************
apicystisnative = native.cond.effects[["apicystis.apicystis_impatiens_abundance"]]
apicystisimpatiens = impatiens.cond.effects[["apicystis.apicystis_impatiens_abundance"]]

#ggplot
apicystis.imp =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = apicystisimpatiens, aes(x = impatiens_abundance, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = apicystisimpatiens, aes(x = impatiens_abundance, y=estimate__), 
            color = dark_gold, linetype = "dashed") +
  geom_ribbon(data = apicystisnative, aes(x = impatiens_abundance, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = apicystisnative, aes(x = impatiens_abundance, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  # labels
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(-0.5,11)) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

apicystis.imp


## ***********************************************************************
## has apicystis ~ doy
## ***********************************************************************
apicystisnative = native.cond.effects[["apicystis.apicystis_julian_date"]]
apicystisimpatiens = impatiens.cond.effects[["apicystis.apicystis_julian_date"]]

#ggplot
apicystis.doy =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = apicystisimpatiens, aes(x = julian_date, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = medium_gold) +
  geom_line(data = apicystisimpatiens, aes(x = julian_date, y=estimate__), 
            color = dark_gold) +
  geom_ribbon(data = apicystisnative, aes(x = julian_date, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "plum") +
  geom_line(data = apicystisnative, aes(x = julian_date, y=estimate__), 
            color = "darkorchid4") +
  
  # labels
  labs(x = "Julian date", y = "*Apicystis spp.* prevalence") +
  scale_x_continuous(breaks = axis.doy,
                     labels =  labs.doy) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

apicystis.doy

## ***********************************************************************
## has apicystis ~ bombus richness
## ***********************************************************************
apicystisnative = native.cond.effects[["apicystis.apicystis_bombus_richness"]]
apicystisimpatiens = impatiens.cond.effects[["apicystis.apicystis_bombus_richness"]]

#ggplot
apicystis.brich =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = apicystisimpatiens, aes(x = bombus_richness, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = medium_gold) +
  geom_line(data = apicystisimpatiens, aes(x = bombus_richness, y=estimate__), 
            color = dark_gold) +
  geom_ribbon(data = apicystisnative, aes(x = bombus_richness, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = apicystisnative, aes(x = bombus_richness, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  # labels
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

apicystis.brich

## ***********************************************************************
## has apicystis ~ prop blueberry
## ***********************************************************************
apicystisnative = native.cond.effects[["apicystis.apicystis_prop_blueberry_500"]]
apicystisimpatiens = impatiens.cond.effects[["apicystis.apicystis_prop_blueberry_500"]]

#ggplot
apicystis.bberry =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = apicystisimpatiens, aes(x = prop_blueberry_500, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = apicystisimpatiens, aes(x = prop_blueberry_500, y=estimate__), 
            color = dark_gold, linetype = "dashed") +
  geom_ribbon(data = apicystisnative, aes(x = prop_blueberry_500, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "plum") +
  geom_line(data = apicystisnative, aes(x = prop_blueberry_500, y=estimate__), 
            color = "darkorchid4") +
  
  # labels
  labs(x = "", y = "*Apicystis spp.* prevalence") +
  scale_x_continuous(breaks = axis.blueberry,
                     labels =  labs.blueberry) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

apicystis.bberry


## ***********************************************************************
## has apicystis ~ floral diversity
## ***********************************************************************
apicystisnative = native.cond.effects[["apicystis.apicystis_floral_diversity"]]
apicystisimpatiens = impatiens.cond.effects[["apicystis.apicystis_floral_diversity"]]

#ggplot
apicystis.fdiv =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = apicystisimpatiens, aes(x = floral_diversity, ymin = lower__, ymax = upper__),
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = apicystisimpatiens, aes(x = floral_diversity, y=estimate__), 
            color = dark_gold, linetype = "dashed") +
  geom_ribbon(data = apicystisnative, aes(x = floral_diversity, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "plum") +
  geom_line(data = apicystisnative, aes(x = floral_diversity, y=estimate__), 
            color = "darkorchid4") +
  
  # labels
  labs(x = "", y = "") +
  scale_x_continuous(breaks = axis.fdiv,
                     labels =  labs.fdiv) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

apicystis.fdiv


## ***********************************************************************
## has nosema ~ impatiens abundance
## ***********************************************************************
hasnosemanative = native.cond.effects[["hasnosema.hasnosema_impatiens_abundance"]]

#ggplot
hasnosema.imp =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hasnosemanative, aes(x = impatiens_abundance, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hasnosemanative, aes(x = impatiens_abundance, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  # labels
  labs(x = "*B. impatiens* abundance", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(-0.5,11)) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.imp


## ***********************************************************************
## has nosema ~ doy
## ***********************************************************************
hasnosemanative = native.cond.effects[["hasnosema.hasnosema_julian_date"]]

#ggplot
hasnosema.doy =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hasnosemanative, aes(x = julian_date, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hasnosemanative, aes(x = julian_date, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  # labels
  labs(x = "Julian date", y = "*Vairimorpha spp.* prevalence") +
  scale_x_continuous(breaks = axis.doy,
                     labels =  labs.doy) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.doy

## ***********************************************************************
## has nosema ~ bombus richness
## ***********************************************************************
hasnosemanative = native.cond.effects[["hasnosema.hasnosema_bombus_richness"]]

#ggplot
hasnosema.brich =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hasnosemanative, aes(x = bombus_richness, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hasnosemanative, aes(x = bombus_richness, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  # labels
  labs(x = "*Bombus* species richness", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.brich

## ***********************************************************************
## has nosema ~ prop blueberry
## ***********************************************************************
hasnosemanative = native.cond.effects[["hasnosema.hasnosema_prop_blueberry_500"]]

#ggplot
hasnosema.bberry =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hasnosemanative, aes(x = prop_blueberry_500, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hasnosemanative, aes(x = prop_blueberry_500, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  # labels
  labs(x = "Proportion blueberry", y = "*Vairimorpha spp.* prevalence") +
  scale_x_continuous(breaks = axis.blueberry,
                     labels =  labs.blueberry) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.bberry


## ***********************************************************************
## has nosema ~ floral diversity
## ***********************************************************************
hasnosemanative = native.cond.effects[["hasnosema.hasnosema_floral_diversity"]]

#ggplot
hasnosema.fdiv =
  ggplot() +
  #plot model prediction with credible interval
  geom_ribbon(data = hasnosemanative, aes(x = floral_diversity, ymin = lower__, ymax = upper__), 
              alpha = 0.4, fill = "lightgrey") +
  geom_line(data = hasnosemanative, aes(x = floral_diversity, y=estimate__), 
            color = "darkorchid4", linetype = "dashed") +
  
  # labels
  labs(x = "Floral diversity", y = "") +
  scale_x_continuous(breaks = axis.fdiv,
                     labels =  labs.fdiv) +
  theme_ms() +
  theme(legend.position = "none") +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.fdiv


## ***********************************************************************
## effect of caste on prevalence of all 3 parasites
## ***********************************************************************
hascrithidianative = native.cond.effects[["hascrithidia.hascrithidia_caste"]]
hasnosema = native.cond.effects[["hasnosema.hasnosema_caste"]]
apicystisnative = native.cond.effects[["apicystis.apicystis_caste"]]

hascrithidianative$parname = "crithidia"

apicystisnative$parname = "apicystis"

hasnosema$parname = "nosema"

colstokeep = c("caste", "estimate__", "lower__", "upper__", "parname")

allpar = rbind(hascrithidianative[colstokeep],
               hasnosema[colstokeep], 
               apicystisnative[colstokeep])

#plot
casteplot = ggplot(allpar, aes(x = parname, y = estimate__, color = caste)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, 
                position = position_dodge(width = 0.3)) +
  scale_colour_manual(values = c("darkorchid4", "plum")) +
  labs(x = "", y = "Parasite Prevalence",
       color = "Caste") +
  scale_x_discrete(labels = 
                     c("apicystis" = expression(italic("Apicystis spp.")), 
                       "crithidia" = expression(italic("Crithidia spp.")), 
                       "nosema" = expression(italic("Vairimorpha spp.")))) +
  scale_y_continuous(limits = c(0,1)) +
  theme_ms() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  theme(legend.position = "bottom") +
  geom_text(data = data.frame(
    x = c(2),
    y = c(0.75),
    label = c("**")
  ), aes(x=x, y=y, label=label),
  color = "black",
  size=10)

casteplot


## ***********************************************************************
## combine individuals plots for manuscript
## ***********************************************************************

#FIGURE 3: floral abundance & landscape effects on bee community
bombusgrid = grid.arrange(brich.fabund, brich.fdiv, brich.edge, brich.bberry,
                          babun.fabun, babun.fdiv, babun.edge, babun.bberry, 
                          iabund.fabund, iabund.fdiv, iabund.edge, iabund.bberry,
                          ncol = 4)
#add labels
bombusgrid <- ggdraw() +
  draw_plot(bombusgrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("(A)", "(D)", "(G)", "(J)", "(B)", "(E)", "(H)", "(K)", "(C)", "(F)", "(I)", "(L)"), 
                  x = c(0, 0.25, 0.5, 0.75, 0, 0.25, 0.5, 0.75, 0, 0.25, 0.5, 0.75), 
                  y = c(1, 1, 1, 1, 0.66, 0.66, 0.66, 0.66, 0.33, 0.33, 0.33, 0.33))

#export and save
ggsave(filename = "figures/manuscript_figures/bombusgrid.jpg", 
       plot = bombusgrid, 
       width = 4000, 
       height = 3500, 
       units = "px")

# #FIGURE FOR BLUEBERRY RESEARCH DAY
# partialbombusgrid = grid.arrange(babun.fabun, babun.bberry, babun.edge,
#                           brich.fabund, brich.bberry, brich.edge,
#                           ncol = 3)
# 
# #export and save
# ggsave(filename = "figures/bcberry.jpg", 
#        plot = partialbombusgrid, 
#        width = 3000, 
#        height = 3000, 
#        units = "px")
# 
# 
# #FIGURE FOR PBESA SLIDES
# partialbombusgrid2 = grid.arrange(babun.bberry, brich.bberry, iabund.bberry,
#                                   babun.edge, brich.edge, iabund.edge,
#                                   ncol = 3)
# 
# #export and save
# ggsave(filename = "figures/pbesa.jpg", 
#        plot = partialbombusgrid2, 
#        width = 3000, 
#        height = 2000, 
#        units = "px")



#FIGURE S3: bee community phenology
bombusdoy = grid.arrange(babun.doy, iabund.doy, brich.doy, ncol = 3)

#add labels
bombusdoy <- ggdraw() +
  draw_plot(bombusdoy, 0.015, 0, 1, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)"), 
                  x = c(0, 0.33, 0.66), 
                  y = c(1, 1, 1))
#save and export
ggsave(filename = "figures/manuscript_figures/bombusdoy.jpg", 
       plot = bombusdoy, 
       width = 4500, 
       height = 1500, 
       units = "px")


#FIGURE 4: parasite prevalence
parasitegrid = grid.arrange(hascrith.bberry, hascrith.imp, hascrith.brich, hascrith.fdiv,
                            apicystis.bberry, apicystis.imp, apicystis.brich, apicystis.fdiv,
                            hasnosema.bberry, hasnosema.imp, hasnosema.brich, hasnosema.fdiv,
                            ncol = 4)
#add labels
parasitegrid <- ggdraw() +
  draw_plot(parasitegrid, 0, 0, 1, 1) +
  draw_plot_label(c("(A)", "(D)", "(G)", "(J)", "(B)", "(E)", "(H)", "(K)", "(C)", "(F)", "(I)", "(L)"), 
                  x = c(0, 0.25, 0.5, 0.75, 0, 0.25, 0.5, 0.75, 0, 0.25, 0.5, 0.75), 
                  y = c(1, 1, 1, 1, 0.66, 0.66, 0.66, 0.66, 0.34, 0.34, 0.34, 0.34))

#export and save
ggsave(filename = "figures/manuscript_figures/parasitegrid.jpg", 
       plot = parasitegrid, 
       width = 4000, 
       height = 3500, 
       units = "px")


#FIGURE S6: parasite phenology
parasitedoy = grid.arrange(hascrith.doy, apicystis.doy, hasnosema.doy, ncol = 3)

#add labels
parasitedoy <- ggdraw() +
  draw_plot(parasitedoy, 0.015, 0, 1, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)"), 
                  x = c(0, 0.33, 0.66), 
                  y = c(1, 1, 1))
#y = c(0.97, 0.73, 0.48)) #if there's a legend
#save and export
ggsave(filename = "figures/manuscript_figures/parasitedoy.jpg", 
       plot = parasitedoy, 
       width = 4500, 
       height = 1500, 
       units = "px")


#Figure S7: parasitism ~ caste
ggsave(filename = "figures/manuscript_figures/casteplot.jpg", 
       plot = casteplot, 
       width = 1600, 
       height = 1200, 
       units = "px")



## ***********************************************************************
## Make some visualizations of landcover data
## ***********************************************************************

landscapemetrics[is.na(landscapemetrics)] = 0
lmet = pivot_longer(landscapemetrics, 
                    cols = -c("X", "sample_pt"),
                    names_to = "metric",
                    values_to = "value")
shdi = filter(lmet, str_detect(metric, "shdi"))
shdi_order = c("landscape_shdi_250", "landscape_shdi_500", "landscape_shdi_750")
edge = filter(lmet, str_detect(metric, "edge"))
edge_order = c("prop_edge_250", "prop_edge_500", "prop_edge_500_750")
blueberry = filter(lmet, str_detect(metric, "blueberry"))
blueberry_order = c("prop_blueberry_250", "prop_blueberry_500", "prop_blueberry_750")

plot_labels = c("250m", "500m", "750m")

panel1 = ggplot(shdi, aes(y = value, x = factor(metric, level = shdi_order))) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2),
              color="darkblue", alpha=0.5) +
  theme_minimal() +
  labs(y = "Shannon Diversity", x = "", ) +
  scale_x_discrete(labels = plot_labels)

panel2 = ggplot(edge, aes(y = value, x = factor(metric, level = edge_order))) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2),
              color="darkblue", alpha=0.5) +
  theme_minimal() +
  labs(y = "Edge density", x = "", ) +
  scale_x_discrete(labels = plot_labels)

panel3 = ggplot(blueberry, aes(y = value, x = factor(metric, level = blueberry_order))) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2),
              color="darkblue", alpha=0.5) +
  theme_minimal() +
  labs(y = "Proportion blueberry", x = "Buffer size", ) +
  scale_x_discrete(labels = plot_labels)

landscapeplots <- grid.arrange(
  arrangeGrob(panel1, panel2, panel3,
              ncol = 1))
ggsave("figures/manuscript_figures/landscapemetrics.jpg", 
       landscapeplots,
       height = 3000, width = 1500,
       units = "px")



## ***********************************************************************
## Plot distributions of landscape metrics
## ***********************************************************************
transect_metrics = landscapemetrics[,colnames(landscapemetrics) %in% c("sample_pt", "site", 
                                                         "prop_blueberry_500", 
                                                         "prop_edge_500", 
                                                         "landscape_shdi_500", 
                                                         "prop_blueberry_1000", 
                                                         "prop_edge_1000", 
                                                         "landscape_shdi_1000")]
transect_metrics$site = gsub("^([A-Za-z]{1,2})[0-9].*", "\\1", transect_metrics$sample_pt)
transect_metrics_long = transect_metrics %>% 
  pivot_longer(-c("sample_pt", "site"),
               names_to = "metric",
               values_to = "values")

split_parts = strsplit(transect_metrics_long$metric, "_")

transect_metrics_long$prefix = sapply(split_parts, function(x) paste(head(x, -1), collapse = "_"))
transect_metrics_long$scale = as.numeric(sapply(split_parts, tail, 1))

plot_list = list()
rows = c("landscape_shdi", "prop_blueberry", "prop_edge")
columns = c(500, 1000)
order = expand.grid(rows, columns)
count = 1

ylims <- list(
  prop_blueberry = c(0, 0.9),
  prop_edge      = c(0, 0.2),
  landscape_shdi = c(0, 2.5)
)

for (r in rows){
  for (c in columns){
    if(c == 500){
      if(r == "prop_edge"){
        ylab = "Edge density"
      } else if(r == "prop_blueberry"){
        ylab = "Proportion blueberry"
      } else if(r== "landscape_shdi"){
        ylab = "Shannon diversity"
      }
    } else{ylab = ""}
    
    if(r == "prop_edge"){
      xlab = "Site"
    } else{xlab = ""}
    
    plot_list[[count]] = ggplot(transect_metrics_long[transect_metrics_long$prefix == r &
                                                        transect_metrics_long$scale == c,], 
                                aes(x = site, y = values)) +
      geom_violin() +
      ylim(ylims[[r]]) +
      ylab(ylab) +
      xlab(xlab) +
      theme_minimal()
    count = count + 1
  }
}

header1 = ggplot() + theme_void() + ggtitle("500 meter buffer") + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
header2 = ggplot() + theme_void() + ggtitle("1000 meter buffer") + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
blank = nullGrob()

final = grid.arrange(blank, header1, blank, header2, 
                     blank, plot_list[[1]], blank, plot_list[[2]],
                     blank, blank, blank, blank,
                     blank, plot_list[[3]], blank, plot_list[[4]],
                     blank, blank, blank, blank,
                     blank, plot_list[[5]], blank, plot_list[[6]],
                     ncol = 4,
                     heights = c(0.15, 1, 0.1, 1, 0.1, 1), widths = c(0.1, 1, 0.1, 1))
site_grid = ggdraw() +
  draw_plot(final, 0.01, 0, 1, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)"), 
                  x = c(0, 0.52, 0, 0.52, 0, 0.52), 
                  y = c(0.98, 0.98, 0.66, 0.66, 0.33, 0.33))

ggsave("figures/manuscript_figures/landscape_metrics.jpg", site_grid,
       units = "px", height = 1500, width = 2500)
