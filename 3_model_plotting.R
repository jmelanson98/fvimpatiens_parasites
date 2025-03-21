## plotting based on tutorial:
## https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#different-kinds-of-averCanopyBin-predictions-with-multilevel-models

setwd('/Users/jenna1/Documents/UBC/Bombus Project/fvimpatiens_parasites')
rm(list=ls())
source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")

## load model results and data
load(file="/Users/jenna1/Documents/UBC/Bombus Project/Rdata_files/fvimpatiens_parasites/AllModels_fv.Rdata")
load(file="/Users/jenna1/Documents/UBC/Bombus Project/Rdata_files/fvimpatiens_parasites/AllModels_fv_inter.Rdata")

# save and/or load conditional effects
save(all.cond.effects, all.cond.effects.interaction,
file="/Users/jenna1/Documents/UBC/Bombus Project/Rdata_files/fvimpatiens_parasites/conditional_effects.Rdata")
load(file="/Users/jenna1/Documents/UBC/Bombus Project/Rdata_files/fvimpatiens_parasites/conditional_effects.Rdata")

## *********************************************************************************
## Prepping conditional effects and axis labels -- veg & native bee abundance models
## *********************************************************************************
new.net <- fvimp_brmsdf[fvimp_brmsdf$Subset == TRUE, ]
new.orig <- orig.spec[fvimp_brmsdf$Subset == TRUE, ]

# calculate all conditional effects
# only run this code if you don't have conditional effects Rdata saved -- it's slow
#all.cond.effects <- conditional_effects(fit.bombus.all)

#create axis values for standardized variables
labs.doy <- (pretty(new.orig$julian_date, n=8))
axis.doy <-  standardize.axis(labs.doy,
                              new.orig$julian_date)

labs.fdiv <- (pretty(new.orig$floral_diversity, n=8))
axis.fdiv <-  standardize.axis(labs.fdiv,
                              new.orig$floral_diversity)

labs.fabun <- (pretty(new.orig$floral_abundance, n=8))
axis.fabun <-  standardize.axis(labs.fabun,
                               new.orig$floral_abundance)

labs.blueberry <- (pretty(new.orig$prop_blueberry, n=9))
axis.blueberry <-  standardize.axis(labs.blueberry,
                                    new.orig$prop_blueberry)

labs.edge <- (pretty(new.orig$prop_edge, n=8))
axis.edge <-  standardize.axis(labs.edge,
                               new.orig$prop_edge)


## ***********************************************************************
## native bombus abundance ~ floral abundance (Fig 3a)
## ***********************************************************************
babun <-
  all.cond.effects[["nativebeeabundance.nativebeeabundance_floral_abundance"]]

#ggplot
babun.fabun <- 
  
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

## ***********************************************************************
## native bombus abundance ~ floral diversity (Fig 3b)
## ***********************************************************************
babun <-
  all.cond.effects[["nativebeeabundance.nativebeeabundance_floral_diversity"]]

#ggplot
babun.fdiv <- 
  
  #plot raw data
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
babun <-
  all.cond.effects[["nativebeeabundance.nativebeeabundance_prop_blueberry"]]

#ggplot
babun.bberry <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_blueberry, y = native_bee_abundance)) +
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
  geom_line(data = babun, aes(x = prop_blueberry, y=estimate__)) +
  geom_ribbon(data = babun, aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "red")

  geom_text(data = data.frame(
    x = c(2.5),
    y = c(40),
    label = c("***p < 0 = 0.98")
  ), aes(x=x, y=y, label=label),
  color = "black",
  size=5)
babun.bberry

## ***********************************************************************
## native bombus abundance ~ edge density (Fig 3d)
## ***********************************************************************
babun <-
  all.cond.effects[["nativebeeabundance.nativebeeabundance_prop_edge"]]

#ggplot
babun.edge <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_edge, y = native_bee_abundance)) +
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
        text = element_text(size=16)) +
  
  #plot model prediction with credible interval
  geom_line(data = babun, aes(x = prop_edge, y=estimate__)) +
  geom_ribbon(data = babun, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "blue") #+

babun.edge

## ***********************************************************************
## native bombus abundance ~ doy (Fig S3a)
## ***********************************************************************
babun <-
  all.cond.effects[["nativebeeabundance.nativebeeabundance_julian_date"]]

#ggplot
babun.doy <- 
  
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
## prepping axis labels -- bombus richness and impatiens abund models
## ***********************************************************************
new.net <- fvimp_brmsdf[fvimp_brmsdf$impSubset == TRUE, ]
new.orig <- orig.spec[fvimp_brmsdf$impSubset == TRUE, ]

#create axis values for standardized variables
labs.doy <- (pretty(new.orig$julian_date, n=8))
axis.doy <-  standardize.axis(labs.doy,
                              new.orig$julian_date)

labs.fdiv <- (pretty(new.orig$floral_diversity, n=8))
axis.fdiv <-  standardize.axis(labs.fdiv,
                               new.orig$floral_diversity)

labs.fabun <- (pretty(new.orig$floral_abundance, n=8))
axis.fabun <-  standardize.axis(labs.fabun,
                                new.orig$floral_abundance)

labs.blueberry <- (pretty(new.orig$prop_blueberry, n=9))
axis.blueberry <-  standardize.axis(labs.blueberry,
                                    new.orig$prop_blueberry)

labs.edge <- (pretty(new.orig$prop_edge, n=8))
axis.edge <-  standardize.axis(labs.edge,
                               new.orig$prop_edge)

## ***********************************************************************
## bombus richness ~ floral abundance (Fig 3i)
## ***********************************************************************
brich <-
  all.cond.effects[["bombusrichness.bombusrichness_floral_abundance"]]

#ggplot
brich.fabund <- 
  
  #plot raw data
  ggplot(new.net, aes(x = floral_abundance, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "Floral abundance (log)", y = "*Bombus* species richness") +
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

## ***********************************************************************
## bombus richness ~ floral diversity (Fig 3j)
## ***********************************************************************
brich <-
  all.cond.effects[["bombusrichness.bombusrichness_floral_diversity"]]

#ggplot
brich.fdiv <- 
  
  #plot raw data
  ggplot(new.net, aes(x = floral_diversity, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "Floral diversity", y = "") +
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

## ***********************************************************************
## bombus richness ~ proportion blueberry (Figure 3k)
## ***********************************************************************
brich <-
  all.cond.effects[["bombusrichness.bombusrichness_prop_blueberry"]]

#ggplot
brich.bberry <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_blueberry, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "Proportion blueberry", y = "") +
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
  geom_line(data = brich, aes(x = prop_blueberry, y=estimate__)) +
  geom_ribbon(data = brich, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red") #+

brich.bberry

## ***********************************************************************
## bombus richness ~ edge density (figure 3l)
## ***********************************************************************
brich <-
  all.cond.effects[["bombusrichness.bombusrichness_prop_edge"]]

#ggplot
brich.edge <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_edge, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "Edge density", y = "") +
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
  geom_line(data = brich, aes(x = prop_edge, y=estimate__)) +
  geom_ribbon(data = brich, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "red") #+

brich.edge

## ***********************************************************************
## bombus richness ~ julian date (Figure S3c)
## ***********************************************************************
brich <-
  all.cond.effects[["bombusrichness.bombusrichness_julian_date"]]

#ggplot
brich.doy <- 
  
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
iabund <-
  all.cond.effects[["impatiensabundance.impatiensabundance_floral_abundance"]]

#ggplot
iabund.fabund <- 
  
  #plot raw data
  ggplot(new.net, aes(x = floral_abundance, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "", y = "*B. impatiens* abundance") +
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
iabund <-
  all.cond.effects[["impatiensabundance.impatiensabundance_floral_diversity"]]

#ggplot
iabund.fdiv <- 
  
  #plot raw data
  ggplot(new.net, aes(x = floral_diversity, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "", y = "") +
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
iabund <-
  all.cond.effects[["impatiensabundance.impatiensabundance_prop_blueberry"]]

#ggplot
iabund.bberry <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_blueberry, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.blueberry,
    labels =  labs.blueberry) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) #+
  geom_text(data = data.frame(
    x = c(2.5),
    y = c(40),
    label = c("n.s.")
  ), aes(x=x, y=y, label=label),
  color = "black",
  size=5)

iabund.bberry

## ***********************************************************************
## impatiens abundance ~ edge density (fig 3h)
## ***********************************************************************
iabund <-
  all.cond.effects[["impatiensabundance.impatiensabundance_prop_edge"]]

#ggplot
iabund.edge <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_edge, y = impatiens_abundance)) +
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
        text = element_text(size=16)) #+
geom_text(data = data.frame(
  x = c(2.5),
  y = c(40),
  label = c("n.s.")
), aes(x=x, y=y, label=label),
color = "black",
size=5)

iabund.edge

## ***********************************************************************
## impatiens abundance ~ julian date (figure S3b)
## ***********************************************************************
iabund <-
  all.cond.effects[["impatiensabundance.impatiensabundance_julian_date"]]

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
data.par <- fvimp_brmsdf[fvimp_brmsdf$subsetPar == TRUE, ]

## ***********************************************************************
## has crithidia ~ impatiens abundance
## ***********************************************************************
hascrith <-
  all.cond.effects[["hascrithidia.hascrithidia_impatiens_abundance"]]

#ggplot
hascrith.imp <- 
  
  #plot model prediction with credible interval
  ggplot(hascrith, aes(x = impatiens_abundance, y = estimate__)) + 
  geom_line(aes(x = impatiens_abundance, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                                   alpha=0.5), fill = "blue") +

  #plot raw data
  geom_jitter(data = data.par, aes(x=impatiens_abundance,y = hascrithidia), height = 0.02, width = 0.3, alpha = 0.4) +
  labs(x = "*B. impatiens* abundance", y = "*Crithidia spp.* prevalence") +
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
hascrith <-
  all.cond.effects[["hascrithidia.hascrithidia_julian_date"]]

#ggplot
hascrith.doy <- 
  
  #plot model prediction with credible interval
  ggplot(hascrith, aes(x = julian_date, y = estimate__)) + 
  geom_line(aes(x = julian_date, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "red") +
  
  #plot raw data
  geom_jitter(data = data.par, aes(x=julian_date,y = hascrithidia), alpha = 0.2, height = 0.02) +
  labs(x = "Julian date", y = "*Crithidia spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = axis.doy,
                     labels =  labs.doy) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hascrith.doy


## ***********************************************************************
## has crithidia ~ status * impatiens abundance
## ***********************************************************************
# calculate all conditional effects
# only run this line of code if you don't have condition effects saved, it's slow:
#all.cond.effects.interaction <- conditional_effects(fit.bombus.inter)

hascrith <-
  all.cond.effects.interaction[["hascrithidia.hascrithidia_impatiens_abundance:status"]]

#ggplot
hascrith.imp.status <- 
  
  #plot model prediction with credible interval
  ggplot(hascrith, aes(x = impatiens_abundance, y = estimate__, color = status)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = status,
                  alpha=1)) +
  geom_line(aes(x = impatiens_abundance, y=estimate__), linewidth = 1.5) +
  scale_fill_manual(values = c("black", "lightblue")) +
  scale_colour_manual(values = c("black", "lightblue")) +
  #scale_color_manual(values = c("black", "blue")) +
  
  #plot raw data
  #geom_jitter(data = data.par, aes(x=impatiens_abundance,y = hascrithidia), height = 0.05) +
  labs(x = "*B. impatiens* abundance", y = "*Crithidia spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,11)) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hascrith.imp.status

## ***********************************************************************
## has crithidia ~ status * native bombus abundance
## ***********************************************************************
hascrith <-
  all.cond.effects.interaction[["hascrithidia.hascrithidia_native_bee_abundance:status"]]
hascrith = hascrith %>% filter((status != "nonnative") | (native_bee_abundance <= max(data.par$native_bee_abundance[data.par$final_id == "Bombus_impatiens"])))
  
#ggplot
hascrith.babund.status <- 
  
  #plot model prediction with credible interval
  ggplot(hascrith, aes(x = native_bee_abundance, y = estimate__, color = status)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = status, alpha = 1)) +
  geom_line(aes(x = native_bee_abundance, y=estimate__), linewidth = 1.5) +
  scale_fill_manual(values = c("black", "lightblue")) +
  scale_color_manual(values = c("black", "lightblue")) +
  
  #plot raw data
  #geom_jitter(data = site_rates, aes(x=native_bee_abundance,y = site_hascrithidia, cex = numbees)) +
  labs(x = "Native *Bombus* abundance", y = "*Crithidia spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "right") +
  scale_alpha(guide = "none") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) + 
  geom_text(data = data.frame(
    x = c(35),
    y = c(1.1),
    label = c("**")
  ), aes(x=x, y=y, label=label),
  color = "black",
  fontface = "bold",
  size=12)

hascrith.babund.status

#save legend for plot grid
g <- ggplotGrob(hascrith.babund.status)
legend_index <- which(g$layout$name == "guide-box-right")
interaction_legend <- g$grobs[[legend_index]]

# remove legend from plot
hascrith.babund.status <- hascrith.babund.status + theme(legend.position = "none")

## ***********************************************************************
## has nosema ~ impatiens abundance
## ***********************************************************************
hasnosema <-
  all.cond.effects[["hasnosema.hasnosema_impatiens_abundance"]]

#ggplot
hasnosema.imp <- 
  
  #plot model prediction with credible interval
  ggplot(hasnosema, aes(x = impatiens_abundance, y = estimate__)) + 
  geom_line(aes(x = impatiens_abundance, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "blue") +
  
  #plot raw data
  geom_jitter(data = data.par, aes(x=impatiens_abundance,y = hasnosema), height = 0.02, width = 0.3, alpha = 0.4) +
  labs(x = "*B. impatiens* abundance", y = "*Vairimorpha spp.* prevalence") +
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
hasnosema <-
  all.cond.effects[["hasnosema.hasnosema_julian_date"]]

#ggplot
hasnosema.doy <- 
  
  #plot raw data
  ggplot(data = data.par, aes(x=julian_date, y = hasnosema)) +
  geom_jitter(alpha = 0.2, height = 0.02) +
  labs(x = "Julian date", y = "*Vairimorpha spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = axis.doy,
                     labels =  labs.doy) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.doy


## ***********************************************************************
## apicystis ~ doy
## ***********************************************************************
apicystis <-
  all.cond.effects[["apicystis.apicystis_julian_date"]]

#ggplot
apicystis.doy <- 
  
  #plot model prediction with credible interval
  ggplot(apicystis, aes(x = julian_date, y = estimate__)) + 
  geom_line(aes(x = julian_date, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "red") +
  
  #plot raw data
  geom_jitter(data = data.par, aes(x=julian_date,y = apicystis), height = 0.02, alpha = 0.2) +
  labs(x = "Julian date", y = "*Apicystis spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = axis.doy,
                     labels =  labs.doy) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

apicystis.doy



## ***********************************************************************
## apicystis ~ floral diversity
## ***********************************************************************
apicystis <-
  all.cond.effects[["apicystis.apicystis_floral_diversity"]]

#ggplot
apicystis.fdiv <- 
  
  #plot model prediction with credible interval
  ggplot(apicystis, aes(x = floral_diversity, y = estimate__)) + 
  geom_line(aes(x = floral_diversity, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "blue") +
  
  #plot raw data
  geom_jitter(data = data.par, aes(x=floral_diversity,y = apicystis), alpha = 0.2, height = 0.02) +
  labs(x = "Floral diversity", y = "*Apicystis spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = axis.fdiv,
                     labels =  labs.fdiv) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

apicystis.fdiv



## ***********************************************************************
## apicystis ~ impatiens abundance
## ***********************************************************************
apicystis <-
  all.cond.effects[["apicystis.apicystis_impatiens_abundance"]]

apicystis.imp <- 
  
  #plot raw data
  ggplot(data = data.par, aes(x=impatiens_abundance,y = apicystis)) +
  geom_jitter(height = 0.02, width = 0.3, alpha = 0.4) +
  labs(x = "*B. impatiens* abundance", y = "*Apicystis spp.* prevalence") +
  theme_ms() +
  #theme(legend.title = element_text(hjust = 0.5), legend.position = "bottom") +
  #guides(size = "none", shape = "none", fill = "none") +
  #guides(
   # size = guide_legend(title.position = "top",
   #                     ncol = 3),
   # color = guide_legend(title.position = "top",
  #    ncol = 2)) +
  scale_x_continuous(limits = c(-0.5,11)) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) 
apicystis.imp

## ***********************************************************************
## effect of caste on prevalence of all 3 parasites
## ***********************************************************************
hascrithidia <-
  all.cond.effects[["hascrithidia.hascrithidia_caste"]]
hasnosema <-
  all.cond.effects[["hasnosema.hasnosema_caste"]]
apicystis <-
  all.cond.effects[["apicystis.apicystis_caste"]]
hascrithidia$parname = "crithidia"
hasnosema$parname = "nosema"
apicystis$parname = "apicystis"
allpar = rbind(hascrithidia, hasnosema, apicystis)

#plot
casteplot = ggplot(allpar, aes(x = parname, y = estimate__, color = caste)) +
  geom_point(position = position_dodge(width = 0.3)) +  # Separate points by model
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2, 
                position = position_dodge(width = 0.3)) +  # Error bars for each model
  labs(x = "", y = "Parasite Prevalence",
       color = "Caste") +
  scale_x_discrete(labels = 
                     c("apicystis" = expression(italic("Apicystis spp.")), 
                       "crithidia" = expression(italic("Crithidia spp.")), 
                       "nosema" = expression(italic("Vairimorpha spp.")))) +
  theme_ms() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  theme(legend.position = "bottom") +
  geom_text(data = data.frame(
    x = c(2),
    y = c(1),
    label = c("**")
  ), aes(x=x, y=y, label=label),
  color = "black",
  size=12)

casteplot


## ***********************************************************************
## combine individuals plots for manuscript
## ***********************************************************************

#FIGURE 3: floral abundance & landscape effects on bee community
bombusgrid = grid.arrange(babun.fabun, babun.fdiv, babun.bberry, babun.edge,
                          iabund.fabund, iabund.fdiv, iabund.bberry, iabund.edge,
                          brich.fabund, brich.fdiv, brich.bberry, brich.edge,
                          ncol = 4)
#add labels
bombusgrid <- ggdraw() +
  draw_plot(bombusgrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), 
                  x = c(0.01, 0.26, 0.51, 0.76, 0.01, 0.26, 0.51, 0.76, 0.01, 0.26, 0.51, 0.76), 
                  y = c(0.97, 0.97, 0.97, 0.97, 0.64, 0.64, 0.64, 0.64, 0.3, 0.3, 0.3, 0.3))

#export and save
ggsave(filename = "figures/manuscript_figures/bombusgrid.jpg", 
       plot = bombusgrid, 
       width = 4000, 
       height = 3500, 
       units = "px")

#FIGURE S3: bee community phenology
bombusdoy = grid.arrange(babun.doy, iabund.doy, brich.doy, ncol = 3)

#add labels
bombusdoy <- ggdraw() +
  draw_plot(doybeegrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c"), 
                  x = c(0.02, 0.35, 0.68), 
                  y = c(0.97, 0.97, 0.97))
#save and export
ggsave(filename = "figures/manuscript_figures/bombusdoy.jpg", 
       plot = bombusdoy, 
       width = 4500, 
       height = 1500, 
       units = "px")


#FIGURE 4: impatiens effects on parasite prevalence
impatiensgrid = grid.arrange(hascrith.imp, hasnosema.imp, apicystis.imp,
                                     ncol = 3)
#add labels
impatiensgrid <- ggdraw() +
  draw_plot(impatiensgrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c"), 
                  x = c(0.02, 0.35, 0.68), 
                  y = c(0.97, 0.97, 0.97))
#save and export
ggsave(filename = "figures/manuscript_figures/impatiensgrid.jpg", 
       plot = impatiensgrid, 
       width = 4500, 
       height = 1500, 
       units = "px")


# #impatiens vs parasites (vertical)
# impatiensparasitegrid = grid.arrange(hascrith.imp, hasnosema.imp, apicystis.imp, #shared_legend,
#                                      ncol = 1)
# #add subplot labels
# final_plot <- ggdraw() +
#   draw_plot(impatiensparasitegrid, 0.015, 0, 1, 1) +
#   draw_plot_label(c("a", "b", "c"), 
#                   x = c(0, 0, 0), 
#                   y = c(0.95, 0.62, 0.3))
# #y = c(0.97, 0.73, 0.48))
# print(final_plot)


#FIGURE S6: parasite phenology
parasitedoy = grid.arrange(hascrith.doy, hasnosema.doy, apicystis.doy, ncol = 3)

#add labels
parasitedoy <- ggdraw() +
  draw_plot(parasitedoy, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c"), 
                  x = c(0, 0.34, 0.68), 
                  y = c(0.97, 0.97, 0.97))
                  #y = c(0.97, 0.73, 0.48)) #if there's a legend
#save and export
ggsave(filename = "figures/manuscript_figures/parasitedoy.jpg", 
       plot = parasitedoy, 
       width = 4500, 
       height = 1500, 
       units = "px")


#FIGURE 5: interaction plots for crithidia
crithinteraction <- grid.arrange(
  arrangeGrob(hascrith.babund.status, hascrith.imp.status,
              ncol = 1),
  interaction_legend, # Add the shared legend
  ncol = 1, # 1 column: grid on top, legend below
  heights = c(10, 1) # Adjust the relative heights
)

#add labels
crithinteraction <- ggdraw() +
  draw_plot(crithinteraction, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b"), 
                  x = c(0.01, 0.01), 
                  y = c(0.97, 0.5))
#add labels
ggsave(filename = "figures/manuscript_figures/crithinteraction.jpg", 
       plot = crithinteraction, 
       width = 1500, 
       height = 3000, 
       units = "px")


