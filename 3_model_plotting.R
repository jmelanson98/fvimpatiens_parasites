
## plotting based on tutorial:
## https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#different-kinds-of-averCanopyBin-predictions-with-multilevel-models

setwd('/Users/jenna1/Documents/UBC/Bombus Project/fvimpatiens_parasites')
rm(list=ls())
source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")
source("src/makeMapGrids.R")

## load model results and data
load(file="saved/AllModels_fv_beerichbeta.Rdata")
load(file="saved/AllModels_fv_beerichbeta_inter.Rdata")

# save and/or load conditional effects
save(all.cond.effects, all.cond.effects.interaction,
file="saved/conditional_effects.Rdata")
load(file="saved/conditional_effects.Rdata")

## ***********************************************************************
## exploratory plots / bar graphs
## ***********************************************************************

##parasite counts
orig.spec$bberry_bin <- cut(orig.spec$prop_blueberry,
                        breaks = c(0, 0.0001, .16, 1), 
                        include.lowest = T, right = T)
parasite.count.table <- orig.spec %>%
  filter(apidae == 1) %>%
  select(barcode_id, apicystis, ascosphaera, cbombii, cexpoeki, 
         crithidiaspp, nbombii, nceranae, bberry_bin, apidae) %>%
  pivot_longer(cols=c(apicystis, ascosphaera, cbombii,
                      crithidiaspp, cexpoeki, nceranae, nbombii, apidae),
               names_to = 'ParasiteName', values_to = 'HasParasite') %>%
  filter(HasParasite == 1)

#using #4C4C5A
parasite.hist <- parasite.count.table %>%
  ggplot(aes(y= fct_infreq(ParasiteName))) + 
  geom_bar(stat = 'count',
           aes(fill = factor(bberry_bin)), position = "dodge") +
  scale_fill_manual(values=c('#004D40', '#FFC107', 'lightblue'),
                    name = "Proportion blueberry", 
                    labels=c('[0]', '(0, 0.16]', '(0.16, 1]')) +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, color = "black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(color = "black"),
        text = element_text(size=16)) +
  labs(y='Parasite', x='Number of \n Infected Individuals') +
  xlim(0,300) +
  scale_y_discrete(labels=c("nceranae"=expression(italic('Nosema ceranae')),
                            "nbombii"=expression(italic('Nosema bombi')),
                            'ascosphaera'=expression(paste(italic('Ascosphaera'), ' spp.')),
                            'apicystis'=expression(paste(italic('Apicystis'), ' spp.')),
                            "cbombii"=expression(italic('Crithidia bombi')),
                            "cexpoeki"=expression(italic('Crithidia expoeki')),
                            "crithidiaspp"=expression(paste(italic('Crithidia'), ' spp.')),
                            parse=TRUE))

parasite.hist


## ***********************************************************************
## prepping for newdata draws -- vegetation & native bee abundance models
## ***********************************************************************
new.net <- fvimp_brmsdf[fvimp_brmsdf$Subset == TRUE, ]
new.orig <- orig.spec[fvimp_brmsdf$Subset == TRUE, ]

# calculate all conditional effects
all.cond.effects <- conditional_effects(fit.bombus.all.beta.base)

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
## veg diversity and doy
## ***********************************************************************
vdiv <-
  all.cond.effects[["floraldiversity.floraldiversity_julian_date"]]

#ggplot
vdiv.doy <- 
  
  #plot model prediction with credible interval
  ggplot(vdiv, aes(x = julian_date,
                               y = estimate__)) +
  geom_line(aes(x = julian_date, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.3, fill = "Red")) +
  
  #plot raw data
  geom_point(data=new.net,
             aes(x=julian_date, y=floral_diversity), cex=2, alpha = 0.2) +
  labs(x = "Day of Year", y = "Floral diversity",
       fill = "Credible Interval") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.doy,
    labels =  labs.doy) +
  scale_y_continuous(
    labels = labs.fdiv,
    breaks = axis.fdiv) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))


  geom_text(data = data.frame(
    x = c(-1.5),
    y = c(2.5),
    label = c("linear term: ***p > 0 = 1
      quadratic term: ***p < 0 = 1")
  ), aes(x=x, y=y, label=label),
            color = "black",
            size=5)

vdiv.doy


## ***********************************************************************
## veg abundance and doy
## ***********************************************************************
vabun <-
  all.cond.effects[["floralabundance.floralabundance_julian_date"]]

#ggplot
vabun.doy <- 
  
  #plot model prediction with credible interval
  ggplot(vabun, aes(x = julian_date,
                   y = estimate__)) +
  geom_line(aes(x = julian_date, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.3, fill = "Reds")) +
  
  #plot raw data
  geom_point(data=new.net,
             aes(x=julian_date, y=floral_abundance), cex=2, alpha = 0.2) +
  labs(x = "Day of Year", y = "Floral abundance",
       fill = "Credible Interval") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.doy,
    labels =  labs.doy) +
  scale_y_continuous(
    labels = labs.fabun,
    breaks = axis.fabun) +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) 
vabun.doy

## ***********************************************************************
## floral abundance and (native) bombus abundance
## ***********************************************************************

#using epred draws
babun.fabun.ci = new.net %>%
  data_grid(floral_abundance = seq_range(floral_abundance, n=101),
            floral_diversity = mean(floral_diversity),
            prop_blueberry = mean(prop_blueberry),
            landscape_shdi = mean(landscape_shdi),
            prop_edge = mean(prop_edge),
            impatiens_abundance = mean(impatiens_abundance),
            bombus_shannon_diversity = mean(bombus_shannon_diversity),
            julian_date = mean(julian_date),
            impSubset = TRUE,
            subsetPar = TRUE,
            native_bee_abundance = mean(native_bee_abundance),
            status = "native",
            caste = "worker",
            Subset = TRUE) %>%
  add_epred_draws(fit.bombus.all.base, resp = "nativebeeabundance", allow_new_levels =TRUE) %>%
  ggplot(aes(x = floral_abundance, y = native_bee_abundance)) +
  stat_lineribbon(aes(y = .epred), .width = c(.99, .95, .8, .5)) +
  scale_fill_brewer(palette = "Reds") +
  #plot raw data
  geom_jitter(data=new.net,
              aes(x=floral_abundance, y=native_bee_abundance), cex=2, alpha = 0.2) +
  labs(y = "Native *Bombus* abundance") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.fabun,
    labels =  labs.fabun) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(axis.title.x = element_blank(), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

babun.fabun.ci


#using conditional effects
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
## floral diversity on (native) bombus abundance 
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
## proportion blueberry on (native) bombus abundance 
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
## proportion edge on (native) bombus abundance 
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
## doy on (native) bombus abundance 
## ***********************************************************************
babun <-
  all.cond.effects[["nativebeeabundance.nativebeeabundance_julian_date"]]

#ggplot
babun.doy <- 
  
  #plot raw data
  ggplot(new.net, aes(x = julian_date, y = native_bee_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Day of Year", y = "Native *Bombus* abundance") +
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
## prepping for newdata draws -- bombus richness and impatiens abund models
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
## veg abundance and bombus richness
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
## veg diversity and bombus richness
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
## prop blueberry and bombus richness
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
## proportion edge and bombus richness
## ***********************************************************************
brich <-
  all.cond.effects[["bombusrichness.bombusrichness_prop_edge"]]

#ggplot
brich.edge <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_edge, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.1, width = 0.1) +
  labs(x = "Proportion edge area", y = "") +
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
## julian date on bombus richness
## ***********************************************************************
brich <-
  all.cond.effects[["bombusrichness.bombusrichness_julian_date"]]

#ggplot
brich.doy <- 
  
  #plot raw data
  ggplot(new.net, aes(x = julian_date, y = bombus_richness)) +
  geom_jitter(cex = 2, alpha = 0.2, height = 0.2, width = 0.2) +
  labs(x = "Day of year", y = "*Bombus* species richness") +
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
## proportion blueberry on bombus diversity
## ***********************************************************************
bdiv <-
  all.cond.effects[["bombusshannondiversity.bombusshannondiversity_prop_blueberry"]]

#ggplot
bdiv.bberry <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_blueberry, y = bombus_shannon_diversity)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Proportion blueberry (500m buffer)", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.blueberry,
    labels =  labs.blueberry) +
  scale_y_continuous() +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) #+
  geom_text(data = data.frame(
    x = c(2.5),
    y = c(1.5),
    label = c("n.s.")
  ), aes(x=x, y=y, label=label),
  color = "black",
  size=5)
bdiv.bberry

## ***********************************************************************
## proportion edge on bombus diversity
## ***********************************************************************
bdiv <-
  all.cond.effects[["bombusshannondiversity.bombusshannondiversity_prop_edge"]]

#ggplot
bdiv.edge <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_edge, y = bombus_shannon_diversity)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Proportion edge area (500m buffer)", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = axis.edge,
    labels =  labs.edge) +
  scale_y_continuous() +
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) #+
  geom_text(data = data.frame(
    x = c(2.5),
    y = c(1.5),
    label = c("n.s.")
  ), aes(x=x, y=y, label=label),
  color = "black",
  size=5)

bdiv.edge

## ***********************************************************************
## impatiens abundance ~ floral abundance
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
                               alpha=0.5), fill = "red") #+
  geom_text(data = data.frame(
    x = c(2.5),
    y = c(40),
    label = c("***p > 0 = 1")
  ), aes(x=x, y=y, label=label),
  color = "black",
  size=5)

iabund.fabund

## ***********************************************************************
## impatiens abundance ~ floral diversity
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
## proportion edge on impatiens abundance 
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
## proportion blueberry on impatiens abundance 
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
## impatiens abundance ~ doy
## ***********************************************************************
iabund <-
  all.cond.effects[["impatiensabundance.impatiensabundance_julian_date"]]

#ggplot
iabund.doy <- 
  
  #plot raw data
  ggplot(new.net, aes(x = julian_date, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Day of year", y = "*B*. *impatiens* abundance") +
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
## prepping for newdata draws -- parasitism models
## ***********************************************************************
data.par <- fvimp_brmsdf[fvimp_brmsdf$subsetPar == TRUE, ]

#for plotting raw data as rates
data.par %>% 
  group_by(sample_id, final_id,status) %>%
  summarize(site_hascrithidia = sum(hascrithidia)/n(),
            site_hasnosema = sum(hasnosema)/n(),
            site_apicystis = sum(apicystis)/n(),
            impatiens_abundance = mean(impatiens_abundance),
            prop_blueberry = mean(prop_blueberry),
            julian_date = mean(julian_date),
            floral_diversity = mean(floral_diversity),
            native_bee_abundance = mean(native_bee_abundance),
            numbees = n()) -> site_rates

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
  labs(x = "Day of year", y = "*Crithidia spp.* prevalence") +
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
all.cond.effects.interaction <- conditional_effects(fit.bombus.all.beta.inter)

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
  labs(x = "", y = "") +
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
  labs(x = "", y = "*Crithidia spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
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
        text = element_text(size=16)) #+
  geom_text(data = data.frame(
    x = c(9),
    y = c(1.1),
    label = c("***p > 0 = 0.99")
  ), aes(x=x, y=y, label=label),
  color = "black",
  fontface = "bold",
  size=4)

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
  labs(x = "Day of year", y = "*Vairimorpha spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = axis.doy,
                     labels =  labs.doy) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))
  
  
  #nosema ~ doy without raw data
  #plot model prediction with credible interval
  hasnosema.doy = ggplot(hasnosema, aes(x = julian_date, y = estimate__)) + 
    geom_line(aes(x = julian_date, y=estimate__), linetype = "dashed") +
    geom_ribbon(aes(ymin = lower__, ymax = upper__,
                    alpha=0.5), fill = "gray") +
    labs(x = "Day of year", y = "*Vairimorpha spp.* prevalence") +
    theme_ms() +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(-0.5,11)) +
    scale_y_continuous(limits = c(0,1)) +
    theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
          axis.title.y = ggtext::element_markdown(size=16),
          text = element_text(size=16))

hasnosema.doy


## ***********************************************************************
## has nosema ~ status * impatiens abundance
## ***********************************************************************
hasnosema <-
  all.cond.effects.interaction[["hasnosema.hasnosema_impatiens_abundance:status"]]

#ggplot
hasnosema.imp.status <- 
  
  #plot model prediction with credible interval
  ggplot(hasnosema, aes(x = impatiens_abundance, y = estimate__, color = status)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = status,
                  alpha=1)) +
  geom_line(aes(x = impatiens_abundance, y=estimate__), linewidth = 1.5) +
  scale_fill_manual(values = c("black", "lightblue")) +
  scale_color_manual(values = c("black", "lightblue")) +
  
  #plot raw data
  #geom_jitter(data = data.par, aes(x=impatiens_abundance,y = hasnosema), height = 0.02, width = 0.2) +
  labs(x = "*B*. *impatiens* abundance", y = "") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,11)) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.imp.status

## ***********************************************************************
## has nosema ~ status * native bee abundance
## ***********************************************************************
hasnosema <-
  all.cond.effects.interaction[["hasnosema.hasnosema_native_bee_abundance:status"]]
hasnosema = hasnosema %>% filter((status != "nonnative") | (native_bee_abundance <= max(data.par$native_bee_abundance[data.par$final_id == "Bombus_impatiens"])))


#ggplot
hasnosema.babund.status <- 
  
  #plot model prediction with credible interval
  ggplot(hasnosema, aes(x = native_bee_abundance, y = estimate__, color = status)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = status,
                  alpha=1)) +
  geom_line(aes(x = native_bee_abundance, y=estimate__), linewidth = 1.5) +
  scale_fill_manual(values = c("black", "lightblue")) +
  scale_color_manual(values = c("black", "lightblue")) +
  
  #plot raw data
  #geom_jitter(data = site_rates, aes(x=native_bee_abundance,y = site_hasnosema, cex = numbees)) +
  labs(x = "Native *Bombus* abundance", y = "*Vairimorpha spp.* prevalence", fill = "Prevalence in", color = "Prevalence in") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) +
  theme(legend.title = element_text(hjust = 0.5), legend.position = "bottom") +
  guides(fill = guide_legend(title = "Prevalence in"), size = "none", shape = "none", alpha = "none") +
  geom_text(data = data.frame(
    x = c(35),
    y = c(1.1),
    label = c("**")
  ), aes(x=x, y=y, label=label),
  color = "black",
  fontface = "bold",
  size=12)

hasnosema.babund.status

g <- ggplotGrob(hasnosema.babund.status)
legend_index <- which(g$layout$name == "guide-box")
shared_legend <- g$grobs[[17]]

# Remove legends from all plots
hasnosema.babund.status <- hasnosema.babund.status + theme(legend.position = "none")



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
  labs(x = "Day of year", y = "*Apicystis spp.* prevalence") +
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

#apicytis without raw data
#plot model prediction with credible interval
apicystis.imp = ggplot(apicystis, aes(x = impatiens_abundance, y = estimate__)) + 
  geom_line(aes(x = impatiens_abundance, y=estimate__), linetype = "dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "gray") +
  labs(x = "", y = "*Apicystis spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(-0.5,11)) +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16)) #+
apicystis.imp


#apicystis with raw data
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
g <- ggplotGrob(apicystis.imp)
legend_index <- which(g$layout$name == "guide-box")
shared_legend <- g$grobs[[17]]

# Remove legends from all plots
apicystis.imp <- apicystis.imp + theme(legend.position = "none")


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
## make grid plots
## ***********************************************************************
# doy vs floral community
doyfloralgrid = grid.arrange(vabun.doy, vdiv.doy, ncol =2)
#add subplot labels
final_plot <- ggdraw() +
  draw_plot(doyfloralgrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b"), 
                  x = c(0, 0.52), 
                  y = c(0.97, 0.97))
print(final_plot)


#doy vs bee community
doybeegrid = grid.arrange(babun.doy, iabund.doy, brich.doy, ncol = 3)
#add subplot labels
final_plot <- ggdraw() +
  draw_plot(doybeegrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c"), 
                  x = c(0.02, 0.35, 0.68), 
                  y = c(0.97, 0.97, 0.97))
print(final_plot)


#floral abundance & landscape effects on bee community
bombusgrid = grid.arrange(babun.fabun, babun.fdiv, babun.bberry, babun.edge,
                                iabund.fabund, iabund.fdiv, iabund.bberry, iabund.edge,
                                brich.fabund, brich.fdiv, brich.bberry, brich.edge,
                                ncol = 4)
#add subplot labels
final_plot <- ggdraw() +
  draw_plot(bombusgrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), 
                  x = c(0.01, 0.26, 0.51, 0.76, 0.01, 0.26, 0.51, 0.76, 0.01, 0.26, 0.51, 0.76), 
                  y = c(0.97, 0.97, 0.97, 0.97, 0.64, 0.64, 0.64, 0.64, 0.3, 0.3, 0.3, 0.3))
print(final_plot)
#export width 1500, height 1000


#doy vs parasites
doyparasitegrid = grid.arrange(hascrith.doy, hasnosema.doy, apicystis.doy, ncol = 3)
#add subplot labels
final_plot <- ggdraw() +
  draw_plot(doyparasitegrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c"), 
                  x = c(0, 0.34, 0.68), 
                  y = c(0.97, 0.97, 0.97))
                  #y = c(0.97, 0.73, 0.48)) #if there's a legend
print(final_plot)



#impatiens vs parasites
impatiensparasitegrid = grid.arrange(hascrith.imp, hasnosema.imp, apicystis.imp, #shared_legend,
                               ncol = 3)
#add subplot labels
final_plot <- ggdraw() +
  draw_plot(impatiensparasitegrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c"), 
                  x = c(0, 0.33, 0.67), 
                  y = c(1, 1, 1))
                  #y = c(0.97, 0.73, 0.48))
print(final_plot)
#export width 1500, height 1000


#impatiens vs parasites (vertical)
impatiensparasitegrid = grid.arrange(hascrith.imp, hasnosema.imp, apicystis.imp, #shared_legend,
                                     ncol = 1)
#add subplot labels
final_plot <- ggdraw() +
  draw_plot(impatiensparasitegrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c"), 
                  x = c(0, 0, 0), 
                  y = c(0.95, 0.62, 0.3))
#y = c(0.97, 0.73, 0.48))
print(final_plot)


#interaction plots for crithidia
interactiongrid <- grid.arrange(
  arrangeGrob(hascrith.babund.status, hascrith.imp.status,
              hasnosema.babund.status, hasnosema.imp.status, ncol = 2), # 2x2 grid of plots
  shared_legend, # Add the shared legend
  ncol = 1, # 1 column: grid on top, legend below
  heights = c(10, 1) # Adjust the relative heights
)

#add subplot labels
final_plot <- ggdraw() +
  draw_plot(interactiongrid, 0.015, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c", "d"), 
                  x = c(0.01, 0.51, 0.01, 0.51), 
                  y = c(0.97, 0.97, 0.5, 0.5))
#y = c(0.97, 0.73, 0.48))
print(final_plot)


## ***********************************************************************
## make map plots
## ***********************************************************************

#Parasite prevalence across sites
#make a df where parasites are pooled at sample point
data.workers = filter(data.par, caste == "worker")
data.par %>% 
  group_by(sample_pt) %>%
  summarize(site_hascrithidia = sum(hascrithidia)/n(),
            site_hasnosema = sum(hasnosema)/n(),
            site_apicystis = sum(apicystis)/n(),
            site_any = sum(any_parasite)/n(),
            site = min(site),
            long = min(long),
            lat = min(lat),
            numbees = n()) -> parbysite

#make a df where surveys are pooled, but keep track of the number of survey events per transect
fvimp_brmsdf %>%
  filter(impSubset == TRUE) %>%
  group_by(sample_pt) %>%
  summarize(
    site = min(site),
    site_any = mean(impatiens_abundance),
    long = min(long),
    lat = min(lat),
    numbees = n() #this is actually the number of survey events 
    #but I'm giving it this name so I can change less code in my function
) -> impatiensabundancepersite

#do this again but with native bee abundance as the "central" var
fvimp_brmsdf %>%
  filter(Subset == TRUE) %>%
  group_by(sample_pt) %>%
  summarize(
    site = min(site),
    site_any = mean(native_bee_abundance),
    long = min(long),
    lat = min(lat),
    numbees = n() #this is actually the number of survey events but I'm giving it this name so I can change less code below
  ) -> beeabundancepersite

makeMapGrids(groupedbysite = parbysite, 
             sampling_effort = "Number of\nspecimens", 
             var_of_interest = "Parasite\nprevalence")


## ***********************************************************************
## variograms testing for spatial autocorrelation
## ***********************************************************************
fvimp_sub = filter(fvimp_brmsdf, impSubset == TRUE)
fvimp_subpar = filter(fvimp_brmsdf, apidae ==1)

#calculate residuals for floral abundance (no predictors)
fabun<- brm(
  bf(floral_abundance ~ 1),
  data = fvimp_sub,
  family = student(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
fabun_resid <- as.data.frame(residuals(fabun))
fvimp_sub$fabun_resid_nopred = fabun_resid$Estimate

#calculate residuals for floral abundance (with predictors)
fabun<- brm(
  bf(floral_abundance ~ julian_date + I(julian_date^2) + (1|sample_pt)),
  data = fvimp_sub,
  family = student(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
fabun_resid <- as.data.frame(residuals(fabun))
fvimp_sub$fabun_resid_pred = fabun_resid$Estimate

#calculate residuals for floral diversity (with no predictors)
fdiv <- brm(
  bf(floral_diversity ~ 1),
  data = fvimp_sub,
  family = student(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
fdiv_resid <- as.data.frame(residuals(fdiv))
fvimp_sub$fdiv_resid_nopred = fdiv_resid$Estimate

#calculate residuals for floral diversity (predictors)
fdiv <- brm(
  bf(floral_diversity ~ julian_date + I(julian_date^2) + (1|sample_pt)),
  data = fvimp_sub,
  family = student(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
fdiv_resid <- as.data.frame(residuals(fdiv))
fvimp_sub$fdiv_resid_pred = fdiv_resid$Estimate

#calculate residuals for bombus diversity (no predictors)
brich <- brm(
  bf(bombus_richness | trials(9) ~ 1),
  data = fvimp_sub,
  family = beta_binomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
brich_resid = as.data.frame(residuals(brich))
fvimp_sub$brich_resid_nopred = brich_resid$Estimate

#calculate residuals for bombus diversity (with predictors)
brich <- brm(
  bf(bombus_richness | trials(9) ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
       (1|sample_pt)),
  data = fvimp_sub,
  family = beta_binomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
brich_resid = as.data.frame(residuals(brich))
fvimp_sub$brich_resid_pred = brich_resid$Estimate

#calculate residuals for wild bee abundance (no predictors)
babun <- brm(
  bf(native_bee_abundance ~ 1),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
babun_resid = as.data.frame(residuals(babun))
fvimp_sub$babun_resid_nopred = babun_resid$Estimate

#calculate residuals for wild bee abundance (with predictors)
babun <- brm(
  bf(native_bee_abundance ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
       (1|sample_pt)),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
babun_resid = as.data.frame(residuals(babun))
fvimp_sub$babun_resid_pred = babun_resid$Estimate

#calculate residuals for impatiens abundance (no predictors)
iabun <- brm(
  bf(impatiens_abundance ~ 1),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
iabun_resid = as.data.frame(residuals(iabun))
fvimp_sub$iabun_resid_nopred = iabun_resid$Estimate

#calculate residuals for impatiens abundance (with predictors)
iabun <- brm(
  bf(impatiens_abundance ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
       (1|sample_pt)),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
iabun_resid = as.data.frame(residuals(iabun))
fvimp_sub$iabun_resid_pred = iabun_resid$Estimate


########################################
### now parasite models
########################################

#calculate residuals for crithidia prevalence (no predictors)
hascrith <- brm(
  bf(hascrithidia ~ 1),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
crith_resid <- as.data.frame(residuals(hascrith))
fvimp_subpar$crith_resid_nopred = crith_resid$Estimate

#calculate residuals for crithidia prevalence (with predictors)
hascrith <- brm(
  bf(hascrithidia ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
crith_resid <- as.data.frame(residuals(hascrith))
fvimp_subpar$crith_resid_pred = crith_resid$Estimate

#calculate residuals for apicystis prevalence (no predictors)
hasapi <- brm(
  bf(apicystis ~ 1 ),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
api_resid <- as.data.frame(residuals(hasapi))
fvimp_subpar$api_resid_nopred = api_resid$Estimate

#calculate residuals for apicystis prevalence (with predictors)
hasapi <- brm(
  bf(apicystis ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
api_resid <- as.data.frame(residuals(hasapi))
fvimp_subpar$api_resid_pred = api_resid$Estimate

#calculate residuals for nosema prevalence (no predictors)
hasnos <- brm(
  bf(hasnosema ~ 1 ),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
nos_resid <- as.data.frame(residuals(hasnos))
fvimp_subpar$nos_resid_nopred = nos_resid$Estimate

#calculate residuals for nosema prevalence (with predictors)
hasnos <- brm(
  bf(hasnosema ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
nos_resid <- as.data.frame(residuals(hasnos))
fvimp_subpar$nos_resid_pred = nos_resid$Estimate

#############################################
### save dataframe 
#############################################
save(fvimp_sub, fvimp_subpar,
     file="saved/data_with_residuals.Rdata")

load(file="saved/data_with_residuals.Rdata")

#############################################
### make data frames spatially explicit
#############################################

# Convert the data frames into sf objects
fvimp_sub <- st_as_sf(fvimp_sub, coords = c("long", "lat"), crs = 4326)
fvimp_subpar <- st_as_sf(fvimp_subpar, coords = c("long", "lat"), crs = 4326)

# Transform the coordinate reference system to EPSG:900913
fvimp_sub900913 <- st_transform(fvimp_sub, crs = 3857)
fvimp_subpar900913 <- st_transform(fvimp_subpar, crs = 3857)


############################################
### variograms with no predictors
############################################

#make variogram for floral abundance
v_fabun.resid.nopred <- variogram(fabun_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_fabun.resid.nopred_plot = ggplot(as.data.frame(v_fabun.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Floral Abundance",
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )
v_fabun.resid.nopred_plot

#make variogram for floral diversity
v_fdiv.resid.nopred <- variogram(fdiv_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_fdiv.resid.nopred_plot = ggplot(as.data.frame(v_fdiv.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Floral Diversity",
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )
v_fdiv.resid.nopred_plot


#make variogram for bombus diversity
v_brich.resid.nopred <- variogram(brich_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_brich.resid.nopred_plot = ggplot(as.data.frame(v_brich.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*Bombus* Species Richness"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_brich.resid.nopred_plot


#make variogram for native bombus abundance
v_babun.resid.nopred <- variogram(babun_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_babun.resid.nopred_plot = ggplot(as.data.frame(v_babun.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("Native *Bombus* Abundance"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_babun.resid.nopred_plot

#make variogram for impatiens abundance
v_iabun.resid.nopred <- variogram(iabun_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_iabun.resid.nopred_plot = ggplot(as.data.frame(v_iabun.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*B*. *impatiens* Abundance"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_iabun.resid.nopred_plot

#make variogram for crithidia prevalence
v_crith.resid.nopred <- variogram(crith_resid_nopred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_crith.resid.nopred_plot = ggplot(as.data.frame(v_crith.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*Crithidia spp* Prevalence"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_crith.resid.nopred_plot

#make variogram for apicystis prevalence
v_api.resid.nopred <- variogram(api_resid_nopred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_api.resid.nopred_plot = ggplot(as.data.frame(v_api.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*Apicystis spp* Prevalence"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_api.resid.nopred_plot

#make variogram for nosema prevalence
v_nos.resid.nopred <- variogram(nos_resid_nopred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_nos.resid.nopred_plot = ggplot(as.data.frame(v_nos.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*Vairimorpha spp* Prevalence"),
    x = "Distance (meters)",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_nos.resid.nopred_plot


#############################################
### variograms with predictors
#############################################
#make variogram for floral abundance
v_fabun.resid.pred <- variogram(fabun_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_fabun.resid.pred_plot = ggplot(as.data.frame(v_fabun.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )
v_fabun.resid.pred_plot

#make variogram for floral diversity
v_fdiv.resid.pred <- variogram(fdiv_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_fdiv.resid.pred_plot = ggplot(as.data.frame(v_fdiv.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )
v_fdiv.resid.pred_plot


#make variogram for bombus diversity
v_brich.resid.pred <- variogram(brich_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_brich.resid.pred_plot = ggplot(as.data.frame(v_brich.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_brich.resid.pred_plot


#make variogram for native bombus abundance
v_babun.resid.pred <- variogram(babun_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_babun.resid.pred_plot = ggplot(as.data.frame(v_babun.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_babun.resid.pred_plot

#make variogram for impatiens abundance
v_iabun.resid.pred <- variogram(iabun_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_iabun.resid.pred_plot = ggplot(as.data.frame(v_iabun.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_iabun.resid.pred_plot

#make variogram for crithidia prevalence
v_crith.resid.pred <- variogram(crith_resid_pred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_crith.resid.pred_plot = ggplot(as.data.frame(v_crith.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_crith.resid.pred_plot

#make variogram for apicystis prevalence
v_api.resid.pred <- variogram(api_resid_pred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_api.resid.pred_plot = ggplot(as.data.frame(v_api.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_api.resid.pred_plot

#make variogram for nosema prevalence
v_nos.resid.pred <- variogram(nos_resid_pred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_nos.resid.pred_plot = ggplot(as.data.frame(v_nos.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "Distance (meters)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_nos.resid.pred_plot




variogramgrid = wrap_plots(v_fdiv.resid.nopred_plot, v_fdiv.resid.pred_plot, 
                           v_fabun.resid.nopred_plot, v_fabun.resid.pred_plot,
                           v_brich.resid.nopred_plot, v_brich.resid.pred_plot, 
                           v_babun.resid.nopred_plot, v_babun.resid.pred_plot,
                           v_iabun.resid.nopred_plot, v_iabun.resid.pred_plot,
                           v_crith.resid.nopred_plot, v_crith.resid.pred_plot,
                           v_api.resid.nopred_plot, v_api.resid.pred_plot,
                           v_nos.resid.nopred_plot, v_nos.resid.pred_plot,
                      ncol = 2)

variogramgrid

#add subplot labels
final_plot <- ggdraw() +
  draw_plot(variogramgrid, 0.02, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p"),
                  x = c(0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52), 
                  y = c(0.99, 0.99, 0.86, 0.86, 0.72, 0.72, 0.6, 0.6, 0.48, 0.48, 0.37, 0.37, 0.25, 0.25, 0.14, 0.14))
print(final_plot)



################################
### Calculate Moran's I
################################
library(sf)
library(spdep)

# Set a distance threshold (e.g., within 1000 meters)
coords <- st_coordinates(fvimp_sub900913)  # Extract coordinates
nb <- dnearneigh(coords, d1 = 0, d2 = 1000)  # Neighbors within 1000 units


# Convert to spatial weights
#ignore groups that don't have neighbors within specified distance
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Perform Moran's I (on raw data)
moran_result <- moran.test(fvimp_sub900913$impatiens_abundance, lw)
print(moran_result)

#okay so I think there is spatial autocorrelation on the raw values which makes sense
#however this spatial autocorrelation is quite weak (I < 0.1 for all relationships)

#test spatial autocorrelation on the model residuals
# Perform Moran's I (on residuals)
moran_result <- moran.test(fvimp_sub900913$fabun_resid_pred, lw)
print(moran_result)

#modelling accounts for spatial autocorrelation!! all p-values are like 0.99 after
#modelling and Moran's I stays low :)
#the only exception to this is the Bombus diversity model -- the pre- and post-model
#residuals show similar patterns of *very low but significant* spatial autocorrelation
#this is not surprising given how absolutely trash the diversity model R2 value is