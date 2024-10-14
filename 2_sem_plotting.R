
## plotting based on tutorial:
## https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#different-kinds-of-averCanopyBin-predictions-with-multilevel-models

rm(list=ls())
source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")

## load model results and data
load(file="saved/AllModels_fv_bernoulli_phylo_prior.Rdata")

## ***********************************************************************
## exploratory plots / bar graphs
## ***********************************************************************

##parasite counts
orig.spec$bberry_bin <- cut(orig.spec$prop_blueberry,
                        breaks = c(0, .25, .75, 1), 
                        include.lowest = T, right = F)
parasite.count.table <- orig.spec %>%
  filter(apidae == 1) %>%
  select(barcode_id, apicystis, ascosphaera, cbombii, cexpoeki, 
         crithidiaspp, nbombii, nceranae, any_parasite, prop_blueberry) %>%
  pivot_longer(cols=c(apicystis, ascosphaera, cbombii,
                      crithidiaspp, cexpoeki, nceranae, nbombii, any_parasite),
               names_to = 'ParasiteName', values_to = 'HasParasite') %>%
  filter(HasParasite == 1)

#using #4C4C5A
parasite.hist <- parasite.count.table %>%
  ggplot(aes(y= fct_infreq(ParasiteName))) + 
  geom_bar(stat = 'count',
           aes(fill = factor(bberry_bin)), position = "dodge") +
  scale_fill_manual(values=c('#004D40', '#FFC107', 'lightblue'),
                    name = "Proportion blueberry", 
                    labels=c('Low', 'Medium', 'High')) +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, color = "black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(color = "black"),
        text = element_text(size=16)) +
  labs(y='Parasite', x='Number of \n Infected Individuals') +
  xlim(0,75) +
  scale_y_discrete(labels=c("nceranae"=expression(italic('Nosema ceranae')),
                            "nbombii"=expression(italic('Nosema bombi')),
                            'ascosphaera'=expression(paste(italic('Ascosphaera'), ' spp.')),
                            'apicystis'=expression(paste(italic('Apicystis'), ' spp.')),
                            "cbombii"=expression(italic('Crithidia bombi')),
                            "cexpoeki"=expression(italic('Crithidia expoeki')),
                            "crithidiaspp"=expression(paste(italic('Crithidia'), ' spp.')),
                            parse=TRUE))  + labs(fill = "Canopy openness")

parasite.hist


## ***********************************************************************
## prepping for newdata draws -- vegetation & native bee abundance models
## ***********************************************************************
new.net <- fvimp_brmsdf[fvimp_brmsdf$Subset == TRUE, ]
new.orig <- orig.spec[fvimp_brmsdf$Subset == TRUE, ]

# calculate all conditional effects
all.cond.effects <- conditional_effects(fit.bombus.all.prior.nointeractions)

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
                  alpha=0.3, fill = "Reds")) +
  
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
babun <-
  all.cond.effects[["nativebeeabundance.nativebeeabundance_floral_abundance"]]

#ggplot
babun.fabun <- 
  
  #plot model prediction with credible interval
  ggplot(babun, aes(x = floral_abundance,
                    y = estimate__)) +
  geom_line(aes(x = floral_abundance, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.3), fill = "red") +
  
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
  theme(axis.title.x = element_text(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

#add annotation with posterior estimates / support level
annotation <- data.frame(
  x = c(0.5),
  y = c(40),
  label = c("***p > 0 = 1")
)

# Add text
babun.fabun + geom_text(data=annotation, aes( x=x, y=y, label=label),                 , 
              color="orange", 
              size=7 , angle=45, fontface="bold" )
babun.fabun

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
  labs(x = "Proportion blueberry \n (500m buffer)", y = "Native *Bombus* abundance") +
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
                  alpha=0.5), fill = "purple")
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
  labs(x = "Proportion edge \n (500m buffer)", y = "Native *Bombus* abundance") +
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
                                alpha=0.5), fill = "purple")
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
## prepping for newdata draws -- bombus div and impatiens abund models
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
## veg abundance and bombus diversity
## ***********************************************************************
bdiv <-
  all.cond.effects[["bombusshannondiversity.bombusshannondiversity_floral_abundance"]]

#ggplot
bdiv.fabund <- 
  
  #plot raw data
  ggplot(new.net, aes(x = floral_abundance, y = bombus_shannon_diversity)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Floral abundance (log)", y = "*Bombus* Shannon diversity") +
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
  geom_line(data = bdiv, aes(x = floral_abundance, y=estimate__)) +
  geom_ribbon(data = bdiv, aes(ymin = lower__, ymax = upper__,
                                alpha=0.5), fill = "purple")
bdiv.fabund


## ***********************************************************************
## julian date on bombus diversity
## ***********************************************************************
bdiv <-
  all.cond.effects[["bombusshannondiversity.bombusshannondiversity_julian_date"]]

#ggplot
bdiv.doy <- 
  
  #plot raw data
  ggplot(new.net, aes(x = julian_date, y = bombus_shannon_diversity)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Day of year", y = "*Bombus* Shannon diversity") +
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
  geom_line(data = bdiv, aes(x = julian_date, y=estimate__)) +
  geom_ribbon(data = bdiv, aes(ymin = lower__, ymax = upper__,
                               alpha=0.5), fill = "red")

bdiv.doy


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
  labs(x = "Proportion blueberry (500m buffer)", y = "*Bombus* Shannon diversity") +
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
  geom_line(data = bdiv, aes(x = prop_blueberry, y=estimate__)) +
  geom_ribbon(data = bdiv, aes(ymin = lower__, ymax = upper__,
                               alpha=0.5), fill = "grey")
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
  labs(x = "Proportion edge area (500m buffer)", y = "*Bombus* Shannon diversity") +
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
  geom_line(data = bdiv, aes(x = prop_edge, y=estimate__)) +
  geom_ribbon(data = bdiv, aes(ymin = lower__, ymax = upper__,
                               alpha=0.5), fill = "grey")

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
  labs(x = "Floral abundance (log)", y = "*Bombus impatiens* abundance") +
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
## proportion edge on impatiens abundance 
## ***********************************************************************
iabund <-
  all.cond.effects[["impatiensabundance.impatiensabundance_prop_edge"]]

#ggplot
iabund.edge <- 
  
  #plot raw data
  ggplot(new.net, aes(x = prop_edge, y = impatiens_abundance)) +
  geom_jitter(cex = 2, alpha = 0.2) +
  labs(x = "Proportion edge area (500m buffer)", y = "*Bombus impatiens* abundance") +
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
  geom_line(data = iabund, aes(x = prop_edge, y=estimate__)) +
  geom_ribbon(data = iabund, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.5), fill = "purple")

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
  labs(x = "Proportion blueberry (500m buffer)", y = "*Bombus impatiens* abundance") +
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
  geom_line(data = iabund, aes(x = prop_blueberry, y=estimate__)) +
  geom_ribbon(data = iabund, aes(ymin = lower__, ymax = upper__,
                                 alpha=0.5), fill = "blue")

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
  labs(x = "Day of year", y = "*Bombus impatiens* abundance") +
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
                                 alpha=0.5), fill = "red")

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
  geom_jitter(data = site_rates, aes(x=impatiens_abundance,y = site_hascrithidia, cex = numbees, color = final_id)) +
  labs(x = "*Bombus impatiens* abundance", y = "*Crithidia spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,11)) +
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
  geom_jitter(data = site_rates, aes(x=julian_date,y = site_hascrithidia, cex = numbees, color = final_id)) +
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
all.cond.effects.interaction <- conditional_effects(fit.bombus.all.prior)

hascrith <-
  all.cond.effects.interaction[["hascrithidia.hascrithidia_impatiens_abundance:status"]]

#ggplot
hascrith.imp.status <- 
  
  #plot model prediction with credible interval
  ggplot(hascrith, aes(x = impatiens_abundance, y = estimate__, color = status)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = status,
                  alpha=0.2)) +
  geom_line(aes(x = impatiens_abundance, y=estimate__)) +
  scale_fill_manual(values = c("pink", "blue")) +
  #scale_color_manual(values = c("black", "blue")) +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=impatiens_abundance,y = site_hascrithidia, cex = numbees)) +
  labs(x = "*Bombus impatiens* abundance", y = "*Crithidia spp.* prevalence") +
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

#ggplot
hascrith.babund.status <- 
  
  #plot model prediction with credible interval
  ggplot(hascrith, aes(x = native_bee_abundance, y = estimate__, color = status)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = status,
                  alpha=0.2)) +
  geom_line(aes(x = native_bee_abundance, y=estimate__)) +
  scale_fill_manual(values = c("pink", "blue")) +
  #scale_color_manual(values = c("black", "blue")) +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=native_bee_abundance,y = site_hascrithidia, cex = numbees)) +
  labs(x = "Native *Bombus* abundance", y = "*Crithidia spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

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
                  alpha=0.5), fill = "purple") +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=impatiens_abundance,y = site_hasnosema, cex = numbees, color = final_id)) +
  labs(x = "*Bombus impatiens* abundance", y = "*Nosema spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,11)) +
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
  
  #plot model prediction with credible interval
  ggplot(hasnosema, aes(x = julian_date, y = estimate__)) + 
  geom_line(aes(x = julian_date, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "blue") +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=julian_date,y = site_hasnosema, cex = numbees, color = final_id)) +
  labs(x = "Day of year", y = "*Nosema spp.* prevalence") +
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
## has nosema ~ blueberry
## ***********************************************************************
hasnosema <-
  all.cond.effects[["hasnosema.hasnosema_prop_blueberry"]]

#ggplot
hasnosema.bberry <- 
  
  #plot model prediction with credible interval
  ggplot(hasnosema, aes(x = prop_blueberry, y = estimate__)) + 
  geom_line(aes(x = prop_blueberry, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "blue") +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=prop_blueberry,y = site_hasnosema, cex = numbees, color = final_id)) +
  labs(x = "Proportion blueberry (500m buffer)", y = "*Nosema spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = axis.blueberry,
                     labels =  labs.blueberry) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.bberry


## ***********************************************************************
## has nosema ~ blueberry (raw data at 0/1)
## ***********************************************************************
hasnosemaeffects <-
  all.cond.effects[["hasnosema.hasnosema_prop_blueberry"]]
data.par$nosemaprev = data.par$hasnosema
#ggplot
hasnosema.bberry <- 
  
  #plot model prediction with credible interval
  ggplot(hasnosemaeffects, aes(x = prop_blueberry, y = estimate__)) + 
  geom_line(aes(x = prop_blueberry, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.4), fill = "blue") +
  
  #plot raw data
  geom_jitter(data = data.par, aes(x=prop_blueberry, y = nosemaprev, alpha = 0.4), height = 0.02, width = 0.05) +
  labs(x = "Proportion blueberry (500m buffer)", y = "*Nosema spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = axis.blueberry,
                     labels =  labs.blueberry) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.bberry



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
                  alpha=0.2)) +
  geom_line(aes(x = impatiens_abundance, y=estimate__)) +
  scale_fill_manual(values = c("pink", "blue")) +
  #scale_color_manual(values = c("black", "blue")) +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=impatiens_abundance,y = site_hasnosema, cex = numbees)) +
  labs(x = "*Bombus impatiens* abundance", y = "*Nosema spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "bottom") +
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

#ggplot
hasnosema.babund.status <- 
  
  #plot model prediction with credible interval
  ggplot(hasnosema, aes(x = native_bee_abundance, y = estimate__, color = status)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = status,
                  alpha=0.2)) +
  geom_line(aes(x = native_bee_abundance, y=estimate__)) +
  scale_fill_manual(values = c("pink", "blue")) +
  #scale_color_manual(values = c("black", "blue")) +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=native_bee_abundance,y = site_hasnosema, cex = numbees)) +
  labs(x = "Native *Bombus* abundance", y = "*Nosema spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

hasnosema.babund.status

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
  geom_jitter(data = site_rates, aes(x=julian_date,y = site_apicystis, cex = numbees, color = final_id)) +
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
                  alpha=0.5), fill = "purple") +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=floral_diversity,y = site_apicystis, cex = numbees, color = final_id)) +
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

#ggplot
apicystis.imp <- 
  
  #plot model prediction with credible interval
  ggplot(apicystis, aes(x = impatiens_abundance, y = estimate__)) + 
  geom_line(aes(x = impatiens_abundance, y=estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.5), fill = "grey") +
  
  #plot raw data
  geom_jitter(data = site_rates, aes(x=impatiens_abundance,y = site_apicystis, cex = numbees, color = final_id)) +
  labs(x = "Floral diversity", y = "*Apicystis spp.* prevalence") +
  theme_ms() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,11)) +
  scale_y_continuous() +
  theme(axis.title.x = ggtext::element_markdown(size = 16), #change to element_blank() for grid plots!
        axis.title.y = ggtext::element_markdown(size=16),
        text = element_text(size=16))

apicystis.imp

## ***********************************************************************
## make grid plots
## ***********************************************************************

#bee / floral community + landscape
#make one big grid with all of them
blank <- grid.rect(gp=gpar(col="white"))
allsignificant = grid.arrange(doy.veg.abund, doy.veg.div, blank, blank,
                              bombus.abund.doy, bombus.abund.fabund, bombus.abund.bberry, bombus.abund.edge,
                              imp.ab.doy, imp.ab.veg, blank, blank,
                              bombus.div.date, bombus.div.plant, blank, blank,
                              hascrith.doy, hascrith.imp, blank, blank,
                              cbombi.doy, cbombi.imp, blank, blank,
                              hasnosema.doy, blank, blank, blank,
                              apicystis.doy, apicystis.fdiv, blank, blank,
                              anyparasite.doy, blank, blank, blank,
                              ncol = 4)
ggsave(allsignificant, file="figures/regressions/allsignificant.pdf",
       height=40, width=30)

#storyboard figures
# doy vs floral/bee community
blank <- grid.rect(gp=gpar(col="white"))
doycommunitygrid = grid.arrange(vabun.doy, vdiv.doy,
                              babun.doy, bdiv.doy,
                              iabund.doy, blank,
                              ncol = 2)
ggsave(doycommunitygrid, file="figures/regressions/doyvscommunity.pdf",
       height=40, width=30)

#other drivers of bee community
blank <- grid.rect(gp=gpar(col="white"))
bombusgrid = grid.arrange(babun.fabun, babun.bberry, babun.edge,
                                iabund.fabund, iabund.bberry, iabund.edge,
                                bdiv.fabund, bdiv.bberry, bdiv.edge,
                                ncol = 3)
ggsave(doycommunitygrid, file="figures/regressions/bombusgrid.pdf",
       height=40, width=30)

#doy vs parasites
blank <- grid.rect(gp=gpar(col="white"))
doyparasitegrid = grid.arrange(hascrith.doy, hasnosema.doy, apicystis.doy,
                                ncol = 3)

#landscape/floral community vs parasites
blank <- grid.rect(gp=gpar(col="white"))
landscapeparasitegrid = grid.arrange(hasnosema.bberry, apicystis.fdiv,
                               ncol = 2)

#impatiens vs parasites
blank <- grid.rect(gp=gpar(col="white"))
impatiensparasitegrid = grid.arrange(hascrith.imp, hasnosema.imp, apicystis.imp,
                               ncol = 3)

#interaction plots for crithidia
blank <- grid.rect(gp=gpar(col='white'))
crithinteractiongrid = grid.arrange(hascrith.babund.status, hascrith.imp.status)




