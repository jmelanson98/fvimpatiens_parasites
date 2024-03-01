
## Prepares the data for model fitting (standardizes continuous
## variables, creates dummy variables to be used as weights to all
## different subsets of data to be used in different model levels),
## builds the models, and fits the models in brms. The model outputs
## are saved as tables, and chain diagnostic plots created.


## setwd('~/Dropbox (University of Oregon)/pnw_survey')

##jesse wd 
setwd("C:/Users/mycol/Dropbox (University of Oregon)/pnw_survey")

library(ggplot2)
library(tidyverse)
library(tidybayes)
library(ggthemes)
library(ggtext)
library(stringr)

## Script for plotting all of the important explanatory variables.

## plotting based on tutorial:
## https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#different-kinds-of-averCanopyBin-predictions-with-multilevel-models

rm(list=ls())
setwd("analyses/parasites")

source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")

## set to the number of cores you would like the models to run on
ncores <- 1

## load model results and data
load(file="saved/CrithidiaFitBombusMod_coast_phylo.Rdata")

## ***********************************************************************
## scatterplot of canopy cover v. dbh, separated by thins 
## ***********************************************************************
#have to do this before centering variables below 
stands <- orig.spec[orig.spec$Weights == 1, ]

#plot w thins highlighted
standsplot <- stands %>% 
  ggplot(aes(x=DomTreeDiam_cm, y=MeanCanopy, color=ThinStatus)) +
  geom_point() + 
  scale_color_manual(values = c("Y" = "gold", "N" = "maroon"),
                     labels=c('Not thinned', 'Thinned'),
                     name="Management history") +
  labs(x="Dominant tree class DBH (cm)", y="Canopy openness") +
  xlim(0,75) +
  theme_classic()+
  theme(
    axis.text.x = element_text(color="black"),
    axis.text.y = element_text(color="black"),
    axis.ticks = element_line(color = "black"))

standsplot

ggsave(standsplot, file="figures/canopyDBH.pdf",
       height=2, width=4)

# #
# thins <- stands[stands$ThinStatus == "Y",]
# range(thins$MeanCanopy, na.rm=TRUE)

## ***********************************************************************
## descriptive bar charts
## ***********************************************************************

#binning by canopy 
orig.spec$CanopyBin <- cut(orig.spec$MeanCanopy,
                     breaks = c(0, 25, 75, 100), 
                     include.lowest = T, right = F)
#reminder 2 self that canopy is unintuitive and 0-25 is actually lowest

orig.spec <- orig.spec[!is.na(orig.spec$CanopyBin),]

####bee genus bars
orig.spec <- orig.spec %>%
  filter(Genus == 'Bombus')
# 
# #need to summarize flowers per n type of stand in each category
# #find stands per stand type 
# standstype <- orig.spec %>% 
#   group_by(age, PlantGenusSpecies) %>% 
#   summarise(count = n())
# #only 61 stands?? 
# 
# #bees per stand type
# beesper <- orig.spec %>% 
#   group_by(age, GenusSpecies) %>% 
#   summarise(count = n())
# 
# 
# #modifying to show prop of total
# bombus.spp.bar <- ggplot(summary, 
#                          aes(y = fct_reorder(GenusSpecies, avg, .desc = TRUE),
#                              x = avg)) + 
#   geom_bar(stat = 'identity',
#            aes(fill = factor(age)), position = "dodge") +
#   scale_fill_manual(values=c('darkgreen', 'gold', 'lightblue'),
#                     name = "Canopy closure", 
#                     labels=c('Closed', 'Intermediate', 'Open')) +
#   theme_classic() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1, 
#                                    face ='italic', color = 'black'),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=14),
#         text = element_text(size=14)) +
#   labs(y=expression(paste('Species')), x='Average individuals 
#        collected per stand type')
# 
# bombus.spp.bar
# 
# 
# ggsave(bombus.spp.bar, file="figures/bombusSppbar.pdf",
#        height=5, width=7.5)


#ok nevermind let's try canopy and not treediam actually. binning... 
orig.spec$CanopyBin <- cut(orig.spec$MeanCanopy,
                     breaks = c(0, 25, 75, 100), 
                     include.lowest = T, right = F)
#reminder 2 self that canopy is unintuitive and 0-25 is closed

orig.spec <- orig.spec[!is.na(orig.spec$CanopyBin),]

####bee genus bars
orig.spec <- orig.spec %>%
  filter(Genus == 'Bombus')


#per canopy type mean and sd of veg abund/diversity
closed <- orig.spec[orig.spec$CanopyBin %in% c("[0,25)"),]
intermed <- orig.spec[orig.spec$CanopyBin %in% c("[25,75)"),]
open <- orig.spec[orig.spec$CanopyBin %in% c("[75,100]"),]




#need to summarize individuals per n type of stand in each category
#find stands per stand type 
standstype <- orig.spec %>% 
  group_by(CanopyBin) %>% 
  summarise(count = n_distinct(Stand))
#bees per stand type
beesper <- orig.spec %>% 
  group_by(CanopyBin, GenusSpecies) %>% 
  summarise(count = n())
#add total stands in each stand type, + divide 
summary <- beesper %>% 
  left_join(standstype, by="CanopyBin") %>% 
  mutate(avg = count.x/count.y) 


#modifying to show prop of total
bombus.spp.bar <- ggplot(summary, 
                         aes(y = fct_reorder(GenusSpecies, avg, .desc = TRUE),
                             x = avg)) + 
  geom_bar(stat = 'identity',
           aes(fill = factor(CanopyBin)), position = "dodge") +
  scale_fill_manual(values=c('#004D40', '#FFC107', 'lightblue'),
                    name = "Canopy closure", 
                    labels=c('Closed', 'Intermediate', 'Open')) +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, 
                                   face ='italic', color = 'black'),
        axis.text.x = element_text(color="black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.ticks = element_line(color = "black"),
        text = element_text(size=14)) +
  labs(y=expression(paste('Species')), x='Average individuals
       collected per canopy type')

bombus.spp.bar

ggsave(bombus.spp.bar, file="figures/bombusSppbar.pdf",
       height=3, width=5)

#plant stats 

#plants per stand type
# plantsper <- orig.spec %>% 
#   group_by(CanopyBin, PlantGenusSpecies) %>% 
#   summarise(count = n())
# 
# #add total stands in each stand type, + divide 
# summary <- plantsper %>% 
#   left_join(standstype, by="CanopyBin") %>% 
#   mutate(avg = count.x/count.y) 
# 
# plant.spp.bar <- ggplot(summary, 
#                          aes(y = fct_reorder(PlantGenusSpecies, avg, .desc = TRUE),
#                              x = avg)) + 
#   geom_bar(stat = 'identity',
#            aes(fill = factor(CanopyBin)), position = "dodge") +
#   theme_classic() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1, 
#                                    face ='italic', color = 'black'),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=14),
#         text = element_text(size=14)) +
#   labs(y=expression(paste('Species')), x='Average individuals 
#        collected per stand type')
# 
# plant.spp.bar
# 

# 
# #this vers shows total collected and not per stand 
# bombus.tot.bar <- ggplot(orig.spec, aes(y = fct_infreq(GenusSpecies))) + 
#   geom_bar(stat = 'count',
#            aes(fill = factor(age)), position = "dodge") +
#   scale_fill_manual(values=c('darkorange4', 'darkorange', 'gold')) +
#   theme_classic() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1, 
#                                    face ='italic', color = 'black'),
#         axis.title.y = element_text(size=16),
#         axis.title.x = element_text(size=16),
#         text = element_text(size=16)) +
#   labs(y=expression(paste('Species')), x='Number of \n Collected Individuals') +
#   geom_text(stat='count', aes(label=..count..), hjust=-0.5) + 
#   xlim(0, 100)
# 
# bombus.tot.bar
# 
# ggsave(bombus.tot.bar, file="figures/bombusTotbar.pdf",
#        height=5, width=7.5)

##parasite counts
parasite.count.table <- orig.spec %>%
  filter(Apidae == 1) %>%
  select(SpecimenID, ApicystisSpp, AscosphaeraSpp, CrithidiaBombi, CrithidiaExpoeki, 
         CrithidiaSpp, NosemaBombi, NosemaCeranae, ParasitePresence, ParasiteRichness, 
         CanopyBin) %>%
  pivot_longer(cols=c(ApicystisSpp, AscosphaeraSpp, CrithidiaBombi,
                      CrithidiaSpp, CrithidiaExpoeki, NosemaCeranae, NosemaBombi),
               names_to = 'ParasiteName', values_to = 'HasParasite') %>%
  filter(HasParasite == 1)

#using #4C4C5A
parasite.hist <- parasite.count.table %>%
  ggplot(aes(y= fct_infreq(ParasiteName))) + 
  geom_bar(stat = 'count',
           aes(fill = factor(CanopyBin)), position = "dodge") +
  scale_fill_manual(values=c('#004D40', '#FFC107', 'lightblue'),
                    name = "Canopy closure", 
                    labels=c('Closed', 'Intermediate', 'Open')) +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, color = "black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(color = "black"),
        text = element_text(size=16)) +
  labs(y='Parasite', x='Number of \n Infected Individuals') +
  xlim(0,75) +
  scale_y_discrete(labels=c("NosemaCeranae"=expression(italic('Nosema ceranae')),
                            "NosemaBombi"=expression(italic('Nosema bombi')),
                            'AscosphaeraSpp'=expression(paste(italic('Ascosphaera'), ' spp.')),
                            'ApicystisSpp'=expression(paste(italic('Apicystis'), ' spp.')),
                            "CrithidiaBombi"=expression(italic('Crithidia bombi')),
                            "CrithidiaExpoeki"=expression(italic('Crithidia expoeki')),
                            'CrithidiaSpp'=expression(paste(italic('Crithidia'), ' spp.')),
                            parse=TRUE))  + labs(fill = "Canopy openness")

parasite.hist

#for geom text count: 1, 7, 18, 32, 54, 55, 149

ggsave(parasite.hist, file="figures/coastparasiteSppHist.pdf",
       height=3, width=5)


## ***********************************************************************
## prepping for newdata draws 
## ***********************************************************************

new.net <- spec.net[spec.net$Weights == 1, ]
new.orig <- orig.spec[orig.spec$Weights == 1, ]

## ***********************************************************************
## forest age ~ floral diversity
## ***********************************************************************

newdata.fl.div <- tidyr::crossing(MeanCanopy =
                                    seq(min(new.net$MeanCanopy, na.rm = TRUE),
                                        max(new.net$MeanCanopy, na.rm = TRUE),
                                        length.out=10),
                                  Stand="100:Camp",
                                  Year = "2021",
                                  ThinStatus = "Y",
                                  DoyStart = mean(new.net$DoyStart, na.rm = TRUE)
)

pred_fldiv <- fit.bombus %>%
  epred_draws(newdata = newdata.fl.div,
              resp = "VegDiversity",
              allow_new_levels = TRUE)

labs.fl.div <- (pretty(c(0, new.orig$MeanCanopy), n=8))
axis.fl.div <-  standardize.axis(labs.fl.div,
                                 new.orig$MeanCanopy)

veg.div.stand <- ggplot(pred_fldiv, aes(x = MeanCanopy,
                                        y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(y = "Flowering plant diversity", x = "Canopy openness",
       fill = "Credible interval") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  geom_point(data=new.net,
             aes(x=MeanCanopy, y=VegDiversity, color = ThinStatus), cex=2) +
  scale_color_manual(values=c("#000000","#999999"))+
  scale_x_continuous(
    labels = labs.fl.div,
    breaks = axis.fl.div) 

veg.div.stand

ggsave(veg.div.stand, file="figures/vegdiv_stand.pdf",
       height=4, width=5)

## ***********************************************************************
## forest age and floral abundance
## ***********************************************************************

labs.fl.ab <- (pretty(c(0, new.orig$MeanCanopy), n=8))
axis.fl.ab <-  standardize.axis(labs.fl.ab,
                                new.orig$MeanCanopy)

newdata.fl.ab <- tidyr::crossing(MeanCanopy =
                                   seq(min(new.net$MeanCanopy, na.rm=TRUE),
                                       max(new.net$MeanCanopy, na.rm=TRUE),
                                       length.out=10),
                                 Stand="100:Camp",
                                 Year = "2021",
                                 DoyStart = mean(new.net$DoyStart, na.rm=TRUE)
)

pred_flab <- fit.bombus %>%
  epred_draws(newdata = newdata.fl.ab,
              resp = "VegAbundance",
              allow_new_levels = TRUE)

flower.ab.stand <- ggplot(pred_flab, aes(x = MeanCanopy,
                                         y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  scale_x_continuous(
    breaks = axis.fl.ab,
    labels = labs.fl.ab) +
  labs(y = str_wrap("Flowering plant abundance (log)", width =20),
       x = "Canopy openness",
       fill = "Credible interval") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  geom_point(data=new.net,
             aes(x=MeanCanopy, y=VegAbundance, color = ThinStatus), cex=2) +
  scale_color_manual(values=c("#000000","#999999")) 

flower.ab.stand

ggsave(flower.ab.stand, file="figures/vegabund_stand.pdf",
       height=4, width=5)

## ***********************************************************************
## veg diversity and doy
## ***********************************************************************

labs.doy <- (pretty(new.orig$DoyStart, n=8))
axis.doy <-  standardize.axis(labs.doy,
                              new.orig$DoyStart)

newdata.vegdiv <- tidyr::crossing(DoyStart =
                                      seq(min(new.net$DoyStart, na.rm=TRUE),
                                          max(new.net$DoyStart, na.rm=TRUE),
                                          length.out=10),
                                    Year = new.net$Year,
                                    ThinStatus = (new.net$ThinStatus),
                                    MeanCanopy = (new.net$MeanCanopy)
)

pred_vegdiv <- fit.bombus %>%
  epred_draws(newdata = newdata.vegdiv,
              resp = "VegDiversity",
              allow_new_levels = TRUE)

doy.veg.div <- ggplot(pred_vegdiv, aes(x = DoyStart, y =
                                           .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(y = str_wrap("Flowering plant diversity", width =20),
       x = "Day of year",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.doy,
    labels =  labs.doy) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  geom_point(data=new.net,
             aes(y=VegDiversity, x=DoyStart), cex=2)

doy.veg.div

# ggsave(doy.veg.div, file="figures/vegDiv_doy.pdf",
#        height=4, width=5)

## ***********************************************************************
## veg abundance and doy
## ***********************************************************************

newdata.vegdiv <- tidyr::crossing(DoyStart =
                                    seq(min(new.net$DoyStart, na.rm=TRUE),
                                        max(new.net$DoyStart, na.rm=TRUE),
                                        length.out=10),
                                  Year = new.net$Year,
                                  ThinStatus = (new.net$ThinStatus),
                                  MeanCanopy = (new.net$MeanCanopy)
)

pred_vegab <- fit.bombus %>%
  epred_draws(newdata = newdata.vegdiv,
              resp = "VegAbundance",
              allow_new_levels = TRUE)

doy.veg.ab <- ggplot(pred_vegab, aes(x = DoyStart, y =
                                         .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(y = str_wrap("Flowering plant abundance (log)", width =20),
       x = "Day of year",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.doy,
    labels =  labs.doy) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  geom_point(data=new.net,
             aes(y=VegAbundance, x=DoyStart), cex=2)

doy.veg.ab

# ggsave(doy.veg.ab, file="figures/vegab_doy.pdf",
#        height=4, width=5)

## ***********************************************************************
## forest age and bombus richness 
## ***********************************************************************

labs.bee.div <- (pretty(new.orig$MeanCanopy, n=8))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  new.orig$MeanCanopy)

newdata.beediv <- tidyr::crossing(MeanCanopy =
                                    seq(min(new.net$MeanCanopy, na.rm=TRUE),
                                        max(new.net$MeanCanopy, na.rm=TRUE),
                                        length.out=10),
                                  FlowerRareRichness=mean(new.net$FlowerRareRichness, na.rm=TRUE),
                                  Stand="100:Camp",
                                  Year = "2020",
                                  ThinStatus = "Y",
                                  DoyStart = mean(new.net$DoyStart, na.rm=TRUE),
                                  VegAbundance = mean(new.net$VegAbundance, na.rm=TRUE),
                                  VegDiversity = mean(new.net$VegDiversity, na.rm=TRUE),
                                  TempCStart = mean(new.net$TempCStart, na.rm=TRUE),
                                  ForageDist_km = mean(new.net$ForageDist_km, na.rm=TRUE)
)

pred_beediv <- fit.bombus %>%
  epred_draws(newdata = newdata.beediv,
              resp = "BombusRareRichness",
              allow_new_levels = TRUE)

bombus.div.canopy <- ggplot(pred_beediv, aes(x = MeanCanopy,
                                             y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = 'Greys') +
  labs(x = "Canopy openness", y = "Rarefied *Bombus* <br>species richness",
       fill = "Credible interval") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_point(data=spec.net,
             aes(x=MeanCanopy, y=BombusRareRichness, color = ThinStatus), cex=2) +
  scale_color_manual(values=c("#000000","#999999")) +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div)

bombus.div.canopy

ggsave(bombus.div.canopy, file="figures/bombusrich_stand.pdf",
       height=4, width=5)

## ***********************************************************************
## forest age on bombus abundance 
## ***********************************************************************

labs.bee.ab <- (pretty(c(0, new.orig$MeanCanopy), n=8))
axis.bee.ab <-  standardize.axis(labs.bee.ab,
                                 new.orig$MeanCanopy)

newdata.beeab <- tidyr::crossing(MeanCanopy =
                                   seq(min(new.net$MeanCanopy, na.rm=TRUE),
                                       max(new.net$MeanCanopy, na.rm=TRUE),
                                       length.out=10),
                                 FlowerRareRichness=mean(new.net$FlowerRareRichness, na.rm=TRUE),
                                 Stand ="100:Camp",
                                 Year = "2020",
                                 ThinStatus = "Y",
                                 DoyStart = mean(new.net$DoyStart, na.rm=TRUE),
                                 VegAbundance = mean(new.net$VegAbundance, na.rm=TRUE),
                                 ForageDist_km = mean(new.net$ForageDist_km, na.rm=TRUE),
                                 TempCStart = mean(new.net$TempCStart, na.rm=TRUE),
                                 VegDiversity = mean(new.net$VegDiversity, na.rm=TRUE)
)

pred_beeab <- fit.bombus %>%
  epred_draws(newdata = newdata.beeab,
              resp = "BombusAbundance",
              allow_new_levels = TRUE)

bombus.ab.canopy <- ggplot(pred_beeab, aes(x = MeanCanopy, 
                                           y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette="Greys") +
  labs(x = "Canopy openness", y = "*Bombus* abundance",
       fill = "Credible interval") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_point(data=new.net,
             aes(x=MeanCanopy, y=BombusAbundance, color = ThinStatus), cex=2) +
  scale_color_manual(values=c("#000000","#999999")) +
  scale_x_continuous(
    breaks = axis.bee.ab,
    labels = labs.bee.ab) 

bombus.ab.canopy

ggsave(bombus.ab.canopy, file="figures/bombusabund_stand.pdf",
       height=4, width=5)

## ***********************************************************************
## floral abundance on bombus abundance 
## ***********************************************************************

#setting axis manually 
labs.bee.ab <- (pretty(c(0, new.orig$VegAbundance), n=8))
axis.bee.ab <- c(0:11)

newdata.beeab <- tidyr::crossing(VegAbundance =
                                   seq(min(new.net$VegAbundance, na.rm=TRUE),
                                       max(new.net$VegAbundance, na.rm=TRUE),
                                       length.out=10),
                                 FlowerRareRichness=mean(new.net$FlowerRareRichness, na.rm=TRUE),
                                 Stand ="100:Camp",
                                 Year = "2020",
                                 ThinStatus = "Y",
                                 DoyStart = mean(new.net$DoyStart, na.rm=TRUE),
                                 ForageDist_km = mean(new.net$ForageDist_km, na.rm=TRUE),
                                 TempCStart = mean(new.net$TempCStart, na.rm=TRUE),
                                 VegDiversity = mean(new.net$VegDiversity, na.rm=TRUE)
)

pred_beeab <- fit.bombus %>%
  epred_draws(newdata = newdata.beeab,
              resp = "BombusAbundance",
              allow_new_levels = TRUE)

bombus.ab.plant <- ggplot(pred_beeab, aes(x = VegAbundance, 
                                           y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette="Blues") +
  labs(x = "Flowering plant abundance (log)", y = "*Bombus* abundance",
       fill = "Credible interval") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_point(data=new.net,
             aes(x=VegAbundance, y=BombusAbundance), cex=2) +
  scale_x_continuous(
    breaks = axis.bee.ab,
    labels = labs.bee.ab) 

bombus.ab.plant

ggsave(bombus.ab.plant, file="figures/bombusabund_plant.pdf",
       height=4, width=5)

## ***********************************************************************
## veg diversity and bombus richness
## ***********************************************************************

labs.bee.div <- (pretty(new.orig$VegDiversity, n=8))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  new.orig$VegDiversity)

newdata.beediv <- tidyr::crossing(VegDiversity =
                                    seq(min(new.net$VegDiversity, na.rm=TRUE),
                                        max(new.net$VegDiversity, na.rm=TRUE),
                                        length.out=10),
                                  FlowerRareRichness=mean(new.net$FlowerRareRichness, na.rm=TRUE),
                                  Stand="100:Camp",
                                  Year = "2020",
                                  ThinStatus = "Y",
                                  DoyStart = mean(new.net$DoyStart, na.rm=TRUE),
                                  VegAbundance = mean(new.net$VegAbundance, na.rm=TRUE),
                                  TempCStart = mean(new.net$TempCStart, na.rm=TRUE),
                                  ForageDist_km = mean(new.net$ForageDist_km, na.rm=TRUE)
)

pred_beediv <- fit.bombus %>%
  epred_draws(newdata = newdata.beediv,
              resp = "BombusRareRichness",
              allow_new_levels = TRUE)

bombus.div.plant <- ggplot(pred_beediv, aes(x = VegDiversity,
                                             y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = 'Greys') +
  labs(x = "Flowering plant diversity", y = "Rarefied *Bombus* <br>species richness",
       fill = "Credible interval") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_point(data=spec.net,
             aes(x=VegDiversity, y=BombusRareRichness), cex=2) +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div)

bombus.div.plant

ggsave(bombus.div.plant, file="figures/bombusrich_veg.pdf",
       height=4, width=5)


## ***********************************************************************
## temp and bombus richness 
## ***********************************************************************

labs.bee.div <- (pretty(new.orig$TempCStart, n=8))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  new.orig$TempCStart)

newdata.beediv <- tidyr::crossing(TempCStart =
                                    seq(min(new.net$TempCStart, na.rm=TRUE),
                                        max(new.net$TempCStart, na.rm=TRUE),
                                        length.out=10),
                                  FlowerRareRichness=mean(new.net$FlowerRareRichness, na.rm=TRUE),
                                  Stand="100:Camp",
                                  Year = "2020",
                                  ThinStatus = "Y",
                                  DoyStart = mean(new.net$DoyStart, na.rm=TRUE),
                                  VegAbundance = mean(new.net$VegAbundance, na.rm=TRUE),
                                  VegDiversity = mean(new.net$VegDiversity, na.rm=TRUE),
                                  ForageDist_km = mean(new.net$ForageDist_km, na.rm=TRUE)
)

pred_beediv <- fit.bombus %>%
  epred_draws(newdata = newdata.beediv,
              resp = "BombusRareRichness",
              allow_new_levels = TRUE)

bombus.div.temp <- ggplot(pred_beediv, aes(x = TempCStart,
                                             y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = 'Blues') +
  labs(fill = "Credible interval") +
  xlab("Survey starting temperature (\u00b0C)") +
  ylab("Rarefied *Bombus* <br> species richness") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_point(data=spec.net,
             aes(x=TempCStart, y=BombusRareRichness), cex=2) +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div)

bombus.div.temp

ggsave(bombus.div.temp, file="figures/bombusrich_temp.pdf",
       height=4, width=5)

## ***********************************************************************
## bombus abundance ~ temp
## ***********************************************************************

labs.bee.ab <- (pretty(c(0, new.orig$TempCStart), n=8))
axis.bee.ab <-  standardize.axis(labs.bee.ab,
                                 new.orig$TempCStart)

newdata.beeab <- tidyr::crossing(TempCStart =
                                   seq(min(new.net$TempCStart, na.rm=TRUE),
                                       max(new.net$TempCStart, na.rm=TRUE),
                                       length.out=10),
                                 FlowerRareRichness=mean(new.net$FlowerRareRichness, na.rm=TRUE),
                                 Stand ="100:Camp",
                                 Year = "2020",
                                 ThinStatus = "Y",
                                 DoyStart = mean(new.net$DoyStart, na.rm=TRUE),
                                 VegAbundance = mean(new.net$VegAbundance, na.rm=TRUE),
                                 ForageDist_km = mean(new.net$ForageDist_km, na.rm=TRUE),
                                 VegDiversity = mean(new.net$VegDiversity, na.rm=TRUE)
)

pred_beeab <- fit.bombus %>%
  epred_draws(newdata = newdata.beeab,
              resp = "BombusAbundance",
              allow_new_levels = TRUE)

bombus.ab.temp <- ggplot(pred_beeab, aes(x = TempCStart, 
                                           y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette="Blues") +
  labs(x = "Survey starting temperature (\u00b0C)", 
       y = "*Bombus* abundance",
       fill = "Credible interval") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_point(data=new.net,
             aes(x=TempCStart, y=BombusAbundance), cex=2) +
  scale_x_continuous(
    breaks = axis.bee.ab,
    labels = labs.bee.ab) 

bombus.ab.temp

ggsave(bombus.ab.temp, file="figures/bombusabund_temp.pdf",
       height=4, width=5)

## ***********************************************************************
## bee community diversity and abundance and parasitism
## ***********************************************************************
############################

## parasitism ~ bombus diversity
# 
# labs.bee.div <- (pretty(c(0, orig.spec$BombusRareRichness), n=8))
# axis.bee.div <-  standardize.axis(labs.bee.div,
#                                   spec.net$BombusRareRichness)

data.par <- spec.net[spec.net$WeightsPar == 1, ]

newdata.beediv <- tidyr::crossing(BombusRareRichness =
                                    seq(min(data.par$BombusRareRichness),
                                        max(data.par$BombusRareRichness),
                                        length.out=10),
                                  FlowerRareRichness=mean(data.par$FlowerRareRichness),
                                  Stand="100:Camp",
                                  Year = "2021",
                                  DoyStart = mean(data.par$DoyStart),
                                  VegAbundance = mean(data.par$VegAbundance),
                                  VegDiversity = mean(data.par$VegDiversity),
                                  BombusAbundance = mean(data.par$BombusAbundance),
                                  ForageDist_km = mean(data.par$ForageDist_km)
)

pred_beediv <- fit.bombus %>%
  epred_draws(newdata = newdata.beediv,
              resp = "HasCrithidia",
              allow_new_levels = TRUE)

bee.div.crithidia <- ggplot(pred_beediv, aes(x = BombusRareRichness, y =
                                               .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "Rarefied *Bombus* <br> species richness", 
       y = "Stand *Crithidia* rate",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous() +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  theme(axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown()) +
  geom_point(data=data.par,
             aes(y=BombusStandCrithidiaRate, x=BombusRareRichness), cex=2)

bee.div.crithidia

##why is there a stand with bombus div = 10 

ggsave(bee.div.crithidia, file="figures/parasite_beeDiv.pdf",
       height=4, width=5)

## parasitism ~ bee abundance
# labs.bee.abund <- (pretty(c(0, unweighted$BombusAbundance)))
# axis.bee.abund <-  standardize.axis(labs.bee.abund,
#                                     unweighted$BombusAbundance)

newdata.beeabund <- tidyr::crossing(BombusAbundance =
                                      seq(min(data.par$BombusAbundance),
                                          max(data.par$BombusAbundance),
                                          length.out=10),
                                    FlowerRareRichness=mean(data.par$FlowerRareRichness),
                                    Stand="100:Camp",
                                    Year = "2020",
                                    DoyStart = data.par$DoyStart,
                                    VegAbundance = mean(data.par$VegAbundance),
                                    TempCStart = mean(data.par$TempCStart),
                                    VegDiversity = mean(data.par$VegDiversity),
                                    ForageDist_km = mean(data.par$ForageDist_km)
)

pred_beeabund <- fit.bombus %>%
  epred_draws(newdata = newdata.beeabund ,
              resp = "HasCrithidia",
              allow_new_levels = TRUE)


bombus.abund.parasite <- ggplot(pred_beeabund, aes(x = BombusAbundance, y =
                                                     .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Greys") +
  labs(x = "*Bombus* abundance", y = "Stand *Crithidia* rate",
       fill = "Credible interval") +
  theme(legend.position = "bottom")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  theme(axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown()) +
  geom_point(data=data.par,
             aes(y=BombusStandCrithidiaRate, x=BombusAbundance), cex=2)

bombus.abund.parasite

ggsave(bombus.abund.parasite, file="figures/parasite_beeAbund.pdf",
       height=4, width=5)

#putting plots together
#canopy and bee/flower plots 
scatter.1 <- grid.arrange(veg.div.stand, 
                            flower.ab.stand,
                            bombus.div.canopy,
                            bombus.ab.canopy,
                            ncol=2)

ggsave(scatter.1, file="figures/CRcanopyplots.pdf",
       height=6, width=9)
                            
scatter.2 <- grid.arrange(bombus.div.plant,
                           bombus.ab.plant,
                          ncol=2)

ggsave(scatter.2, file="figures/CRplantplots.pdf",
       height=3, width=9)

scatter.3 <- grid.arrange(bombus.div.temp,
                          bombus.ab.temp,
                          ncol=2)

ggsave(scatter.3, file="figures/CRtempplots.pdf",
       height=3, width=9)

scatter.4 <- grid.arrange(bee.div.crithidia,
                          bombus.abund.parasite,
                          ncol=2)

ggsave(scatter.4, file="figures/CRdiseaseplots.pdf",
       height=3, width=9)

scatter.5 <- grid.arrange(doy.veg.div,
                          doy.veg.ab,
                          ncol=2)

ggsave(scatter.5, file="figures/CRdoyplots.pdf",
       height=3, width=9)

## ggsave(parasite.all, file="figures/presentation_dark/all_parasite.pdf",
##        height=4, width=10)

# #############
# 
# labs.veg.div <- (pretty(spec.net$FlowerRareRichness, n=8))
# axis.veg.div <-  standardize.axis(labs.veg.div,
#                                   spec.net$FlowerRareRichness)
# 
# ## parasitism ~ floral div
# newdata.beeabund <- tidyr::crossing(FlowerRareRichness =
#                                       seq(min(data.par$FlowerRareRichness),
#                                           max(data.par$FlowerRareRichness),
#                                           length.out=10),
#                                     BombusDiversity=mean(data.par$BombusDiversity),
#                                     Stand="100:Camp",
#                                     Year = "2021",
#                                     DoyStart = mean(data.par$DoyStart),
#                                     VegAbundance = mean(data.par$VegAbundance),
#                                     VegDiversity = mean(data.par$VegDiversity),
#                                     BombusAbundance = mean(data.par$BombusAbundance),
#                                     BombusRareRichness = mean(data.par$BombusRareRichness),
#                                     ForageDist_km = mean(data.par$ForageDist_km)
# )
# 
# pred_beeabund <- fit.bombus %>%
#   epred_draws(newdata = newdata.beeabund ,
#               resp = "HasCrithidia",
#               allow_new_levels = TRUE)
# 
# 
# veg.div.parasite <- ggplot(pred_beeabund, aes(x = FlowerRareRichness, y =
#                                                 .epred)) +
#   stat_lineribbon() +
#   scale_fill_brewer(palette = "Reds") +
#   labs(x = "Veg diversity", y = "Bombus Stand Parasitism Rate",
#        fill = "Credible interval") +
#   theme(legend.position = "bottom") +
#   scale_x_continuous(
#     breaks = axis.veg.div,
#     labels =  labs.veg.div) +
#   theme(axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16),
#         text = element_text(size=16)) +
#   theme_ms() +
#   geom_point(data=data.par,
#              aes(y=BombusStandParasitismRate, x=FlowerRareRichness), cex=2)
# 
# veg.div.parasite
# 
# ggsave(veg.div.parasite, file="figures/parasite_vegDiv.pdf",
#        height=4, width=5)

############################

## parasitism ~ doy
# 
# labs.doy <- (pretty(spec.net.orig$DoyStart, n=8))
# axis.doy <-  standardize.axis(labs.doy,
#                               spec.net.orig$DoyStart)
# 
# 
# newdata.beeabund <- tidyr::crossing(DoyStart =
#                                       seq(min(data.par$DoyStart),
#                                           max(data.par$DoyStart),
#                                           length.out=10),
#                                     BombusDiversity=mean(data.par$BombusDiversity),
#                                     Stand="100:Camp",
#                                     FireName = data.par$FireName,
#                                     Year = data.par$Year,
#                                     VegDiversity = mean(data.par$VegDiversity),
#                                     VegAbundance = mean(data.par$VegAbundance),
#                                     BombusAbundance = mean(data.par$BombusAbundance),
#                                     LargestHighSevPatch = mean(data.par$LargestHighSevPatch)
# )
# 
# pred_beeabund <- fit.bombus %>%
#   epred_draws(newdata = newdata.beeabund ,
#               resp = "ParasitePresence",
#               allow_new_levels = TRUE)
# 
# 
# doy.parasite <- ggplot(pred_beeabund, aes(x = DoyStart, y =
#                                             .epred)) +
#   stat_lineribbon() +
#   scale_fill_brewer(palette = "Reds") +
#   labs(x = "DoyStart", y = "Bombus Stand Parasitism Rate",
#        fill = "Credible interval") +
#   theme(legend.position = "bottom") +
#   scale_x_continuous(
#     breaks = axis.doy,
#     labels =  labs.doy) +
#   theme(axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16),
#         text = element_text(size=16)) +
#   theme_ms() +
#   geom_point(data=data.par,
#              aes(y=BombusStandParasitismRate, x=DoyStart), cex=2)
# 
# doy.parasite
# 
# ggsave(doy.parasite, file="figures/parasite_DOY.pdf",
#        height=4, width=5)

# ########################

# parasite.all <- grid.arrange(p1.parasite, p2.parasite, ncol=2)

## ggsave(parasite.all, file="figures/presentation_dark/all_parasite.pdf",
##        height=4, width=10)

