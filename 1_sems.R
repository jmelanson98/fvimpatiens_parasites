setwd("~/Documents/UBC/Bombus Project/fvimpatiens_parasites")

## Prepares the data for model fitting (standardizes continuous
## variables, creates dummy variables to be used as weights to all
## different subsets of data to be used in different model levels),
## builds the models, and fits the models in brms. The model outputs
## are saved as tables, and chain diagnostic plots created.

rm(list=ls())

## set to the number of cores you would like the models to run on
ncores <- 1

fvimp_brmsdf <- read.csv("data/fvimp_brmsdf.csv", row.names = 1)
negbin_brmsdf <- read.csv("data/negativebinomial_brmsdf.csv", row.names = 1)
source("src/init.R")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")

## **********************************************************
## formula for site effects on the bee community
## **********************************************************

## all of the variables that are explanatory variables and thus need
## to be centered. Because brms needs a single dataset, which must be
## at the indivudal-level for the parasite model, we need to carefully
## standardize the data at the correct level. 

## standardize by transect -- landscape variables which are the same during all rounds, vary only 
## by transect location
vars_transect <- c("prop_blueberry",
               "prop_edge",
               "landscape_shdi")

## standardize by transect and sample round -- day of year and vegetation data
vars_transect_sr<- c("julian_date",
               "floral_diversity",
               "floral_abundance" #this is already log-scaled so don't do it again
               )


## variables to log but add 1 first (due to zeros)
variables.to.log.p1 <- c()


## Create a dummy variable "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms
## 
## ********************************************************
## log variables here
## ********************************************************

orig.spec <- negbin_brmsdf

## Make SEM weights and standardize data.
negbin_brmsdf <- prepDataSEM(spec.data = negbin_brmsdf,
                            vars_transect = vars_transect,
                            vars_transect_sr = vars_transect_sr)




## **********************************************************
## Model 1.1: formula for landscape effects on floral community
## **********************************************************
## define all the formulas for the different parts of the models

formula.flower.div <- formula(floral_diversity | weights(Weights) ~
                                  poly(julian_date, 2) +
                                  (1|sample_pt)
                              )

## flower abund with simpson div
formula.flower.abund <- formula(floral_abundance | weights(Weights) ~
                                    poly(julian_date, 2) +
                                    (1|sample_pt)
                                )

## **********************************************************
## Model 1.2: formula for landscape & floral effects on bee community
## **********************************************************

formula.bee.div <- formula(bombus_shannon_diversity | weights(Weights)~
                             floral_diversity + floral_abundance +
                              poly(julian_date, 2) + landscape_shdi +
                               prop_blueberry + prop_edge +
                               (1|sample_pt)
                           )

formula.bee.abund <- formula(native_bee_abundance | weights(Weights)~
                               floral_abundance + floral_diversity +
                               prop_blueberry + prop_edge + landscape_shdi + 
                               poly(julian_date, 2) +
                                 (1|sample_pt)  
                             )

formula.imp.abund <- formula(impatiens_abundance | weights(Weights)~
                               floral_abundance + floral_diversity +
                               prop_blueberry + prop_edge + landscape_shdi + 
                               poly(julian_date, 2) +
                               (1|sample_pt)  
)

## **********************************************************
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************

xvars.fv <- c("bombus_shannon_diversity",
              "native_bee_abundance",
              "impatiens_abundance",
              "floral_abundance",
              "floral_diversity",
              "prop_blueberry",
              "poly(julian_date, 2)",
              "final_id",
              "(1|sample_pt)",
              "(1|subsite)",
              "offset(number_screened)"
                 )
## **********************************************************
## Without parasite model
## **********************************************************
bf.fdiv <- bf(formula.flower.div, family="student")
bf.fabund <- bf(formula.flower.abund, family = "student")
bf.bdiv <- bf(formula.bee.div, family="hurdle_lognormal")
bf.babund <- bf(formula.bee.abund, family = "hurdle_poisson")
bf.iabund <- bf(formula.imp.abund, family = "hurdle_poisson")

## convert to brms format
bform <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund +
  set_rescor(FALSE)

## run model
fit.bombus.nopar <- brm(bform, fvimp_brmsdf,
                  cores=ncores,
                  iter = (10^4),
                  chains =1,
                  thin=1,
                  init=0,
                  open_progress = FALSE,
                  control = list(adapt_delta = 0.999,
                                 stepsize = 0.001,
                                 max_treedepth = 20)
)

write.ms.table(fit.bombus.nopar, "NoParasiteModel_fv")
save(fit.bombus.nopar, fvimp_brmsdf, orig.spec,
     file="saved/NoParasiteModel_fv.Rdata")

load(file="saved/NoParasiteModel_fv.Rdata")

plot.res(fit.bombus.nopar, "NoParasiteModel_fv")

summary(fit.bombus.nopar)

bayes_R2(fit.bombus.nopar)

plot(pp_check(fit.bombus.nopar, resp="floraldiversity"))
plot(pp_check(fit.bombus.nopar, resp="floralabundance"))
plot(pp_check(fit.bombus.nopar, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.nopar, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus.nopar, resp = "impatiensabundance"))


## **********************************************************
## Crithidia 
## **********************************************************

formula.crithidia <-  runParasiteModels(fvimp_brmsdf,
                                       "hascrithidia", 
                                       xvars.fv)

bf.fdiv <- bf(formula.flower.div, family="student")
bf.fabund <- bf(formula.flower.abund, family = "student")
bf.bdiv <- bf(formula.bee.div, family="hurdle_lognormal")
bf.babund <- bf(formula.bee.abund, family = "hurdle_poisson")
bf.iabund <- bf(formula.imp.abund, family = "hurdle_poisson")

## convert to brms format
bf.par <- bf(formula.crithidia, family="bernoulli")
bform <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund + bf.par +
    set_rescor(FALSE)

## run model
fit.bombus <- brm(bform, fvimp_brmsdf,
                  cores=ncores,
                  iter = (10^4),
                  chains =1,
                  thin=1,
                  init=0,
                  open_progress = FALSE,
                  control = list(adapt_delta = 0.999,
                                 stepsize = 0.001,
                                 max_treedepth = 20)
                  )

write.ms.table(fit.bombus, "Crithidia_allbee_fv")
save(fit.bombus, fvimp_brmsdf, orig.spec,
     file="saved/CrithidiaFitAllBee_fv.Rdata")

load(file="saved/CrithidiaFitAllBee_fv.Rdata")

plot.res(fit.bombus, "Crithidia_allbee_fv")

summary(fit.bombus)

bayes_R2(fit.bombus)

plot(pp_check(fit.bombus, resp="floraldiversity"))
plot(pp_check(fit.bombus, resp="floralabundance"))
plot(pp_check(fit.bombus, resp="nativebeeabundance"))
plot(pp_check(fit.bombus, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus, resp = "impatiensabundance"))
plot(pp_check(fit.bombus, resp = "hascrithidia"))


## **********************************************************
## Apicystis
## **********************************************************
#ran this without "subsite" as random effect in parasite model, and without "site" as fixed effect
#in bombus abundance/diversity models. Will rerun!

formula.api <- runParasiteModels(negbin_brmsdf,
                                 "spp_apicystis", 
                                 xvars.fv)

## convert to brms format
bf.par.api <- bf(formula.api, family="negbinomial")
bform.api <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund + bf.par.api +
  set_rescor(FALSE)

## run model
fit.bombus.api <- brm(bform.api, negbin_brmsdf,
                      cores=ncores,
                      iter = (10^4),
                      chains =1,
                      thin=1,
                      init=0,
                      open_progress = FALSE,
                      control = list(adapt_delta = 0.999,
                                     stepsize = 0.001,
                                     max_treedepth = 20)
)


write.ms.table(fit.bombus.api, "Apicystis_allbee_fv_negbinomial")
save(fit.bombus.api, negbin_brmsdf, orig.spec,
     file = "saved/ApicystisAllBee_fv_negbinomial.Rdata")

load(file = "saved/ApicystisAllBee_fv_negbinomial.Rdata")

plot.res(fit.bombus.api, "Apicystis_allbee_fv_negbinomial")

summary(fit.bombus.api)

bayes_R2(fit.bombus.api)

plot(pp_check(fit.bombus.api, resp="floraldiversity"))
plot(pp_check(fit.bombus.api, resp="floralabundance"))
plot(pp_check(fit.bombus.api, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.api, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus.api, resp = "impatiensabundance"))
plot(pp_check(fit.bombus.api, resp = "sppapicystis"))

## **********************************************************
## Nosema bombi
## **********************************************************
#ran this without "subsite" as random effect in parasite model, and without "site" as fixed effect
#in bombus abundance/diversity models. Will rerun!

formula.nbom <- runParasiteModels(fvimp_brmsdf,
                                 "nbombii", 
                                 xvars.fv)

## convert to brms format
bf.par.nbom <- bf(formula.nbom, family="bernoulli")
bform.nbom <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund + bf.par.nbom +
  set_rescor(FALSE)

## run model
fit.bombus.nbom <- brm(bform.nbom, fvimp_brmsdf,
                      cores=ncores,
                      iter = (10^4),
                      chains =2,
                      thin=1,
                      init=0,
                      open_progress = FALSE,
                      control = list(adapt_delta = 0.999,
                                     stepsize = 0.001,
                                     max_treedepth = 20)
)


write.ms.table(fit.bombus.nbom, "Nbombi_allbee_fv")
save(fit.bombus.nbom, fvimp_brmsdf, orig.spec,
     file = "saved/NbombiAllBee_fv.Rdata")

load(file = "saved/NbombiAllBee_fv.Rdata")

plot.res(fit.bombus.nbom, "Nbombi_allbee_fv")

summary(fit.bombus.nbom)

bayes_R2(fit.bombus.nbom)

plot(pp_check(fit.bombus.nbom, resp="floraldiversity"))
plot(pp_check(fit.bombus.nbom, resp="floralabundance"))
plot(pp_check(fit.bombus.nbom, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.nbom, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus.nbom, resp = "impatiensabundance"))
plot(pp_check(fit.bombus.nbom, resp = "nbombii"))

## **********************************************************
## Any parasite
## **********************************************************
#ran this without "subsite" as random effect in parasite model, and without "site" as fixed effect
#in bombus abundance/diversity models: saved as AnyParasiteAllBee_fv.Rdata and Anyparasite_allbee_fv
# also ran WITH both of the above modifications: saved as AnyParasiteAllBeeExtraEffects_fv.Rdata and 
# Anyparasite_allbee_extraeffects_fv

formula.any <- runParasiteModels(fvimp_brmsdf,
                                  "any_parasite", 
                                  xvars.fv)

## convert to brms format
bf.par.any <- bf(formula.any, family="bernoulli")
bform.any <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund + bf.par.any +
  set_rescor(FALSE)

## run model
fit.bombus.any <- brm(bform.any, fvimp_brmsdf,
                       cores=ncores,
                       iter = (10^4),
                       chains =1,
                       thin=1,
                       init=0,
                       open_progress = FALSE,
                       control = list(adapt_delta = 0.999,
                                      stepsize = 0.001,
                                      max_treedepth = 20)
)


write.ms.table(fit.bombus.any, "Anyparasite_allbee_extraeffects_fv")
save(fit.bombus.any, fvimp_brmsdf, orig.spec,
     file = "saved/AnyParasiteAllBeeExtraEffects_fv.Rdata")

load(file = "saved/AnyParasiteAllBeeExtraEffects_fv.Rdata")

plot.res(fit.bombus.any, "AnyParasiteAllBeeExtraEffects_allbee_fv")

summary(fit.bombus.any)

bayes_R2(fit.bombus.any)

plot(pp_check(fit.bombus.any, resp="floraldiversity"))
plot(pp_check(fit.bombus.any, resp="floralabundance"))
plot(pp_check(fit.bombus.any, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.any, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus.any, resp = "impatiensabundance"))
plot(pp_check(fit.bombus.any, resp = "any_parasite"))
