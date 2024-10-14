setwd("~/Documents/UBC/Bombus Project/fvimpatiens_parasites")

## Prepares the data for model fitting (standardizes continuous
## variables, creates dummy variables to be used as weights to all
## different subsets of data to be used in different model levels),
## builds the models, and fits the models in brms. The model outputs
## are saved as tables, and chain diagnostic plots created.

rm(list=ls())

## set to the number of cores you would like the models to run on
ncores <- 1

negbin_brmsdf <- read.csv("data/negativebinomial_brmsdf.csv", row.names = 1)
source("src/init.R")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")
source("src/runPlotFreqModelDiagnostics.R")

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
negbin_brmsdf <- prepDataSEM_negbin(spec.data = negbin_brmsdf,
                             vars_transect = vars_transect,
                             vars_transect_sr = vars_transect_sr)




## **********************************************************
## Model 1.1: formula for landscape effects on floral community
## **********************************************************
## define all the formulas for the different parts of the models

formula.flower.div <- formula(floral_diversity | subset(Subset) ~
                                julian_date + I(julian_date^2) +
                                (1|sample_pt)
)

## flower abund with simpson div
formula.flower.abund <- formula(floral_abundance | subset(Subset) ~
                                  julian_date + I(julian_date^2) +
                                  (1|sample_pt)
)

## **********************************************************
## Model 1.2: formula for landscape & floral effects on bee community
## **********************************************************

formula.bee.div <- formula(bombus_shannon_diversity | subset(Subset)~
                             floral_diversity + floral_abundance +
                             julian_date + I(julian_date^2) + landscape_shdi +
                             prop_blueberry + prop_edge +
                             (1|sample_pt)
)

formula.bee.abund <- formula(native_bee_abundance | subset(Subset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry + prop_edge + landscape_shdi + 
                               julian_date, I(julian_date^2) +
                               (1|sample_pt)  
)

formula.imp.abund <- formula(impatiens_abundance | subset(Subset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry + prop_edge + landscape_shdi + 
                               julian_date + I(julian_date^2) +
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
              "julian_date",
              "I(julian_date^2)",
              "final_id",
              "(1|sample_pt)",
              "(1|subsite)",
              "offset(number_screened)"
)

## **********************************************************
## Crithidia 
## **********************************************************

formula.crithidia <-  runParasiteModels(negbin_brmsdf,
                                        "spp_anycrithidia", 
                                        xvars.fv)

bf.fdiv <- bf(formula.flower.div, family="student")
bf.fabund <- bf(formula.flower.abund, family = "student")
bf.bdiv <- bf(formula.bee.div, family="hurdle_lognormal")
bf.babund <- bf(formula.bee.abund, family = "negbinomial")
bf.iabund <- bf(formula.imp.abund, family = "negbinomial")

## convert to brms format
bf.par <- bf(formula.crithidia, family="negbinomial")
bform <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund + bf.par +
  set_rescor(FALSE)


##set priors
prior<-c(set_prior("normal(log(0.5), 1)",class = "Intercept", coef = "", resp = "sppanycrithidia"),
         set_prior("normal(0,1)", class = "b", coef = "", resp = "sppanycrithidia"))
## run model
fit.bombus.crith.priors <- brm(bform, negbin_brmsdf,
                  prior = prior,
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

write.ms.table(fit.bombus.crith.priors, "Crithidia_fv_negbinomial_priors")
save(fit.bombus.crith.priors, negbin_brmsdf, orig.spec,
     file="saved/CrithidiaFit_fv_negbinomial_priors.Rdata")

load(file="saved/CrithidiaFit_fv_negbinomial_priors.Rdata")

plot.res(fit.bombus.crith.priors, "Crithidia_fv_negbinomial_priors")

summary(fit.bombus.crith.priors)

bayes_R2(fit.bombus.crith)

plot(pp_check(fit.bombus.crith, resp="floraldiversity"))
plot(pp_check(fit.bombus.crith, resp="floralabundance"))
plot(pp_check(fit.bombus.crith, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.crith, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus.crith, resp = "impatiensabundance"))
plot(pp_check(fit.bombus.crith, resp = "sppanycrithidia"))


## **********************************************************
## Apicystis
## **********************************************************

formula.api <- runParasiteModels(negbin_brmsdf,
                                 "spp_apicystis", 
                                 xvars.fv)

## convert to brms format
bf.par.api <- bf(formula.api, family="hurdle_negbinomial")
bform.api <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund + bf.par.api +
  set_rescor(FALSE)

##set priors
prior<-c(set_prior("normal(log(0.5), 1)",class = "Intercept", coef = "", resp = "sppapicystis"),
         set_prior("normal(0,1)", class = "b", coef = "", resp = "sppapicystis"))

## run model
fit.bombus.api <- brm(bform.api, negbin_brmsdf,
                              cores=ncores,
                              iter = (10^4),
                              chains =1,
                              prior = prior,
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

formula.nbom <- runParasiteModels(negbin_brmsdf,
                                  "spp_nbombi", 
                                  xvars.fv)

## convert to brms format
bf.par.nbom <- bf(formula.nbom, family="negbinomial")
bform.nbom <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund + bf.par.nbom +
  set_rescor(FALSE)

#set priors
prior<-c(set_prior("normal(log(0.5), 1)",class = "Intercept", coef = "", resp = "sppnbombi"),
         set_prior("normal(0,1)", class = "b", coef = "", resp = "sppnbombi"))


## run model
fit.bombus.nbom <- brm(bform.nbom, negbin_brmsdf,
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


write.ms.table(fit.bombus.nbom, "Nbombi_allbee_fv_negbinomial")
save(fit.bombus.nbom, negbin_brmsdf, orig.spec,
     file = "saved/NbombiAllBee_fv_negbinomial.Rdata")

load(file = "saved/NbombiAllBee_fv_negbinomial.Rdata")

plot.res(fit.bombus.nbom, "Nbombi_allbee_fv_negbinomial")

summary(fit.bombus.nbom)

bayes_R2(fit.bombus.nbom)

plot(pp_check(fit.bombus.nbom, resp="floraldiversity"))
plot(pp_check(fit.bombus.nbom, resp="floralabundance"))
plot(pp_check(fit.bombus.nbom, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.nbom, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus.nbom, resp = "impatiensabundance"))
plot(pp_check(fit.bombus.nbom, resp = "sppnbombi"))

## **********************************************************
## Any parasite
## **********************************************************

formula.any <- runParasiteModels(negbin_brmsdf,
                                 "spp_anyparasite", 
                                 xvars.fv)

## convert to brms format
bf.par.any <- bf(formula.any, family="negbinomial")
bform.any <-  bf.fdiv + bf.fabund + bf.babund + bf.bdiv + bf.iabund + bf.par.any +
  set_rescor(FALSE)

#set priors
prior<-c(set_prior("normal(log(0.5), 1)",class = "Intercept", coef = "", resp = "sppnanyparasite"),
         set_prior("normal(0,1)", class = "b", coef = "", resp = "sppanyparasite"))

## run model
fit.bombus.any <- brm(bform.any, negbin_brmsdf,
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


write.ms.table(fit.bombus.any, "Anyparasite_allbee_fv_negbinomial")
save(fit.bombus.any, negbin_brmsdf, orig.spec,
     file = "saved/AnyParasiteAllBee_fv_negbinomial.Rdata")

load(file = "saved/AnyParasiteAllBee_fv_negbinomial.Rdata")

plot.res(fit.bombus.any, "AnyParasiteAllBee_allbee_fv_negbinomial")

summary(fit.bombus.any)

bayes_R2(fit.bombus.any)

plot(pp_check(fit.bombus.any, resp="floraldiversity"))
plot(pp_check(fit.bombus.any, resp="floralabundance"))
plot(pp_check(fit.bombus.any, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.any, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus.any, resp = "impatiensabundance"))
plot(pp_check(fit.bombus.any, resp = "sppanyparasite"))




## **********************************************************
## Frequentist Model Checks
## **********************************************************

#formulas to remove subset from parasite & bee models
remove_subset_parasite_formula <- function(form){
  char.form <- as.character(form)
  no.sub <-
    gsub("\\| subset\\(subsetPar}*\\)",
         "", char.form[2])
  form.out <- formula(paste(no.sub, "~", char.form[3]))
  return(form.out)
}

remove_subset_other_formula <- function(form){
  char.form <- as.character(form)
  no.sub <-
    gsub("\\| subset\\(Subset}*\\)",
         "", char.form[2])
  form.out <- formula(paste(no.sub, "~", char.form[3]))
  return(form.out)
}


#run frequentist checks
run_plot_freq_model_diagnostics(remove_subset_parasite_formula(formula.crithidia),
                                this_data=negbin_brmsdf[negbin_brmsdf$subsetPar == TRUE,],
                                this_family="negbinomial", site.lat="site")
