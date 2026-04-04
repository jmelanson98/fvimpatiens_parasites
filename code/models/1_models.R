setwd("~/Documents/UBC/bombus_project/fvimpatiens_parasites")

## Prepares the data for model fitting (standardizes continuous
## variables, creates dummy variables to be used as weights to all
## different subsets of data to be used in different model levels),
## builds the models, and fits the models in brms. The model outputs
## are saved as tables, and chain diagnostic plots created.

rm(list=ls())

fvimp_brmsdf <- read.csv("data/fvimp_brmsdf.csv", sep = ",", header = T, row.names = 1)
source("code/src/init.R")
source("code/src/misc.R")
source("code/src/writeResultsTable.R")
source("code/src/runParasiteModels.R")
source("code/src/standardize_weights.R")
source("code/src/getPhyloMatrix.R")
source("code/src/runPlotFreqModelDiagnostics.R")

## **********************************************************
## Standardize data at correct level(s)
## **********************************************************

## all of the variables that are explanatory variables and thus need
## to be centered. Because brms needs a single dataset, which must be
## at the individual-level for the parasite model, we need to carefully
## standardize the data at the correct level. 

## standardize by transect -- landscape variables which are the same during all rounds, vary only 
## by transect location
vars_transect <- c("prop_blueberry_500",
               "prop_edge_500",
               "landscape_shdi_500",
               "prop_blueberry_250",
               "prop_edge_250",
               "landscape_shdi_250",
               "prop_blueberry_750",
               "prop_edge_750",
               "landscape_shdi_750")

## standardize by transect and sample round -- day of year and vegetation data
vars_transect_sr<- c("julian_date",
               "floral_diversity",
               "floral_abundance" #this is already log-scaled so don't do it again
               )


## variables to log but add 1 first (due to zeros)
variables.to.log.p1 <- c()

## save original data for later
orig.spec <- fvimp_brmsdf

## Make SEM subsets and standardize data.
#insert Bombus medius as filler species for NA rows
fvimp_brmsdf <- prepDataSEM_bernoulli(spec.data = fvimp_brmsdf,
                            vars_transect = vars_transect,
                            vars_transect_sr = vars_transect_sr)


## **********************************************************
## Random effects of phylogeny
## **********************************************************
#create phylo variance-covariance matrix for group of study species
studyspecies = c("Bombus_mixtus", "Bombus_flavifrons", "Bombus_rufocinctus", "Bombus_californicus", 
                 "Bombus_impatiens", "Bombus_vosnesenskii", "Bombus_sitkensis", "Bombus_melanopygus", 
                 "Bombus_nevadensis", "Bombus_medius")
studycov = getPhyloMatrix(studyspecies)


## ************************************************************************************
## Model 1.1: formula for landscape & floral effects on bee community -- 500m buffer
## ************************************************************************************
# use different subsets because there were several days at the beginning of the season when we did not collect impatiens;
# these sampling efforts can be included in native bombus abundance models, 
# but they must be subseted out of all other models; luckily we did not screen any bees for parasites from
# those days, so this does not effect subsetPar

formula.bee.rich <- formula(bombus_richness | trials(9) + subset(impSubset) ~
                              floral_abundance + floral_diversity +
                              prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 + 
                              julian_date + I(julian_date^2) +
                              (1|sample_pt))
formula.bee.abund <- formula(native_bee_abundance | subset(Subset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 +
                               julian_date + I(julian_date^2) +
                                 (1|sample_pt)  )
formula.imp.abund <- formula(impatiens_abundance | subset(impSubset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 +
                               julian_date + I(julian_date^2) +
                               (1|sample_pt))

## ******************************************************************************
## Model 1.2: formula for bee community effects on parasitism -- 500m buffer
## ******************************************************************************

xvars.base.native <- c("floral_abundance",
                    "floral_diversity",
                    "bombus_richness",
                    "native_bee_abundance",
                    "impatiens_abundance",
                    "prop_blueberry_500",
                    "prop_edge_500",
                    "landscape_shdi_500",
                    "caste",
                    "julian_date",
                    "I(julian_date^2)",
                    "(1|subsite)",
                    "(1|sample_pt)",
                   "(1|gr(final_id,cov = studycov))")

xvars.base.impatiens <- c("floral_abundance",
                          "floral_diversity",
                          "bombus_richness",
                          "native_bee_abundance",
                          "impatiens_abundance",
                          "prop_blueberry_500",
                          "prop_edge_500",
                          "landscape_shdi_500",
                          "julian_date",
                          "I(julian_date^2)",
                          "(1|subsite)",
                          "(1|sample_pt)")

## **********************************************************
## Run 500 m buffer models
## **********************************************************

formula.allcrith.native <-  runParasiteModels(fvimp_brmsdf,
                                        "hascrithidia", 
                                        xvars.base.native,
                                        subsetvar = "subsetNativePar")
formula.allnos.native <-  runParasiteModels(fvimp_brmsdf,
                                   "hasnosema", 
                                   xvars.base.native,
                                   subsetvar = "subsetNativePar")
formula.apicystis.native <-  runParasiteModels(fvimp_brmsdf,
                                        "apicystis", 
                                        xvars.base.native,
                                        subsetvar = "subsetNativePar")
formula.allcrith.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                              "hascrithidia", 
                                              xvars.base.impatiens,
                                              subsetvar = "subsetImpPar")
formula.allnos.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                            "hasnosema", 
                                            xvars.base.impatiens,
                                            subsetvar = "subsetImpPar")
formula.apicystis.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                               "apicystis", 
                                               xvars.base.impatiens,
                                               subsetvar = "subsetImpPar")

#convert to brms format
bf.brich <- bf(formula.bee.rich, family="beta_binomial")
bf.babund <- bf(formula.bee.abund, family = "negbinomial")
bf.iabund <- bf(formula.imp.abund, family = "negbinomial")
bf.crith.native <- bf(formula.allcrith.native, family="bernoulli")
bf.nos.native <- bf(formula.allnos.native, family="bernoulli")
bf.api.native <- bf(formula.apicystis.native, family="bernoulli")
bf.crith.impatiens <- bf(formula.allcrith.impatiens, family="bernoulli")
bf.nos.impatiens <- bf(formula.allnos.impatiens, family="bernoulli")
bf.api.impatiens <- bf(formula.apicystis.impatiens, family="bernoulli")

bform.par.native <- bf.crith.native + bf.nos.native + bf.api.native +
  set_rescor(FALSE)
bform.par.impatiens <- bf.crith.impatiens + bf.api.impatiens +
  set_rescor(FALSE)
bform.base <-  bf.brich + bf.babund + bf.iabund +
  set_rescor(FALSE)

#this will likely not run on brms 2.22.0
#not entirely sure why but gradient evaluation on the beta binomial
#becomes EXTREMELY slow
#run on brms version 2.20.4 (and make sure the c++ file is 
#properly recompiled when you switch over...otherwise it will still
#be slow. this makes me think it's a difference in how the model
#gets specified in stan...?)
fit.bombus.base <- brm(bform.base, fvimp_brmsdf,
                  cores=4, chains = 4,
                  iter = 10^4, init=0,
                  control = list(adapt_delta = 0.999,
                                 stepsize = 0.001,
                                 max_treedepth = 20))

write.ms.table(fit.bombus.base, "Base500m_32610")
save(fit.bombus.base, fvimp_brmsdf, orig.spec,
     file="saved/Base500m_32610.Rdata")
load(file="saved/Base500m_32610.Rdata")

# Some quick model checks
summary(fit.bombus.base)
bayes_R2(fit.bombus.base)
plot(pp_check(fit.bombus.base, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.base, resp="bombusrichness"))
plot(pp_check(fit.bombus.base, resp = "impatiensabundance"))


fit.par.native <- brm(bform.par.native, fvimp_brmsdf,
                       cores=4,
                       iter = (10^4),
                       chains = 4,
                       thin=1,
                       init=0,
                       control = list(adapt_delta = 0.999,
                                      stepsize = 0.001,
                                      max_treedepth = 20),
                       data2 = list(studycov = studycov))

write.ms.table(fit.par.native, "NativePar500m_32610")
save(fit.par.native, fvimp_brmsdf, orig.spec,
     file="saved/NativePar500m_32610.Rdata")
load(file="saved/NativePar500m_32610.Rdata")

# Some quick model checks
summary(fit.par.native)
bayes_R2(fit.par.native)
plot(pp_check(fit.par.native, resp = "hascrithidia"))
plot(pp_check(fit.par.native, resp = "hasnosema"))
plot(pp_check(fit.par.native, resp = "apicystis"))


fit.par.impatiens <- brm(bform.par.impatiens, fvimp_brmsdf,
                      cores=4,
                      iter = (10^4),
                      chains = 4,
                      thin=1,
                      init=0,
                      control = list(adapt_delta = 0.999,
                                     stepsize = 0.001,
                                     max_treedepth = 20),
                      data2 = list(studycov = studycov))

write.ms.table(fit.par.impatiens, "ImpatiensPar500m_32610")
save(fit.par.impatiens, fvimp_brmsdf, orig.spec,
     file="saved/ImpatiensPar500m_32610.Rdata")
load(file="saved/ImpatiensPar500m_32610.Rdata")

# Some quick model checks
summary(fit.par.impatiens)
bayes_R2(fit.par.impatiens)
plot(pp_check(fit.par.impatiens, resp = "hascrithidia"))
plot(pp_check(fit.par.impatiens, resp = "apicystis"))





## ************************************************************************************
## Model 2.1: formula for landscape & floral effects on bee community -- 250m buffer
## ************************************************************************************
# use different subsets because there were several days at the beginning of the season when we did not collect impatiens;
# these sampling efforts can be included in native bombus abundance models, 
# but they must be subseted out of all other models; luckily we did not screen any bees for parasites from
# those days, so this does not effect subsetPar

formula.bee.rich <- formula(bombus_richness | trials(9) + subset(impSubset) ~
                              floral_abundance + floral_diversity +
                              prop_blueberry_250 + prop_edge_250 + landscape_shdi_250 + 
                              julian_date + I(julian_date^2) +
                              (1|sample_pt))
formula.bee.abund <- formula(native_bee_abundance | subset(Subset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry_250 + prop_edge_250 + landscape_shdi_250 + 
                               julian_date + I(julian_date^2) +
                               (1|sample_pt)  )
formula.imp.abund <- formula(impatiens_abundance | subset(impSubset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry_250 + prop_edge_250 + landscape_shdi_250 + 
                               julian_date + I(julian_date^2) +
                               (1|sample_pt))

## ******************************************************************************
## Model 2.2: formula for bee community effects on parasitism -- 250m buffer
## ******************************************************************************

xvars.base.native <- c("floral_abundance",
                       "floral_diversity",
                       "bombus_richness",
                       "native_bee_abundance",
                       "impatiens_abundance",
                       "prop_blueberry_250",
                       "prop_edge_250",
                       "landscape_shdi_250",
                       "caste",
                       "julian_date",
                       "I(julian_date^2)",
                       "(1|subsite)",
                       "(1|sample_pt)",
                       "(1|gr(final_id,cov = studycov))")

xvars.base.impatiens <- c("floral_abundance",
                          "floral_diversity",
                          "bombus_richness",
                          "native_bee_abundance",
                          "impatiens_abundance",
                          "prop_blueberry_250",
                          "prop_edge_250",
                          "landscape_shdi_250",
                          "julian_date",
                          "I(julian_date^2)",
                          "(1|subsite)",
                          "(1|sample_pt)")

## **********************************************************
## Run 250 m buffer models
## **********************************************************

formula.allcrith.native <-  runParasiteModels(fvimp_brmsdf,
                                              "hascrithidia", 
                                              xvars.base.native,
                                              subsetvar = "subsetNativePar")
formula.allnos.native <-  runParasiteModels(fvimp_brmsdf,
                                            "hasnosema", 
                                            xvars.base.native,
                                            subsetvar = "subsetNativePar")
formula.apicystis.native <-  runParasiteModels(fvimp_brmsdf,
                                               "apicystis", 
                                               xvars.base.native,
                                               subsetvar = "subsetNativePar")
formula.allcrith.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                                 "hascrithidia", 
                                                 xvars.base.impatiens,
                                                 subsetvar = "subsetImpPar")
formula.allnos.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                               "hasnosema", 
                                               xvars.base.impatiens,
                                               subsetvar = "subsetImpPar")
formula.apicystis.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                                  "apicystis", 
                                                  xvars.base.impatiens,
                                                  subsetvar = "subsetImpPar")

#convert to brms format
bf.brich <- bf(formula.bee.rich, family="beta_binomial")
bf.babund <- bf(formula.bee.abund, family = "negbinomial")
bf.iabund <- bf(formula.imp.abund, family = "negbinomial")
bf.crith.native <- bf(formula.allcrith.native, family="bernoulli")
bf.nos.native <- bf(formula.allnos.native, family="bernoulli")
bf.api.native <- bf(formula.apicystis.native, family="bernoulli")
bf.crith.impatiens <- bf(formula.allcrith.impatiens, family="bernoulli")
bf.nos.impatiens <- bf(formula.allnos.impatiens, family="bernoulli")
bf.api.impatiens <- bf(formula.apicystis.impatiens, family="bernoulli")

bform.par.native <- bf.crith.native + bf.nos.native + bf.api.native +
  set_rescor(FALSE)
bform.par.impatiens <- bf.crith.impatiens + bf.api.impatiens +
  set_rescor(FALSE)
bform.base <-  bf.brich + bf.babund + bf.iabund +
  set_rescor(FALSE)

#this will likely not run on brms 2.22.0
#not entirely sure why but gradient evaluation on the beta binomial
#becomes EXTREMELY slow
#run on brms version 2.20.4 (and make sure the c++ file is 
#properly recompiled when you switch over...otherwise it will still
#be slow. this makes me think it's a difference in how the model
#gets specified in stan...?)
fit.bombus.base <- brm(bform.base, fvimp_brmsdf,
                       cores=4, chains = 4,
                       iter = 10^4, init=0,
                       control = list(adapt_delta = 0.999,
                                      stepsize = 0.001,
                                      max_treedepth = 20))

write.ms.table(fit.bombus.base, "Base250m_32610")
save(fit.bombus.base, fvimp_brmsdf, orig.spec,
     file="saved/Base250m_32610.Rdata")
load(file="saved/Base250m_32610.Rdata")

# Some quick model checks
summary(fit.bombus.base)
bayes_R2(fit.bombus.base)
plot(pp_check(fit.bombus.base, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.base, resp="bombusrichness"))
plot(pp_check(fit.bombus.base, resp = "impatiensabundance"))


fit.par.native <- brm(bform.par.native, fvimp_brmsdf,
                      cores=4,
                      iter = (10^4),
                      chains = 4,
                      thin=1,
                      init=0,
                      control = list(adapt_delta = 0.999,
                                     stepsize = 0.001,
                                     max_treedepth = 20),
                      data2 = list(studycov = studycov))

write.ms.table(fit.par.native, "NativePar250m_32610")
save(fit.par.native, fvimp_brmsdf, orig.spec,
     file="saved/NativePar250m_32610.Rdata")
load(file="saved/NativePar250m_32610.Rdata")

# Some quick model checks
summary(fit.par.native)
bayes_R2(fit.par.native)
plot(pp_check(fit.par.native, resp = "hascrithidia"))
plot(pp_check(fit.par.native, resp = "hasnosema"))
plot(pp_check(fit.par.native, resp = "apicystis"))


fit.par.impatiens <- brm(bform.par.impatiens, fvimp_brmsdf,
                         cores=4,
                         iter = (10^4),
                         chains = 4,
                         thin=1,
                         init=0,
                         control = list(adapt_delta = 0.999,
                                        stepsize = 0.001,
                                        max_treedepth = 20),
                         data2 = list(studycov = studycov))

write.ms.table(fit.par.impatiens, "ImpatiensPar250m_32610")
save(fit.par.impatiens, fvimp_brmsdf, orig.spec,
     file="saved/ImpatiensPar250m_32610.Rdata")
load(file="saved/ImpatiensPar250m_32610.Rdata")

# Some quick model checks
summary(fit.par.impatiens)
bayes_R2(fit.par.impatiens)
plot(pp_check(fit.par.impatiens, resp = "hascrithidia"))
plot(pp_check(fit.par.impatiens, resp = "apicystis"))



## ************************************************************************************
## Model 2.1: formula for landscape & floral effects on bee community -- 750m buffer
## ************************************************************************************
# use different subsets because there were several days at the beginning of the season when we did not collect impatiens;
# these sampling efforts can be included in native bombus abundance models, 
# but they must be subseted out of all other models; luckily we did not screen any bees for parasites from
# those days, so this does not effect subsetPar

formula.bee.rich <- formula(bombus_richness | trials(9) + subset(impSubset) ~
                              floral_abundance + floral_diversity +
                              prop_blueberry_750 + prop_edge_750 + landscape_shdi_750 + 
                              julian_date + I(julian_date^2) +
                              (1|sample_pt))
formula.bee.abund <- formula(native_bee_abundance | subset(Subset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry_750 + prop_edge_750 + landscape_shdi_750 + 
                               julian_date + I(julian_date^2) +
                               (1|sample_pt)  )
formula.imp.abund <- formula(impatiens_abundance | subset(impSubset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry_750 + prop_edge_750 + landscape_shdi_750 + 
                               julian_date + I(julian_date^2) +
                               (1|sample_pt))

## ******************************************************************************
## Model 2.2: formula for bee community effects on parasitism -- 750m buffer
## ******************************************************************************

xvars.base.native <- c("floral_abundance",
                       "floral_diversity",
                       "bombus_richness",
                       "native_bee_abundance",
                       "impatiens_abundance",
                       "prop_blueberry_750",
                       "prop_edge_750",
                       "landscape_shdi_750",
                       "caste",
                       "julian_date",
                       "I(julian_date^2)",
                       "(1|subsite)",
                       "(1|sample_pt)",
                       "(1|gr(final_id,cov = studycov))")

xvars.base.impatiens <- c("floral_abundance",
                          "floral_diversity",
                          "bombus_richness",
                          "native_bee_abundance",
                          "impatiens_abundance",
                          "prop_blueberry_750",
                          "prop_edge_750",
                          "landscape_shdi_750",
                          "julian_date",
                          "I(julian_date^2)",
                          "(1|subsite)",
                          "(1|sample_pt)")

## **********************************************************
## Run 750 m buffer models
## **********************************************************

formula.allcrith.native <-  runParasiteModels(fvimp_brmsdf,
                                              "hascrithidia", 
                                              xvars.base.native,
                                              subsetvar = "subsetNativePar")
formula.allnos.native <-  runParasiteModels(fvimp_brmsdf,
                                            "hasnosema", 
                                            xvars.base.native,
                                            subsetvar = "subsetNativePar")
formula.apicystis.native <-  runParasiteModels(fvimp_brmsdf,
                                               "apicystis", 
                                               xvars.base.native,
                                               subsetvar = "subsetNativePar")
formula.allcrith.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                                 "hascrithidia", 
                                                 xvars.base.impatiens,
                                                 subsetvar = "subsetImpPar")
formula.allnos.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                               "hasnosema", 
                                               xvars.base.impatiens,
                                               subsetvar = "subsetImpPar")
formula.apicystis.impatiens <-  runParasiteModels(fvimp_brmsdf,
                                                  "apicystis", 
                                                  xvars.base.impatiens,
                                                  subsetvar = "subsetImpPar")

#convert to brms format
bf.brich <- bf(formula.bee.rich, family="beta_binomial")
bf.babund <- bf(formula.bee.abund, family = "negbinomial")
bf.iabund <- bf(formula.imp.abund, family = "negbinomial")
bf.crith.native <- bf(formula.allcrith.native, family="bernoulli")
bf.nos.native <- bf(formula.allnos.native, family="bernoulli")
bf.api.native <- bf(formula.apicystis.native, family="bernoulli")
bf.crith.impatiens <- bf(formula.allcrith.impatiens, family="bernoulli")
bf.nos.impatiens <- bf(formula.allnos.impatiens, family="bernoulli")
bf.api.impatiens <- bf(formula.apicystis.impatiens, family="bernoulli")

bform.par.native <- bf.crith.native + bf.nos.native + bf.api.native +
  set_rescor(FALSE)
bform.par.impatiens <- bf.crith.impatiens +bf.api.impatiens +
  set_rescor(FALSE)
bform.base <-  bf.brich + bf.babund + bf.iabund +
  set_rescor(FALSE)

#this will likely not run on brms 2.22.0
#not entirely sure why but gradient evaluation on the beta binomial
#becomes EXTREMELY slow
#run on brms version 2.20.4 (and make sure the c++ file is 
#properly recompiled when you switch over...otherwise it will still
#be slow. this makes me think it's a difference in how the model
#gets specified in stan...?)
fit.bombus.base <- brm(bform.base, fvimp_brmsdf,
                       cores=4, chains = 4,
                       iter = 10^4, init=0,
                       control = list(adapt_delta = 0.999,
                                      stepsize = 0.001,
                                      max_treedepth = 20))

write.ms.table(fit.bombus.base, "Base750m_32610")
save(fit.bombus.base, fvimp_brmsdf, orig.spec,
     file="saved/Base750m_32610.Rdata")
load(file="saved/Base750m_32610.Rdata")

# Some quick model checks
summary(fit.bombus.base)
bayes_R2(fit.bombus.base)
plot(pp_check(fit.bombus.base, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.base, resp="bombusrichness"))
plot(pp_check(fit.bombus.base, resp = "impatiensabundance"))


fit.par.native <- brm(bform.par.native, fvimp_brmsdf,
                      cores=4,
                      iter = (10^4),
                      chains = 4,
                      thin=1,
                      init=0,
                      control = list(adapt_delta = 0.999,
                                     stepsize = 0.001,
                                     max_treedepth = 20),
                      data2 = list(studycov = studycov))

write.ms.table(fit.par.native, "NativePar750m_32610")
save(fit.par.native, fvimp_brmsdf, orig.spec,
     file="saved/NativePar750m_32610.Rdata")
load(file="saved/NativePar750m_32610.Rdata")

# Some quick model checks
summary(fit.par.native)
bayes_R2(fit.par.native)
plot(pp_check(fit.par.native, resp = "hascrithidia"))
plot(pp_check(fit.par.native, resp = "hasnosema"))
plot(pp_check(fit.par.native, resp = "apicystis"))


fit.par.impatiens <- brm(bform.par.impatiens, fvimp_brmsdf,
                         cores=4,
                         iter = (10^4),
                         chains = 4,
                         thin=1,
                         init=0,
                         control = list(adapt_delta = 0.999,
                                        stepsize = 0.001,
                                        max_treedepth = 20),
                         data2 = list(studycov = studycov))

write.ms.table(fit.par.impatiens, "ImpatiensPar750m_32610")
save(fit.par.impatiens, fvimp_brmsdf, orig.spec,
     file="saved/ImpatiensPar750m_32610.Rdata")
load(file="saved/ImpatiensPar750m_32610.Rdata")

# Some quick model checks
summary(fit.par.impatiens)
bayes_R2(fit.par.impatiens)
plot(pp_check(fit.par.impatiens, resp = "hascrithidia"))
plot(pp_check(fit.par.impatiens, resp = "apicystis"))



## **********************************************************
## Model Checks with DHARMa
## **********************************************************
check_brms <- function(model,             # brms model
                       integer = TRUE,   # integer response? (TRUE/FALSE)
                       plot = TRUE,       # make plot?
                       ...                # further arguments for DHARMa::plotResiduals
) {
  mdata <- brms::standata(model)
  if (!"Y" %in% names(mdata))
    stop("Cannot extract the required information from this brms model")
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model, ndraws = 1000)),
    observedResponse = mdata$Y,
    fittedPredictedResponse = apply(
      t(brms::posterior_epred(model, ndraws = 1000, re.form = NA)),
      1,
      mean),
    integerResponse = integer)
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  invisible(dharma.obj)
}

checked.crith.base <- check_brms(fit.bombus.crith)
checked.nos.base <- check_brms(fit.bombus.nos)
checked.api.base <- check_brms(fit.bombus.api.inter)
checked.brich <- check_brms(fit.brich)
checked.babund <- check_brms(fit.babund)
checked.iabund <- check_brms(fit.iabund)

testZeroInflation(checked.babund)
testQuantiles(checked.babund)



## **********************************************************
## Model Comparisons
## **********************************************************
#comparison of simple vs two-way interaction models
loo_crith = loo::loo(
  fit.bombus.par,
  fit.bombus.par.inter,
  compare = TRUE,
  resp = "hascrithidia",
  pointwise = FALSE,
  k_threshold = 0.7,
  model_names = NULL
)

loo_apicystis = loo::loo(
  fit.bombus.par,
  fit.bombus.par.inter,
  compare = TRUE,
  resp = "apicystis",
  pointwise = FALSE,
  k_threshold = 0.7,
  model_names = NULL
)

crith_diff = as.data.frame(loo_crith$diffs)
crith_diff$parasite = "crithidia"
api_diff = as.data.frame(loo_apicystis$diffs)
api_diff$parasite = "apicystis"
loo_df = rbind(crith_diff, api_diff)

write.csv(loo_df, "saved/tables/two_way_loo.csv")

## **********************************************************
## Models for pairwise comparison of parasitism between species
## **********************************************************
bf.cbombi <- bf(formula(cbombii ~ final_id), family="bernoulli")
bf.cexpoeki <- bf(formula(cexpoeki ~ final_id), family="bernoulli")
bf.crithidiaspp <- bf(formula(crithidiaspp ~ final_id), family = "bernoulli")
bf.apicystis <- bf(formula(apicystis ~ final_id), family="bernoulli")
bf.nceranae <- bf(formula(nceranae ~ final_id), family = "bernoulli")
bf.nbombi <- bf(formula(nbombii ~ final_id), family = "bernoulli")
bf.asco <- bf(formula(ascosphaera ~ final_id), family = "bernoulli")


prior <- c(
  prior(normal(0, 5), class = "b"), # Prior for the coefficients (fixed effects)
  prior(normal(0, 5), class = "Intercept") # Prior for the intercept
)

fvimp_subpar$final_id = as.factor(fvimp_subpar$final_id)

fit.asco.species <- brm(bf.asco, fvimp_subpar,
                                cores=3,
                                iter = (10^4),
                                chains = 3,
                                prior = prior,
                                thin=1,
                                init=0,
                                save_pars = save_pars(all = TRUE),
                                open_progress = FALSE,
                                control = list(adapt_delta = 0.999,
                                               stepsize = 0.001,
                                               max_treedepth = 20)
)


write.ms.table.differences(fit.cbombi.species, "cbombispeciestable")
write.ms.table.differences(fit.cexpoeki.species, "cexpoekispeciestable")
write.ms.table.differences(fit.crithidiaspp.species, "crithidiasppspeciestable")
write.ms.table.differences(fit.nbombi.species, "nbombispeciestable")
write.ms.table.differences(fit.nceranae.species, "nceranaespeciestable")
write.ms.table.differences(fit.apicystis.species, "apicystisspeciestable")
write.ms.table.differences(fit.asco.species, "ascospeciestable")



