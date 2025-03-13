setwd("~/Documents/UBC/Bombus Project/fvimpatiens_parasites")

## Prepares the data for model fitting (standardizes continuous
## variables, creates dummy variables to be used as weights to all
## different subsets of data to be used in different model levels),
## builds the models, and fits the models in brms. The model outputs
## are saved as tables, and chain diagnostic plots created.

rm(list=ls())

## set to the number of cores you would like the models to run on
ncores <- 3

fvimp_brmsdf <- read.csv("data/fvimp_brmsdf.csv", sep = ",", header = T, row.names = 1)
source("src/init.R")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")
source("src/getPhyloMatrix.R")
source("src/runPlotFreqModelDiagnostics.R")

## **********************************************************
## formula for site effects on the bee community
## **********************************************************

## all of the variables that are explanatory variables and thus need
## to be centered. Because brms needs a single dataset, which must be
## at the individual-level for the parasite model, we need to carefully
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

## save original data for later
orig.spec <- fvimp_brmsdf

## Make SEM subsets and standardize data.
#insert Bombus medius as filler species for NA rows
fvimp_brmsdf <- prepDataSEM_bernoulli(spec.data = fvimp_brmsdf,
                            vars_transect = vars_transect,
                            vars_transect_sr = vars_transect_sr)



## **********************************************************
## Model 1.1: formula for landscape & floral effects on bee community
## **********************************************************
#use different subsets because there were several days at the beginning of the season when we did not collect impatiens;
#these sampling efforts can be included in native bombus abundance models, because it does
#not affect either. But they must be subseted out of all other models; luckily we did not screen any bees for parasites from
#those days, so this does not effect subsetPar


formula.bee.rich <- formula(bombus_richness | trials(9) + subset(impSubset) ~
                             prop_blueberry + prop_edge + landscape_shdi +
                              floral_abundance + floral_diversity +
                             julian_date + I(julian_date^2) +
                             (1|sample_pt)
)


formula.bee.abund <- formula(native_bee_abundance | subset(Subset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry + prop_edge + landscape_shdi + 
                               julian_date + I(julian_date^2) +
                                 (1|sample_pt)  
                             )

formula.imp.abund <- formula(impatiens_abundance | subset(impSubset)~
                               floral_abundance + floral_diversity +
                               prop_blueberry + prop_edge + landscape_shdi + 
                               julian_date + I(julian_date^2) +
                               (1|sample_pt)  
)

## **********************************************************
## Model 1.2: formula for bee community effects on parasitism
## **********************************************************

xvars.fv.base <- c("floral_abundance",
                    "floral_diversity",
                    "bombus_richness",
                    "native_bee_abundance",
                    "impatiens_abundance",
                    "prop_blueberry",
                    "prop_edge",
                    "landscape_shdi",
                    "caste",
                    "julian_date",
                    "I(julian_date^2)",
                    "(1|sample_pt)",
                   "(1|subsite)",
                   "(1|gr(final_id,cov = studycov))"
                 )

## **********************************************************
## Random effects of phylogeny
## **********************************************************
#create phylo variance-covariance matrix for group of study species
studyspecies = c("Bombus_mixtus", "Bombus_flavifrons", "Bombus_rufocinctus", "Bombus_californicus", 
                 "Bombus_impatiens", "Bombus_vosnesenskii", "Bombus_sitkensis", "Bombus_melanopygus", 
                 "Bombus_nevadensis", "Bombus_medius")
studycov = getPhyloMatrix(studyspecies)


## **********************************************************
## Run all models
## **********************************************************

#no interactions
formula.allcrith.base <-  runParasiteModels(fvimp_brmsdf,
                                        "hascrithidia", 
                                        xvars.fv.base)
formula.allnos.base <-  runParasiteModels(fvimp_brmsdf,
                                   "hasnosema", 
                                   xvars.fv.base)
formula.apicystis.base <-  runParasiteModels(fvimp_brmsdf,
                                        "apicystis", 
                                        xvars.fv.base)

bf.brich <- bf(formula.bee.rich, family="beta_binomial")
bf.babund <- bf(formula.bee.abund, family = "negbinomial")
bf.iabund <- bf(formula.imp.abund, family = "negbinomial")

## convert to brms format
bf.allcrith.base <- bf(formula.allcrith.base, family="bernoulli")
bf.allnos.base <- bf(formula.allnos.base, family="bernoulli")
bf.api.base <- bf(formula.apicystis.base, family="bernoulli")
bform.base <-  bf.brich + bf.babund + bf.iabund + 
  bf.allcrith.base + bf.allnos.base + bf.api.base +
  set_rescor(FALSE)


##set priors for pesky interaction term
prior <- c(set_prior("normal(0, 1)", class = "b", 
                      coef = "native_bee_abundance:statusnonnative", 
                      resp = "hasnosema"))


## run model without interactions
fit.bombus <- brm(bform.base, fvimp_brmsdf,
                              cores=4,
                              iter = (10^4),
                              chains = 4,
                              thin=1,
                              init=0,
                              save_pars = save_pars(all = TRUE),
                              open_progress = FALSE,
                              control = list(adapt_delta = 0.999,
                                             stepsize = 0.001,
                                             max_treedepth = 20),
                              data2 = list(studycov = studycov)
)



write.ms.table(fit.bombus.all.beta.base, "AllModels_fv_beerichbeta")
save(fit.bombus.all.beta.base, fvimp_brmsdf, orig.spec,
     file="saved/AllModels_fv_beerichbeta.Rdata")

load(file="saved/AllModels_fv_beerichbeta.Rdata")

plot.res(fit.bombus.all.beta.base, "AllModels_fv_beerichbeta")

summary(fit.all.base.brichnb)

bayes_R2(fit.bombus.all.beta)

plot(pp_check(fit.bombus.all.base, resp="floraldiversity"))
plot(pp_check(fit.bombus.all.base, resp="floralabundance"))
plot(pp_check(fit.bombus.all.base, resp="nativebeeabundance"))
plot(pp_check(fit.all.base.brichnb, resp="bombusrichness"))
plot(pp_check(fit.bombus.all.base, resp = "impatiensabundance"))
plot(pp_check(fit.bombus.all.base, resp = "hascrithidia"))
plot(pp_check(fit.bombus.all.base, resp = "hasnosema"))
plot(pp_check(fit.bombus.all.base, resp = "apicystis"))



## **********************************************************
## Frequentist Model Checks
## **********************************************************

#formulas to remove subset from parasite & bee models
remove_subset_parasite_formula <- function(form){
  char.form <- as.character(form)
  no.sub <-
    gsub("\\| subset\\(subsetPar}*\\)",
         "", char.form[2])
  no.group <-
    gsub("\\(1 \\| gr\\(final_id, cov = studycov\\)\\)", "\\(1 \\| final_id\\)", char.form[3])
  form.out <- formula(paste(no.sub, "~", no.group))
  return(form.out)
}

remove_subset_imp_formula <- function(form){
  char.form <- as.character(form)
  no.sub <-
    gsub("\\| subset\\(impSubset}*\\)",
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


#run frequentist checks for parasite models
run_plot_freq_model_diagnostics(remove_subset_parasite_formula(formula.allcrith.base),
                                this_data=fvimp_brmsdf[fvimp_brmsdf$subsetPar == TRUE,],
                                this_family="bernoulli", site.lat="site")
#impatiens_abundance VIF becomes intermediate for apicystis interaction model (which we already know is not great)
#still big error bars on nosema interaction VIF still (and binned residuals are worse?)

#run frequentist checks for models including impatiens
run_plot_freq_model_diagnostics(remove_subset_imp_formula(formula.bee.div),
                                  this_data=fvimp_brmsdf[fvimp_brmsdf$impSubset == TRUE,],
                                  this_family="negbinomial", site.lat="site")


#run frequentist checks for floral/native bee abundance models
run_plot_freq_model_diagnostics(remove_subset_other_formula(formula.bee.abund),
                                this_data=fvimp_brmsdf[fvimp_brmsdf$Subset == TRUE,],
                                this_family="negbinomial", site.lat="site")



## **********************************************************
## Model Comparisons
## **********************************************************
#comparison of two best models
waic(
  fit.bombus.all.prior.natinteraction,
  fit.bombus.all.prior.impinteraction,
  fit.bombus.all.prior.nointeraction,
  compare = TRUE,
  resp = "hasnosema",
  pointwise = FALSE,
  model_names = NULL
)


loofit = loo::loo(
  fit.bombus.all.inter,
  compare = TRUE,
  resp = "hasnosema",
  pointwise = FALSE,
  k_threshold = 0.7,
  model_names = NULL
)

loofit = loo::loo(
  fit.nosema.nointeractions.momentmatching,
  fit.nosema.natinteraction.momentmatching,
  fit.nosema.impinteraction.momentmatching,
  fit.nosema.bothinteractions.momentmatching,
  compare = TRUE,
  resp = "hasnosema",
  pointwise = FALSE,
  k_threshold = 0.7,
  moment_match = TRUE,
  model_names = NULL
)

loo_moment_match(
  fit.nosema.nointeractions.momentmatching,
  fit.nosema.natinteraction.momentmatching,
  fit.nosema.impinteraction.momentmatching,
  fit.nosema.bothinteractions.momentmatching,
  loofit,
  k_threshold = 0.7,
  check = FALSE
)


par = fvimp_brmsdf[fvimp_brmsdf$subsetPar == TRUE,]
pareto_k <- loofit$diagnostics$pareto_k
high_pareto_k_indices <- which(pareto_k > 0.7)
high_pareto_k_data_new <- par[high_pareto_k_indices, ]

# View or save the high k data points
print(high_pareto_k_data)


## **********************************************************
## Models for parasitism between species
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



