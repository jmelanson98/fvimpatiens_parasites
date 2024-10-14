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

orig.spec <- fvimp_brmsdf

## Make SEM weights and standardize data.
#insert Bombus medius as filler species for NA rows
fvimp_brmsdf <- prepDataSEM_bernoulli(spec.data = fvimp_brmsdf,
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
#use different subsets because there were several days at the beginning of the season when we did not collect impatiens;
#these sampling efforts can be included in i) the vegetation models, ii) native bombus abundance models, because it does
#not effect either. But they must be subseted out of all other models; luckily we did not screen any bees for parasites from
#those days, so this does not effect subsetPar

formula.bee.div <- formula(bombus_shannon_diversity | subset(impSubset)~
                             floral_abundance + floral_diversity +
                              julian_date + I(julian_date^2) + 
                             landscape_shdi + prop_edge + prop_blueberry +
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
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************

xvars.fv <- c("bombus_shannon_diversity",
              "native_bee_abundance",
              "impatiens_abundance",
              "status",
              "floral_abundance",
              "floral_diversity",
              "prop_blueberry",
              "julian_date",
              "I(julian_date^2)",
              "(1|sample_pt)",
              "(1|subsite)",
              "(1|gr(final_id, cov = studycov))"
                 )

## **********************************************************
## Random effects of phylogeny
## **********************************************************
#create phylo variance-covariance matrix for group of study species
studyspecies = c("Bombus_mixtus", "Bombus_flavifrons", "Bombus_rufocinctus", "Bombus_californicus", 
                 "Bombus_impatiens", "Bombus_vosnesenskii", "Bombus_sitkensis", "Bombus_melanopygus", 
                 "Bombus_medius")
studycov = getPhyloMatrix(studyspecies)


## **********************************************************
## Run all models
## **********************************************************

formula.allcrith <-  runParasiteModels(fvimp_brmsdf,
                                        "hascrithidia", 
                                        xvars.fv)
formula.allnos <-  runParasiteModels(fvimp_brmsdf,
                                   "hasnosema", 
                                   xvars.fv)
formula.apicystis <-  runParasiteModels(fvimp_brmsdf,
                                        "apicystis", 
                                        xvars.fv)



bf.fdiv <- bf(formula.flower.div, family="student")
bf.fabund <- bf(formula.flower.abund, family = "student")
bf.bdiv <- bf(formula.bee.div, family="hurdle_lognormal")
bf.babund <- bf(formula.bee.abund, family = "negbinomial")
bf.iabund <- bf(formula.imp.abund, family = "negbinomial")

## convert to brms format
bf.allcrith <- bf(formula.allcrith, family="bernoulli")
bf.allnos <- bf(formula.allnos, family="bernoulli")
bf.api <- bf(formula.apicystis, family="bernoulli")
bform <-  bf.fdiv + bf.fabund + 
  bf.babund + bf.bdiv + bf.iabund + 
  bf.allcrith + bf.allnos + bf.api +
  set_rescor(FALSE)

##set priors for pesky interaction term
prior<-c(set_prior("normal(0,1)", class = "b", coef = "", resp = "hasnosema"))
prior<-c(set_prior("normal(0,1)", class = "b", coef = ""))

## run model
fit.bombus.all.prior.nointeractions <- brm(bform, fvimp_brmsdf,
                              cores=ncores,
                              iter = (10^4),
                              chains =3,
                              prior = prior,
                              thin=1,
                              init=0,
                              save_pars = save_pars(all = TRUE),
                              open_progress = FALSE,
                              control = list(adapt_delta = 0.999,
                                             stepsize = 0.001,
                                             max_treedepth = 20),
                              data2 = list(studycov = studycov)
)

write.ms.table(fit.bombus.all.prior.nointeractions, "AllModels_fv_bernoulli_phylo_prior_nointeractions")
save(fit.bombus.all.prior.natinteraction.reduced, fvimp_brmsdf, orig.spec,
     file="saved/AllModels_fv_bernoulli_phylo_prior_natinteraction_reduced.Rdata")

load(file="saved/AllModels_fv_bernoulli_phylo_prior_natinteraction_reduced.Rdata")

plot.res(fit.bombus.all.prior.natinteraction.reduced, "AllModels_fv_bernoulli_phylo_prior_natinteraction_reduced")

summary(fit.bombus.all.prior.natinteraction.reduced)

bayes_R2(fit.bombus.all.prior.natinteraction.reduced)

plot(pp_check(fit.bombus.all.prior.impinteraction, resp="floraldiversity"))
plot(pp_check(fit.bombus.all.prior.impinteraction, resp="floralabundance"))
plot(pp_check(fit.bombus.all.prior.impinteraction, resp="nativebeeabundance"))
plot(pp_check(fit.bombus.all.prior.impinteraction, resp="bombusshannondiversity"))
plot(pp_check(fit.bombus.all.prior.impinteraction, resp = "impatiensabundance"))
plot(pp_check(fit.bombus.all.prior.impinteraction, resp = "hascrithidia"))
plot(pp_check(fit.bombus.all.prior.impinteraction, resp = "hasnosema"))
plot(pp_check(fit.bombus.all.prior.impinteraction, resp = "apicystis"))


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
run_plot_freq_model_diagnostics(remove_subset_parasite_formula(formula.allnos),
                                this_data=fvimp_brmsdf[fvimp_brmsdf$subsetPar == TRUE,],
                                this_family="bernoulli", site.lat="site")

#run frequentist checks for models including impatiens
  run_plot_freq_model_diagnostics(remove_subset_imp_formula(formula.imp.abund),
                                  this_data=fvimp_brmsdf[fvimp_brmsdf$impSubset == TRUE,],
                                  this_family="negbinomial", site.lat="site")

#run frequentist checks for floral/native bee abundance models
run_plot_freq_model_diagnostics(remove_subset_other_formula(formula.flower.abund),
                                this_data=fvimp_brmsdf[fvimp_brmsdf$Subset == TRUE,],
                                this_family="students", site.lat="site")



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


par = fvimp_reduced[fvimp_reduced$subsetPar == TRUE,]
pareto_k <- loo_result$diagnostics$pareto_k
high_pareto_k_indices <- which(pareto_k > 0.7)
high_pareto_k_data_new <- par[high_pareto_k_indices, ]

# View or save the high k data points
print(high_pareto_k_data)