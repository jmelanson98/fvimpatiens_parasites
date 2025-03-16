## Script started: March 15, 2025

# Took Lizzie Wolkovich's Bayesian hierarchical modeling class and if I learned one thing
# it's that there's no point in being Bayesian if I'm not gonna simulate, simulate, simulate.
# Given how far I was into analyzing this dataset, I was like "eh, I'll do it next time."
# But then I had a conversation with one of Lizzie's students, and it occurred to me that
# I would never show this monstrosity of a model output to Lizzie without checking that it actually works.

# And then I was like....but you would show it to the rest of the scientific community?? with your name on it??? for forever????

#So here we are...


# The first thing I will do (most pressing) is simulate data to see how often its
# DHARMa residuals fail a KS test (or test of quantile deviations). This will be
# especially useful because my parasite models DO fail this test, although the magnitude
# of the deviation is (visually) minor.

# use rstanarm for precompiled code--so that I can run multiple simple models quickly
library(rstanarm)
library(fdrtool)

# simulate data for a bernoulli distribution with one predictor (impatiens abundance)

N = 750
sim_imp = rhalfnorm(N, sd2theta(5))

b_imp = -0.5
intercept = 0

eta = intercept + b_imp*sim_imp
prob = 1 / (1 + exp(-eta))

sim_crith = rbinom(N, 1, prob)


#fit simple (correct) model to data
simple.fit = rstanarm::stan_glm(sim_crith ~ sim_imp, family = "binomial")

#change dharma helper function to work with stanarm object
check_rstanarm <- function(model,             # brms model
                       integer = TRUE,   # integer response? (TRUE/FALSE)
                       plot = TRUE,       # make plot?
                       ...                # further arguments for DHARMa::plotResiduals
) {
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(rstanarm::posterior_predict(model, ndraws = 1000)),
    observedResponse = model$y,
    fittedPredictedResponse = apply(
      t(rstanarm::posterior_epred(model, ndraws = 1000, re.form = NA)),
      1,
      mean),
    integerResponse = integer)
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  invisible(dharma.obj)
}

#run DHARMa on simulated data
checked.simple.sim = check_rstanarm(simple.fit)


#okayyyy now loop 'er
pvals = c()

for (i in 1:100) {
  #simulate data
  N = 750
  sim_imp = rhalfnorm(N, sd2theta(5))
  
  b_imp = -0.5
  intercept = 0
  
  eta = intercept + b_imp*sim_imp
  prob = 1 / (1 + exp(-eta))
  
  sim_crith = rbinom(N, 1, prob)
  
  #fit model
  loop.fit = rstanarm::stan_glm(sim_crith ~ sim_imp, family = "binomial")
  
  #make DHARMa object
  loop.dharma = check_rstanarm(loop.fit, plot = FALSE)
  KS = testUniformity(loop.dharma, plot = FALSE)
  
  #add KS test p value to the list pvalues
  pvals = c(pvals, KS$p.value)
  
}


# alright so that's pretty spot on. what if we add some noise variables...ones that
# don't matter for outcome?
pvals_noise = c()

for (i in 1:100) {
  #simulate data
  N = 750
  sim_imp = rhalfnorm(N, sd2theta(5))
  noise_1 = runif(N, 0, 10)
  noise_2 = as.factor(rep(1:5, N/5))
  
  b_imp = -0.5
  intercept = 0
  
  eta = intercept + b_imp*sim_imp
  prob = 1 / (1 + exp(-eta))
  
  sim_crith = rbinom(N, 1, prob)
  
  #fit model
  loop.fit = rstanarm::stan_glm(sim_crith ~ sim_imp + noise_1 + noise_2, family = "binomial")
  
  #make DHARMa object
  loop.dharma = check_rstanarm(loop.fit, plot = FALSE)
  KS = testUniformity(loop.dharma, plot = FALSE)
  
  #add KS test p value to the list pvalues
  pvals_noise = c(pvals_noise, KS$p.value)
  
}



# what about random effects? epred draw zeroes out random effects and i'm not
# sure how that will impact the KS test.

# first simulate with a model MISSING the group effect. the residuals should show 
# overdispersion...I think?
pvals_random = c()
disp_random = c()

for (i in 1:10) {
  #simulate data
  N = 750
  numgroups = 10
  group = as.factor(rep(1:numgroups, N/numgroups))
  group_sds = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  group_effects = sin(seq(0, pi, length.out = numgroups)) * 5
  sim_imp = rhalfnorm(N, sd2theta(5))
  noise_1 = runif(N, 0, 10)
  data = as.data.frame(cbind(group, sim_imp, noise_1))
  
  b_imp = -1.5
  #group_effects = rnorm(numgroups, 5, 5)
  intercept = 0
  
  data$eta = intercept + b_imp*data$sim_imp + group_effects[data$group]*sim_imp
  data$prob = 1 / (1 + exp(-data$eta))
  
  data$sim_crith = rbinom(N, 1, data$prob)
  
  #fit model
  loop.fit = rstanarm::stan_glm(sim_crith ~ sim_imp, data = data, family = "binomial")
  
  #make DHARMa object
  loop.dharma = check_rstanarm(loop.fit, plot = FALSE)
  KS = testUniformity(loop.dharma, plot = FALSE)
  disp = testDispersion(loop.dharma, plot = FALSE)
  
  #add KS test p value to the list pvalues
  pvals_random = c(pvals_random, KS$p.value)
  disp_random = c(disp_random, disp$p.value)
  
}
