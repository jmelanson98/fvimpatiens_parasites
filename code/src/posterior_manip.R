get_ame = function(fit,
                   subset_data, 
                   variable, # as string
                   increment){
  library(tidybayes)
  
  # Build counterfactual datasets
  newdata_low = subset_data
  newdata_high = subset_data
  
  newdata_low[variable]  = 0
  newdata_high[variable] = increment
  
  # Get posterior draws on outcome scale
  draws_low  = posterior_epred(fit, newdata = newdata_low,  allow_new_levels = TRUE)
  draws_high = posterior_epred(fit, newdata = newdata_high, allow_new_levels = TRUE)
  
  # Contrast within each draw, then marginalize
  contrast_draws = rowMeans(draws_high - draws_low)
  
  # Summarize
  mean(contrast_draws)
  quantile(contrast_draws, c(0.025, 0.975))
}