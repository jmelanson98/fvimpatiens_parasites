
pstars <- function(x){
  if(x >= 0.975){
    out <- "***"
  } else if(x < 0.975 & x >= 0.95){
    out <- "**"
  } else if(x < 0.95 & x >= 0.9){
    out <- "*"
  } else{
    out <- ""
  }
  return(out)
}


write.ms.table <- function(mod.output, mod.name){
    sum.mod <- as.data.frame(round(summary(mod.output)$fixed,2))

    coeffs <- c(paste0("b_",
                       rownames(sum.mod)),
                paste0("bs_",
                       rownames(sum.mod)))

    samps.mod <- posterior_samples(mod.output)

    coeffs <- coeffs[coeffs %in% colnames(samps.mod)]

    samps.mod <- samps.mod[, coeffs]

    coeff.samps <- colnames(samps.mod)
    coeff.samps.sum <- sub("[a-z]*_", "", coeff.samps)

   samps.mod <- samps.mod[order(match( coeff.samps.sum, rownames(sum.mod)))]

    sum.mod$Pgt0  <- round(apply(samps.mod, 2, function(x)
        sum(x > 0)/length(x)), 2)

    sum.mod$Plt0  <- round(apply(samps.mod, 2, function(x)
        sum(x < 0)/length(x)),2)

    sum.mod$Pgt0Stars  <- sapply(sum.mod$Pgt0, pstars)
    sum.mod$Plt0Stars  <- sapply(sum.mod$Plt0, pstars)

    write.table(sum.mod,
                file=sprintf("saved/tables/%s.txt", mod.name),
                sep="&")

    write.csv(sum.mod,
              file=sprintf("saved/tables/%s.csv", mod.name))
}


write.ms.table.differences <- function(mod.output, mod.name){
  # Get summary of fixed effects
  sum.mod <- as.data.frame(round(summary(mod.output)$fixed, 2))
  
  # Define coefficients for the model
  coeffs <- c(paste0("b_", rownames(sum.mod)), paste0("bs_", rownames(sum.mod)))
  
  # Extract posterior samples
  samps.mod <- posterior_samples(mod.output)
  
  # Ensure only the relevant coefficients are included
  coeffs <- coeffs[coeffs %in% colnames(samps.mod)]
  samps.mod <- samps.mod[, coeffs]
  
  coeff.samps <- colnames(samps.mod)
  coeff.samps.sum <- sub("[a-z]*_", "", coeff.samps)
  
  # Sort the posterior samples by matching with the summary
  samps.mod <- samps.mod[order(match(coeff.samps.sum, rownames(sum.mod)))]
  
  # Initialize a data frame for results
  diff.results <- data.frame()
  
  # Get the levels of the predictor (final_id)
  levels <- unique(rownames(sum.mod))  # Modify this line if necessary
  
  # Loop through the levels to compute differences
  for(i in 1:(length(levels)-1)){
    for(j in (i+1):length(levels)){
      # Create the names of the coefficients for the difference (e.g., speciesB - speciesA)
      level1 <- levels[i]
      level2 <- levels[j]
      coeff_name1 <- paste0("b_", level1)
      coeff_name2 <- paste0("b_", level2)
      
      # Check if the coefficients are in the posterior samples
      if(coeff_name1 %in% colnames(samps.mod) & coeff_name2 %in% colnames(samps.mod)){
        # Calculate the difference between the two coefficients
        diff.samps <- samps.mod[[coeff_name1]] - samps.mod[[coeff_name2]]
        
        # Calculate the proportion of posterior > or < 0
        p_gt_0 <- sum(diff.samps > 0) / length(diff.samps)
        p_lt_0 <- sum(diff.samps < 0) / length(diff.samps)
        
        #calculate a bayesian p value (1-posterior probability)
        bayes_p = 1-max(c(p_gt_0, p_lt_0))
        
        # Add results to the data frame
        diff.results <- rbind(diff.results, data.frame(Comparison = paste(level1, "-", level2),
                                                       Estimate = mean(diff.samps),
                                                       EstError = sd(diff.samps),
                                                       Pgt0 = round(p_gt_0, 4),
                                                       Plt0 = round(p_lt_0, 4),
                                                       BayesP = bayes_p,
                                                       Pgt0Stars = pstars(p_gt_0),
                                                       Plt0Stars = pstars(p_lt_0)))
      }
    }
  }
  diff.results$adjusted_p = p.adjust(diff.results$BayesP, method = "bonferroni")
  
  # Add results to the summary table
  sum.mod$Pgt0  <- round(apply(samps.mod, 2, function(x) sum(x > 0) / length(x)), 2)
  sum.mod$Plt0  <- round(apply(samps.mod, 2, function(x) sum(x < 0) / length(x)), 2)
  sum.mod$Pgt0Stars  <- sapply(sum.mod$Pgt0, pstars)
  sum.mod$Plt0Stars  <- sapply(sum.mod$Plt0, pstars)
  
  # Write results to files
  write.table(sum.mod, file = sprintf("saved/tables/%s.txt", mod.name), sep = "&")
  write.csv(sum.mod, file = sprintf("saved/tables/%s.csv", mod.name))
  
  # Write the pairwise comparisons (differences between levels)
  write.table(diff.results, file = sprintf("saved/tables/%s_pairwise_comparisons.txt", mod.name), sep = "&")
  write.csv(diff.results, file = sprintf("saved/tables/%s_pairwise_comparisons.csv", mod.name))
}

# pstars function for significance stars
pstars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}
