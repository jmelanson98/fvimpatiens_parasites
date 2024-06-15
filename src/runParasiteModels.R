
runParasiteModels <- function(spec.data,
                              parasite,
                              xvars){

  formula.parasite  <- as.formula(paste(
    paste(parasite, "| weights(weightsPar)"),  #+ trials(1)
    paste(xvars,
          collapse=" + "),
    sep=" ~ "))

  return(formula.parasite)
}
