standardize <- function(x)
  ## subtract mean and devide by sd to standardize data
(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)


makeDataMultiLevel <- function(indiv.data, site.col, year.col="Year"){
  ## Sets weights to 0s and 1s so that a single dataset can be passed
  ## to brms but site-level data isn't duplicated by the number of
  ## specimens

  ## split data by year
  indiv.data.split <- split(indiv.data, indiv.data[, year.col])

  out.indiv.data <- lapply(indiv.data.split, addWeightCol,
                           site.col=site.col)

  out.indiv.data <- do.call(rbind, out.indiv.data)

  return(out.indiv.data)

}


addWeightCol <- function(each.year.dat, site.col){
  ## for all individuals from the same site, 
  site.ids <- unlist(tapply(each.year.dat[, site.col],
                            each.year.dat[, site.col],
                            function(x) 1:length(x)))


  names(site.ids) <- NULL
  each.year.dat$SiteIDs <- site.ids
  each.year.dat$Weights <- each.year.dat$SiteIDs
  each.year.dat$Weights[each.year.dat$Weights > 1] <- 0
  return(each.year.dat)
}

standardizeVars <- function(spec.data, vars, key, by.stand=TRUE){
  ##  center all of the x variables, need to use unique values to avoid
  ##  repetition by the number of specimens
  if(by.stand){
    unique.site.vals <-  unique(spec.data[,c("Stand", key, vars)])
  } else {
    unique.site.vals <-  unique(spec.data[,c(key, vars)])
  }
  unique.site.vals[, vars] <- apply(as.matrix(unique.site.vals[, vars]), 2, standardize)
  print("Dimensions of the data before merging the standardize data")
  print(dim(spec.data))
  spec.data[, vars] <- NULL
  spec.data <- merge(spec.data, unique.site.vals, all.x=TRUE)
  print("Dimensions of the data after merging the standardize data")
  print(dim(spec.data))
  layout(matrix(1:(2*round(length(vars)/2)), nrow=2))
  for(var in vars){
    hist(unique.site.vals[, var], main=var)
  }
  
  return(spec.data)
}

prepParasiteWeights <- function(spec.data){
  ## create a dummy varaible "WeightPar" for the parasite data. The
  ## intention is to keep stan from dropping data for site-level models,
  ## but weight is 0 for parasite models.
  spec.data$WeightsPar <- 1
  spec.data$WeightsPar[spec.data$Apidae == 0 |
                       is.na(spec.data$Apidae)] <- 0
  ## stan drops all NA data, so can set AnyParasite to 0 with WeightsPar
  ## to keep it in the models
  spec.data$ParasitePresence[is.na(spec.data$ParasitePresence)] <- 0
  spec.data$CrithidiaBombi[is.na(spec.data$CrithidiaBombi)] <- 0
  spec.data$HasCrithidia[is.na(spec.data$HasCrithidia)] <- 0
  spec.data$ApicystisSpp[is.na(spec.data$ApicystisSpp)] <- 0
  spec.data$Year <- as.factor(spec.data$Year)
  return(spec.data)
}


prepDataSEM <-
  function(spec.data,#individual level specimen data
           variables.to.log = NULL, #variables to be logged
           variables.to.log.1 = NULL, #variables to be logged + 1
           vars_yearsr = NULL,#variables to standardize at year site
                                        #sampling round level
           vars_year = NULL,#variables to standardize at year site sampling round level 
           vars_sp = NULL,  #variables to standardize at the species
                                        #level
           vars_sp_yearsr= NULL) #variables to standardize at the species level
{
  ## Function for making the SEM weights and standarizing variables.
  spec.data <- spec.data[order(spec.data$Stand), ]
  
  ## create a dummy variable "Weight" to deal with the data sets being at
  ## different levels to get around the issue of having to pass in one
  ## data set into brms
  spec.data$YearSR <-
    paste(spec.data$Year, spec.data$DoyStart, sep = ";")

  ## species, year, SR key
  spec.data$YearSRGenusSpecies <-
    paste(spec.data$Year, spec.data$DoyStart, spec.data$GenusSpecies, sep = ";")
  
  print("Number of unique stand, year, sampling round combinations")
  print(length(unique(paste(spec.data$Stand, spec.data$YearSR))))
  spec.data <- makeDataMultiLevel(spec.data, "Stand", "YearSR")
  print("Number of individuals with Weights == 1, should be the same as above")
  print(sum(spec.data$Weights))

  ## log variables before standardizing
  if(!is.null(variables.to.log)){
    spec.data[, variables.to.log] <-
      log(spec.data[, variables.to.log])
  }
  if(!is.null(variables.to.log.1)){
    spec.data[, variables.to.log.1] <-
      log(spec.data[, variables.to.log.1]+ 1)
  }
  
  ##  center all of the x variables, need to use unique values to avoid
  ##  repetition by the number of specimens
  
  if(!is.null(vars_yearsr)){
    print("Standardizing variables with year, sampling round, stand combinations")
    spec.data <- standardizeVars(spec.data, vars_yearsr, "YearSR")
  }
  if(!is.null(vars_sp)){
    print("Standardizing variables by species")
    spec.data <-
      standardizeVars(spec.data, vars_sp, "GenusSpecies", by.stand=FALSE)
  }
  if(!is.null(vars_year)){
    print("Standardizing variables with year, sampling round, stand combinations")
    spec.data <- standardizeVars(spec.data, vars_yearsr, "YearSR")
  }

  if(!is.null(vars_sp_yearsr)){
    print("Standardizing variables with year, sampling round, stand, species combinations")
    spec.data <- standardizeVars(spec.data, vars_sp_yearsr, "YearSRGenusSpecies")
  }
  
  ## create a dumby varaible "WeightPar" for the parasite data. The
  ## original intention was to keep stan from dropping data for
  ## site-level models, but weight is 0 for parasite models.
  
  print("Number of successful parasite screenings")
  print(sum(spec.data$Apidae, na.rm = TRUE))
  spec.data <- prepParasiteWeights(spec.data)
  print("Number of of individuals with WeightsPar == 1, should be the same as above")
  print(sum(spec.data$WeightsPar))
  print("Final dim of data after adding WeightsPar")
  print(dim(spec.data))
  
  return(spec.data)
}
