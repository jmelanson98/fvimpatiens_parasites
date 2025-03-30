## FVBombus impatiens parasite dataset
# code for creating a data frame for input into brms
# every row is a specimen, but create "empty" rows for sampling points which were surveyed 
#but specimens not found

rm(list=ls())
setwd('/Users/jenna1/Documents/UBC/Bombus Project/FVBombus_code/parasitecode')
#load packages
source('src/init.R')
load.packages()
source('src/prepDF.R')


## **********************************************************
## Load in necessary data
## **********************************************************
#load bombus survey data
specimenData = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")
sampleEffort = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022sampledata.csv", sep = ","))
vegData = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022vegetationdata.csv", sep = ","))

#load parasite data
parasiteScores = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022parasitedata.csv", sep = ","))
parasiteScores$barcode_id[parasiteScores$barcode_id == "HR1_14_01"] = "HR1_14A_01"

# load landscape raster/shapefiles
landscape1 <- raster(paste0(here(),"/landscape/FValley_lc_1res.tif"))
fv_points <- read_sf(paste0(here(),"/landscape/fvbombus/fvbombus_points.shp"))
fv_points2022 = fv_points[fv_points$site_id %in% sampleEffort$sample_point,]
landcover = read_csv(paste0(here(), "/landscape/landcover.csv"))


#create a list of all the flowers visited by bumblebees in 2022 + 2023 (used for veg abundance/diversity)
beeflowers = unique(c(specimenData$active_flower, specimenData2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]



## **********************************************************
## Prepare specimen dataframe
## **********************************************************
#this function will:
#### - add parasite scores to the specimen dataframe and clean unwanted individuals

specimenTable = prepSpecimenTable(specimenData, parasiteScores)



## **********************************************************
## Calculate landscape metrics
## **********************************************************
#this function will:
#### - calculate area of different land classes in 500m buffer around each sample_pt
#### - calculate landscape shdi in 500m buffer around each sample_pt
#### - return a dataframe containing prop_blueberry, prop_edge, and shdi for each sample_pt

landscapeMetrics = calculateLandscapeMetrics(landscape1, fv_points2022, landcover)

#OR
setwd('/Users/jenna1/Documents/UBC/Bombus Project/FVBombus_code/pollencode')
load(file="saved/landscapeMetrics.Rdata")
setwd('/Users/jenna1/Documents/UBC/Bombus Project/FVBombus_code/parasitecode')


## **********************************************************
## Prepare single dataframe for brms
## **********************************************************
#this function will add the following to the sample effort dataframe:
#### - add GPS coords
#### - julian date
#### - floral/bee diversity/abundance data
#### - landscape metrics (prop_blueberry, prop_edge, shdi)
# it will then merge the sample effort and specimen data frames, and remove columns that are not
# used by brms
# The output is a dataframe where every row is an observation of a single bee, or an "empty row" for a sampling
#event where no bees were captured. Suitable for bernoulli family parasite models.

brmsdf = prepBRMSdf(sampleEffort, samplePoints, specimenTable, beeflowers, vegData, landscapeMetrics)
write.csv(brmsdf, "fvimp_brmsdf.csv")

## **********************************************************
## Create "adjusted julian date" to account for colony phase
## **********************************************************
#calculate sample effort per day
#add julian date to sample effort data frame
sampleEffort$date = paste(sampleEffort$day, sampleEffort$month, sampleEffort$year)
sampleEffort$date = gsub(" ", "", sampleEffort$date, fixed = TRUE)
sampleEffort$date <- as.POSIXlt(sampleEffort$date, format = "%d%b%y")
sampleEffort$julian_date = sampleEffort$date$yday

sampleEffort %>% group_by(julian_date) %>%
  summarize(effort = 5*n()) -> effortperday

brmsdf %>% group_by(julian_date) %>%
  summarize(rufosum = sum(rufocinctus_abundance),
            mixsum = sum(mixtus_abundance),
            impsum = sum(impatiens_abundance),
            calsum = sum(californicus_abundance),
            vossum = sum(vos_abundance),
            flavsum = sum(flavifrons_abundance)) -> beesperday

calendar = left_join(beesperday, effortperday, by ="julian_date")
calendar$rufopermin = calendar$rufosum/calendar$effort
calendar$mixpermin = calendar$mixsum/calendar$effort
calendar$imppermin = calendar$impsum/calendar$effort
calendar$flavpermin = calendar$flavsum/calendar$effort
calendar$calpermin = calendar$calsum/calendar$effort
calendar$vospermin = calendar$vossum/calendar$effort

# calculate the weighted mean of the julian dates for each species
weighted_mean_mix <- sum(calendar$julian_date * calendar$mixpermin) / sum(calendar$mixpermin)
weighted_mean_imp <- sum(calendar$julian_date * calendar$imppermin) / sum(calendar$imppermin)
weighted_mean_flav <- sum(calendar$julian_date * calendar$flavpermin) / sum(calendar$flavpermin)
weighted_mean_rufo <- sum(calendar$julian_date * calendar$rufopermin) / sum(calendar$rufopermin)
weighted_mean_cal <- sum(calendar$julian_date * calendar$calpermin) / sum(calendar$calpermin)
weighted_mean_vos <- sum(calendar$julian_date * calendar$vospermin) / sum(calendar$vospermin)

# calculate the weighted variance of the julian dates, for each species
weighted_sd_mix <- sqrt(sum(calendar$mixpermin * (calendar$julian_date - weighted_mean_mix)^2) / sum(calendar$mixpermin))
weighted_sd_imp <- sqrt(sum(calendar$imppermin * (calendar$julian_date - weighted_mean_imp)^2) / sum(calendar$imppermin))
weighted_sd_flav <- sqrt(sum(calendar$flavpermin * (calendar$julian_date - weighted_mean_flav)^2) / sum(calendar$flavpermin))
weighted_sd_rufo <- sqrt(sum(calendar$rufopermin * (calendar$julian_date - weighted_mean_rufo)^2) / sum(calendar$rufopermin))
weighted_sd_cal <- sqrt(sum(calendar$calpermin * (calendar$julian_date - weighted_mean_cal)^2) / sum(calendar$calpermin))
weighted_sd_vos <- sqrt(sum(calendar$vospermin * (calendar$julian_date - weighted_mean_vos)^2) / sum(calendar$vospermin))

# calculate the weighted z-scores for each Julian date
calendar$weighted_z_mix <- (calendar$julian_date - weighted_mean_mix) / weighted_sd_mix
calendar$weighted_z_imp <- (calendar$julian_date - weighted_mean_imp) / weighted_sd_imp
calendar$weighted_z_flav <- (calendar$julian_date - weighted_mean_flav) / weighted_sd_flav
calendar$weighted_z_rufo <- (calendar$julian_date - weighted_mean_rufo) / weighted_sd_rufo
calendar$weighted_z_cal <- (calendar$julian_date - weighted_mean_cal) / weighted_sd_cal
calendar$weighted_z_vos <- (calendar$julian_date - weighted_mean_vos) / weighted_sd_vos

ggplot(calendar, aes(x = julian_date, y = mixpermin)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_point(data = calendar, aes(x = julian_date, y = weighted_z_mix, color = "mix")) +
  geom_point(data = calendar, aes(x = julian_date, y = weighted_z_imp, color = "imp")) +
  geom_point(data = calendar, aes(x = julian_date, y = weighted_z_flav, color = "flav")) +
  geom_point(data = calendar, aes(x = julian_date, y = weighted_z_rufo, color = "rufo")) +
  geom_point(data = calendar, aes(x = julian_date, y = weighted_z_cal, color = "cal")) +
  geom_point(data = calendar, aes(x = julian_date, y = weighted_z_vos, color = "vos")) +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Value",
    title = "Bar Plot of Values Over Time"
  )

colstokeep = c("julian_date", "weighted_z_mix","weighted_z_imp","weighted_z_flav",
               "weighted_z_rufo","weighted_z_cal","weighted_z_vos")
cal_z = calendar[,colnames(calendar) %in% colstokeep]
brmsdf = left_join(brmsdf, cal_z, by = "julian_date")

# create the adjusted julian date column using ifelse statements
brmsdf$adjusted_julian <- ifelse(
  brmsdf$final_id == "B. mixtus", brmsdf$weighted_z_mix,
  ifelse(
    brmsdf$final_id == "B. impatiens", brmsdf$weighted_z_imp,
    ifelse(
      brmsdf$final_id == "B. flavifrons", brmsdf$weighted_z_flav,
      ifelse(
        brmsdf$final_id == "B. rufocinctus", brmsdf$weighted_z_rufo,
        ifelse(
          brmsdf$final_id == "B. californicus", brmsdf$weighted_z_cal,
          ifelse(
            brmsdf$final_id == "B. vosnesenskii", brmsdf$weighted_z_vos, 0)
          )
        )
      )
    )
)
brmsdf$adjusted_julian[is.na(brmsdf$adjusted_julian)] <- 0
write.csv(brmsdf, "fvimp_brmsdf.csv")


## **********************************************************
## Update data frame for parasite model family = negbinomial
## **********************************************************
#####

# Modify this dataframe for use with a negative binomial for parasite prevalence at transect level
brmsdf %>%
  group_by(sample_pt, round, final_id, sample_id, subsite) %>%
  summarize(
    julian_date = mean(julian_date),
    impatiens_abundance = mean(impatiens_abundance),
    native_bee_abundance = mean(native_bee_abundance),
    bombus_shannon_diversity = mean(bombus_shannon_diversity),
    floral_diversity = mean(floral_diversity),
    floral_abundance = mean(floral_abundance),
    prop_blueberry = mean(prop_blueberry),
    prop_edge = mean(prop_edge),
    landscape_shdi = mean(landscape_shdi),
    number_screened = sum(apidae),
    spp_apicystis = sum(apicystis),
    spp_nbombi = sum(nbombii),
    spp_ascosphaera = sum(ascosphaera),
    spp_nceranae = sum(nceranae),
    spp_crithidia = sum(crithidiaspp),
    spp_cexpoeki = sum(cexpoeki),
    spp_cbombi = sum(cbombii),
    spp_anycrithidia = sum(hascrithidia),
    spp_anynosema = sum(hasnosema),
    spp_anyparasite = sum(any_parasite)
  ) -> negbin_brmsdf

negbin_brmsdf$weightsPar = negbin_brmsdf$number_screened
negbin_brmsdf$weightsPar[negbin_brmsdf$weightsPar > 0] = 1
write.csv(negbin_brmsdf, "negativebinomial_brmsdf.csv")




## **********************************************************
## Covariance checks
## **********************************************************
#check for covariance--largest covariance 0.07
samp_data = data.frame(floral_abundance = sampleEffort$floral_abundance,
                       floral_diversity = sampleEffort$floral_diversity,
                       prop_blue = sampleEffort$prop_blueberry,
                       prop_edge = sampleEffort$prop_edge,
                       shdi = sampleEffort$landscape_shdi)
samp_data = na.omit(samp_data)
covariance_mat = cov(samp_data)
correlation_mat = cor(samp_data) #correlation a bit higher but all below .52

#check covariance again but with no repeats in landscape values -- largest covariance 0.05
sampleEffort %>% group_by(sample_pt) %>%
  summarize(prop_blue = mean(prop_blueberry),
            prop_edge = mean(prop_edge),
            shdi = mean(landscape_shdi)) -> reduced
reduced = na.omit(reduced[,-1])
covariance_mat2 = cov(reduced)
correlation_mat2 = cor(reduced) #again....higher but below 0.55




## **********************************************************
## Exploratory plotting
## **********************************************************

#scatterplot of landscape metrics vs parasite prev (parasites grouped at transect level)
library(ggplot2)
library(cowplot)

brmsdf %>%
  filter(apidae == 1) %>%
  group_by(sample_pt, final_id) %>%
  summarize(blueberry = mean(prop_blueberry, na.rm = TRUE),
            shdi = mean(landscape_shdi, na.rm = TRUE),
            prop_edge = mean(prop_edge, na.rm = TRUE),
            parasites = sum(any_parasite, na.rm = TRUE),
            total_par = n()) -> by_transect
by_transect$parasites = by_transect$parasites/by_transect$total_par

p1 <- ggplot(by_transect, aes(blueberry, parasites, color = final_id)) + 
  geom_jitter() +
  xlab("Proportion blueberry (500m buffer)") + 
  ylab("Parasite prevalence") +
  theme_cowplot(12)
p2 <- ggplot(by_transect, aes(shdi, parasites, color = final_id)) + 
  geom_jitter() +
  xlab("Landscape Shannon Diversity (500m buffer") + 
  ylab("Parasite prevalence") +
  theme_cowplot(12)
p3 <- ggplot(by_transect, aes(prop_edge, parasites, color = final_id)) + 
  geom_jitter() +
  xlab("Proportion edge area (500m buffer") + 
  ylab("Parasite prevalence") +
  theme_cowplot(12)
plot_grid(p1, p2, p3)

#scatterplots of community indices (parasites grouped at sampling event level)
brmsdf %>%
  filter(apidae == 1) %>%
  group_by(sample_id) %>%
  summarize(impatiens = mean(impatiens_abundance, na.rm = TRUE),
            native_abundance = mean(native_bee_abundance, na.rm = TRUE),
            prop_edge = mean(prop_edge, na.rm = TRUE),
            prop_blueberry = mean(prop_blueberry, na.rm = TRUE),
            shdi = mean(landscape_shdi, na.rm = TRUE),
            fdiv = mean(floral_diversity, na.rm = TRUE),
            fabun = mean(floral_abundance, na.rm = TRUE),
            julian_date = mean(julian_date, na.rm = TRUE),
            parasites = sum(apicystis, na.rm = TRUE),
            total_par = n()) -> by_sampling_event
by_sampling_event$parasites = by_sampling_event$parasites/by_sampling_event$total_par

#plots of parasites vs community
ggplot(by_sampling_event, aes(poly(julian_date,2)[2], parasites, color = log(native_abundance))) + 
  geom_jitter() +
  xlab("Julian Date") + 
  ylab("Apicystis prevalence") +
  theme_cowplot(12)
p5 <- ggplot(by_sampling_event, aes(native_abundance, parasites, color = final_id)) + 
  geom_jitter() +
  xlab("Native bombus abundance") + 
  ylab("Parasite prevalence") +
  theme_cowplot(12)
p6 <- ggplot(by_sampling_event, aes(fdiv, parasites, color = final_id)) + 
  geom_jitter() +
  xlab("Floral diversity") + 
  ylab("Parasite prevalence") +
  theme_cowplot(12)
p7 <- ggplot(by_sampling_event, aes(fabun, parasites, color = final_id)) + 
  geom_jitter() +
  xlab("Floral abundance") + 
  ylab("Parasite prevalence") +
  theme_cowplot(12)
plot_grid(p4, p5, p6, p7)

#plots of landscape vs community
p8 <- ggplot(by_sampling_event, aes(prop_blueberry, native_abundance, color = julian_date)) + 
  geom_point() +
  xlab("Proportion blueberry (500m buffer)") + 
  ylab("Native bombus abundance") +
  theme_cowplot(12)
p9 <- ggplot(by_sampling_event, aes(prop_edge, native_abundance, color = julian_date)) + 
  geom_point() +
  xlab("Proportion edge area (500m buffer)") + 
  ylab("Native bombus abundance") +
  theme_cowplot(12)
p10 <- ggplot(by_sampling_event, aes(shdi, native_abundance, color = julian_date)) + 
  geom_point() +
  xlab("Landscape SHDI (500m buffer)") + 
  ylab("Native bombus abundance") +
  theme_cowplot(12)
plot_grid(p8, p9, p10)


#scatterplot of julian date
brmsdf %>%
  filter(apidae == 1) %>%
  group_by(julian_date, final_id) %>%
  summarize(julian_date = mean(julian_date, na.rm = TRUE),
            parasites = sum(any_parasite, na.rm = TRUE),
            total_par = n()) -> by_date
by_date$parasites = by_date$parasites/by_date$total_par

ggplot(by_date, aes(julian_date, parasites, color = final_id)) + 
  geom_point() +
  xlab("Julian Date") + 
  ylab("Parasite prevalence") +
  theme_cowplot(12)
