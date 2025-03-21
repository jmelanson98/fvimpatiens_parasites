## Create map grids of parasite intensity; plot semivariograms of raw data and model residuals

setwd('/Users/jenna1/Documents/UBC/Bombus Project/fvimpatiens_parasites')
rm(list=ls())
source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")
source("src/makeMapGrids.R")

## ***********************************************************************
## make map plots
## ***********************************************************************

#Parasite prevalence across sites
#make a df where parasites are pooled at sample point
data.workers = filter(data.par, caste == "worker")
data.par %>% 
  group_by(sample_pt) %>%
  summarize(site_hascrithidia = sum(hascrithidia)/n(),
            site_hasnosema = sum(hasnosema)/n(),
            site_apicystis = sum(apicystis)/n(),
            site_any = sum(any_parasite)/n(),
            site = min(site),
            long = min(long),
            lat = min(lat),
            numbees = n()) -> parbysite

#make a df where surveys are pooled, but keep track of the number of survey events per transect
fvimp_brmsdf %>%
  filter(impSubset == TRUE) %>%
  group_by(sample_pt) %>%
  summarize(
    site = min(site),
    site_any = mean(impatiens_abundance),
    long = min(long),
    lat = min(lat),
    numbees = n() #this is actually the number of survey events 
    #but I'm giving it this name so I can change less code in my function
  ) -> impatiensabundancepersite

#do this again but with native bee abundance as the "central" var
fvimp_brmsdf %>%
  filter(Subset == TRUE) %>%
  group_by(sample_pt) %>%
  summarize(
    site = min(site),
    site_any = mean(native_bee_abundance),
    long = min(long),
    lat = min(lat),
    numbees = n() #this is actually the number of survey events but I'm giving it this name so I can change less code below
  ) -> beeabundancepersite

makeMapGrids(groupedbysite = parbysite, 
             sampling_effort = "Number of\nspecimens", 
             var_of_interest = "Parasite\nprevalence")


## ***********************************************************************
## variograms testing for spatial autocorrelation
## ***********************************************************************
fvimp_sub = filter(fvimp_brmsdf, impSubset == TRUE)
fvimp_subpar = filter(fvimp_brmsdf, apidae ==1)

#calculate residuals for floral abundance (no predictors)
fabun<- brm(
  bf(floral_abundance ~ 1),
  data = fvimp_sub,
  family = student(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
fabun_resid <- as.data.frame(residuals(fabun))
fvimp_sub$fabun_resid_nopred = fabun_resid$Estimate

#calculate residuals for floral abundance (with predictors)
fabun<- brm(
  bf(floral_abundance ~ julian_date + I(julian_date^2) + (1|sample_pt)),
  data = fvimp_sub,
  family = student(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
fabun_resid <- as.data.frame(residuals(fabun))
fvimp_sub$fabun_resid_pred = fabun_resid$Estimate

#calculate residuals for floral diversity (with no predictors)
fdiv <- brm(
  bf(floral_diversity ~ 1),
  data = fvimp_sub,
  family = student(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
fdiv_resid <- as.data.frame(residuals(fdiv))
fvimp_sub$fdiv_resid_nopred = fdiv_resid$Estimate

#calculate residuals for floral diversity (predictors)
fdiv <- brm(
  bf(floral_diversity ~ julian_date + I(julian_date^2) + (1|sample_pt)),
  data = fvimp_sub,
  family = student(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
fdiv_resid <- as.data.frame(residuals(fdiv))
fvimp_sub$fdiv_resid_pred = fdiv_resid$Estimate

#calculate residuals for bombus diversity (no predictors)
brich <- brm(
  bf(bombus_richness | trials(9) ~ 1),
  data = fvimp_sub,
  family = beta_binomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
brich_resid = as.data.frame(residuals(brich))
fvimp_sub$brich_resid_nopred = brich_resid$Estimate

#calculate residuals for bombus diversity (with predictors)
brich <- brm(
  bf(bombus_richness | trials(9) ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
       (1|sample_pt)),
  data = fvimp_sub,
  family = beta_binomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
brich_resid = as.data.frame(residuals(brich))
fvimp_sub$brich_resid_pred = brich_resid$Estimate

#calculate residuals for wild bee abundance (no predictors)
babun <- brm(
  bf(native_bee_abundance ~ 1),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
babun_resid = as.data.frame(residuals(babun))
fvimp_sub$babun_resid_nopred = babun_resid$Estimate

#calculate residuals for wild bee abundance (with predictors)
babun <- brm(
  bf(native_bee_abundance ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
       (1|sample_pt)),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
babun_resid = as.data.frame(residuals(babun))
fvimp_sub$babun_resid_pred = babun_resid$Estimate

#calculate residuals for impatiens abundance (no predictors)
iabun <- brm(
  bf(impatiens_abundance ~ 1),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
iabun_resid = as.data.frame(residuals(iabun))
fvimp_sub$iabun_resid_nopred = iabun_resid$Estimate

#calculate residuals for impatiens abundance (with predictors)
iabun <- brm(
  bf(impatiens_abundance ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
       (1|sample_pt)),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
iabun_resid = as.data.frame(residuals(iabun))
fvimp_sub$iabun_resid_pred = iabun_resid$Estimate


########################################
### now parasite models
########################################

#calculate residuals for crithidia prevalence (no predictors)
hascrith <- brm(
  bf(hascrithidia ~ 1),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
crith_resid <- as.data.frame(residuals(hascrith))
fvimp_subpar$crith_resid_nopred = crith_resid$Estimate

#calculate residuals for crithidia prevalence (with predictors)
hascrith <- brm(
  bf(hascrithidia ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
crith_resid <- as.data.frame(residuals(hascrith))
fvimp_subpar$crith_resid_pred = crith_resid$Estimate

#calculate residuals for apicystis prevalence (no predictors)
hasapi <- brm(
  bf(apicystis ~ 1 ),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
api_resid <- as.data.frame(residuals(hasapi))
fvimp_subpar$api_resid_nopred = api_resid$Estimate

#calculate residuals for apicystis prevalence (with predictors)
hasapi <- brm(
  bf(apicystis ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
api_resid <- as.data.frame(residuals(hasapi))
fvimp_subpar$api_resid_pred = api_resid$Estimate

#calculate residuals for nosema prevalence (no predictors)
hasnos <- brm(
  bf(hasnosema ~ 1 ),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
nos_resid <- as.data.frame(residuals(hasnos))
fvimp_subpar$nos_resid_nopred = nos_resid$Estimate

#calculate residuals for nosema prevalence (with predictors)
hasnos <- brm(
  bf(hasnosema ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 1,
  thin = 1,
  init = 0,
  cores = 1,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
nos_resid <- as.data.frame(residuals(hasnos))
fvimp_subpar$nos_resid_pred = nos_resid$Estimate

#############################################
### save dataframe 
#############################################
save(fvimp_sub, fvimp_subpar,
     file="saved/data_with_residuals.Rdata")

load(file="saved/data_with_residuals.Rdata")

#############################################
### make data frames spatially explicit
#############################################

# Convert the data frames into sf objects
fvimp_sub <- st_as_sf(fvimp_sub, coords = c("long", "lat"), crs = 4326)
fvimp_subpar <- st_as_sf(fvimp_subpar, coords = c("long", "lat"), crs = 4326)

# Transform the coordinate reference system to EPSG:900913
fvimp_sub900913 <- st_transform(fvimp_sub, crs = 3857)
fvimp_subpar900913 <- st_transform(fvimp_subpar, crs = 3857)


############################################
### variograms with no predictors
############################################

#make variogram for floral abundance
v_fabun.resid.nopred <- variogram(fabun_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_fabun.resid.nopred_plot = ggplot(as.data.frame(v_fabun.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Floral Abundance",
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )
v_fabun.resid.nopred_plot

#make variogram for floral diversity
v_fdiv.resid.nopred <- variogram(fdiv_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_fdiv.resid.nopred_plot = ggplot(as.data.frame(v_fdiv.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Floral Diversity",
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )
v_fdiv.resid.nopred_plot


#make variogram for bombus diversity
v_brich.resid.nopred <- variogram(brich_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_brich.resid.nopred_plot = ggplot(as.data.frame(v_brich.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*Bombus* Species Richness"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_brich.resid.nopred_plot


#make variogram for native bombus abundance
v_babun.resid.nopred <- variogram(babun_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_babun.resid.nopred_plot = ggplot(as.data.frame(v_babun.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("Native *Bombus* Abundance"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_babun.resid.nopred_plot

#make variogram for impatiens abundance
v_iabun.resid.nopred <- variogram(iabun_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_iabun.resid.nopred_plot = ggplot(as.data.frame(v_iabun.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*B*. *impatiens* Abundance"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_iabun.resid.nopred_plot

#make variogram for crithidia prevalence
v_crith.resid.nopred <- variogram(crith_resid_nopred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_crith.resid.nopred_plot = ggplot(as.data.frame(v_crith.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*Crithidia spp* Prevalence"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_crith.resid.nopred_plot

#make variogram for apicystis prevalence
v_api.resid.nopred <- variogram(api_resid_nopred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_api.resid.nopred_plot = ggplot(as.data.frame(v_api.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*Apicystis spp* Prevalence"),
    x = "",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_api.resid.nopred_plot

#make variogram for nosema prevalence
v_nos.resid.nopred <- variogram(nos_resid_nopred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_nos.resid.nopred_plot = ggplot(as.data.frame(v_nos.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression("*Vairimorpha spp* Prevalence"),
    x = "Distance (meters)",
    y = "Semivariance"
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_nos.resid.nopred_plot


#############################################
### variograms with predictors
#############################################
#make variogram for floral abundance
v_fabun.resid.pred <- variogram(fabun_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_fabun.resid.pred_plot = ggplot(as.data.frame(v_fabun.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )
v_fabun.resid.pred_plot

#make variogram for floral diversity
v_fdiv.resid.pred <- variogram(fdiv_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_fdiv.resid.pred_plot = ggplot(as.data.frame(v_fdiv.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )
v_fdiv.resid.pred_plot


#make variogram for bombus diversity
v_brich.resid.pred <- variogram(brich_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_brich.resid.pred_plot = ggplot(as.data.frame(v_brich.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_brich.resid.pred_plot


#make variogram for native bombus abundance
v_babun.resid.pred <- variogram(babun_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_babun.resid.pred_plot = ggplot(as.data.frame(v_babun.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_babun.resid.pred_plot

#make variogram for impatiens abundance
v_iabun.resid.pred <- variogram(iabun_resid_pred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_iabun.resid.pred_plot = ggplot(as.data.frame(v_iabun.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_iabun.resid.pred_plot

#make variogram for crithidia prevalence
v_crith.resid.pred <- variogram(crith_resid_pred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_crith.resid.pred_plot = ggplot(as.data.frame(v_crith.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_crith.resid.pred_plot

#make variogram for apicystis prevalence
v_api.resid.pred <- variogram(api_resid_pred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_api.resid.pred_plot = ggplot(as.data.frame(v_api.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_api.resid.pred_plot

#make variogram for nosema prevalence
v_nos.resid.pred <- variogram(nos_resid_pred ~ 1, data = fvimp_subpar900913, cutoff = 2000, width = 250)

v_nos.resid.pred_plot = ggplot(as.data.frame(v_nos.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
    x = "Distance (meters)",
    y = ""
  ) +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()  
  )
v_nos.resid.pred_plot




variogramgrid = wrap_plots(v_fdiv.resid.nopred_plot, v_fdiv.resid.pred_plot, 
                           v_fabun.resid.nopred_plot, v_fabun.resid.pred_plot,
                           v_brich.resid.nopred_plot, v_brich.resid.pred_plot, 
                           v_babun.resid.nopred_plot, v_babun.resid.pred_plot,
                           v_iabun.resid.nopred_plot, v_iabun.resid.pred_plot,
                           v_crith.resid.nopred_plot, v_crith.resid.pred_plot,
                           v_api.resid.nopred_plot, v_api.resid.pred_plot,
                           v_nos.resid.nopred_plot, v_nos.resid.pred_plot,
                           ncol = 2)

variogramgrid

#add subplot labels
final_plot <- ggdraw() +
  draw_plot(variogramgrid, 0.02, 0, 1, 1) +
  draw_plot_label(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p"),
                  x = c(0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52, 0, 0.52), 
                  y = c(0.99, 0.99, 0.86, 0.86, 0.72, 0.72, 0.6, 0.6, 0.48, 0.48, 0.37, 0.37, 0.25, 0.25, 0.14, 0.14))
print(final_plot)



################################
### Calculate Moran's I
################################
library(sf)
library(spdep)

# Set a distance threshold (e.g., within 1000 meters)
coords <- st_coordinates(fvimp_sub900913)  # Extract coordinates
nb <- dnearneigh(coords, d1 = 0, d2 = 1000)  # Neighbors within 1000 units


# Convert to spatial weights
#ignore groups that don't have neighbors within specified distance
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Perform Moran's I (on raw data)
moran_result <- moran.test(fvimp_sub900913$impatiens_abundance, lw)
print(moran_result)

#okay so I think there is spatial autocorrelation on the raw values which makes sense
#however this spatial autocorrelation is quite weak (I < 0.1 for all relationships)

#test spatial autocorrelation on the model residuals
# Perform Moran's I (on residuals)
moran_result <- moran.test(fvimp_sub900913$fabun_resid_pred, lw)
print(moran_result)

#modelling accounts for spatial autocorrelation!! all p-values are like 0.99 after
#modelling and Moran's I stays low :)
#the only exception to this is the Bombus diversity model -- the pre- and post-model
#residuals show similar patterns of *very low but significant* spatial autocorrelation
#this is not surprising given how absolutely trash the diversity model R2 value is