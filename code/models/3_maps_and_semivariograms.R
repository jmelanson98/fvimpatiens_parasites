## Create map grids of parasite intensity; plot semivariograms of raw data and model residuals

setwd('/Users/jenna1/Documents/UBC/bombus_project/fvimpatiens_parasites')
rm(list=ls())
source("code/src/ggplotThemes.R")
source("code/src/init.R")
source("code/src/misc.R")
source("code/src/makeMapGrids.R")
source("code/src/getPhyloMatrix.R")
source("code/src/writeResultsTable.R")

# load r data
fvimp_brmsdf <- read.csv("data/fvimp_brmsdf.csv", sep = ",", header = T, row.names = 1)
load(file="saved/Base500m_32610.Rdata")
load(file="saved/NativePar500m_32610.Rdata")
load(file="saved/ImpatiensPar500m_32610.Rdata")

  
## ***********************************************************************
## make map plots
## ***********************************************************************

#Parasite prevalence across sites
#make a df where parasites are pooled at sample point
data.workers = filter(data.par, caste == "worker")
data.workers %>% 
  group_by(sample_pt) %>%
  summarize(site_hascrithidia = sum(hascrithidia)/n(),
            site_hasnosema = sum(hasnosema)/n(),
            site_apicystis = sum(apicystis)/n(),
            site_any = sum(any_parasite)/n(),
            site = min(site),
            long = min(long),
            lat = min(lat),
            numbees = n()) -> parbysite

# #make a df where surveys are pooled, but keep track of the number of survey events per transect
# fvimp_brmsdf %>%
#   filter(impSubset == TRUE) %>%
#   group_by(sample_pt) %>%
#   summarize(
#     site = min(site),
#     site_any = mean(impatiens_abundance),
#     long = min(long),
#     lat = min(lat),
#     numbees = n() #this is actually the number of survey events 
#     #but I'm giving it this name so I can change less code in my function
#   ) -> impatiensabundancepersite
# 
# #do this again but with native bee abundance as the "central" var
fvimp_brmsdf %>%
  filter(!is.na(barcode_id)) %>%
  filter(round == 1 | round ==2) %>%
  group_by(sample_pt) %>%
  summarize(
    site = min(site),
    site_any = 1,
    long = min(long),
    lat = min(lat),
    numbees = n() #this is actually the number of survey events but I'm giving it this name so I can change less code below
  ) -> beeabundancepersite

#note: if you're planning to run this code you need to UNZIP the file "zipped_clipped_rasters.zip"
#i've compressed it to save space locally and so that it can be pushed to github
#you will receive a warning that some rows are outside the scale range; don't worry about it
abundancemap_period1 = makeMapGrids(groupedbysite = beeabundancepersite, 
                      sampling_effort = "Number of\nspecimens", 
                      var_of_interest = "NA")

parmap = makeMapGrids(groupedbysite = parbysite, 
             sampling_effort = "Number of\nspecimens", 
             var_of_interest = "Parasite\nprevalence")

ggsave(filename = "figures/manuscript_figures/parasitemap.jpg", 
       plot = parmap, 
       width = 3000, 
       height = 2000, 
       units = "px")


## ***********************************************************************
## Variograms testing for spatial autocorrelation
## ***********************************************************************
fvimp_sub = filter(fvimp_brmsdf, impSubset == TRUE)
fvimp_par = filter(fvimp_brmsdf, subsetImpPar ==1)
fvnative_par = filter(fvimp_brmsdf, subsetNativePar ==1)
studyspecies = c("Bombus_mixtus", "Bombus_flavifrons", "Bombus_rufocinctus", "Bombus_californicus", 
                 "Bombus_impatiens", "Bombus_vosnesenskii", "Bombus_sitkensis", "Bombus_melanopygus", 
                 "Bombus_nevadensis", "Bombus_medius")
studycov = getPhyloMatrix(studyspecies)

runSVmodels = function(fvimp_sub,
                       fvimp_par,
                       fvnative_par,
                       studycov){
#calculate residuals for bombus richness (with predictors)
brich <- brm(
  bf(bombus_richness | trials(9) ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 +
       (1|sample_pt)),
  data = fvimp_sub,
  family = beta_binomial(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
brich_resid = as.data.frame(residuals(brich))
fvimp_sub$brich_resid_pred = brich_resid$Estimate

#calculate residuals for wild bee abundance (with predictors)
babun <- brm(
  bf(native_bee_abundance ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 +
       (1|sample_pt)),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4))
babun_resid = as.data.frame(residuals(babun))
fvimp_sub$babun_resid_pred = babun_resid$Estimate

#calculate residuals for impatiens abundance (with predictors)
iabun <- brm(
  bf(impatiens_abundance ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 +
       (1|sample_pt)),
  data = fvimp_sub,
  family = negbinomial(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
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

#calculate residuals for crithidia prevalence (with predictors)
hascrith <- brm(
  bf(hascrithidia ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 + caste +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvnative_par,
  family = bernoulli(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
crith_resid <- as.data.frame(residuals(hascrith))
fvnative_par$crith_resid_pred = crith_resid$Estimate

#calculate residuals for apicystis prevalence (with predictors)
hasapi <- brm(
  bf(apicystis ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 + caste +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvnative_par,
  family = bernoulli(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
api_resid <- as.data.frame(residuals(hasapi))
fvnative_par$api_resid_pred = api_resid$Estimate

#calculate residuals for nosema prevalence (with predictors)
hasnos <- brm(
  bf(hasnosema ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 + caste +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvnative_par,
  family = bernoulli(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
nos_resid <- as.data.frame(residuals(hasnos))
fvnative_par$nos_resid_pred = nos_resid$Estimate


#calculate residuals for crithidia prevalence (with predictors)
hascrith <- brm(
  bf(hascrithidia ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite)),
  data = fvimp_par,
  family = bernoulli(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
crith_resid <- as.data.frame(residuals(hascrith))
fvimp_par$crith_resid_pred = crith_resid$Estimate

#calculate residuals for apicystis prevalence (with predictors)
hasapi <- brm(
  bf(apicystis ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry_500 + prop_edge_500 + landscape_shdi_500 +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite)),
  data = fvimp_par,
  family = bernoulli(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
  open_progress = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001,
                 max_treedepth = 20),
  iter = (10^4),
  data2 = list(studycov = studycov))
api_resid <- as.data.frame(residuals(hasapi))
fvimp_par$api_resid_pred = api_resid$Estimate


return(list(fvimp_sub, fvimp_par, fvnative_par))}

#only run this code if models have changed -- otherwise load from .Rdata file
dfs = runSVmodels(fvimp_sub = fvimp_sub,
                  fvnative_par = fvnative_par,
                  fvimp_par = fvimp_par,
                  studycov = studycov)
fvimp_withresid = dfs[[1]]
fvimp_par_withresid = dfs[[2]]
fvimp_native_withresid = dfs[[3]]

#############################################
### save dataframe 
#############################################
save(fvimp_withresid, fvimp_par_withresid, fvimp_native_withresid,
     file="saved/data_with_residuals.Rdata")

load(file="saved/data_with_residuals.Rdata")

#############################################
### make data frames spatially explicit
#############################################

# Convert the data frames into sf objects
fvimp_sub = st_as_sf(fvimp_withresid, coords = c("long", "lat"), crs = 4326)
fvimp_par = st_as_sf(fvimp_par_withresid, coords = c("long", "lat"), crs = 4326)
fvnative_par = st_as_sf(fvimp_native_withresid, coords = c("long", "lat"), crs = 4326)

# Transform the coordinate reference system to EPSG:900913
fvimp_sub32610 = st_transform(fvimp_sub, crs = 32610)
fvimp_par32610 = st_transform(fvimp_par, crs = 32610)
fvnative_par32610 = st_transform(fvnative_par, crs = 32610)


#############################################
### variograms with predictors
#############################################

#make variogram for bombus diversity
v_brich.resid.pred <- variogram(brich_resid_pred ~ 1, data = fvimp_sub32610, cutoff = 2000, width = 250)

v_brich.resid.pred_plot = ggplot(as.data.frame(v_brich.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = "*Bombus* species richness",
    x = "",
    y = "Semivariance") +
  theme_minimal() +
  theme(title = ggtext::element_markdown(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20))
v_brich.resid.pred_plot


#make variogram for native bombus abundance
v_babun.resid.pred <- variogram(babun_resid_pred ~ 1, data = fvimp_sub32610, cutoff = 2000, width = 250)

v_babun.resid.pred_plot = ggplot(as.data.frame(v_babun.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = "Native *Bombus* abundance",
    x = "",
    y = "") +
  theme_minimal() +
  theme(title = ggtext::element_markdown(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20))
v_babun.resid.pred_plot

#make variogram for impatiens abundance
v_iabun.resid.pred <- variogram(iabun_resid_pred ~ 1, data = fvimp_sub32610, cutoff = 2000, width = 250)

v_iabun.resid.pred_plot = ggplot(as.data.frame(v_iabun.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = "*B. impatiens* abundance",
    x = "",
    y = "") +
  theme_minimal() +
  theme(title = ggtext::element_markdown(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20))
v_iabun.resid.pred_plot

#make variogram for crithidia prevalence -- native Bombus
v_crith.resid.pred <- variogram(crith_resid_pred ~ 1, data = fvnative_par32610, cutoff = 2000, width = 250)

v_crith.resid.pred_plot_native = ggplot(as.data.frame(v_crith.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = "*Crithidia* prev. -- native *Bombus*",
    x = "",
    y = "Semivariance") +
  theme_minimal() +
  theme(
    title = ggtext::element_markdown(),  
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20))
v_crith.resid.pred_plot_native

#make variogram for apicystis prevalence -- native Bombus
v_api.resid.pred <- variogram(api_resid_pred ~ 1, data = fvnative_par32610, cutoff = 2000, width = 250)

v_api.resid.pred_plot_native = ggplot(as.data.frame(v_api.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = "*Apicystis* prev. -- native *Bombus*",
    x = "",
    y = "") +
  theme_minimal() +
  theme(title = ggtext::element_markdown(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20))
v_api.resid.pred_plot_native

#make variogram for nosema prevalence -- native Bombus
v_nos.resid.pred <- variogram(nos_resid_pred ~ 1, data = fvnative_par32610, cutoff = 2000, width = 250)

v_nos.resid.pred_plot_native = ggplot(as.data.frame(v_nos.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = "*Vairimorpha* prev. -- native *Bombus*",
    x = "Distance (meters)",
    y = "") +
  theme_minimal() +
  theme(title = ggtext::element_markdown(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20))
v_nos.resid.pred_plot_native



#make variogram for crithidia prevalence -- Bombus impatiens
v_crith.resid.pred <- variogram(crith_resid_pred ~ 1, data = fvimp_par32610, cutoff = 2000, width = 250)

v_crith.resid.pred_plot_impatiens = ggplot(as.data.frame(v_crith.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = "*Crithidia* prev. -- *B. impatiens*",
    x = "Distance (meters)",
    y = "Semivariance") +
  theme_minimal() +
  theme(title = ggtext::element_markdown(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20))
v_crith.resid.pred_plot_impatiens

#make variogram for apicystis prevalence -- Bombus impatiens
v_api.resid.pred <- variogram(api_resid_pred ~ 1, data = fvimp_par32610, cutoff = 2000, width = 250)

v_api.resid.pred_plot_impatiens = ggplot(as.data.frame(v_api.resid.pred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = "*Apicystis* prev. -- *B. impatiens*",
    x = "",
    y = "Semivariance") +
  theme_minimal() +
  theme(title = ggtext::element_markdown(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20))
v_api.resid.pred_plot_impatiens

variogramgrid = grid.arrange(v_brich.resid.pred_plot, v_babun.resid.pred_plot, v_iabun.resid.pred_plot,
                             v_crith.resid.pred_plot_native, v_api.resid.pred_plot_native, v_nos.resid.pred_plot_native,
                             v_crith.resid.pred_plot_impatiens, v_api.resid.pred_plot_impatiens, ncol = 3)
  
#add labels
variogramgrid <- ggdraw() +
  draw_plot(variogramgrid, 0.01, 0, 0.93, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(G)", "(H)"),
                  x = c(0, 0.32, 0.63, 0, 0.32, 0.63, 0, 0.32), 
                  y = c(0.98, 0.98, 0.98, 0.67, 0.67, 0.67, 0.33, 0.33))
#export and save
ggsave(filename = "figures/manuscript_figures/variograms.jpg", 
       plot = variogramgrid, 
       width = 3000, 
       height = 2000, 
       units = "px")



################################
### Calculate Moran's I
################################
library(sf)
library(spdep)

# extract coordinates and get neighbors within 1000 units
coords <- st_coordinates(fvimp_sub32610) 
nb <- dnearneigh(coords, d1 = 0, d2 = 1000)

coords_par <- st_coordinates(fvimp_par32610) 
nb_par_imp <- dnearneigh(coords_par, d1 = 0, d2 = 1000)

coords_par <- st_coordinates(fvnative_par32610) 
nb_par_native <- dnearneigh(coords_par, d1 = 0, d2 = 1000)

# convert to spatial weights
#ignore groups that don't have neighbors within specified distance
lw = nb2listw(nb, style = "W", zero.policy = TRUE)
lwparimp = nb2listw(nb_par_imp, style = "W", zero.policy = TRUE)
lwparnative = nb2listw(nb_par_native, style = "W", zero.policy = TRUE)

# Perform Moran's I on the model residuals
moran_nbees_resid <- moran.test(fvimp_sub32610$babun_resid_pred, lw)
moran_richness_resid <- moran.test(fvimp_sub32610$brich_resid_pred, lw)
moran_impatiens_resid <- moran.test(fvimp_sub32610$iabun_resid_pred, lw)
moran_crith_resid_nat <- moran.test(fvnative_par32610$crith_resid_pred, lwparnative)
moran_api_resid_nat <- moran.test(fvnative_par32610$api_resid_pred, lwparnative)
moran_nos_resid_nat <- moran.test(fvnative_par32610$nos_resid_pred, lwparnative)
moran_crith_resid_imp <- moran.test(fvimp_par32610$crith_resid_pred, lwparimp)
moran_api_resid_imp <- moran.test(fvimp_par32610$api_resid_pred, lwparimp)



# extract coordinates and get neighbors within 1000 units
coords <- st_coordinates(fvimp_sub32610) 
nb <- dnearneigh(coords, d1 = 0, d2 = 500)

coords_par <- st_coordinates(fvimp_par32610) 
nb_par_imp <- dnearneigh(coords_par, d1 = 0, d2 = 500)

coords_par <- st_coordinates(fvnative_par32610) 
nb_par_native <- dnearneigh(coords_par, d1 = 0, d2 = 500)

# convert to spatial weights
#ignore groups that don't have neighbors within specified distance
lw = nb2listw(nb, style = "W", zero.policy = TRUE)
lwparimp = nb2listw(nb_par_imp, style = "W", zero.policy = TRUE)
lwparnative = nb2listw(nb_par_native, style = "W", zero.policy = TRUE)

# Perform Moran's I on the model residuals
moran_nbees_resid_500 <- moran.test(fvimp_sub32610$babun_resid_pred, lw)
moran_richness_resid_500 <- moran.test(fvimp_sub32610$brich_resid_pred, lw)
moran_impatiens_resid_500 <- moran.test(fvimp_sub32610$iabun_resid_pred, lw)
moran_crith_resid_nat_500 <- moran.test(fvnative_par32610$crith_resid_pred, lwparnative)
moran_api_resid_nat_500 <- moran.test(fvnative_par32610$api_resid_pred, lwparnative)
moran_nos_resid_nat_500 <- moran.test(fvnative_par32610$nos_resid_pred, lwparnative)
moran_crith_resid_imp_500 <- moran.test(fvimp_par32610$crith_resid_pred, lwparimp)
moran_api_resid_imp_500 <- moran.test(fvimp_par32610$api_resid_pred, lwparimp)

# Make list of results
results = list(
  "\\emph{Bombus} species richness (500m)" = list(test = moran_richness_resid_500, scale = "500m"),
  "\\emph{Bombus} species richness (1000m)" = list(test = moran_richness_resid, scale = "1000m"),
  "Native \\emph{Bombus} abundance (500m)" = list(test = moran_nbees_resid_500, scale = "500m"),
  "Native \\emph{Bombus} abundance (1000m)" = list(test = moran_nbees_resid, scale = "1000m"),
  "\\emph{B. impatiens} abundance (500m)" = list(test = moran_impatiens_resid_500, scale = "500m"),
  "\\emph{B. impatiens} abundance (1000m)" = list(test = moran_impatiens_resid, scale = "1000m"),
  "\\emph{Crithidia} prevalence -- native \\emph{Bombus} (500m)" = list(test = moran_crith_resid_nat_500, scale = "500m"),
  "\\emph{Crithidia} prevalence -- native \\emph{Bombus} (1000m)" = list(test = moran_crith_resid_nat, scale = "1000m"),
  "\\emph{Apicystis} prevalence -- native \\emph{Bombus} (500m)" = list(test = moran_api_resid_nat_500, scale = "500m"),
  "\\emph{Apicystis} prevalence -- native \\emph{Bombus} (1000m)" = list(test = moran_api_resid_nat, scale = "1000m"),
  "\\emph{Vairimorpha} prevalence -- native \\emph{Bombus} (500m)" = list(test = moran_nos_resid_nat_500, scale = "500m"),
  "\\emph{Vairimorpha} prevalence -- native \\emph{Bombus} (1000m)" = list(test = moran_nos_resid_nat, scale = "1000m"),
  "\\emph{Crithidia} prevalence -- \\emph{B. impatiens} (500m)" = list(test = moran_crith_resid_imp_500, scale = "500m"),
  "\\emph{Crithidia} prevalence -- \\emph{B. impatiens} (1000m)" = list(test = moran_crith_resid_imp, scale = "1000m"),
  "\\emph{Apicystis} prevalence -- \\emph{B. impatiens} (500m)" = list(test = moran_api_resid_imp_500, scale = "500m"),
  "\\emph{Apicystis} prevalence -- \\emph{B. impatiens} (1000m)" = list(test = moran_api_resid_imp, scale = "1000m")
)

morans_to_latex(results, file = "docs/tables/morans_table.tex",
                caption = "Results of Moran's I test at two spatial scales for all model residuals. We tested for positive spatial autocorrelation in model residuals, i.e. whether residuals of spatially proximate points were more similar than expected under spatial randomness. All p-values are $\\geq$ 0.05, indicating a lack of positive spatial autocorrelation, except in the case of \\emph{B. impatiens} abundance at the 1000m spatial scale. In this case, because Moran's I < 0.1, the degree of spatial autocorrelation is likely too small to hold much biological relevance.  E[I]: expectation of Moran's I.")

