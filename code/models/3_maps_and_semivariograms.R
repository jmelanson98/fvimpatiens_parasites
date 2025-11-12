## Create map grids of parasite intensity; plot semivariograms of raw data and model residuals

setwd('/Users/jenna1/Documents/UBC/bombus_project/fvimpatiens_parasites')
rm(list=ls())
source("code/src/ggplotThemes.R")
source("code/src/init.R")
source("code/src/misc.R")
source("code/src/makeMapGrids.R")
source("code/src/getPhyloMatrix.R")

# load r data
load(file="saved/AllModels_fv.Rdata")


## ***********************************************************************
## Plot distributions of landscape metrics
## ***********************************************************************
transects = orig.spec %>% 
  group_by(sample_pt) %>%
  sample_n(1)
transect_metrics = transects[,colnames(transects) %in% c("sample_pt", "site", 
                                                         "prop_blueberry_500", 
                                                         "prop_edge_500", 
                                                         "landscape_shdi_500", 
                                                         "prop_blueberry_1000", 
                                                         "prop_edge_1000", 
                                                         "landscape_shdi_1000")]
transect_metrics_long = transect_metrics %>% 
  pivot_longer(-c("sample_pt", "site"),
               names_to = "metric",
               values_to = "values")

split_parts = strsplit(transect_metrics_long$metric, "_")

transect_metrics_long$prefix = sapply(split_parts, function(x) paste(head(x, -1), collapse = "_"))
transect_metrics_long$scale = as.numeric(sapply(split_parts, tail, 1))

plot_list = list()
rows = c("landscape_shdi", "prop_blueberry", "prop_edge")
columns = c(500, 1000)
order = expand.grid(rows, columns)
count = 1

ylims <- list(
  prop_blueberry = c(0, 1),
  prop_edge      = c(0, 0.3),
  landscape_shdi = c(0, 2.5)
)

for (r in rows){
  for (c in columns){
    if(c == 500){
      if(r == "prop_edge"){
        ylab = "Edge density"
      } else if(r == "prop_blueberry"){
        ylab = "Proportion blueberry"
      } else if(r== "landscape_shdi"){
        ylab = "Shannon diversity"
      }
    } else{ylab = ""}
    
    if(r == "prop_edge"){
      xlab = "Site"
    } else{xlab = ""}
    
    plot_list[[count]] = ggplot(transect_metrics_long[transect_metrics_long$prefix == r &
                                                        transect_metrics_long$scale == c,], 
                                aes(x = site, y = values)) +
      geom_violin() +
      ylim(ylims[[r]]) +
      ylab(ylab) +
      xlab(xlab) +
      theme_minimal()
    count = count + 1
  }
}

header1 = ggplot() + theme_void() + ggtitle("500 meter buffer") + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
header2 = ggplot() + theme_void() + ggtitle("1000 meter buffer") + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
blank = nullGrob()

final = grid.arrange(blank, header1, blank, header2, 
  blank, plot_list[[1]], blank, plot_list[[2]],
  blank, blank, blank, blank,
  blank, plot_list[[3]], blank, plot_list[[4]],
  blank, blank, blank, blank,
  blank, plot_list[[5]], blank, plot_list[[6]],
  ncol = 4,
  heights = c(0.15, 1, 0.1, 1, 0.1, 1), widths = c(0.1, 1, 0.1, 1))
site_grid = ggdraw() +
  draw_plot(final, 0.01, 0, 1, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)"), 
                  x = c(0, 0.52, 0, 0.52, 0, 0.52), 
                  y = c(0.98, 0.98, 0.66, 0.66, 0.33, 0.33))

ggsave("figures/manuscript_figures/landscape_metric_dist.jpg", site_grid,
       units = "px", height = 1500, width = 2500)
  


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
fvimp_subpar = filter(fvimp_brmsdf, apidae ==1)
studyspecies = c("Bombus_mixtus", "Bombus_flavifrons", "Bombus_rufocinctus", "Bombus_californicus", 
                 "Bombus_impatiens", "Bombus_vosnesenskii", "Bombus_sitkensis", "Bombus_melanopygus", 
                 "Bombus_nevadensis", "Bombus_medius")
studycov = getPhyloMatrix(studyspecies)

runSVmodels = function(fvimp_sub,
                       fvimp_subpar,
                       studycov){
#calculate residuals for bombus richness (no predictors)
brich <- brm(
  bf(bombus_richness | trials(9) ~ 1),
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
fvimp_sub$brich_resid_nopred = brich_resid$Estimate

#calculate residuals for bombus richness (with predictors)
brich <- brm(
  bf(bombus_richness | trials(9) ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
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

#calculate residuals for wild bee abundance (no predictors)
babun <- brm(
  bf(native_bee_abundance ~ 1),
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
fvimp_sub$babun_resid_nopred = babun_resid$Estimate

#calculate residuals for wild bee abundance (with predictors)
babun <- brm(
  bf(native_bee_abundance ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
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

#calculate residuals for impatiens abundance (no predictors)
iabun <- brm(
  bf(impatiens_abundance ~ 1),
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
fvimp_sub$iabun_resid_nopred = iabun_resid$Estimate

#calculate residuals for impatiens abundance (with predictors)
iabun <- brm(
  bf(impatiens_abundance ~ julian_date + I(julian_date^2) + 
       floral_abundance + floral_diversity +
       prop_blueberry + prop_edge + landscape_shdi +
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

#calculate residuals for crithidia prevalence (no predictors)
hascrith <- brm(
  bf(hascrithidia ~ 1),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
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
       prop_blueberry + prop_edge + landscape_shdi +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
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
fvimp_subpar$crith_resid_pred = crith_resid$Estimate

#calculate residuals for apicystis prevalence (no predictors)
hasapi <- brm(
  bf(apicystis ~ 1 ),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
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
       prop_blueberry + prop_edge + landscape_shdi +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
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
fvimp_subpar$api_resid_pred = api_resid$Estimate

#calculate residuals for nosema prevalence (no predictors)
hasnos <- brm(
  bf(hasnosema ~ 1 ),
  data = fvimp_subpar,
  family = bernoulli(),
  chains = 4,
  thin = 1,
  init = 0,
  cores = 4,
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
       prop_blueberry + prop_edge + landscape_shdi +
       native_bee_abundance + impatiens_abundance + bombus_richness +
       (1|sample_pt) + (1|subsite) + (1|gr(final_id, cov = studycov))),
  data = fvimp_subpar,
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
fvimp_subpar$nos_resid_pred = nos_resid$Estimate


return(list(fvimp_sub, fvimp_subpar))}

#only run this code if models have changed -- otherwise load from .Rdata file
# dfs = runSVmodels(fvimp_sub = fvimp_sub,
#                   fvimp_subpar = fvimp_subpar,
#                   studycov = studycov)
# fvimp_sub_withresid = dfs[[1]]
# fvimp_subpar_withresid = dfs[[2]]

#############################################
### save dataframe 
#############################################
#save(fvimp_sub_withresid, fvimp_subpar_withresid,
#     file="saved/data_with_residuals.Rdata")

load(file="saved/data_with_residuals.Rdata")

#############################################
### make data frames spatially explicit
#############################################

# Convert the data frames into sf objects
fvimp_sub = st_as_sf(fvimp_sub, coords = c("long", "lat"), crs = 4326)
fvimp_subpar = st_as_sf(fvimp_subpar, coords = c("long", "lat"), crs = 4326)

# Transform the coordinate reference system to EPSG:900913
fvimp_sub900913 <- st_transform(fvimp_sub, crs = 3857)
fvimp_subpar900913 <- st_transform(fvimp_subpar, crs = 3857)


############################################
### variograms with no predictors
############################################

#make variogram for bombus richness
v_brich.resid.nopred <- variogram(brich_resid_nopred ~ 1, data = fvimp_sub900913, cutoff = 2000, width = 250)

v_brich.resid.nopred_plot = ggplot(as.data.frame(v_brich.resid.nopred), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(
    title = expression(""),
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
    title = expression(""),
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
    title = expression(""),
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
    title = expression(""),
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
    title = expression(""),
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
    title = expression(""),
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


#### make grid plots
blank <- grid.rect(gp=gpar(col="white"))
brich = ggplot() +
  theme_void() +
  annotate("richtext", x = 0, y = 0, 
           label = "<i>Bombus</i> species richness", 
           fill = NA, label.color = NA,
           size = 5)
babun = ggplot() +
  theme_void() +
  annotate("richtext", x = 0, y = 0, 
           label = "Native <i>Bombus</i> abundance", 
           fill = NA, label.color = NA,
           size = 5)
iabun = ggplot() +
  theme_void() +
  annotate("richtext", x = 0, y = 0, 
           label = "<i>Bombus impatiens</i> abundance", 
           fill = NA, label.color = NA,
           size = 5)
crith = ggplot() +
  theme_void() +
  annotate("richtext", x = 0, y = 0, 
           label = "<i>Crithidia spp.</i> prevalence", 
           fill = NA, label.color = NA,
           size = 5)
api = ggplot() +
  theme_void() +
  annotate("richtext", x = 0, y = 0, 
           label = "<i>Apicystis spp.</i> prevalence", 
           fill = NA, label.color = NA,
           size = 5)
nos = ggplot() +
  theme_void() +
  annotate("richtext", x = 0, y = 0, 
           label = "<i>Vairimorpha spp.</i> prevalence", 
           fill = NA, label.color = NA,
           size = 5)




row1 = grid.arrange(
  v_brich.resid.nopred_plot, blank, v_brich.resid.pred_plot,
  ncol = 3,
  widths = c(4.5, 1, 4.5))
row2 = grid.arrange(
  v_babun.resid.nopred_plot, blank, v_babun.resid.pred_plot,
  ncol = 3,
  widths = c(4.5, 1, 4.5))
row3 = grid.arrange(
  v_iabun.resid.nopred_plot, blank, v_iabun.resid.pred_plot,
  ncol = 3,
  widths = c(4.5, 1, 4.5))
row4 = grid.arrange(
  v_crith.resid.nopred_plot, blank, v_crith.resid.pred_plot,
  ncol = 3,
  widths = c(4.5, 1, 4.5))
row5 = grid.arrange(
  v_api.resid.nopred_plot, blank, v_api.resid.pred_plot,
  ncol = 3,
  widths = c(4.5, 1, 4.5))
row6 = grid.arrange(
  v_nos.resid.nopred_plot, blank, v_nos.resid.pred_plot,
  ncol = 3,
  widths = c(4.5, 1, 4.5))

variogramgrid = grid.arrange(
  blank, brich, row1, 
  blank, babun, row2, 
  blank, iabun, row3, 
  blank, crith, row4, 
  blank, api, row5, 
  blank, nos, row6,
  ncol = 1,
  heights = c(rep(c(0.5, 0.5, 2), 6))
) 
  
#add labels
variogramgrid <- ggdraw() +
  draw_plot(variogramgrid, 0.07, 0, 0.93, 1) +
  draw_plot_label(c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(G)", "(H)", "(I)", "(J)", "(K)", "(L)"),
                  x = c(0, 0.51, 0, 0.51, 0, 0.51, 0, 0.51, 0, 0.51, 0, 0.51), 
                  y = c(0.95, 0.95, 0.78, 0.78, 0.61, 0.61, 0.45, 0.45, 0.28, 0.28, 0.11, 0.11))
#export and save
ggsave(filename = "figures/manuscript_figures/variograms.jpg", 
       plot = variogramgrid, 
       width = 2000, 
       height = 3000, 
       units = "px")



################################
### Calculate Moran's I
################################
library(sf)
library(spdep)

# extract coordinates and get neighbors within 1000 units
coords <- st_coordinates(fvimp_sub900913) 
nb <- dnearneigh(coords, d1 = 0, d2 = 1000)

coords_par <- st_coordinates(fvimp_subpar900913) 
nb_par <- dnearneigh(coords_par, d1 = 0, d2 = 1000)

# convert to spatial weights
#ignore groups that don't have neighbors within specified distance
lw = nb2listw(nb, style = "W", zero.policy = TRUE)
lwpar = nb2listw(nb_par, style = "W", zero.policy = TRUE)

# Perform Moran's I (on raw data)
moran_nbees_raw <- moran.test(fvimp_sub900913$native_bee_abundance, lw)
moran_richness_raw <- moran.test(fvimp_sub900913$bombus_richness, lw)
moran_impatiens_raw <- moran.test(fvimp_sub900913$impatiens_abundance, lw)
moran_crith_raw <- moran.test(fvimp_subpar900913$hascrithidia, lwpar)
moran_api_raw <- moran.test(fvimp_subpar900913$apicystis, lwpar)
moran_nos_raw <- moran.test(fvimp_subpar900913$hasnosema, lwpar)


#okay so I think there is spatial autocorrelation on the raw values which makes sense
#however this spatial autocorrelation is quite weak (I < 0.1 for all relationships)

#test spatial autocorrelation on the model residuals
moran_nbees_resid <- moran.test(fvimp_sub900913$babun_resid_pred, lw)
moran_richness_resid <- moran.test(fvimp_sub900913$brich_resid_pred, lw)
moran_impatiens_resid <- moran.test(fvimp_sub900913$iabun_resid_pred, lw)
moran_crith_resid <- moran.test(fvimp_subpar900913$crith_resid_pred, lwpar)
moran_api_resid <- moran.test(fvimp_subpar900913$api_resid_pred, lwpar)
moran_nos_resid <- moran.test(fvimp_subpar900913$nos_resid_pred, lwpar)



#  check for residuals, with 500 m buffer this time
coords <- st_coordinates(fvimp_sub900913) 
nb <- dnearneigh(coords, d1 = 0, d2 = 500)

coords_par <- st_coordinates(fvimp_subpar900913) 
nb_par <- dnearneigh(coords_par, d1 = 0, d2 = 500)

# convert to spatial weights
#ignore groups that don't have neighbors within specified distance
lw = nb2listw(nb, style = "W", zero.policy = TRUE)
lwpar = nb2listw(nb_par, style = "W", zero.policy = TRUE)

moran_nbees_resid <- moran.test(fvimp_sub900913$babun_resid_pred, lw)
moran_richness_resid <- moran.test(fvimp_sub900913$brich_resid_pred, lw)
moran_impatiens_resid <- moran.test(fvimp_sub900913$iabun_resid_pred, lw)
moran_crith_resid <- moran.test(fvimp_subpar900913$crith_resid_pred, lwpar)
moran_api_resid <- moran.test(fvimp_subpar900913$api_resid_pred, lwpar)
moran_nos_resid <- moran.test(fvimp_subpar900913$nos_resid_pred, lwpar)


plot(lwpar, coords = coords, add=T, col = "red")
