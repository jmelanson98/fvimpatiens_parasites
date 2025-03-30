prepSpecimenTable <- function(specimen.data, parasite.scores){
  ## add parasite scores to the specimen dataframe
  ## clean unwanted individuals
  
  #join dataframes to get specimenData and parasite scores into a single dataframe
  parasiteScores = mutate(parasiteScores, across(-c(barcode_id), as.numeric))
  allspecimens = left_join(specimenData, parasiteScores, by = "barcode_id")
  
  #clean up: remove males / bees without final_id
  allspecimens %>% filter(final_id != "") %>%
    filter(final_id != "B. huntii") -> cleanedspecimens
  cleanedspecimens = cleanedspecimens[str_detect(cleanedspecimens$notes, "released|male|lost", negate = T),]
  
  return(cleanedspecimens)
  
}



# internal data needs to be read
landscape1 <- raster(paste0(here(),"/FValley_lc_1res.tif"))

#load as shp file vs csv
fv_points <- read_sf(paste0(here(),"/fvbombus/fvbombus_points.shp"))
landcover = read_csv(paste0(here(), "/landcover.csv"))

#set and change CRS
crs(landscape1) <- 900913
st_crs(fv_points) <- 900913

#summary of landscape1
print(landscape1)
mapview(landscape1)

#calculate landscape metrics at 500m buffer
lsm <- list_lsm()

#area of each class in each buffer zone
res1_lsm_c_ca <- sample_lsm(landscape1, 
                            y = fv_points, 
                            plot_id = fv_points$site_id, 
                            shape = "circle",
                            size = 500,
                            return_raster = TRUE,
                            what = "lsm_c_ca")

#set IDs
res1_lsm_c_ca = merge(res1_lsm_c_ca, landcover[,c("lc_type", "class")], by = "class")
colnames(res1_lsm_c_ca)[colnames(res1_lsm_c_ca) == "plot_id"] <- "sample_pt"


#shannon diversity
res1_lsm_l_shdi <- sample_lsm(landscape1, 
                              y = fv_points, 
                              plot_id = fv_points$site_id, 
                              shape = "circle",
                              size = 500,
                              return_raster = TRUE,
                              what = "lsm_l_shdi")
colnames(res1_lsm_l_shdi)[colnames(res1_lsm_l_shdi) == "plot_id"] <- "sample_pt"
colnames(res1_lsm_l_shdi)[colnames(res1_lsm_l_shdi) == "value"] <- "landscape_shdi"
