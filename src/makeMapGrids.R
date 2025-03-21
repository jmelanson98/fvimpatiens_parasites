makeMapGrids <- function(groupedbysite, 
                         #a dataframe which includes columns for: 
                         #site ID (e.g., "W", "SD", etc.)
                         #numbees (sampling effort, e.g., number of specimens screened or number of survey efforts),
                         #site_any (variable of interest, e.g., parasite prev or bee abundance),
                         #lat (latitude),
                         #long (longitude)
                         sampling_effort, #a string describing the relevant sample effort variable, e.g., "Number of specimens"
                         var_of_interest #a strong describing the variable being plotted, e.g., "Parasite Prevalence"
                         ){
  
  
  # Load raster files
  westham_basemap <- raster("data/zipped_clipped_rasters/westham_masked.tif")
  sd_basemap <- raster("data/zipped_clipped_rasters/sd_masked.tif")
  ed_basemap <- raster("data/zipped_clipped_rasters/ed_masked.tif")
  nr_basemap <- raster("data/zipped_clipped_rasters/nr_masked.tif")
  hr_basemap <- raster("data/zipped_clipped_rasters/hr_masked.tif")
  pm_basemap <- raster("data/zipped_clipped_rasters/pm_masked.tif")
  
  #make input dataframe into an sf object
  site_sf <- st_as_sf(groupedbysite, coords = c("long", "lat"), crs = 4326)
  
  #make sure the projections match
  site_sf <- st_transform(site_sf, crs = crs(westham_basemap))
  
  #just bees from the site in question
  westham_sf = filter(site_sf, site == "W")
  sd_sf = filter(site_sf, site == "SD")
  ed_sf = filter(site_sf, site == "ED")
  nr_sf = filter(site_sf, site == "NR")
  hr_sf = filter(site_sf, site == "HR")
  pm_sf = filter(site_sf, site == "PM")
  
  #downsample the rasters because they're too hefty
  westham_downsampled = aggregate(westham_basemap, fact = 10)
  sd_downsampled = aggregate(sd_basemap, fact = 10)
  ed_downsampled = aggregate(ed_basemap, fact = 10)
  nr_downsampled = aggregate(nr_basemap, fact = 10)
  hr_downsampled = aggregate(hr_basemap, fact = 10)
  pm_downsampled = aggregate(pm_basemap, fact = 10)
  
  
  #remove outliers from site_any to make the color scaling better
  # Calculate IQR
  iqr_value <- IQR(groupedbysite$site_any)
  
  # Calculate Q1 and Q3
  q1 <- quantile(groupedbysite$site_any, 0.25)
  q3 <- quantile(groupedbysite$site_any, 0.75)
  
  # Define the lower and upper bounds for non-outlier data
  lower_bound <- q1 - 3 * iqr_value
  upper_bound <- q3 + 3 * iqr_value
  
  # Remove outliers
  values_no_outliers <- groupedbysite$site_any[groupedbysite$site_any >= lower_bound & groupedbysite$site_any <= upper_bound]
  
  #calculate range
  var_range = c(min(values_no_outliers), max(values_no_outliers))
  print(var_range)
  
  #make plots
  westham_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(westham_downsampled, xy = TRUE), aes(x = x, y = y, fill = westham_masked)) +
    #scale_fill_gradient2() +  # Optional: Using the 'viridis' color scale
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    scale_y_continuous(limits = c(6287171, 6292570)) + 
    scale_x_continuous(limits = c(-13711570, -13707050)) +
    # Add GPS points
    geom_sf(data = westham_sf, aes(color = site_any, size = numbees)) +
    scale_color_viridis(option = "C", limits = var_range, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), limits = range(c(1, max(groupedbysite$numbees)))) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", title = "Westham Island") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.4)
  
  
  
  sd_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(sd_downsampled, xy = TRUE), aes(x = x, y = y, fill = sd_masked)) +
    #scale_fill_gradient2() +  # Optional: Using the 'viridis' color scale
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    scale_y_continuous(limits = c(6283496, 6287350)) + 
    scale_x_continuous(limits = c(-13702481, -13697791)) +
    # Add GPS points
    geom_sf(data = sd_sf, aes(color = site_any, size = numbees)) +
    scale_color_viridis(option = "C", limits = var_range, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), limits = range(c(1, max(groupedbysite$numbees)))) +
    theme_minimal() +
    labs(x = "Longitude", 
         y = "Latitude", 
         title = "South Delta") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.3)
  
  ed_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(ed_downsampled, xy = TRUE), aes(x = x, y = y, fill = ed_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +  # White background for NA
    scale_y_continuous(limits = c(6289400, 6293500)) + 
    scale_x_continuous(limits = c(-13688350, -13683650)) +
    # Add GPS points
    geom_sf(data = ed_sf, aes(color = site_any, size = numbees)) +
    scale_color_viridis(option = "C", limits = var_range, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), limits = range(c(1, max(groupedbysite$numbees)))) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", title = "East Delta") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.3)
  
  nr_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(nr_downsampled, xy = TRUE), aes(x = x, y = y, fill = nr_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    scale_y_continuous(limits = c(6286543, 6290900)) + 
    scale_x_continuous(limits = c(-13668600, -13664000)) +
    # Add parasite data points
    geom_sf(data = nr_sf, aes(color = site_any, size = numbees)) +
    scale_color_viridis(option = "C", limits = var_range, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), limits = range(c(1, max(groupedbysite$numbees)))) +
    theme_minimal() +
    labs(x = "Longitude", 
         y = "Latitude", 
         title = "Nicomekl River") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.4)
  
  hr_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(hr_downsampled, xy = TRUE), aes(x = x, y = y, fill = hr_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    scale_y_continuous(limits = c(6299000, 6303300)) + 
    scale_x_continuous(limits = c(-13663050, -13658300)) +
    # Add GPS points
    geom_sf(data = hr_sf, aes(color = site_any, size = numbees)) +
    scale_color_viridis(option = "C", limits = var_range, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), limits = range(c(1, max(groupedbysite$numbees)))) +
    theme_minimal() +
    labs(x = "Longitude", 
         y = "Latitude", 
         title = "Harvie Road") +
    theme(legend.position = "none",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.3)
  
  #plot PM twice to extract relevant legends
  #first, plot with value legend
  pm_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(pm_downsampled, xy = TRUE), aes(x = x, y = y, fill = pm_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    scale_y_continuous(limits = c(6315100, 6319400)) + 
    scale_x_continuous(limits = c(-13656000, -13651250)) +
    # Add GPS points
    geom_sf(data = pm_sf, aes(color = site_any, size = numbees)) +
    scale_color_viridis(option = "C", limits = var_range, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), limits = range(c(1, max(groupedbysite$numbees)))) +
    theme_minimal() +
    labs(color = var_of_interest, 
         size = sampling_effort, 
         x = "Longitude", 
         y = "Latitude", 
         title = "Pitt Meadows") +
    theme(legend.position = "right",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none", size = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.3)
  
  #snag the value legend!
  g <- ggplotGrob(pm_plot)
  value_legend <- g$grobs[[15]]
  
  #now replot with size legend
  pm_plot = ggplot() +
    # Add raster
    geom_tile(data = as.data.frame(pm_downsampled, xy = TRUE), aes(x = x, y = y, fill = pm_masked)) +
    scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = "white") +
    scale_y_continuous(limits = c(6315100, 6319400)) + 
    scale_x_continuous(limits = c(-13656000, -13651250)) +
    # Add GPS points
    geom_sf(data = pm_sf, aes(color = site_any, size = numbees)) +
    scale_color_viridis(option = "C", limits = var_range, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), limits = range(c(1, max(groupedbysite$numbees)))) +
    theme_minimal() +
    labs(color = var_of_interest, 
         size = sampling_effort, 
         x = "Longitude", 
         y = "Latitude", 
         title = "Pitt Meadows") +
    theme(legend.position = "right",
          panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6)
    ) +
    guides(fill = "none", color = "none") +
    # Add map features (scale bar, north arrow, etc.)
    annotation_scale(location = "bl", width_hint = 0.3)
  
  
  #snag the size legend
  g <- ggplotGrob(pm_plot)
  size_legend <- g$grobs[[15]]
  
  # Remove legends from plot
  pm_plot <- pm_plot + theme(legend.position = "none")
  
  #make each plot into a grob
  g1 = ggplotGrob(westham_plot)
  g2 = ggplotGrob(sd_plot)
  g3 = ggplotGrob(ed_plot)
  g4 = ggplotGrob(nr_plot)
  g5 = ggplotGrob(hr_plot)
  g6 = ggplotGrob(pm_plot)
  
  #make spacers
  blank_plot <- ggplot() + 
    theme_void() + 
    theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
  b = ggplotGrob(blank_plot)
  
  #make grid plot of maps
  plot_grid(
    g1, b, g2, b, g3, b, value_legend,
    b, b, b, b, b, b, b,
    g5, b, g4, b, g6, b, size_legend,
    ncol = 7,
    nrow = 3,
    rel_heights = c(1, 0.1, 1),
    rel_widths = c(1, 0.1, 1, 0.1, 1, 0.1, 0.25),
    align = "hv",
    axis = "tb"
  )
  
  
  #final_plot <- ggdraw() +
  #  draw_plot(interactiongrid, 0.015, 0, 1, 1) +
  #  draw_plot_label(c("a", "b", "c"), 
  #                  x = c(0, 0, 0), 
   #                 y = c(0.97, 0.63, 0.3))
  #y = c(0.97, 0.73, 0.48))
  #print(final_plot)
  
  
}