
# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)
library(sf)
library(fasterize)

# Studyregion and grid ----------------------------------------------------

# The following code needs disturbance maps downloaded from 10.5281/zenodo.3924381 (both version 1.0 and 1.1) and saved into a folder "disturbances"
path_dist_maps <- "disturbances/version1.1/"
path_forest <- "disturbances/version1.0/"

countries <- list.files(path_dist_maps)
countries <- countries[!countries %in% c("andorra", "liechtenstein", "luxembourg")]

# Find bounding box
min_lng <- c()
min_lat <- c()
max_lng <- c()
max_lat <- c()

for (i in 1:length(countries)) {
  
  cntr <- countries[i]
  
  print(cntr)
  
  disturbance <- raster(paste0(path_dist_maps, cntr, "/disturbance_year_", cntr, ".tif"))
  ext <- as(extent(disturbance), 'SpatialPolygons')
  ext <- st_as_sf(ext)
  st_crs(ext) <- projection(disturbance)
  ext_latlng <- ext %>% st_transform(., crs = "epsg:4326")
  min_lng <- c(min_lng ,min(st_coordinates(ext_latlng)[,1]))
  min_lat <- c(min_lat ,min(st_coordinates(ext_latlng)[,2]))
  max_lng <- c(max_lng ,max(st_coordinates(ext_latlng)[,1]))
  max_lat <- c(max_lat ,max(st_coordinates(ext_latlng)[,2]))
  
  }

# Corner-points for grid (currently ERA5 0.5 degree)
min_lng <- floor(min(min_lng)) - 0.25
min_lat <- floor(min(min_lat)) - 0.25
max_lng <- ceiling(max(max_lng)) + 0.25
max_lat <- ceiling(max(max_lat)) + 0.25

# Create grid
grid <- raster(xmn = min_lng, 
               xmx = max_lng, 
               ymn = min_lat, 
               ymx = max_lat, 
               nrows = (max_lat - min_lat) / 0.5,
               ncols = (max_lng - min_lng) / 0.5,
               crs = "+proj=longlat +datum=WGS84")

values(grid) <- 1:ncell(grid)

writeRaster(grid, "grid.tif", overwrite = TRUE)

# Create sf object of grid and transform to epsg:3035 for intersection
grid_sf <- grid %>%
  rasterToPolygons() %>% 
  st_as_sf(grid) %>%
  rename(gridindex = layer)

grid_sf_epsg3035 <- grid_sf %>%
  st_transform(., "epsg:3035")

# Aggregate to climate grid -----------------------------------------------

for (i in 1:length(countries)) {
  
  cntr <- countries[i]
  
  dir.create("temp", showWarnings = FALSE)
  
  if (!file.exists(paste0("temp/disturbances_aggregated_to_grid_", cntr, ".csv"))) {
  
    print(cntr)
    
    disturbance <- raster(paste0(path_dist_maps, cntr, "/disturbance_year_", cntr, ".tif"))
    forest <- raster(paste0(path_forest, cntr, "/prediction_forestcover_", cntr, ".tif"))
    
    ext <- as(extent(disturbance), 'SpatialPolygons')
    ext <- st_as_sf(ext)
    st_crs(ext) <- st_crs("epsg:3035")
    
    grid_sel <- st_intersection(st_as_sf(grid_sf_epsg3035), st_as_sf(ext))
    grid_sel_ras <- fasterize(grid_sel, disturbance, field = "gridindex")
    grid_values <- values(grid_sel_ras)
    
    disturbance <- data.frame(gridindex = grid_values,
                              disturbance = values(disturbance),
                              country = cntr) %>%
      na.omit(.) %>%
      group_by(gridindex, year = disturbance) %>%
      summarize(disturbance_ha = n() * 0.09,
                country = unique(country)) %>%
      ungroup(.)
    
    forest <- data.frame(gridindex = grid_values,
                         forest = values(forest),
                         country = cntr) %>%
      filter(!is.na(forest)) %>%
      group_by(gridindex) %>%
      summarize(forest_ha = sum(forest == 1, na.rm = TRUE) * 0.09,
                land_ha = n() * 0.09,
                country = unique(country)) %>%
      ungroup(.)
    
    write_csv(disturbance, paste0("temp/disturbances_aggregated_to_grid_", cntr, ".csv"))
    write_csv(forest, paste0("temp/forest_aggregated_to_grid_", cntr, ".csv"))
    
    gc()
    
  }
  
}

# Combine everything ------------------------------------------------------

disturbance <- list.files("temp/", pattern = "disturbances_aggregated_to_grid*", full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()

countries_grid <- disturbance %>% 
  group_by(gridindex) %>%
  summarise(country = unique(country)) 

disturbance <- disturbance %>%
  filter(!(country == "norway" & year == 2020))

forest <- list.files("temp/", pattern = "forest_aggregated_to_grid*", full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()

dat <- disturbance %>%
  filter(!is.na(gridindex)) %>%
  group_by(gridindex, year) %>%
  summarize(disturbance_ha = sum(disturbance_ha, na.rm = TRUE)) %>%
  ungroup() %>%
  split(.$gridindex) %>%
  map(~ right_join(., data.frame(gridindex = unique(.$gridindex),
                                 year = 1986:2020), 
                   by = c("gridindex", "year"))) %>%
  map(~ mutate_at(., .vars = vars(disturbance_ha), .funs = function(x) ifelse(is.na(x), 0, x))) %>%
  bind_rows()

# Drop anomaly for Norway in 2020, as there is no data available
dat <- dat %>% 
  left_join(countries_grid, by = "gridindex") %>%
  filter(!(country == "Norway" & year == 2020))

forest <- forest %>%
  group_by(gridindex) %>%
  summarise(forest_ha = sum(forest_ha),
            land_ha = sum(land_ha)) %>%
  mutate(forestcover = forest_ha / land_ha) %>%
  ungroup()

forest_grid <- grid_sf %>%
  right_join(forest)

reference_period <- 2000:2020

dat <- dat %>%
  left_join(forest) %>%
  filter(sum(disturbance_ha) > 35) %>% # Exclude areas with less than 1 ha/yr of disturbances on average
  filter(sum(disturbance_ha[year %in% reference_period]) > length(reference_period)) %>% # Exclude areas with less than 1 ha/yr of disturbances on average within the reference period
  mutate(rate = disturbance_ha / forest_ha) %>%
  group_by(gridindex) %>%
  mutate(anomaly = disturbance_ha / mean(disturbance_ha[year %in% reference_period], na.rm = TRUE) - 1) %>% 
  ungroup()

dat <- dat %>% filter(!is.na(forest_ha))

head(dat)

save(dat, file = "temp/dat.RData")

dat_grid <- grid_sf_epsg3035 %>% right_join(dat, by = "gridindex")

dat_grid <- dat_grid %>%
  mutate(rate_capped = ifelse(rate > 0.05, 0.05, rate),
         anomaly_capped = ifelse(anomaly > 5, 5, anomaly))

ggplot() +
  geom_sf(data = dat_grid %>% filter(year == 2000),
          aes(fill = anomaly_capped * 100), col = NA) +
  scale_fill_gradient2(low = "#2166ac", mid = "#FFFFFF", high = "#b2182b",
                       breaks = c(-100, 0, 100, 200, 300, 400, 500),
                       labels = c("-100%", "0%", "100%", "200%", "300%", "400%", ">500%")) +
  scale_color_gradient2(low = "#2166ac", mid = "#FFFFFF", high = "#b2182b",
                        breaks = c(-100, 0, 100, 200, 300, 400, 500),
                        labels = c("-100%", "0%", "100%", "200%", "300%", "400%", ">500%")) +
  theme_linedraw() +
  theme(panel.spacing = unit(0, "cm"),
        #panel.background = element_rect(fill = "#d1e5f0", color = "black", size = 1.125),
        panel.background = element_rect(fill = "white", color = "black", size = 1.125),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.125, "cm"),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 9, color = "black"),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(col = NULL, fill = NULL, 
       title = "Forest disturbance anomalies")

ggplot() +
  geom_sf(data = dat_grid %>% filter(year == 2019),
          aes(fill = rate_capped * 100), col = NA) +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  theme_linedraw() +
  theme(panel.spacing = unit(0, "cm"),
        #panel.background = element_rect(fill = "#d1e5f0", color = "black", size = 1.125),
        panel.background = element_rect(fill = "white", color = "black", size = 1.125),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.125, "cm"),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 9, color = "black"),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(col = NULL, fill = NULL, 
       title = "Forest disturbance rate")

anomalies <- vector("list", length(1986:2020))
rates <- vector("list", length(1986:2020))

k <- 0

for (y in 1986:2020) {
  
  k <- k + 1
  
  reclass_anomalies <- expand_grid(gridindex = values(grid)) %>%
    left_join(dat %>%
                filter(year == y) %>%
                dplyr::select(gridindex, anomaly))
  
  reclass_rates <- expand_grid(gridindex = values(grid)) %>%
    left_join(dat %>%
                filter(year == y) %>%
                dplyr::select(gridindex, rate))
  
  anomalies[[k]] <- reclassify(grid, as.matrix(reclass_anomalies))
  rates[[k]] <- reclassify(grid, as.matrix(reclass_rates))
  
}

dir.create("results", showWarnings = FALSE)

anomalies <- stack(anomalies)
names(anomalies) <- 1986:2020
writeRaster(anomalies, "results/anomalies.tif", overwrite = TRUE)

rates <- stack(rates)
names(rates) <- 1986:2020
writeRaster(rates, "results/rates.tif", overwrite = TRUE)
