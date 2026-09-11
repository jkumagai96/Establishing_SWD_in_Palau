# Date: June 23rd 2026
# Purpose: Create a map of seagrass wasting disease pathogenicity studies around the globe
# Seagrass in Palau

##### Load packages ############################################################
library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(geodata)
library(ggrepel)
library(ggnewscale)

##### Load data ################################################################
temp <- geodata::worldclim_global(var = "tavg", res = 10, path = tempdir())
studies <- read.csv("Data/Locations_SWD_Global.csv")
seagrass <- read_sf("Data/014_001_WCMC013-014_SeagrassPtPy2021_v7_1/01_Data/WCMC_013_014_SeagrassesPt_v7_1.shp")

# Convert points to sf object
points <- st_as_sf(studies, coords = c("Long", "Lat"), crs = 4326)
my_study <- points %>% filter(This_Study == 1)

##### Making Map ###############################################################
## Setting up map
crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
graticules <- st_graticule(lat = c(-89.99, -66.5, -23.5, 0, 23.5, 66.5, 89.99), lon = c(-180, 180), crs = 4326) 

map_with_temp <- ggplot() +
  geom_sf(data = graticules, color = "gray50", linewidth = 0.3, linetype = "dashed") +
  geom_spatraster(data = temp, aes(fill = wc2.1_10m_tavg_01), na.rm = TRUE, maxcell = 5e7) +
  scale_fill_viridis_c(name = "Temp (degrees C)", na.value = "transparent") +
  new_scale_fill() +  # resets fill scale for layers below
  geom_sf(data = points, shape = 21, color = "black", aes(fill = Year_sampled), size = 3) +
  scale_fill_gradient(name = "Year Sampled", low = "black", high = "orange") +
  new_scale_fill() +
  geom_sf(data = my_study, shape = 21, color = "black", fill = "red", size = 5) +
  geom_label_repel(stat = "sf_coordinates",
                   data = points,
                   aes(label = Last.name.and.Year, geometry = geometry),
                   size = 3, 
                   fill = "#FFFFFF80",
                   label.size = NA,
                   fontface = "bold",
                   max.overlaps = 15) +
  coord_sf(crs = crs) +
  theme(panel.background = element_rect(fill = NA), 
        legend.position = "bottom") +
  scale_y_continuous(breaks = c(-89.99, -66.5, -23.5, 0, 23.5, 66.5, 89.99),
                     labels = c("90°S", "66.5°S", "23.5°S", "0°", "23.5°N", "66.5°N", "90°N")) +
  scale_x_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL)
map_with_temp

map_with_temp_grey <- ggplot() +
  geom_sf(data = graticules, color = "gray50", linewidth = 0.3, linetype = "dashed") +
  geom_spatraster(data = temp, aes(fill = wc2.1_10m_tavg_01), na.rm = TRUE, maxcell = 5e7) +
  scale_fill_gradient(name = "Temp (degrees C)", high = "grey90", low = "grey20", na.value = "transparent") +
  new_scale_fill() +  # resets fill scale for layers below
  geom_sf(data = points, shape = 21, color = "black", aes(fill = Year_sampled), size = 3) +
  scale_fill_gradient(name = "Year Sampled", low = "#063b00", high = "#14E505") +
  new_scale_fill() +
  geom_sf(data = my_study, shape = 21, color = "black", fill = "yellow", size = 5) +
  geom_label_repel(stat = "sf_coordinates",
                   data = points,
                   aes(label = Last.name.and.Year, geometry = geometry),
                   size = 3, 
                   fill = "#FFFFFF80",
                   label.size = NA,
                   fontface = "bold",
                   max.overlaps = 15) +
  coord_sf(crs = crs) +
  theme(panel.background = element_rect(fill = NA), 
        legend.position = "bottom") +
  scale_y_continuous(breaks = c(-89.99, -66.5, -23.5, 0, 23.5, 66.5, 89.99),
                     labels = c("90°S", "66.5°S", "23.5°S", "0°", "23.5°N", "66.5°N", "90°N")) +
  scale_x_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL)
map_with_temp_grey

map_with_temp_seagrass <- ggplot() +
  geom_sf(data = graticules, color = "gray50", linewidth = 0.3, linetype = "dashed") +
  geom_sf(data = seagrass, color = "#389f0a") +
  new_scale_fill() +
  geom_spatraster(data = temp, aes(fill = wc2.1_10m_tavg_01), na.rm = TRUE, maxcell = 5e7) +
  scale_fill_gradient(name = "Temp (degrees C)", high = "grey90", low = "grey20", na.value = "transparent") +
  new_scale_fill() +  # resets fill scale for layers below
  geom_sf(data = points, shape = 21, color = "black", aes(fill = Year_sampled), size = 3) +
  scale_fill_gradient(name = "Year Sampled", low = "#fff4df", high = "#ff6700") +
  new_scale_fill() +
  geom_sf(data = my_study, shape = 21, color = "black", fill = "yellow", size = 5) +
  geom_label_repel(stat = "sf_coordinates",
                   data = points,
                   aes(label = Last.name.and.Year, geometry = geometry),
                   size = 3, 
                   fill = "#FFFFFF80",
                   label.size = NA,
                   fontface = "bold",
                   max.overlaps = 15) +
  coord_sf(crs = crs) +
  theme(panel.background = element_rect(fill = NA), 
        legend.position = "bottom") +
  scale_y_continuous(breaks = c(-89.99, -66.5, -23.5, 0, 23.5, 66.5, 89.99),
                     labels = c("90°S", "66.5°S", "23.5°S", "0°", "23.5°N", "66.5°N", "90°N")) +
  scale_x_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL)
map_with_temp_seagrass

##### Export ###################################################################
# Map with temp (color)
png("Figures/SWD_pathogencity_Map.png", width = 9, height = 6, units = "in", res = 600)
map_with_temp
dev.off()

# Map with temp (grey)

png("Figures/SWD_pathogencity_Map_greyscale.png", width = 9, height = 6, units = "in", res = 600)
map_with_temp_grey
dev.off()

# Map with seagrass

png("Figures/SWD_pathogencity_Map_seagrass.png", width = 9, height = 6, units = "in", res = 600)
map_with_temp_seagrass
dev.off()
