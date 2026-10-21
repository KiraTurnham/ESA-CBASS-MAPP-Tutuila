rm(list=ls())
library(here)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scatterpie)

# 1. Source prep script and load shapefile
source(here("Scripts", "00_visualization_prep.R"))
tutuila_shape <- st_read(here("Data", "Tut_shapefiles", "TUT.shp"))

# Ensure shapefile uses WGS84 (lat/lon in degrees)
tutuila_shape <- st_transform(tutuila_shape, crs = 4326)

# 2. Prepare pie chart data
pie_data <- meta_filteredK %>%
  # Force Longitude to be negative for Western Hemisphere (-170.x)
  mutate(Longitude = -abs(as.numeric(Longitude)),
         Latitude  = as.numeric(Latitude)) %>%
  group_by(Site, Longitude, Latitude, Species) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(
    names_from  = Species, 
    values_from = N, 
    values_fill = 0
  )

spp_cols <- names(species_colors) # c("AABR", "AGLO", "AHYA", "ICRA")

# Get exact spatial bounds of Tutuila
bbox <- st_bbox(tutuila_shape)

# 3. Render Spatial Map
ggplot() +
  # Island shapefile base
  geom_sf(data = tutuila_shape, fill = "grey90", color = "grey40", linewidth = 0.3) +
  
  # Pie charts (r = 0.0015 creates small, distinct pies at ~150m scale)
  geom_scatterpie(
    data      = pie_data,
    aes(x = Longitude, y = Latitude, group = Site),
    cols      = spp_cols,
    r         = 0.001,     # Adjust between 0.001 and 0.003 if needed
    color     = "black",
    linewidth = 0.2
  ) +
  
  # Site labels placed slightly above each pie chart
  geom_text(
    data     = pie_data,
    aes(x = Longitude, y = Latitude - 0.01, label = Site),
    fontface = "bold",
    size     = 5
  ) +
  
  # Color palette & italicized species names in legend
  scale_fill_manual(
    values = species_colors, 
    labels = species_labels,
    name   = species_title
  ) +
  
  # Lock coordinate limits to Tutuila bounding box
  coord_sf(
    xlim = c(bbox["xmin"] - 0.01, bbox["xmax"] + 0.01),
    ylim = c(bbox["ymin"] - 0.01, bbox["ymax"] + 0.01),
    expand = FALSE
  ) +
  theme_coral_pub()+
  theme_coral_pub() +
  theme(
    
    # Place legend inside lower-right corner
    legend.position        = "inside",
    legend.position.inside = c(0.82, 0.18),
    legend.justification   = c(0.5, 0.5),
    
    # Add optional background box & margin padding around legend
    legend.background      = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.margin          = ggplot2::margin(t = 5, r = 8, b = 5, l = 8, unit = "pt") # Namespaced fix
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  )+
  plot_annotation(tag_levels = list(c("A")))
  
#####################################################