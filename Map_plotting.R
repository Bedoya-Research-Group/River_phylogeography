library(sf)
library(ggplot2)
library(ggspatial)

panama <- st_read("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/pa_shp/pa.shp") #Downloaded from https://data.humdata.org/dataset/hotosm_pan_waterways
rivers <- st_read("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/hotosm_pan_waterways_lines_shp/hotosm_pan_waterways_lines_shp.shp") #downloaded from https://stridata-si.opendata.arcgis.com/datasets/d8f99ed58d88463c87577a9e329413a4_6/explore?location=8.547599%2C-80.341862%2C9.00

#Making sure both layers match
rivers <- st_transform(rivers, st_crs(panama))

# Clipping rivers to Panama bounding box
bbox <- st_bbox(panama)
rivers_panama <- st_crop(rivers, bbox)

# Plot
ggplot() +
  geom_sf(data = panama, fill = "gray90", color = "black") +
  geom_sf(data = rivers_panama, color = "blue") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_fancy_orienteering)
