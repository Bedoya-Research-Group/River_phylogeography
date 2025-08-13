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


##Panama and Colombia and all points

library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(readr)
library(dplyr)


countries <- ne_countries(scale = "medium", returnclass = "sf") |>
  filter(admin %in% c("Panama", "Colombia")) |>
  st_transform(4326)   # keep in WGS84 lon/lat



points_raw <- read_csv("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/all_samples_coords.csv", show_col_types = FALSE)

points_sf <- points_raw |>
  rename(lon = !!rlang::sym(names(points_raw)[grep("^lon$", names(points_raw), ignore.case=TRUE)]),
         lat = !!rlang::sym(names(points_raw)[grep("^lat$", names(points_raw), ignore.case=TRUE)])) |>
  st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE)


bb_all <- st_bbox(st_union(st_geometry(countries), st_geometry(points_sf)))
pad <- 0.5  # degrees of padding on each side
xlim <- c(bb_all["xmin"] - pad, bb_all["xmax"] + pad)
ylim <- c(bb_all["ymin"] - pad, bb_all["ymax"] + pad)


 xlim <- c(-84, -71)  # longitudes covering Panama->Colombia
 ylim <- c(5, 13)    # latitudes


base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  )

#Plot

aes_color <- if ("species" %in% names(points_raw)) aes(color = species) else NULL

p <- ggplot() +
  geom_sf(data = countries, fill = "grey92", color = "grey35", linewidth = 0.4) +
  geom_point(data = points_sf, mapping = aes(x = lon, y = lat), size = 2.5, alpha = 0.9, inherit.aes = FALSE) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7, line_width = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.6, "cm"), pad_y = unit(1.0, "cm"),
                         style = north_arrow_fancy_orienteering) +
  labs(title = "Sampling Map: Panama & Colombia") +
  base_theme

if (!is.null(aes_color)) {
  p <- ggplot() +
    geom_sf(data = countries, fill = "grey92", color = "grey35", linewidth = 0.4) +
    geom_point(data = points_sf, mapping = aes(x = lon, y = lat, color = species),
               size = 1.5, alpha = 0.95, inherit.aes = FALSE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7, line_width = 0.3) +
    annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.6, "cm"), pad_y = unit(1.0, "cm"),
                           style = north_arrow_fancy_orienteering) +
    labs(title = "Sampling Map: Panama & Colombia", color = "Species") +
    base_theme
}

print(p)


ggsave("panama_colombia_sampling_map.png", p, width = 7.5, height = 6, dpi = 300)
