# Code to create publication study area map of Starkey
# Kenneth Loonam
# February 2023

#Environment====================================================================

require(ggplot2); require(rnaturalearth); require(maps); require(elevatr)
require(rnaturalearthdata); require(ggspatial); require(dplyr); require(sf);
require(terra); require(ggnewscale); require(whitebox); require(purrr)
require(cowplot); require(tidyterra); require(giscoR); require(ggblend)
require(ggsci); require(grDevices); require(grid); require(raster)

#Data prep variables============================================================

study_area_file <- "data//Spatial//Starkey_main"
sa_rect_buffer <- c(x = 0.015, y = 0.015) # size of area around SA to show
elevation_zoom_level <- 12 # resolution for elevation, 1-14, 1 is low detail
projection <- 4326 # What projection do you want to work in


#Data Prep======================================================================

sarb <- sa_rect_buffer # short object names make me happier

starkey <- st_read(study_area_file) %>% # read in
  st_transform(projection) %>% # reproject
  mutate(place = "starkey_main") %>% # rename
  select(place, geometry) # and pare down the study area shapefile

world <- ne_countries(scale = "medium", returnclass = "sf") %>% # pull world dat
  st_transform(projection) %>% # reporoject it
  mutate(place = sovereignt) %>% # get consistent column names
  select(place, geometry) %>% # pare it down
  filter(place != "United States of America") # we get US data from states

states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  st_transform(projection) %>% # reporoject it
  mutate(place = ID) %>% # repeat all those steps for the state data
  rename(geometry = geom) %>%
  select(place, geometry)

map_dat <- rbind(states, world, starkey) # pull all your boundary data together
map_dat2 <- rbind(states, world) %>%
  st_transform("+proj=laea +lat_0=30 +lon_0=-95") 

# Build the study area rectangle for the inset map
sa_rect <- st_bbox(starkey) + c(-sarb[1], -sarb[2], sarb[1], sarb[2])
# turn that rectangle into an sf object
geometry <- st_sfc(st_polygon(x = list(matrix(c(
  sa_rect[1], sa_rect[2],
  sa_rect[1], sa_rect[4],
  sa_rect[3], sa_rect[4],
  sa_rect[3], sa_rect[2],
  sa_rect[1], sa_rect[2]
), ncol = 2, byrow = T
))))
# Make consistent names and define the projection
sa_rect_sf <- st_sf(name = "name", geometry = geometry)
st_crs(sa_rect_sf) <- projection
sa_star <- sa_rect_sf %>%
  st_transform("+proj=laea +lat_0=30 +lon_0=-95") %>%
  # pull(geometry) %>%
  st_centroid()
# pull the elevation data
sa_elev <- get_elev_raster(sa_rect_sf, z = elevation_zoom_level) %>%
  rast() %>%
  mask(vect(sa_rect_sf))
# make it a data frame and give the data column a nice name
sa_ele_df <- as.data.frame(sa_elev, xy = T)
names(sa_ele_df)[3] <- "ele"

# pull slope and aspect data
slope <- terrain(sa_elev, "slope", unit = "radians")
aspect <- terrain(sa_elev, "aspect", unit = "radians")

# calculate the "shadows" and put it into a data frame
hill_single <- shade(slope, aspect, 
                     angle = 30, 
                     direction = 300,
                     normalize= TRUE)
hill1_df <- as.data.frame(hill_single, xy = T)

hill_multi <- map(c(270, 15, 60, 330), function(dir){ 
  shade(slope, aspect, 
        angle = 45, 
        direction = dir,
        normalize= TRUE)}
) %>% rast() %>% sum()
hill4_df <- as.data.frame(hill_multi, xy = T) %>%
  mutate(hillshade = sum)

hills_sf <- st_as_sf(hill4_df, coords = c("x", "y"), crs = projection)

# plot(hill_multi, col = grey(1:100/100))
# plot(hill_single, col = grey(1:100/100))

#Plotting=======================================================================

loc_fig <- ggplot(map_dat2) +
  theme_bw() +
  geom_sf(fill = "#f5f5dc") +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 5)) +
  # Set map bounds
  scale_y_continuous(limits = c(400000, 2800000)) +
  scale_x_continuous(limits = c(-2900000, -600000)) +
  theme(
    panel.grid.major = element_line(
      color = "black",
      linetype = "dashed",
      size = 0.5),
    panel.background = element_rect(fill = "#afeeee")) +
  # Draw the inset rectangle
  geom_point(
    aes(geometry = geometry),
    data = sa_star, 
    shape = 21, 
    size = 5, 
    stat = "sf_coordinates",
    color = "black",
    fill = "#ff0000"
  ) +
  labs(x = NULL, y = NULL)

stky_fig <- ggplot() +
  geom_raster( # add the hill shade layer
    data = hill_multi,
    aes(fill = hillshade),
    show.legend = F) + 
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_raster(
    data = sa_ele_df,
    aes(x, y, fill = ele),
    alpha = 0.7
  ) +
  scale_fill_gradientn(
    colors = hcl.colors(n = 50, palette = "Mako"), 
    limits = c(1100, 1525)
  ) +
  # redraw study area shape so it is on top of the hill shade layer
  geom_sf(data = starkey, alpha = 0, color = "black", size = 2.0) +
  theme_bw() +
  guides(fill = guide_colorbar(title.position = "top")) +
  theme(legend.position = "right") +
  labs(fill = "Elevation (meters)") + 
  annotation_scale( # add the scale
    location = "br", 
    text_cex = .85) +
  annotation_north_arrow( # add the north arrow
    location = "bl",
    style = north_arrow_orienteering(),
    height = unit(1.25, "cm"),
    width = unit(1, "cm")) +
  labs(y = "Latitude", x = "Longitude") +
  theme(
    panel.ontop = T, 
    panel.background = element_blank(),
    panel.grid = element_line(color = "#606060", linetype = "dashed")
  )

stky_fig %>%
  ggdraw() +
  draw_plot(
    loc_fig,
    x = 0.47,
    y = 0.375,
    width = 0.6,
    height = 0.6
  )

# This save file will look different than the plot viewer in RStudio!!!
# I like ggsave for publication figures because it is easy and consistent to 
# remake figures and makes saving at high resolutions easy. The downside is that 
# it adds a step when you are making the fine adjustments for the figure
# (saving the file and viewing the result)
ggsave( # save at high resolution
  "figures//study_area_figure.png", 
  width = 4, 
  height = 4, 
  units = "in", 
  dpi = 500)

