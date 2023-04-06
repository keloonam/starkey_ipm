# Code to create publication study area map of Starkey
# Kenneth Loonam
# February 2023

#Environment====================================================================

require(ggplot2); require(rnaturalearth); require(maps); require(elevatr)
require(rnaturalearthdata); require(ggspatial); require(dplyr); require(sf);
require(terra); require(ggnewscale); require(whitebox); require(purrr)
require(cowplot)

#Data prep variables============================================================

study_area_file <- "data//Spatial//Starkey_main"
sa_rect_buffer <- c(x = 0.1, y = 0.1) # size of area around SA to show
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
  mutate(place = ID) %>% # repeat all those steps for the state data
  rename(geometry = geom) %>%
  select(place, geometry)

map_dat <- rbind(states, world, starkey) # pull all your boundary data together

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

# pull the elevation data
sa_elev <- get_elev_raster(sa_rect_sf, z = elevation_zoom_level) %>%
  rast() %>%
  mask(vect(sa_rect_sf))
# make it a data frame and give the data column a nice name
sa_ele_df <- as.data.frame(sa_elev, xy = T)
names(sa_ele_df)[3] <- "ele"

# pull slope and aspect data
s1 <- terrain(sa_elev, "slope", unit = "radians")
aspect <- terrain(sa_elev, "aspect", unit = "radians")

# calculate the "shadows" and put it into a data frame
hill_single <- shade(s1, aspect, 
                     angle = 30, 
                     direction = 300,
                     normalize= TRUE)
hill_df <- as.data.frame(hill_single, xy = T)

#Plotting=======================================================================

sa_location <- ggplot(data = map_dat) +
  theme_bw() +
  geom_sf() +
  # Set map bounds
  coord_sf(xlim = c(-110, -125.5), ylim = c(42,53)) +
  theme(
    panel.grid.major = element_line(
      color = "gray",
      linetype = "dashed",
      size = 0.5)) +
  # Draw the inset rectangle
  geom_rect(
    aes( # Adjust exact position of rectangle to make the inset border look good
      xmin = sa_rect[1] - 0.0008, 
      ymin = sa_rect[2] - 0.002, 
      xmax = sa_rect[3] + 0.002, 
      ymax = sa_rect[4] - 0.000), 
    fill = NA, 
    color = "black",
    size = 1.0) +
  annotation_scale( # add a scale bar and adjust its position
    location = "br", 
    text_cex = .85,
    pad_y = unit(0.09, "cm"),
    pad_x = unit(0.32, "cm")
    ) 

sa_location %>% # Start with the big map
  ggdraw() +
  draw_plot( # Draw another map on it
    {
      sa_location +
        geom_raster( # add the hill shade layer
          data = hill_df,
          aes(x, y, fill = lyr1)) + 
        scale_fill_distiller(palette = "Greys") +
        # redraw study area shape so it is on top of the hill shade layer
        geom_sf(data = starkey, alpha = 0, color = "black", size = 1.0) +
        coord_sf( # zoom in on a section of the big map for the inset map
          # adjust the exact zoom location so that it matches the rectangle
          xlim = c(sa_rect[1] + 0.0115, sa_rect[3] - 0.0115), 
          ylim = c(sa_rect[2] + 0.0125, sa_rect[4] - 0.0125)) +
        theme_void() +
        theme(legend.position = "none") +
        annotation_scale( # add the scale
          location = "br", 
          text_cex = .85) +
        annotation_north_arrow( # add the north arrow
            location = "tr",
            style = north_arrow_orienteering(),
            height = unit(1.25, "cm"),
            width = unit(1, "cm"))
          
    }, # adjust the size and position of the inset map
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

