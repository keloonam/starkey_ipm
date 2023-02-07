# Code to create publication study area map of Starkey
# Kenneth Loonam
# February 2023

#Environment====================================================================

require(ggplot2); require(sf); require(rnaturalearth); require(maps)
require(rnaturalearthdata); require(ggspatial)
theme_set(theme_bw())

starkey <- st_read("data//Spatial//Starkey_main")
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))
starkey <- st_transform(starkey, 4326)
names(states) <- names(starkey)
st_geometry(states) <- "geometry"
world <- st_transform(world, 4326)[, c(4,ncol(world))]
world <- world[-which(world$sovereignt == "United States of America"),]
names(world) <- names(starkey)
states <- rbind(starkey, states, world)

ggplot(data = states) +
  geom_sf() +
  coord_sf(xlim = c(-115, -126), ylim = c(42,50)) +
  theme(
    # panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_line(
      color = "gray", 
      linetype = "dashed", 
      size = 0.5)) +
  # annotation_scale(location = "bl", width_hint = 0.5) +
  # annotation_north_arrow(
  #   location = "tr", 
  #   which_north = "T", 
  #   pad_x = unit(0.75, "in"), 
  #   pad_y = unit(0.5, "in"), 
  #   style = north_arrow_fancy_orienteering) +
  geom_rect()
