# Code to create publication study area map of Starkey
# Kenneth Loonam
# February 2023

#Environment====================================================================

require(ggplot2); require(rnaturalearth); require(maps); require(elevatr)
require(rnaturalearthdata); require(ggspatial); require(dplyr); require(sf);
require(terra)
theme_set(theme_bw())

#Variables=====================================================================

study_area_file <- "data//Spatial//Starkey_main"
sa_rect_buffer <- c(x = 0.1, y = 0.06)
elevation_zoom_level <- 10

#Data Prep======================================================================

sarb <- sa_rect_buffer

starkey <- st_read(study_area_file) %>%
  st_transform(4326) %>%
  mutate(place = "starkey_main") %>%
  select(place, geometry)

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(4326) %>%
  mutate(place = sovereignt) %>%
  select(place, geometry) %>%
  filter(place != "United States of America")

states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) %>%
  mutate(place = ID) %>%
  rename(geometry = geom) %>%
  select(place, geometry)

map_dat <- rbind(states, world, starkey)

sa_rect <- st_bbox(starkey) + c(-sarb[1], -sarb[2], sarb[1], sarb[2])
geometry <- st_sfc(st_polygon(x = list(matrix(c(
  sa_rect[1], sa_rect[2],
  sa_rect[1], sa_rect[4],
  sa_rect[3], sa_rect[4],
  sa_rect[3], sa_rect[2],
  sa_rect[1], sa_rect[2]
), ncol = 2, byrow = T
))))
sa_rect_sf <- st_sf(name = "name", geometry = geometry)
st_crs(sa_rect_sf) <- 4326

sa_elev <- get_elev_raster(sa_rect_sf, z = 10) %>%
  rast() %>%
  mask(vect(sa_rect_sf))
sa_ele_df <- as.data.frame(sa_elev, xy = T)
names(sa_ele_df)[3] <- "ele"

ggplot()+
  geom_raster(
    data = sa_ele_df,
    aes(x, y, fill = ele)
  ) +
  geom_sf(data = starkey, alpha = 0, color = "black")


#https://dominicroye.github.io/en/2022/hillshade-effects/

#Plotting=======================================================================

ggplot(data = map_dat) +
  geom_sf() +
  coord_sf(xlim = c(-115, -126), ylim = c(42,50)) +
  theme(
    panel.grid.major = element_line(
      color = "gray", 
      linetype = "dashed", 
      size = 0.5)) +
  geom_rect(
    aes(
      xmin = sa_rect[1], 
      ymin = sa_rect[2], 
      xmax = sa_rect[3], 
      ymax = sa_rect[4]), 
    fill = NA, 
    color = "black")

ggplot(data = map_dat) + 
  geom_sf() + 
  coord_sf(
    xlim = c(sa_rect[1], sa_rect[3]), 
    ylim = c(sa_rect[2], sa_rect[4]))
