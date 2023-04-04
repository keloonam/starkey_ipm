# Code to create publication study area map of Starkey
# Kenneth Loonam
# February 2023

#Environment====================================================================

require(ggplot2); require(rnaturalearth); require(maps); require(elevatr)
require(rnaturalearthdata); require(ggspatial); require(dplyr); require(sf);
require(terra); require(ggnewscale); require(whitebox); require(purrr)
require(cowplot)
# theme_set(theme_bw())

#Variables=====================================================================

study_area_file <- "data//Spatial//Starkey_main"
sa_rect_buffer <- c(x = 0.1, y = 0.1)
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

states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
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

sa_elev <- get_elev_raster(sa_rect_sf, z = elevation_zoom_level) %>%
  rast() %>%
  mask(vect(sa_rect_sf))
sa_ele_df <- as.data.frame(sa_elev, xy = T)
names(sa_ele_df)[3] <- "ele"

s1 <- terrain(sa_elev, "slope", unit = "radians")
aspect <- terrain(sa_elev, "aspect", unit = "radians")
hill_single <- shade(s1, aspect, 
                     angle = 30, 
                     direction = 300,
                     normalize= TRUE)
hill_df <- as.data.frame(hill_single, xy = T)

# hillmulti <- map(c(270, 15, 60, 330), function(dir){ 
#   shade(s1, aspect, 
#         angle = 45, 
#         direction = dir,
#         normalize= TRUE)}
# )
# 
# # create a multidimensional raster and reduce it by summing up
# hillmulti <- rast(hillmulti) %>% sum()
# 
# # multidirectional
# plot(hillmulti, col = grey(1:100/100))
# hill_multi_df <- as.data.frame(hillmulti, xy = T)

#https://dominicroye.github.io/en/2022/hillshade-effects/

#Plotting=======================================================================

sa_location <- ggplot(data = map_dat) +
  theme_bw() +
  geom_sf() +
  coord_sf(xlim = c(-114, -125.5), ylim = c(41,51)) +
  # theme(
  #   panel.grid.major = element_line(
  #     color = "gray", 
  #     linetype = "dashed", 
  #     size = 0.5)) +
  geom_rect(
    aes(
      xmin = sa_rect[1], 
      ymin = sa_rect[2], 
      xmax = sa_rect[3], 
      ymax = sa_rect[4]), 
    fill = NA, 
    color = "black") +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text = element_blank())
# ggsave(
#   "figures//study_area_setting_fig.png", 
#   width = 5, 
#   height = 5, 
#   units = "in", 
#   dpi = 500)

# sa_hillshade <- ggplot() +
#   geom_raster(
#     data = hill_df,
#     aes(x, y, fill = lyr1)
#   ) +
#   scale_fill_distiller(palette = "Greys") +
#   geom_sf(data = starkey, alpha = 0, color = "black") +
#   coord_sf(xlim = c(-118.625, -118.45), ylim = c(45.18, 45.40)) +
#   theme(
#     panel.grid.major = element_line(
#       color = "gray", 
#       linetype = "dashed", 
#       size = 0.5)) +
#   theme(legend.position = "none") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
#   theme(axis.text.y = element_text(angle = 0, hjust = 0)) +
#   scale_x_continuous(breaks = c(-118.70, -118.60, -118.50, -118.40)) +
#   labs(x = "", y = "") 
sa_hillshade <- ggplot() +
  geom_raster(
    data = hill_df,
    aes(x, y, fill = lyr1)
  ) +
  scale_fill_distiller(palette = "Greys") +
  geom_sf(data = starkey, alpha = 0, color = "black") +
  coord_sf(
    xlim = c(sa_rect[1] + .02, sa_rect[3] - .02), 
    ylim = c(sa_rect[2] + .02, sa_rect[4] - .02)) +
  theme(
    panel.grid.major = element_line(
      color = "gray", 
      linetype = "dashed", 
      size = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "", y = "") 

sa_location %>%
  ggdraw() +
  draw_plot(
    {
      sa_hillshade +
    },
    x = 0.45,
    y = 0.6,
    width = 0.25,
    height = 0.25
  )

ggsave(
  "figures//study_area_hillshade_fig.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 500)


plot_grid(
  setting, starkey_hillshade, 
  align = "v", 
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//study_area_setting_fig.png", 
  width = 5, 
  height = 5, 
  units = "in", 
  dpi = 500)
