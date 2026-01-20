# Code to create publication study area map of Starkey
# Kenneth Loonam
# February 2023

#Environment====================================================================

require(raster)
require(ggplot2); require(rnaturalearth); require(maps); require(elevatr)
require(rnaturalearthdata); require(ggspatial); require(dplyr); require(sf);
require(terra); require(ggnewscale); require(purrr)
require(cowplot); require(tidyterra); require(giscoR); require(ggblend)
require(ggsci); require(grDevices); require(grid); 

#Data prep variables============================================================

study_area_file <- "data//Spatial//Starkey_main"
sa_rect_buffer <- c(x = 0.015, y = 0.015) # size of area around SA to show
elevation_zoom_level <- 12 # resolution for elevation, 1-14, 1 is low detail
projection <- 4326 # What projection do you want to work in
ins_lim_buff <- 4e5

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
starkey_2 <- st_transform(starkey, "+proj=laea +lat_0=30 +lon_0=-95")
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
sa_centroid <- sa_star$geometry %>% unlist()

x_ins_lims <- c(sa_centroid[1] - ins_lim_buff, sa_centroid[1] + ins_lim_buff)
y_ins_lims <- c(sa_centroid[2] - ins_lim_buff, sa_centroid[2] + ins_lim_buff)

geometry2 <- st_sfc(st_polygon(x = list(matrix(c(
  sa_rect[1]-.0001, sa_rect[2]-.0002,
  sa_rect[1]-.0001, sa_rect[4]-.0002,
  sa_rect[3]+.0001, sa_rect[4]-.0002,
  sa_rect[3]+.0001, sa_rect[2]-.0002,
  sa_rect[1]-.0001, sa_rect[2]-.0002
), ncol = 2, byrow = T
))))
# Make consistent names and define the projection
sa_rect_sf2 <- st_sf(name = "name", geometry = geometry2)
st_crs(sa_rect_sf2) <- projection

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

#Plotting=======================================================================

loc_fig_1 <- ggplot(map_dat2) +
  theme_bw() +
  geom_sf(fill = "#f5f5dc") +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  # Set map bounds
  scale_y_continuous(limits = c(-600000, 3100000)) +
  scale_x_continuous(limits = c(-2700000, -1200000)) +
  theme(
    panel.grid.major = element_line(
      color = "#909090",
      linetype = "dashed",
      linewidth = 0.3),
    panel.background = element_rect(fill = "#afeeee")) +
  # Draw the inset rectangle
  annotate("rect",
           xmin = x_ins_lims[1], 
           xmax = x_ins_lims[2],
           ymin = y_ins_lims[1],
           ymax = y_ins_lims[2],
           alpha = 0,
           color = "black") + 
  labs(x = NULL, y = NULL)

loc_fig_2 <- ggplot(map_dat2) +
  theme_bw() +
  geom_sf(fill = "#f5f5dc", linewidth = 0.5, color = "#909090") +
  geom_point(
    x = sa_centroid[1], 
    y = sa_centroid[2],
    shape = 21,
    color = "black",
    fill = "red"
    )+
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(
    panel.grid.major = element_line(
      color = "#909090",
      linetype = "dashed",
      linewidth = 0.3),
    panel.background = element_rect(fill = "#afeeee")) +
  # Set map bounds
  scale_y_continuous(limits = y_ins_lims) +
  scale_x_continuous(limits = x_ins_lims) +
  labs(x = NULL, y = NULL)

stky_fig <- ggplot() +
  geom_sf(data = sa_rect_sf2, fill = NA, color = "black", lwd = 1) +
  geom_raster( # add the hill shade layer
    data = hill1_df,
    aes(x = x, y = y, fill = hillshade),
    show.legend = F) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_raster(
    data = sa_ele_df,
    aes(x, y, fill = ele),
    alpha = 0.7
  ) +
  scale_fill_gradientn(
    colors = hcl.colors(n = 50, palette = "Terrain"), 
    limits = c(1100, 1525)
  ) +
  # redraw study area shape so it is on top of the hill shade layer
  geom_sf(data = starkey, alpha = 0, color = "black", linewidth = .750) +
  theme_bw() +
  guides(fill = guide_colorbar(
    title.theme = element_text(size = 5.75),
    label.theme = element_text(size = 4.75),
    title.position = "top",
    title.hjust = 0,
    label.hjust = 0.25,
    direction = "horizontal",
    barheight = 0.15,
    barwidth = 4
  )) +
  theme(legend.position = c(0.37, 0.135)) +
  labs(fill = "Elevation (meters)") + 
  annotation_scale( # add the scale
    location = "br", 
    text_cex = .6,
    pad_x = unit(0.4, "cm"),
    pad_y = unit(1.03, "cm"),
    height = unit(0.12, "cm")
  ) +
  annotation_north_arrow( # add the north arrow
    location = "br",
    style = north_arrow_orienteering(text_size = 5),
    height = unit(0.65, "cm"),
    width = unit(0.45, "cm"),
    pad_x = unit(0.89, "cm"),
    pad_y = unit(1.3, "cm")
  ) +
  labs(y = NULL, x = NULL) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 5.75),
    axis.text.x.bottom = element_text(vjust = 15),
    legend.background = element_rect(
      fill = alpha('#999999', .3),
      linewidth = 0.2,
      color = "#404040"),
    legend.margin = margin(t = 0.03, r = 0.2, b = 0.02, l = 0.13, unit = 'cm')
  ) +
  scale_x_continuous(limits = c(-118.65, -118.489)) +
  scale_y_continuous(limits = c(45.165, 45.329)) +
  theme(panel.border = element_blank(),
        axis.text.y.left = element_blank(),
        axis.text.x.bottom = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank()) +
  annotate(
    geom = "text", 
    label = "45.30°N", 
    x = -118.646, 
    y = 45.30, 
    size = 2,
    hjust = 0) +
  annotate(
    x = -118.62529, xend = -118.489, 
    y = 45.30, yend = 45.30, 
    colour = "#666666", lwd=0.4, geom="segment", linetype = "dashed") +
  annotate(
    geom = "text", 
    label = "45.25°N", 
    x = -118.646, 
    y = 45.25, 
    size = 2,
    hjust = 0) +
  annotate(
    x = -118.62529, xend = -118.489, 
    y = 45.25, yend = 45.25, 
    colour = "#666666", lwd=0.4, geom="segment", linetype = "dashed") +
  annotate(
    geom = "text", 
    label = "45.20°N", 
    x = -118.646, 
    y = 45.20, 
    size = 2,
    hjust = 0) +
  annotate(
    x = -118.62529, xend = -118.489, 
    y = 45.20, yend = 45.20, 
    colour = "#666666", lwd=0.4, geom="segment", linetype = "dashed") +
  annotate(
    geom = "text", 
    label = "118.60°W", 
    x = -118.60, 
    y = 45.1725, 
    size = 2,
    hjust = 0.5,
    vjust = 1) +
  annotate(
    x = -118.60, xend = -118.60, 
    y = 45.17392, yend = 45.3290, 
    colour = "#666666", lwd=0.4, geom="segment", linetype = "dashed") +
  annotate(
    geom = "text", 
    label = "118.55°W", 
    x = -118.55, 
    y = 45.1725, 
    size = 2,
    hjust = 0.5,
    vjust = 1) +
  annotate(
    x = -118.55, xend = -118.55, 
    y = 45.17392, yend = 45.3290, 
    colour = "#666666", lwd=0.4, geom="segment", linetype = "dashed") +
  annotate(
    geom = "text", 
    label = "118.50°W", 
    x = -118.50, 
    y = 45.1725, 
    size = 2,
    hjust = 0.5, 
    vjust = 1) +
  annotate(
    x = -118.50, xend = -118.50, 
    y = 45.17392,  yend = 45.329,
    colour = "#666666", lwd=0.4, geom="segment", linetype = "dashed")

tmp <- stky_fig %>%
  ggdraw() +
  draw_plot(
    loc_fig_1,
    x = 0.621 - 0.0254 + 0.001,
    y = 0.681 - 0.042 + 0.002,
    width = 0.3,
    height = 0.3
  ) 
tmp %>%
  ggdraw() +
  draw_plot(
    loc_fig_2,
    x = 0.564 - 0.026 +.001,
    y = 0.8305 - 0.0415 + 0.002,
    width = 0.15,
    height = 0.15
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

