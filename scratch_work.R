ggplot() +
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

ggsave( # save at high resolution
  "figures//study_area_figure.png", 
  width = 4, 
  height = 4, 
  units = "in", 
  dpi = 500)
