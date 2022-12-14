require(readr); require(dplyr); require(lubridate); require(ggplot2)
avhrr <- read_csv("data//climate//ndvi_avhrr.csv") %>%
  mutate(source = "avhrr") %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(NDVI))
# avhrr$ndvi <- as.numeric(scale(avhrr$ndvi))
modis <- read_csv("data//climate//ndvi_modis.csv") %>%
  mutate(source = "modis") %>%
  filter(SummaryQA <= 1) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(NDVI))
# modis$ndvi <- as.numeric(scale(modis$ndvi))

combo <- inner_join(avhrr, modis, by = c("mn", "yr"), suffix = c("_a", "_m")) %>%
  arrange(yr, mn) %>%
  filter(yr < 2015) %>%
  filter(mn > 4, mn < 10)

m1 <- with(combo, lm(ndvi_m ~ ndvi_a))
m2 <- with(combo, lm(ndvi_a ~ ndvi_m))
diff <- mean(combo$ndvi_a - combo$ndvi_m)

avhrr_m1 <- avhrr %>%
  mutate(
    ndvi = ndvi*m1$coefficients[2] + m1$coefficients[1]
  ) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(ndvi))

avhrr_m2 <- avhrr %>%
  mutate(
    ndvi = (ndvi - m2$coefficients[1]) / m2$coefficients[2]
  ) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(ndvi))

ndvi_m1 <- full_join(avhrr_m1, modis, by = c("mn", "yr"), suffix = c("_a", "_m")) %>%
  mutate(ndvi = case_when(
    !is.na(ndvi_m) ~ ndvi_m,
    T ~ ndvi_a
  )) %>%
  filter(mn > 4, mn < 10) %>%
  group_by(yr) %>%
  summarise(summer_ndvi = mean(ndvi))
ndvi_m2 <- full_join(avhrr_m2, modis, by = c("mn", "yr"), suffix = c("_a", "_m")) %>%
  mutate(ndvi = case_when(
    !is.na(ndvi_m) ~ ndvi_m,
    T ~ ndvi_a
  )) %>%
  filter(mn > 4, mn < 10) %>%
  group_by(yr) %>%
  summarise(summer_ndvi = mean(ndvi))
plot(ndvi_m1)
plot(ndvi_m2)

####################
avhrr %>%
  filter(mn > 4, mn < 10) %>%
  group_by(yr) %>%
  summarise(mean_ndvi = mean(ndvi)) %>%
  filter(yr < 2015) %>%
  plot()
  