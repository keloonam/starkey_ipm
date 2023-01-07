require(readr); require(dplyr); require(lubridate); require(ggplot2)
avhrr <- read_csv("data//climate//ndvi_avhrr.csv") %>%
  mutate(source = "avhrr") %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(NDVI)) %>%
  filter(yr < 2015)

modis <- read_csv("data//climate//ndvi_modis.csv") %>%
  mutate(source = "modis") %>%
  filter(SummaryQA <= 1) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(NDVI))

combo <- inner_join(avhrr, modis, by = c("mn", "yr"), suffix = c("_a", "_m")) %>%
  arrange(yr, mn) %>%
  filter(yr < 2015) %>%
  filter(mn > 4, mn < 10)

m1 <- with(combo, lm(ndvi_m ~ ndvi_a))
avhrr_m1 <- avhrr %>%
  mutate(
    ndvi = ndvi*m1$coefficients[2] + m1$coefficients[1]
  ) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(ndvi))
ndvi_modis <- full_join(avhrr_m1, modis, by = c("mn", "yr"), suffix = c("_a", "_m")) %>%
  mutate(ndvi = case_when(
    !is.na(ndvi_m) ~ ndvi_m,
    T ~ ndvi_a
  )) %>%
  filter(mn > 4, mn < 10) %>%
  group_by(yr) %>%
  summarise(summer_ndvi = mean(ndvi)) %>%
  filter(yr %in% 1988:2021) %>%
  select(value) %>%
  scale() %>%
  as.vector()

m2 <- with(combo, lm(ndvi_a ~ ndvi_m))
modis_m2 <- modis %>%
  mutate(
    ndvi = ndvi*m2$coefficients[2] + m2$coefficients[1]
  ) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(ndvi))
ndvi_avhrr <- full_join(modis_m2, avhrr, by = c("mn", "yr"), suffix = c("_m", "_a")) %>%
  mutate(ndvi = case_when(
    !is.na(ndvi_a) ~ ndvi_a,
    T ~ ndvi_m
  )) %>%
  filter(mn > 4, mn < 10) %>%
  group_by(yr) %>%
  summarise(summer_ndvi = mean(ndvi)) %>%
  filter(yr %in% 1988:2021) %>%
  select(value) %>%
  scale() %>%
  as.vector()

  