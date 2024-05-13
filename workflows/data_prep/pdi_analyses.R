require(dplyr); require(tidyr); require(ggplot2); require(readr)

pdi <- read_csv("data//climate//pdi_3508_ne_or.csv") %>%
  pivot_longer(cols = 2:13, names_to = "month", values_to = "pdi") %>%
  filter(year > 1986) %>% 
  filter(month %in% c("may", "june", "july", "august", "september"))

pdi_comparison <- pdi %>% group_by(year) %>%
  summarise(mean_season = mean(pdi)) %>%
  mutate(september_pdi = pdi$pdi[pdi$month == "september"]) 

ggplot(pdi_comparison, aes(x = mean_season, y = september_pdi, color = year)) +
  geom_point() +
  theme_classic() +
  xlab("Mean May-September PDI") + ylab("September PDI") +
  labs(title = "1987 - 2022 PDI Comparison", color = "Year")

lm(pdi_comparison$september_pdi ~ pdi_comparison$mean_season) |> summary()  

full_sep_pdi <- read_csv("data//climate//pdi_3508_ne_or.csv") %>%
  pivot_longer(cols = 2:13, names_to = "month", values_to = "pdi") %>%
  filter(month == "september")

full_pdi <- read_csv("data//climate//pdi_3508_ne_or.csv") %>%
  pivot_longer(cols = 2:13, names_to = "month", values_to = "pdi") %>%
  filter(month %in% c("may", "june", "july", "august", "september")) %>%
  group_by(year) %>%
  summarise(mean_season = mean(pdi)) %>%
  mutate(september_pdi = full_sep_pdi$pdi) 

ggplot(full_pdi, aes(x = mean_season, y = september_pdi, color = year)) +
  geom_point() +
  theme_classic() +
  xlab("Mean May-September PDI") + ylab("September PDI") +
  labs(title = "1895 - 2022 PDI Comparison", color = "Year")

lm(full_pdi$september_pdi ~ full_pdi$mean_season) |> summary()  





avhrr <- read_csv("data//climate//ndvi_avhrr.csv") %>% 
  group_by(mn, yr) %>%
  summarise(avhrr = mean(NDVI)) %>%
  ungroup() %>%
  filter(yr > 2000 & yr < 2004)
modis <- read_csv("data//climate//ndvi_modis.csv") %>% 
  group_by(mn, yr) %>%
  summarise(modis = mean(NDVI)) %>%
  ungroup() %>%
  filter(yr > 2000 & yr < 2004)

ndvi <- full_join(avhrr, modis)

ggplot(ndvi, aes(x = modis, y = avhrr, color = as.factor(yr))) +
  geom_point() +
  xlab("Modis")+
  ylab("AVHRR") +
  labs(title = "Comparision of NDVI sources (2001-2003)", color = "Year") +
  theme_classic() 

logit <- function(x){log(x/(1-x))}
ilogit <- function(x){1/(1+exp(-x))}
rnorm(10000, 0, 4) |> ilogit() |> hist()
