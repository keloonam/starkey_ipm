# Prep summary stats and graphs for decision on age structure in the model
# Kenneth Loona
# August 2024

#Variables======================================================================

calf <- 1
young <- 2:3
prime <- 4:13
old <- 14
yr_a <- 1988
yr_z <- 2021

#Packages=======================================================================

require(dplyr); require(lubridate); require(ggplot2); require(purrr)
require(readr); require(tidyr)

#Build Data=====================================================================

load("data//elk_data.Rdata")

calf <- elk_data$age_tib %>%
  arrange(id) %>%
  select(as.character(yr_a:yr_z)) %>%
  as.matrix()

herd <- elk_data$hrd_tib %>%
  arrange(id) %>%
  select(as.character(yr_a:yr_z)) %>%
  mutate(across(all_of(as.character(yr_a:yr_z)), ~ (.x == "main")*1)) %>%
  as.matrix()

live <- elk_data$liv_tib %>%
  arrange(id) %>% 
  select(as.character(yr_a:yr_z)) %>%
  as.matrix()

sex <- elk_data$sex_tib %>%
  arrange(id)

age_matrix <- calf * herd * live 

adt <- matrix(0, nrow = max(age_matrix, na.rm = T), ncol = ncol(age_matrix))
for(t in 1:ncol(adt)){
  for(k in 1:nrow(adt)){
    adt[k,t] <- sum((age_matrix[,t]*(sex$Sex=="F")) == k, na.rm = T)
  }
}

tot <- colSums(adt[-1,])
yng <- colSums(adt[young,], na.rm = T) / tot
prm <- colSums(adt[prime,], na.rm = T) / tot
old <- colSums(adt[old:nrow(adt),], na.rm = T) / tot
rd <- tibble(
  year = 1988:2021,
  total = tot,
  young = yng,
  prime = colSums(adt[prime,], na.rm = T) / total,
  old = old
) %>%
  filter(year >= 1990) %>%
  pivot_longer(cols = c(young, prime, old))

ggplot(rd, aes(x = year, y = value, color = name)) +
  geom_line() +
  theme_classic() +
  xlab("Year") + ylab("Proportion") + 
  labs(title = "Elk age classes through time") +
  scale_color_manual(
    labels = c("Old (>13)", "Prime", "Young (<4)"),
    name = NULL,
    values = c(
      "#ffa800",
      "#005cff", 
      "#09be00"
    )
  ) +
  theme(panel.grid.minor.y = element_line(color = "#dddddd", linetype = "12")) + 
  theme(panel.grid.major.y = element_line(color = "#dddddd", linetype = "12"))
ggsave(
  "figures//age_distribution.tif", 
  dpi = 600, 
  units = "cm", 
  width = 18, 
  height = 10
)
