# Source some functions:
source("functions//cjs_data_prep_functions.R")

# Manually pull data from capture/handling summaries:
# I've hard coded this data to avoid sending 10 million .csvs with the IPM
summary_count_data <- tibble(
  yr = 1989:2021,
  fe = c(
    306, 169,   3,  48,  48, 102,  76, 112,  84, 122, 113,  83,  84,  56,  89, 
     45,  60,  49, 105, 101, 124, 150, 107, 165, 155, 190, 234, 276, 149, 131,
    122,  46,  13
    ),
  ma = c(
    58, 46,  1,  2,  2,  8, 12, 17,  9, 23, 37, 28, 22, 16, 46,  9, 16, 12, 39,
    27, 28, 52,  3, 30, 21, 33, 46, 70, 64, 51, 48,  8,  0
    ),
  ca = c(
    128,  91,   0,   8,   0,   7,   2,  84,  10,  62,  45,  35,  39,  16,  39,
     11,  16,  10,  23,  28,  48,  44,  31,  62,  49,  75,  89,  89,  26,  36,
     35,  12,  14
    )
)

# Pull the capture and handling database data
fd <- readRDS(capture_handling_data)
# Pull capture history (1 = counted that year)
ch <- fd$capture_history %>% to_matrix_wo_filter(yr_range)
# Pull sex information (1 = female)
sx <- ((fd$sex %>% arrange(id) %>% pull(sex)) == "F") * 1
# Pull herd information (1 = in main)
he <- ((fd$herd_assignment %>% to_matrix_wo_filter(yr_range)) == "main") * 1
# Pull age information (1 = is a calf)
ca <- ((fd$annual_age %>% to_matrix_wo_filter(yr_range)) == 1) * 1
# Count the cows and calves
nf <- rep(NA, ncol(ch))
nc <- rep(NA, ncol(ch))
for(i in 1:ncol(ch)){
  nf[i] <- sum(ch[,i] * sx * he[,i] * (1-ca[,i]), na.rm = T)
  nc[i] <- sum(ch[,i] * sx * he[,i] *    ca[,i],  na.rm = T)
}

# Combine the two data sources
ratio_data <- tibble(
  yr = as.numeric(dimnames(ch)[[2]]),
  fy = nf,
  cy = nc
) %>% 
  full_join(summary_count_data) %>%
  mutate(females = case_when(
    is.na(fe) ~ fy,
    is.na(fy) ~ fe,
    fe > fy ~ fe,
    fy > fe ~ fy,
    T ~ fy
  )) %>%
  mutate(calves = case_when(
    is.na(fe) ~ cy,
    is.na(fy) ~ ca,
    fe > fy ~ ca,
    fy > fe ~ cy,
    T ~ cy
  )) %>%
  mutate(females = case_when(
    females < calves ~ calves,
    T ~ females
  )) %>%
  select(yr, females, calves) %>%
  filter(females > 0) %>%
  as.matrix()

dtlist <- list(
  data = list(
    n_calf = ratio_data[,3],
    n_cow  = ratio_data[,2]
  ),
  constants = list(
    n_years = nrow(ratio_data),
    years = ratio_data[,1]
  ),
  initial_values = list(R = rep(0.5, nrow(ratio_data)))
)

saveRDS(dtlist, "data//recruitment_data.rds")  
rm(list = ls())
