# Prepare the miscellaneous data that the IPM uses
# The data hard coded here is hand compiled from capture and handling summaries
# That is compared to the capture/handling database data where appropriate
# The larger value of the two is used in all cases 
  # (Missing elk is more likely than counting extra elk)

#Functions======================================================================
source("functions//cjs_data_prep_functions.R")

ilogit <- function(x){
  1 / (1 + exp(-x))
}
#Data===========================================================================
chd <- readRDS(capture_handling_data)

fg_counts <- tibble(
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
    11,   16,  10,  23,  28,  48,  44,  31,  62,  49,  75,  89,  89,  26,  36,
    35,   12,  14
  )
)

fg_management <- tibble(
  yr  = rep(1989:2021, 4),
  sex = c(rep("f", 33), rep("m", 33), rep("f", 33), rep("m", 33)),
  age = c(rep("ad", 66), rep("ca", 66)),
  n_moved = c(
    # Female Adults
    0, 0, 0, -5, 0, 0, 0, 0, -32, 0, -54, -71, -29, 0, 13, 7, 13, 0, 9, 0, 1, 1,
    0, 0, 6, 13, -17, 0, -111, 0, 0, -95, -47,
    # Male Adults
    0, 0, 0, -3, 0, 0, 0, -4, -7, 0, -16, -3, 5, 30, 15, 7, 0, 0, 7, 0, 3, 1, 0,
    0, 0, 3, 0, 0, 0, 0, 0, -46, -5,
    # Female Calves
    0, 0, 0, -5, 0, 0, 25, 0, 0, 0, 0, 0, 0, 0, 2, 8, 0, 0, 1, 8, 1, 13, 0, 1,
    3, 3, -19, 0, -20, 0, 0, -24, -7,
    # Male Calves
    0, 0, 0, -5, 0, 0, 25, -17, -35, 0, -26, -32, -19, -13, 2, 8, 14, 0, 2, 16,
    9, 13, 9, 4, 6, 4, -15, 0, -22, 0, 0, -11, -1
  )
)

#Clean the capture and handling data============================================

ch <- chd$capture_history %>% to_matrix_wo_filter(yr_range)
herd <- chd$herd_assignment %>% to_matrix_wo_filter(yr_range)
age <- chd$annual_age %>% to_matrix_wo_filter(yr_range)
known_alive <- chd$known_alive %>% 
  mutate(across(2:ncol(.), as.numeric)) %>%
  to_matrix_wo_filter(yr_range)
harvest <- chd$harvest_year %>% to_matrix_wo_filter(yr_range)
sex <- matrix(NA, nrow = nrow(ch), ncol = ncol(ch))
for(i in 1:ncol(sex)){
  sex[,i] <- chd$sex %>% arrange(id) %>% pull(sex)
}  

#Minimum known alive============================================================

nf_min_ch <- (known_alive * (sex == "F") * (herd == "main") * (age != 1)) %>%
  colSums(na.rm = T)
nm_min_ch <- (known_alive * (sex == "M") * (herd == "main") * (age != 1)) %>%
  colSums(na.rm = T)
ncf_min_ch <- (known_alive * (sex == "F") * (herd == "main") * (age == 1)) %>%
  colSums(na.rm = T)
ncm_min_ch <- (known_alive * (sex == "M") * (herd == "main") * (age == 1)) %>%
  colSums(na.rm = T)
nc_min_ch <- ncm_min_ch + ncf_min_ch
n_min <- tibble(
  yr = names(nf_min_ch),
  nf = nf_min_ch,
  nm = nm_min_ch,
  nc = nc_min_ch
) %>%
  mutate(yr = as.numeric(yr)) %>%
  full_join(fg_counts) %>%
  group_by(yr) %>%
  summarise(
    female = max(nf, fe, na.rm = T),
    male   = max(nm, ma, na.rm = T),
    calf   = max(nc, ca, na.rm = T)
  ) %>% ungroup()

#N Harvested====================================================================
nf_har <- (harvest * (sex == "F") * (herd == "main") * (age != 1)) %>%
  colSums(na.rm = T)
nm_har <- (harvest * (sex == "M") * (herd == "main") * (age != 1)) %>%
  colSums(na.rm = T)
fc_har <- (harvest * (sex == "F") * (herd == "main") * (age == 1)) %>%
  colSums(na.rm = T)
mc_har <- (harvest * (sex == "M") * (herd == "main") * (age == 1)) %>%
  colSums(na.rm = T)
nr <- length(nf_har)
n_har <- tibble(
  yr = rep(as.numeric(names(nf_har)), 4),
  sex = c(rep("f", nr), rep("m", nr), rep("f", nr), rep("m", nr)),
  age = c(rep("ad", nr*2), rep("ca", nr*2)),
  harvest = c(nf_har, nm_har, fc_har, mc_har)
)

#Counts=========================================================================
# This pulls the number observed on the feed-grounds (CJS data) and divides by
# the estimated detection probability (CJS) to give a probability correct count
cjs_rslt <- readRDS(cjs_results_summary)

nf_ch <- (ch * (sex == "F") * (herd == "main") * (age != 1)) %>%
  colSums(na.rm = T)
nm_ch <- (ch * (sex == "M") * (herd == "main") * (age != 1)) %>%
  colSums(na.rm = T)
ncf_ch <- (ch * (sex == "F") * (herd == "main") * (age == 1)) %>%
  colSums(na.rm = T)
ncm_ch <- (ch * (sex == "M") * (herd == "main") * (age == 1)) %>%
  colSums(na.rm = T)

index_available <- which(names(nf_ch) %in% cjs_rslt$sf_logit$yr)
fac <- nf_ch[index_available]
mac <- nm_ch[index_available] 
cc <- ncf_ch[index_available] + ncm_ch[index_available]
fp <- ilogit(cjs_rslt$pf_logit$mn)
mp <- ilogit(cjs_rslt$pm_logit$mn)
nr <- length(fac)
counts_cjs <- tibble(
  yr = rep(as.numeric(names(fac)), 3),
  sex = c(rep("f", nr), rep("m", nr), rep("c", nr)),
  age = c(rep("ad", nr*2), rep("ca", nr)),
  count = c(fac, mac, cc),
  prob = c(fp, mp, fp)
) %>%
  filter(yr != 1991) %>%
  filter(count > 0)

#Clean-up=======================================================================

dtlist <- list(
  n_min = n_min,
  n_har = n_har,
  count = counts_cjs,
  n_mov = fg_management
)
saveRDS(dtlist, "data//misc_data.rds")
rm(list = ls())
