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

fg_counts <- read_csv(feedground_counts) %>% 
  mutate(yr = year) %>% 
  select(-year)

fg_management <- read_csv(management_moves)
  

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
    female = max(nf, fe_ad, na.rm = T),
    male   = max(nm, ma_ad, na.rm = T),
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
