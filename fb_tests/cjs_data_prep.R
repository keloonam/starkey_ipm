# Prepare data for the cjs run
# Kenneth Loonam
# August 2024

#Environment====================================================================
require(tidyverse)
source("functions//cjs_data_prep_functions.R")
dtlist <- readRDS("data//capture_handling_data.rds")

ids_to_remove <- find_removals(
  data_list = dtlist, 
  yr_range = yr_range)

y <- to_matrix(
  x = dtlist$capture_history, 
  yr_range = yr_range, 
  removals = ids_to_remove
)

hy <- to_matrix(
  x = dtlist$harvest_year, 
  yr_range = yr_range, 
  removals = ids_to_remove
)

y <- ((y + hy) > 0)*1

f <- find_first_value_position(y)

z <- build_z(y, f)

male <- as.numeric(to_vector(
  dtlist$sex, 
  removals = ids_to_remove, 
  target_column = "sex"
) == "M")

herd <- (to_matrix(
  dtlist$herd_assignment,
  yr_range = yr_range,
  removals = ids_to_remove
) != "main") * 1

calf <- (to_matrix(
  dtlist$annual_age,
  yr_range = yr_range,
  removals = ids_to_remove
) == 1) * 1

nocc <- ncol(y)
nind <- nrow(y)
l <- dtlist$harvest_year %>%
  to_matrix(yr_range = yr_range, removals = ids_to_remove) %>%
  find_first_value_position()

full_cjs_data <- list(
  constants = list(
    nocc = nocc,
    nind = nind,
    l = l,
    f = f
  ),
  data = list(
    y = y,
    z = z,
    m = male,
    h = herd,
    c = calf
  ),
  initial_values = list(
    sf = runif(nocc, 0, 1),
    sm = runif(nocc, 0, 1),
    sc = runif(nocc, 0, 1),
    pf = runif(nocc, 0, 1),
    pm = runif(nocc, 0, 1),
    sh = runif(1, 0, 1),
    ph = runif(1, 0, 1)
  )
)

saveRDS(full_cjs_data, "fb_tests//cjs_data.rds")
rm(list = ls()[-which(ls() == "yr_range")])
