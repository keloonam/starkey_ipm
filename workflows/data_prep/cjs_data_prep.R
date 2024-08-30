# Prepare data for the cjs run
# Kenneth Loonam
# August 2024

#Environment====================================================================

source("functions//cjs_data_prep_functions.R")
dtlist <- readRDS("data//capture_handling_data.rds")


full_cjs_data <- list(
  constants = list(
    nocc = nocc,
    nind = nind,
    l = last_occ,
    f = first_occ
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