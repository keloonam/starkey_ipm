# Population simulation
# Kenneth Loonam
# January 2020

#Variables======================================================================

# Model
model <- "cjs_phiyear_pyear"

# Sampler variables
ni <- 50000
nt <- 1
nb <- 5000
nc <- 3

# Starting population
n_af <- 50
n_am <- 10
n_jf <- 20
n_jm <- 20

# Mean demographic rates
af_s <- 0.9
am_s <- 0.7
jf_s <- 0.7
jm_s <- 0.7
recruit <- 0.5

# Capture rates
cap_p <- 0.5
trap_hap <- 0.2
n_years <- 5

#Environment====================================================================

require(tidyverse); require(R2jags); require(mcmcplots)
source("r//simulation_functions.R"); source("r//data_prep_functions.R")

#Simulate_population============================================================

pop <- sim_pop_fn(
  n_af = n_af,
  n_am = n_am,
  n_jf = n_jf,
  n_jm = n_jm,
  af_s = af_s,
  am_s = am_s,
  jf_s = jf_s,
  jm_s = jm_s,
  recruit = recruit,
  n_years = n_years
)

eh <- sim_cap_fn(
  population = pop,
  cap_p = cap_p
)

#Prep_data======================================================================

# data
f <- apply(eh, 1, function(x) min(which(x != 0)) )
model_data <- list(
  n_occ = ncol(eh),
  n_ind = nrow(eh),
  y = eh,
  f = f,
  pr_p = 0.5
)

# initial values
z <- ch_init_fn(eh, f)

model_inits <- function(){
  list(
    z = z,
    b0.phi = runif(ncol(eh), -4, 4),
    b0.p = runif(ncol(eh), -4, 4)
  )
}

#Fit_model======================================================================

model_file <- paste0("models//", model, ".txt")

params <- c("mean_phi", "mean_p")


# run the MCMC chain in JAGS
rslt <- jags(
  data = model_data, 
  inits = model_inits,
  parameters.to.save = params,
  model.file = model_file,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.thin = nt
)
 mcmcplot(rslt)
 