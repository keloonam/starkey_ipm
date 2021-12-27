# Starkey elk CJS model fit workflow
# Kenneth Loonam

#Variables======================================================================

start_year <- 1988
end_year <- 2020

# Parameters to track
params <- c(
  "phi", 
  "p",
  "sd_ind"
)

# File names/paths
model_file <- "models//survival//cjs_model_elk.R"
result_file <- "results//survival//cjs_rslt_24dec2021.Rdata"


# Sampler variables
n_i <- 20
n_t <- 2
n_b <- 10
n_c <- 3

params <- c("survival_af", "survival_am", "survival_ca")

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots); require(nimble)
load("data//elk_data.Rdata")

ch_init_fn <- function(ch, f, l){
  for(i in 1:nrow(ch)){
    ch[i,1:f[i]] <- NA
  }
  out <- ifelse(!is.na(ch), 1, ch)
  for(i in 1:nrow(out)){
    if(l[i] < ncol(out)){
      out[i,(l[i]+1):ncol(out)] <- NA
    }
  }
  return(out)
}

#Prep_data======================================================================

y <- elk_data$cap_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()

f <- apply(y, 1, function(x) min(which(x != 0)))

l <- elk_data$hnt_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  mutate('2020' = 1) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  apply(., 1, function(x) min(which(x != 0)))

w <- elk_data$hnt_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  magrittr::subtract(1) %>%
  abs()

for(i in 1:length(l)){
  y[i,l[i]] <- 1
}

male <- elk_data$sex_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male)

calf <- elk_data$age_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  replace_na(0)

herd <- elk_data$hrd_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  mutate_all(funs(case_when(
    . == "main" ~ 0,
    T ~ 1
  ))) %>%
  as.matrix()

gone_elk <- apply(herd - 1, 1, sum) == ncol(herd)
weird_elk <- f == l
bad_elk <- which((gone_elk + weird_elk) != 0)
y <- y[-bad_elk,]
f <- f[-bad_elk]
l <- l[-bad_elk]
w <- w[-bad_elk,]
male <- male[-bad_elk]
calf <- calf[-bad_elk,]
herd <- herd[-bad_elk,] - 1

#Fit_model======================================================================

source("models//survival//cjs_model_elk.R")
cjs_constants <- list(
  f = f,
  l = l,
  male = male - 1,
  calf = calf,
  herd = herd,
  n_occ = ncol(y),
  n_ind = nrow(y)
)

z_data <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
for(i in 1:nrow(y)){
  z_data[i, f[i]] <- 1
}
cjs_data <- list(
  y = y,
  z = z_data
)

cjs_inits <- list(
  s0_ps  = runif(1, 0.5, 0.9),
  p0_ps  = runif(1, 0, 1),
  mu_ms  = rnorm(1),
  mu_cs  = rnorm(1),
  mu_mp  = rnorm(1),
  sd_ms  = runif(1,1,5),
  sd_cs  = runif(1,1,5),
  sd_mp  = runif(1,1,5),
  sd_s0  = runif(1,1,5),
  sd_p0  = runif(1,1,5),
  s0     = rnorm(nrow(y)),
  sm     = rnorm(nrow(y)),
  sc     = rnorm(nrow(y)),
  p0     = rnorm(nrow(y)),
  pm     = rnorm(nrow(y)),
  sh     = rnorm(1),
  ph     = rnorm(1),
  sd_ind = rnorm(1),
  b_ind  = rnorm(1),
  z      = ch_init_fn(y, f, l)
)

cjs_rslt <- nimbleMCMC(
  code = cjs_code,
  data = cjs_data,
  monitors = params,
  thin = 1,
  niter = 5000,
  nburnin = 2000,
  nchains = 3,
  # inits = cjs_inits,
  constants = cjs_constants
)

# cjs_nimble_model <- nimbleModel(
#   code = cjs_code,
#   name = "elk_cjs",
#   constants = cjs_constants,
#   data = cjs_data,
#   inits = cjs_inits
# )
# 
# MCMC_cjs <- configureMCMC(cjs_nimble_model, monitors = params, print = T)
# cjs_MCMC <- buildMCMC(MCMC_cjs)
# matt_type <- compileNimble(cjs_nimble_model, showCompilerOutput = TRUE)
# comp_cjs <- compileNimble(cjs_MCMC, project = matt_type)
# cjs_rslt <- runMCMC(
#   comp_cjs,
#   mcmc = comp_cjs,
#   niter = 1000,
#   nburnin = 100,
#   thin = 1,
#   nchains = 3
# )



























