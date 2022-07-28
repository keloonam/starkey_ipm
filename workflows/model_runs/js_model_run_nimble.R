# Starkey elk CJS model fit workflow
# Kenneth Loonam

#Variables======================================================================

# amount to MULTIPLY the population by for augmentation
aug <- 2

start_year <- 1988
end_year <- 2021

# Parameters to track
params <- c(
  "survival_af", 
  "survival_am", 
  "survival_ca",
  "detection_f",
  "detection_m",
  "detection_c",
  "N_super",
  "N",
  "lambda"
)

# File names/paths
result_file <- "results//survival//js_rslt_28jul2022.Rdata"

# Sampler variables
n_i <- 750
n_t <- 1
n_b <- 250
n_c <- 1
n_a <- 10



#Environment====================================================================

require(tidyverse); require(mcmcplots); require(nimble)
load("data//elk_data.Rdata")

#Prep_data======================================================================

y <- elk_data$cap_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()

# f <- apply(y, 1, function(x) min(which(x != 0)))

l <- elk_data$hnt_tib %>%
  # filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  mutate('2021' = 1) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  apply(., 1, function(x) min(which(x != 0)))

w <- elk_data$hnt_tib %>%
  # filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() 

y <- w + y # harvests count as observed alive that year (did not die naturally)

male <- elk_data$sex_tib %>%
  # filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male)

calf <- elk_data$age_tib %>%
  # filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()
for(i in 1:nrow(calf)){
  if(sum(calf[i,], na.rm = T) > 0){
    calf_sess <- min(which(calf[i,] == 1))
    calf[i,calf_sess:ncol(calf)] <- 1
  }
  if(all(calf[i,] == 0)){
    calf[i,1] <- 2
    calf[i,2:ncol(calf)] <- 1
  }
  if(any(!is.na(calf[i,]))){
    start_age <- min(which(!is.na(calf[i])))
    calf[i,start_age:ncol(calf)] <- cumsum(calf[i,start_age:ncol(calf)])
  }
}

herd <- elk_data$hrd_tib %>%
  # filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  mutate_all(funs(case_when(
    . == "main" ~ 0,
    T ~ 1
  ))) %>%
  as.matrix()

# We don't care about elk that were never in main study
gone_elk <- which(apply(herd, 1, sum) == ncol(herd))
# weird_elk <- f == l
# bad_elk <- which((gone_elk + weird_elk) != 0)
y <- y[-gone_elk,]
f <- f[-gone_elk]
l <- l[-gone_elk]
w <- w[-gone_elk,]
male <- male[-gone_elk]
calf <- calf[-gone_elk,]
herd <- herd[-gone_elk,]

# y <- y[1:50,]
# f <- f[1:50]
# l <- l[1:50]
# w <- w[1:50,]
# male <- male[1:50]
# calf <- calf[1:50,]
# herd <- herd[1:50,]

# z <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
# l_k <- rep(NA, nrow(y))
# for(i in 1:nrow(y)){
#   l_k[i] <- max(which(y[i,] == 1))
#   z[i, (f[i]):l_k[i]] <- 1
# }

z <- elk_data$liv_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()

#Augment========================================================================

n_aug <- nrow(y) * aug
y <- rbind(y, matrix(0,  nrow = n_aug, ncol = ncol(y)))
z <- rbind(z, matrix(NA, nrow = n_aug, ncol = ncol(y)))
l <- c(l, rep(34, n_aug))
male <- c(male, rbinom(n_aug, 1, 0.5))
herd <- rbind(herd, matrix(0,  nrow = n_aug, ncol = ncol(y)))
calf <- rbind(calf, matrix(NA, nrow = n_aug, ncol = ncol(y)))

#Fit_model======================================================================



js_data <- list(
  y = y,
  z = z,
  L = l,
  M = male,
  H = herd,
  K = ncol(y),
  M = nrow(y),
  AGE = 
)

js_inits <- list(
  s0_ps  = runif(ncol(y), 0.5, 0.9),
  p0_ps  = runif(ncol(y), 0, 1),
  a0_ps  = runif(ncol(y), 0, 1),
  sm     = rnorm(ncol(y)),
  sc     = rnorm(ncol(y)),
  pm     = rnorm(ncol(y)),
  sh     = rnorm(1),
  ph     = rnorm(1),
  sd_ind = rnorm(1),
  b_ind  = rnorm(nrow(y))
)

# Load model as js_superpop_model
source("models/survival/js_model_nimble.R")

js_conf <- configureMCMC(
  model = js_superpop_model,
  monitors = params,
  control = list(adaptInterval = n_a), 
  thin = n_t, 
  useConjugacy = TRUE
)
js_mcmc <- buildMCMC(js_conf)
js_comp <- compileNimble(js_superpop_model)
js_mcmc_comp <- compileNimble(
  js_mcmc,
  project = js_superpop_model
)
rslt <- runMCMC(
  js_mcmc_comp,
  niter = n_i,
  nburnin = n_b,
  nchains = n_c,
  inits = inits,
  setSeed = F,
  progressBar = T,
  samplesAsCodaMCMC = T
)

save(rslt, file = result_file)
