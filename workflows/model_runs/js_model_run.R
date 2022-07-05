# Starkey elk CJS model fit workflow
# Kenneth Loonam

#Variables======================================================================

# amount to MULTIPLY the population by for augmentation
aug <- 2

start_year <- 1988
end_year <- 2021

# Parameters to track
params <- c(
  # "sd_ind",
  "survival_af", 
  "survival_am", 
  "survival_ca",
  "s0_ps",
  "p0_ps",
  "N"
)

# File names/paths
model_file <- "models//survival//js_model_elk.txt"
result_file <- "results//survival//js_rslt_28jun2022.Rdata"

# Sampler variables
n_i <- 750
n_t <- 1
n_b <- 250
n_c <- 3
n_a <- 10



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
  mutate('2021' = 1) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  apply(., 1, function(x) min(which(x != 0)))

w <- elk_data$hnt_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() 

y <- w + y # harvests count as observed alive that year (did not die naturally)

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
herd <- herd[-bad_elk,]

# y <- y[1:50,]
# f <- f[1:50]
# l <- l[1:50]
# w <- w[1:50,]
# male <- male[1:50]
# calf <- calf[1:50,]
# herd <- herd[1:50,]

z <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
l_k <- rep(NA, nrow(y))
for(i in 1:nrow(y)){
  l_k[i] <- max(which(y[i,] == 1))
  z[i, (f[i]):l_k[i]] <- 1
}

#Augment========================================================================

n_aug <- nrow(y) * aug
y <- rbind(y, matrix(0,  nrow = n_aug, ncol = ncol(y)))
z <- rbind(z, matrix(NA, nrow = n_aug, ncol = ncol(y)))
l <- c(l, rep(34, n_aug))
male <- c(male, rbinom(n_aug, 1, 0.5))
herd <- rbind(herd, matrix(0, nrow = n_aug, ncol = ncol(y)))

#Fit_model======================================================================



js_data <- list(
  y = y,
  z = z,
  # f = f,
  l = l,
  m = male,
  # c = calf,
  h = herd,
  n_occ = ncol(y),
  M = nrow(y)
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

require(rjags)
js_model <- jags.model(
  file = model_file,
  data = js_data,
  n.chains = 1,
  n.adapt = 100
)

update(
  object = js_model,
  n.iter = n_b
)

coda.samples(
  model = js_model,
  variable.names = params,
  n.iter = n_i,
  n.thin = n_t
)

# js_nimble_model <- readBUGSmodel(
#   model = model_file,
#   data = js_data,
#   inits = js_inits,
#   calculate = F,
#   check = F
# )
# 
# js_rslt <- nimbleMCMC(
#   model = cjs_nimble_model,
#   monitors = params,
#   niter = n_i,
#   nburnin = n_b,
#   nchains = n_c,
#   thin = n_t
# )
# 
# mcmcplots::mcmcplot(js_rslt)
# 
# save(js_rslt, file = result_file)