# Starkey elk CJS model fit workflow
# Kenneth Loonam

#Variables======================================================================

start_year <- 1988
end_year <- 2020

# Parameters to track
params <- c(
  # "sd_ind",
  "survival_af", 
  "survival_am", 
  "survival_ca",
  "s0_ps",
  "p0_ps"
)

# File names/paths
model_file <- "models//survival//cjs_model_elk_cip.txt"
result_file <- "results//survival//cjs_rslt_cip_31dec2021.Rdata"

# Sampler variables
n_i <- 25000
n_t <- 1
n_b <- 15000
n_c <- 3
n_a <- 1000



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

# y <- y[1:100,]
# f <- f[1:100]
# l <- l[1:100]
# w <- w[1:100,]
# male <- male[1:100]
# calf <- calf[1:100,]
# herd <- herd[1:100,]

#Fit_model======================================================================

z_data <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
l_k <- rep(NA, nrow(y))
for(i in 1:nrow(y)){
  l_k[i] <- max(which(y[i,] == 1))
  z_data[i, (f[i]):l_k[i]] <- 1
}

cjs_data <- list(
  y = y,
  z = z_data,
  f = f,
  l = l,
  m = male,
  c = calf,
  h = herd,
  n_occ = ncol(y),
  n_ind = nrow(y)
)

cjs_inits <- list(
  s0_ps  = runif(ncol(y), 0.5, 0.9),
  p0_ps  = runif(ncol(y), 0, 1),
  sm     = rnorm(ncol(y)),
  sc     = rnorm(ncol(y)),
  pm     = rnorm(ncol(y)),
  sh     = rnorm(1),
  ph     = rnorm(1),
  sd_ind = rnorm(1),
  b_ind  = rnorm(nrow(y))
)

cjs_nimble_model <- readBUGSmodel(
  model = model_file,
  data = cjs_data,
  inits = cjs_inits,
  calculate = F,
  check = F
)

cjs_rslt <- nimbleMCMC(
  model = cjs_nimble_model,
  monitors = params,
  niter = n_i,
  nburnin = n_b,
  nchains = n_c,
  thin = n_t
)

mcmcplots::mcmcplot(cjs_rslt)

save(cjs_rslt, file = result_file)