# Starkey elk CJS model fit workflow
# Kenneth Loonam

#Variables======================================================================

start_year <- 1988
end_year <- 2021

# Parameters to track
params <- c(
  "sf", 
  "sm", 
  "sc",
  "pf",
  "pm"
)

# File names/paths
result_file <- "results//survival//cjs_rslt_02sep2022.Rdata"

# Sampler variables
ni <- 5000
nt <- 1
nb <- 1000
nc <- 1
na <- 1000

#Environment====================================================================

require(tidyverse); require(mcmcplots); require(nimble)
load("data//elk_data.Rdata")

#Prep_data======================================================================

y <- elk_data$cap_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()

l <- elk_data$hnt_tib %>%
  arrange(id) %>%
  mutate('2021' = 1) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  apply(., 1, function(x) min(which(x != 0)))

w <- elk_data$hnt_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() 

y <- w + y # harvests count as observed alive that year (did not die naturally)

male <- elk_data$sex_tib %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male)

calf <- elk_data$age_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()

herd <- elk_data$hrd_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  mutate_all(funs(case_when(
    . == "main" ~ 0,
    T ~ 1
  ))) %>%
  as.matrix()

z <- elk_data$liv_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()

# We don't care about elk that were never in main study
gone_elk <- apply(herd, 1, sum) == ncol(herd)
unseen_elk <- apply(y, 1, sum) == 0
bad_elk <- which((gone_elk + unseen_elk) != 0)

y <- y[-bad_elk,]
l <- l[-bad_elk]
w <- w[-bad_elk,]
male <- male[-bad_elk]
calf <- calf[-bad_elk,]
herd <- herd[-bad_elk,]
calf <- (calf == 1) * 1
calf[is.na(calf)] <- 0

f <- apply(y, 1, function(x) min(which(x != 0)))
weird_elk <- f == l

y <- y[-weird_elk,]
l <- l[-weird_elk]
w <- w[-weird_elk,]
male <- male[-weird_elk]
calf <- calf[-weird_elk,]
herd <- herd[-weird_elk,]

calf <- (calf == 1) * 1
calf[is.na(calf)] <- 0

#Fit_model======================================================================

data <- list(
  y = y,
  z = z,
  m = male,
  h = herd,
  c = calf
)
constants <- list(
  l     = l,
  f     = f,
  nocc  = ncol(y),
  nind  = nrow(y)
)

nocc <- ncol(y)
inits <- list(
  sf = runif(nocc, 0, 1),
  sm = runif(nocc, 0, 1),
  sc = runif(nocc, 0, 1),
  pf = runif(nocc, .9, 1),
  pm = runif(nocc, .9, 1),
  sh = runif(1, 0, 1),
  ph = runif(1, 0, 1)
)

source("models/survival/cjs_elk.R")

rslt <- nimbleMCMC(
  code              = code,
  constants         = constants,
  data              = data,
  monitors          = params,
  inits             = inits,
  niter             = ni,
  nburnin           = nb,
  nchains           = nc,
  progressBar       = T,
  samplesAsCodaMCMC = T,
  summary           = F,
  check             = T
)  

# save(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
