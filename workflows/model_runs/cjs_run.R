# Starkey elk CJS model fit workflow
# Kenneth Loonam

#Variables======================================================================

start_year <- 1988
end_year <- 2021

# Parameters to track
params <- c(
  "survival_af", 
  "survival_am", 
  "survival_ca",
  "prob_af",
  "prob_am"
)

# File names/paths
result_file <- "results//survival//cjs_rslt_13sep2022.Rdata"

# Sampler variables
ni <- 25000
nt <- 1
nb <- 10000
nc <- 3
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
# had a replace_na(0) line appended

herd <- elk_data$hrd_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  mutate_all(funs(case_when(
    . == "main" ~ 0,
    T ~ 1
  ))) %>%
  as.matrix()

# z <- elk_data$liv_tib %>%
#   arrange(id) %>%
#   select(as.character(start_year:end_year)) %>%
#   as.matrix()
# was built manually lower down

# We don't care about elk that were never in main study
gone_elk <- apply(herd, 1, sum) == ncol(herd)
# this was set up to not filter anything?
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

z <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
l_k <- rep(NA, nrow(y))
for(i in 1:nrow(y)){
  l_k[i] <- max(which(y[i,] == 1))
  z[i, (f[i]):l_k[i]] <- 1
}

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
  pf = runif(nocc, 0, 1),
  pm = runif(nocc, 0, 1),
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
  summary           = F,
  check             = F
)  

save(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
