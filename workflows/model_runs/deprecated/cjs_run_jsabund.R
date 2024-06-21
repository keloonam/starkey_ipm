# Starkey elk CJS model fit workflow
# Kenneth Loonam

#Variables======================================================================

start_year <- 1988
end_year <- 2021

# Parameters to track
params <- c(
  # "Nf",
  # "Nm",
  # "m.a",
  # "f.a",
  # "ntot.a",
  "s.f",
  "s.m",
  "s.c",
  "p.f",
  "p.m"
)

# File names/paths
model_file <- "models//survival//cjs_elk.R"
result_file <- "results//survival//cjs_rslt_24aug2022.Rdata"

# Sampler variables
ni <- 35000
nt <- 1
nb <- 10000
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
rm <- which(apply(y, 1, sum) == 0)
y  <- y[-rm,]

f <- apply(y, 1, function(x) min(which(x != 0)))

l <- elk_data$hnt_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  mutate('2021' = 1) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  apply(., 1, function(x) min(which(x != 0)))
l <- l[-rm]

w <- elk_data$hnt_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() 
w <- w[-rm,]

y <- w + y # harvests count as observed alive that year (did not die naturally)

male <- elk_data$sex_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male)
male <- male[-rm]

age <- elk_data$age_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  replace_na(0) 
calf <- age == 1
calf <- calf[-rm,]

herd <- elk_data$hrd_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  mutate_all(funs(case_when(
    . == "main" ~ 0,
    T ~ 1
  ))) %>%
  as.matrix()
herd <- herd[-rm,]

gone_elk <- apply(1 - herd, 1, sum) == ncol(herd)
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

z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
l_k <- rep(NA, nrow(y))
for(i in 1:nrow(y)){
  l_k[i] <- max(which(y[i,] == 1))
  z[i, (f[i]):l_k[i]] <- 1
  z[i, (l_k[i]):l[i]] <- NA
}

z_con <- matrix(1, nrow = nrow(z), ncol = ncol(z))

Nf.oa <- apply(y * (1 - calf) * (1 - male), 2 , sum)[-1]
Nm.oa <- apply(y * (1 - calf) *      male,  2 , sum)[-1]

cjs_data <- list(
  y     = y,
  z     = z,
  z_con = z_con
  # yf    = Nf.oa,
  # ym    = Nm.oa
)

cjs_constants <- list(
  f     = f,
  l     = l,
  nocc  = ncol(y),
  nind  = nrow(y),
  m     = male,
  c     = 1 * calf,
  h     = herd
)

cjs_inits <- list(
  s0  = rnorm(ncol(y)),
  p0  = rnorm(ncol(y)),
  sm  = rnorm(ncol(y)),
  sc  = rnorm(ncol(y)),
  pm  = rnorm(ncol(y)),
  sh  = rnorm(1),
  ph  = rnorm(1),
  z   = z
)

source(model_file)

rslt <- nimbleMCMC(
  code      = nimble_code,
  constants = cjs_constants,
  data      = cjs_data,
  inits     = cjs_inits,
  monitors  = params,
  niter     = ni,
  nburnin   = nb,
  nchains   = nc,
  thin      = nt,
  check     = F,
  samplesAsCodaMCMC = T
)

mcmcplots::mcmcplot(rslt)

# save(cjs_rslt, file = result_file)
cnames <- dimnames(rslt)[[2]]
pm.mat <- rslt[,grep("p.m", cnames)]
pf.mat <- rslt[,grep("p.f", cnames)]



Nf.exp <- apply(Nf.oa/pf.mat, 2, mean)
Nm.exp <- apply(Nm.oa/pf.mat, 2, mean)

out_file <- list(
  samples = rslt,
  exp.f   = Nf.exp,
  exp.n   = Nm.exp
)

save(cjs_rslt, file = result_file)