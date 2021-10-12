# Starkey elk CJS model fit workflow
# Kenneth Loonam
# March 2020

#Variables======================================================================

start_year <- 1988
end_year <- 2020

# Parameters to track
params <- c(
  "phi", 
  "p",
  "hnt",
  "tau_ind"
)

# Tau for logistic transformed priors
pr_p <- 0.5

# File names/paths
model_file <- "models//survival//survival_multi_state_model.txt"
result_file <- "results//survival//survival_multistate_result.Rdata"


# Sampler variables
n_i <- 20000
n_t <- 1
n_b <- 10000
n_c <- 3

#Environment====================================================================

require(tidyverse); require(nimble); require(mcmcplots)
load("data//elk_data.Rdata")

ch_init_fn <- function(ch, f){
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
  y[i,l[i]] <- 0
}

male <- elk_data$sex_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male) %>%
  magrittr::add(1)

calf <- elk_data$age_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  replace_na(0) %>%
  magrittr::subtract(2) %>%
  abs()

herd <- elk_data$hrd_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  mutate_all(funs(case_when(
    . == "main" ~ 1,
    . == "ne_e" ~ 2,
    . == "ne_w" ~ 3,
    . == "ne_o" ~ 4,
    . == "camp" ~ 5,
    . == "otsd" ~ 6,
    . == "hand" ~ 7
  ))) %>%
  as.matrix()

bad_elk <- which(f == l)
y <- y[-bad_elk,]
f <- f[-bad_elk]
l <- l[-bad_elk]
w <- w[-bad_elk,]
male <- male[-bad_elk]
calf <- calf[-bad_elk,]
herd <- herd[-bad_elk,]

jags_data <- list(
  y = y,
  f = f,
  l = l,
  w = w,
  male = male,
  calf = calf,
  herd = herd,
  n_ind = nrow(y),
  n_occ = ncol(y),
  pr_p = pr_p
)

z <- ch_init_fn(jags_data$y, jags_data$f)

inits <- function(){
  list(
    z = z
  )
}

#Fit_model======================================================================

# run the MCMC chain in nimble
model_object <- readBUGSmodel(
  model = model_file,
  data = jags_data,
  inits = list(z = z)
)

rslt <- nimbleMCMC(
  model = model_object,
  monitors = params,
  thin = n_t,
  niter = n_i,
  nburnin = n_b,
  nchains = n_c
)

summary(rslt)