# Starkey elk CJS model fit workflow
# Kenneth Loonam

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
model_file <- "models//survival//cjs_harvest_rate_model.txt"
result_file <- "results//survival//survival_harvest_rate_result.Rdata"


# Sampler variables
n_i <- 20
n_t <- 1
n_b <- 10
n_c <- 3

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
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

n_h <- array(data = 0, dim = c(2,2,7,ncol(y)))
for(i in 1:nrow(y)){
  for(t in 1:ncol(y)){
    n_h[calf[i,t],male[i],herd[i,t],t]<-n_h[calf[i,t],male[i],herd[i,t], t]+(w[i,t]==0)
  }
}

jags_data <- list(
  y = y,
  f = f,
  l = l,
  n_h = n_h,
  male = male,
  calf = calf,
  herd = herd,
  n_ind = nrow(y),
  n_occ = ncol(y),
  pr_p = pr_p
)

z <- ch_init_fn(jags_data$y, jags_data$f)

inits <- list(
    z = z,
    n_l = array(data = 0, dim = c(2,2,7,ncol(y),nrow(y)))
  )


#Fit_model======================================================================

mdl <- nimble::readBUGSmodel(
  model = model_file,
  data = jags_data,
  inits = inits
)

rslt <- nimble::nimbleMCMC(
  model = mdl,
  thin = n_t,
  niter = n_i,
  nburnin = n_b,
  nchains = n_c
)

# # run the MCMC chain in JAGS
# jgs_mdl <- jags.model(
#   file = model_file,
#   data = jags_data,
#   inits = inits,
#   n.chains = n_c
# )
# 
# # save(jgs_mdl, file = "temp_model_file.Rdata")
# 
# update(jgs_mdl, n.iter = n_b)
# 
# # save(jgs_mdl, file = "temp_model_file.Rdata")
# # load("temp_model_file.Rdata")
# 
# rslt <- coda.samples(
#   jgs_mdl,
#   variable.names = params,
#   n.iter = n_i,
#   thin = n_t
# )
# 
# save(rslt, file = result_file)