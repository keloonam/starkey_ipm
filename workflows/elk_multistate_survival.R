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

# Model
model_file <- "modells//survival//survival_multi_state_model.txt"

# Sampler variables
ni <- 100
nt <- 1
nb <- 20
nc <- 3

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
male <- elk_data$sex_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male)
calf <- elk_data$age_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()
herds <- elk_data$hrd_tib %>%
  filter(id %in% elk_data$cap_tib$id) %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()
bad_elk <- which(f == l)
y <- y[-bad_elk,]
f <- f[-bad_elk]
l <- l[-bad_elk]
male <- male[-bad_elk]
calf <- calf[-bad_elk,]
herds <- herds[-bad_elk,]
camp <- herds == "camp"
hand <- herds == "hand"
ne_e <- herds == "ne_e"
ne_o <- herds == "ne_o"
ne_w <- herds == "ne_w"
otsd <- herds == "otsd"


elk_cjs_data <- list(
  y = y,
  f = f,
  l = l,
  male = male,
  calf = calf,
  camp = camp,
  hand = hand,
  ne_e = ne_e,
  ne_o = ne_o,
  ne_w = ne_w,
  otsd = otsd,
  n_ind = nrow(y),
  n_occ = ncol(y),
  pr_p = pr_p
)

z <- ch_init_fn(elk_cjs_data$y, elk_cjs_data$f)

model_inits <- function(){
  list(
    z = z
  )
}

#Fit_model======================================================================

# run the MCMC chain in JAGS
rslt <- jags(
  data = elk_cjs_data, 
  inits = model_inits,
  parameters.to.save = params,
  model.file = model_file,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.thin = nt
)
mcmcplot(rslt)
save(rslt, file = "rand_effect.Rdata")
