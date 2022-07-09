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

js_superpop_model <- nimbleCode({
  # Nimble version of jolly-seber super population formulation with age data
  
  # Priors
  psi  ~ dbeta(1, 1)    # probability animals is in superpopulation
  p0   ~ dbeta(1, 1)    # mean capture probability
  phi0 ~ dbeta(1, 1)    # prior for survival at theoretical age 0
  alpha0 <- logit(phi0) 
  
  # starting age distribution: not yet recruited (age = 0), or age
  piAGE[1:max.age] ~ ddirch(a[1:max.age])
  # age distribution conditioned on alive at t=1
  piAGEuncond[1:(max.age + 1)] <- c((1 - eta[1]), eta[1] * piAGE[1:max.age])
  
  # recruitment rate from 'not entered population' at t 
  beta[1:K] ~ ddirch(b[1:K])
  eta[1] <- beta[1]
  for(k in 2:K){
    eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
  }
  
  # Likelihoods 
  for (i in 1:M){
    # is individual i real?
    w[i] ~ dbern(psi)
    
    # initial ages
    agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) # where agePlusOne are data
    age[i,1] <- (agePlusOne[i]-1) 
    # I think age+1 is solving the age = zero indexing issue??
    
    # state process
    u[i,1] <- step(age[i,1]-.1) # sets u to zero if age is zero
    z[i,1] <- u[i,1]*w[i] # z is the "real" state
    
    # Observation process
    y[i,1] ~ dbern(z[i,1]*p)
    
    # derived stuff
    avail[i,1] <- 1 - u[i,1] # still available -- i.e. not yet recruited
    
    # for occasions > 1     
    for (t in 2:K){
      
      # State process
      u[i,t] ~ dbern(u[i,t-1]*phi[i,t] + avail[i,t-1]*eta[t])   
      logit(phi[i,t]) <- alpha0 
      z[i,t] <- u[i,t]*w[i]
      
      # Age process
      age[i,t] <- age[i,t-1] + max(u[i,1:t]) 
      # ages by one year after recruitment (NIMBLE allows this syntax)
      
      # Observation process
      y[i,t] ~ dbern(z[i,t]*p)
      
      # derived stuff
      avail[i,t] <- 1- max(u[i,1:t]) # still available -- i.e. not yet recruited
    } #t
  } #i
  
  ## Derived population level stuff
  # Annual abundance
  for (t in 1:K){
    N[t] <- sum(z[1:M,t])               
  } #t
  
  # Annual growth rate
  for (t in 1:(K-1)){
    lambda[t] <- N[t+1]/N[t]               
  } #t
  
  # Superpopulation size
  Nsuper <- sum(w[1:M])       
  
})# end model

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
