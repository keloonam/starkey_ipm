# Starkey elk CJS model fit workflow
# Kenneth Loonam

#Variables======================================================================

# amount to MULTIPLY the population by for augmentation
aug <- 2

start_year <- 1988
end_year <- 2021

# Parameters to track
params <- c(
  "sf", 
  "sm", 
  "sc",
  "pf",
  "pm",
  "pc",
  "Nf",
  "Nm",
  "Nt",
  "Ns"
)

# File names/paths
result_file <- "results//survival//js_rslt_25aug2022.Rdata"

# Sampler variables
ni <- 10000
nt <- 1
nb <- 5000
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

z <- elk_data$liv_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()

# We don't care about elk that were never in main study
gone_elk <- which(apply(herd, 1, sum) == ncol(herd))
# weird_elk <- f == l
# bad_elk <- which((gone_elk + weird_elk) != 0)
y <- y[-gone_elk,]
l <- l[-gone_elk]
w <- w[-gone_elk,]
male <- male[-gone_elk]
calf <- calf[-gone_elk,]
herd <- herd[-gone_elk,]
calf <- (calf == 1) * 1
#Switch_to_Multistate_Model_Data================================================

y[y == 0] <- 2
y <- cbind(rep(2, nrow(y)), y)

l <- l + 1



herd <- cbind(herd, rep(0, nrow(herd)))

z[z == 1] <- 2
calf_sess <- apply(calf, 1, function(x)min(which(x == 1)))
calf_sess[calf_sess == Inf] <- NA
last_sess <- apply(z,    1, function(x) max(which(x==2)))
last_sess[last_sess == Inf] <- NA
for(i in 1:nrow(z)){
  if(!is.na(calf_sess[i])){
    z[i,1:calf_sess[i]] <- 1
  }
  if(!is.na(last_sess[i])){
    for(j in last_sess[i]:ncol(z)){
      if(!is.na(z[i,j])){
        if(z[i,j] == 0){
          z[i,j] <- 3
        }
      }
    }
  }
}

calf <- cbind(calf, rep(0, nrow(calf)))

#Augment========================================================================

n_aug <- nrow(y) * aug
y <- rbind(y, matrix(0,  nrow = n_aug, ncol = ncol(y)))
z <- rbind(z, matrix(NA, nrow = n_aug, ncol = ncol(y)))
l <- c(l, rep(ncol(y), n_aug))
male <- c(male, rbinom(n_aug, 1, 0.5))
herd <- rbind(herd, matrix(0,  nrow = n_aug, ncol = ncol(y)))
calf <- rbind(calf, matrix(NA, nrow = n_aug, ncol = ncol(y)))

#Fit_model======================================================================

js_data <- list(
  y    = y,
  z    = z,
  AGE  = calf,
  A_P1 = calf[,1] + 1,
  M    = male,
  H    = herd,
  c    = calf
)
js_constants <- list(
  L     = l,
  K     = ncol(y),
  NAUG  = nrow(y),
  MAX_A = max(calf[,1], na.rm = T),
  A     = rep(1, max(calf[,1], na.rm = T))
)

#####FAKE DATA#####
y <- matrix(0,  nrow = 100, ncol = 6)
z <- matrix(NA, nrow = 100, ncol = 6)
f <- sample(1:5, nrow(y), replace = T)
age <- y

for(i in 1:nrow(y)){
  z[i,f[i]] <- 1
  for(j in (f[i] + 1):ncol(y)){
    z[i,j] <- rbinom(1, 1, z[i,j-1] * 0.8)
    y[i,j] <- rbinom(1, 1, 0.6)
    age[i,j] <- sum(z[i,1:(j-1)], na.rm = T)
  }
}

male <- rbinom(nrow(y), 1, 0.5)
herd <- matrix(1, nrow = nrow(y)*3, ncol = ncol(y))

add_n <- nrow(y)*2
y <- rbind(y, matrix(0, nrow = add_n, ncol = ncol(y)))
z <- rbind(z, matrix(NA, nrow = add_n, ncol = ncol(y)))
age <- rbind(age, matrix(NA, nrow = add_n, ncol = ncol(y)))
male <- c(male, rep(NA, add_n))


fake_js_data <- list(
  y = y,
  z = z,
  AGE = age,
  A_P1 = age + 1,
  M = male,
  H = herd,
  c = age ==1
)

fake_js_constants <- list(
  L     = rep(ncol(y), nrow(y)),
  K     = ncol(y),
  NAUG  = nrow(y),
  MAX_A = 3,
  A     = rep(1, 2)
)

#####END FAKE DATA#####

js_inits <- list(
  s0   = runif(ncol(y), 0.5, 0.9),
  p0   = runif(ncol(y), 0, 1),
  e0   = runif(1, 0, 1),
  a0   = rep(1/fake_js_constants$MAX_A, fake_js_constants$MAX_A),
  s.ma = rnorm(ncol(y)),
  s.ca = rnorm(ncol(y)),
  s.he = rnorm(ncol(y)),
  p.ma = rnorm(ncol(y)),
  p.he = rnorm(ncol(y))
)

jsnm <- list(
  constants = js_constants,
  data = js_data,
  monitors = params
)

# save(jsnm, file = "cqls//js//js_nimble_data_cqls.Rdata")

# Load model as js_superpop_model
source("models/survival/js_model_nimble.R")

rslt <- nimbleMCMC(
  code = js_model_code,
  constants = fake_js_constants,
  data = fake_js_data,
  monitors = params,
  # inits = js_inits,
  niter = ni,
  nburnin = nb,
  nchains = nc,
  progressBar = T,
  samplesAsCodaMCMC = T,
  summary = F,
  check = F
)  

# save(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
