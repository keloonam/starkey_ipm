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
# for(i in 1:nrow(calf)){
#   if(sum(calf[i,], na.rm = T) > 0){
#     calf_sess <- min(which(calf[i,] == 1))
#     calf[i,calf_sess:ncol(calf)] <- 1
#   }
#   if(all(calf[i,] == 0)){
#     calf[i,1] <- 2
#     calf[i,2:ncol(calf)] <- 1
#   }
#   if(any(!is.na(calf[i,]))){
#     start_age <- min(which(!is.na(calf[i])))
#     calf[i,start_age:ncol(calf)] <- cumsum(calf[i,start_age:ncol(calf)])
#   }
# }

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
calf[is.na(calf)] <- 0
#Switch_to_Multistate_Model_Data================================================

y[y == 0] <- 2
y <- cbind(rep(2, nrow(y)), y)

l <- l + 1



herd <- cbind(herd, rep(0, nrow(herd)))

z[z == 1] <- 2
calf_sess <- apply(calf, 1, function(x)min(which(x == 1)))
calf_sess[calf_sess == Inf] <- NA
last_sess <- apply(z,    1, function(x) max(which(x==2)))
last_sess[last_sess == Inf | last_sess == -Inf] <- NA
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
z[z == 0] <- NA
z <- cbind(rep(NA, nrow(z)), z)
calf <- cbind(calf, rep(0, nrow(calf)))

#Augment========================================================================

n_aug <- nrow(y) * aug
y <- rbind(y, matrix(2,  nrow = n_aug, ncol = ncol(y)))
z <- rbind(z, matrix(NA, nrow = n_aug, ncol = ncol(y)))
l <- c(l, rep(ncol(y), n_aug))
male <- c(male, rep(NA, n_aug))
herd <- rbind(herd, matrix(0, nrow = n_aug, ncol = ncol(y)))
calf <- rbind(calf, matrix(0, nrow = n_aug, ncol = ncol(y)))

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
  nocc  = ncol(y),
  nind  = nrow(y)
)

nocc <- ncol(y)
inits <- list(
  sf = runif(nocc-1, 0, 1),
  sm = runif(nocc-1, 0, 1),
  sc = runif(nocc-1, 0, 1),
  pf = runif(nocc-1, 0, 1),
  pm = runif(nocc-1, 0, 1),
  sh = runif(1, 0, 1),
  ph = runif(1, 0, 1)
)

source("models/survival/js_multistate_elk_nocalves.R")

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
