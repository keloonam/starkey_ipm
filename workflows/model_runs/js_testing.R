require(nimble)

# Sampler variables
ni <- 100000
nt <- 1
nb <- 15000
nc <- 3

#####FAKE DATA#####
y <- matrix(0,  nrow = 100, ncol = 6)
z <- matrix(NA, nrow = 100, ncol = 6)
f <- sample(1:5, nrow(y), replace = T)
age <- y

for(i in 1:nrow(y)){
  z[i,f[i]] <- 1
  age[i,f[i]] <- 1
  if(f[i] == 1){
    age[i,1] <- age[i,1] + rbinom(1, 6, 0.1)
  }
  for(j in (f[i] + 1):ncol(y)){
    z[i,j] <- rbinom(1, 1, z[i,j-1] * 0.8)
    y[i,j] <- rbinom(1, 1, 0.6)
    if(age[i,j-1] > 0){
      age[i,j] <- age[i,j-1] + 1
    } 
  }
}

male <- rbinom(nrow(y), 1, 0.5)
herd <- matrix(1, nrow = nrow(y)*3, ncol = ncol(y))

add_n <- nrow(y)*2
y <- rbind(y, matrix(0, nrow = add_n, ncol = ncol(y)))
z <- rbind(z, matrix(NA, nrow = add_n, ncol = ncol(y)))
age <- rbind(age, matrix(NA, nrow = add_n, ncol = ncol(y)))
male <- c(male, rep(NA, add_n))

for(i in 1:nrow(z)){
  for(j in 1:ncol(z)){
    if(z[i,j] == 0){
      z[i,j] <- NA
    }
  }
}


fake_js_data <- list(
  y = y,
  z = z
  # AGE = age,
  # A_P1 = age[,1] + 1,
  # H = herd,
  # c = age ==1
)

MAX_A <- max(age[,1], na.rm = T)
A <- rep(1, MAX_A)

fake_js_constants <- list(
  L     = rep(ncol(y), nrow(y)),
  nocc  = ncol(y),
  naug  = nrow(y),
  R     = rep(1, ncol(y))
  # MAX_A = MAX_A,
  # A     = A
)

#####END FAKE DATA#####

js_inits <- list(
  s0   = runif(1, 0.5, 0.9),
  p0   = runif(1, 0, 1),
  e0   = runif(1, 0, 1)
  # a0   = rep(1/fake_js_constants$MAX_A, fake_js_constants$MAX_A),
  # s.ma = rnorm(ncol(y)),
  # s.ca = rnorm(ncol(y)),
  # s.he = rnorm(ncol(y)),
  # p.ma = rnorm(ncol(y)),
  # p.he = rnorm(ncol(y))
)

# jsnm <- list(
#   constants = js_constants,
#   data = js_data,
#   monitors = params
# )

# save(jsnm, file = "cqls//js//js_nimble_data_cqls.Rdata")

# Load model as js_superpop_model
source("models/survival/js_model_null.R")

params <- c(
  "N_super",
  "N",
  "survival_af",
  "detection_f"
)

rslt <- nimbleMCMC(
  code = js_nimble_code,
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
