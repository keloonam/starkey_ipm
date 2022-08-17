require(nimble)

n <- 100
s <- 0.9
p <- 0.5
nocc <- 5
naug <- 200

# Sampler variables
ni <- 25000
nt <- 1
nb <- 15000
nc <- 1

params <- c(
  "N_super",
  "N",
  "survival_af",
  "detection_f",
  "s.b0",
  "p.b0",
  "e"
)

alive <- matrix( 0, nrow = n, ncol = nocc)
y     <- matrix( 0, nrow = n, ncol = nocc)
z     <- matrix(NA, nrow = n, ncol = nocc)
f     <- sample(1:(nocc-1), replace = T, size = n)
for(i in 1:nrow(alive)){
  alive[i, f[i]] <- 1
  for(j in 2:ncol(alive)){
    if(alive[i,j-1] == 1){
      alive[i,j] <- rbinom(1, 1, s)
    }
    y[i,j] <- rbinom(1, 1, alive[i,j] * p)
    if(y[i,j] == 1){
      z[i,j] <- 1
    }
  }
}

y <- rbind(y, matrix( 0, nrow = naug, ncol = nocc))
z <- rbind(z, matrix(NA, nrow = naug, ncol = nocc))

fake_js_data <- list(
  y = y,
  z = z
)

fake_js_constants <- list(
  nocc  = ncol(y),
  naug  = nrow(y),
  R     = rep(1, ncol(y))
)

fake_js_inits <- list(
  p.b0 = -4,
  s.b0 = 4,
  e    = .5
)

source("models/survival/js_model_null.R")

rslt <- nimbleMCMC(
  code = js_nimble_code,
  constants = fake_js_constants,
  data = fake_js_data,
  inits = fake_js_inits,
  monitors = params,
  niter = ni,
  nburnin = nb,
  nchains = nc,
  progressBar = T,
  samplesAsCodaMCMC = T,
  summary = F,
  check = F
)  

# y > z

# save(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
