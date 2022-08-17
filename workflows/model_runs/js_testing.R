require(nimble)

n <- 100
s <- 0.9
p <- 0.5
nocc <- 5
naug <- 200

# Sampler variables
ni <- 75000
nt <- 1
nb <- 50000
nc <- 1
na <- 1000

params <- c(
  "N.seen",
  "N",
  # "N.old",
  "mean.s",
  "mean.p",
  "s.b0",
  "p.b0"
)

alive <- matrix( 0, nrow = n, ncol = nocc)
y     <- matrix( 0, nrow = n, ncol = nocc)
z     <- matrix(NA, nrow = n, ncol = nocc)
f     <- sample(1:(nocc-1), replace = T, size = n)
drop  <- c()
for(i in 1:nrow(alive)){
  alive[i, f[i]] <- 1
  for(j in 2:ncol(alive)){
    if(alive[i,j-1] == 1){
      alive[i,j] <- rbinom(1, 1, s)
    }
  }
  for(j in 1:ncol(alive)){
    y[i,j] <- rbinom(1, 1, alive[i,j] * p)
  }
  if(sum(y[i,]) == 0){
    drop <- c(drop, i)
  }
}

y <- y[-drop,]
z <- replace(y, z == y, NA)
# y <- rbind(y, matrix( 0, nrow = naug, ncol = nocc))
# z <- rbind(z, matrix(NA, nrow = naug, ncol = nocc))


f <- apply(y, 1, function(x){min(which(x == 1))})


fake_js_data <- list(
  y = y,
  z = z
)

fake_js_constants <- list(
  nocc = ncol(y),
  nind = nrow(y),
  f    = f
)
jags_data <- list(
  y = y,
  z = z,
  nocc = ncol(y),
  nind = nrow(y),
  f = f
)

fake_js_inits <- list(
  z = y
)

source("models/survival/cjs_model_null.R")

rslt <- nimbleMCMC(
  code = js_nimble_code,
  constants = fake_js_constants,
  data = fake_js_data,
  # inits = fake_js_inits,
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
jags_model <- rjags::jags.model(
  file = "models/survival/cjs_model_null.txt",
  data = jags_data,
  n.chains = nc,
  n.adapt = na
  
)

# save(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
apply(alive, 2, sum)
