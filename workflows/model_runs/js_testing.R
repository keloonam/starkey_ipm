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
  "N",
  # "z",
  # "N.old",
  "N.new",
  "expN.new",
  "sdN.new",
  "N.alive",
  "mean.s",
  "mean.p",
  "s.b0",
  "p.b0",
  "N.bin"
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
z <- replace(y, y == 0, NA)
# w <- rep(1, nrow(y))
# y <- rbind(y, matrix( 0, nrow = naug, ncol = nocc))
# z <- rbind(z, matrix(NA, nrow = naug, ncol = nocc))
# 
# w <- c(w, rep(NA, naug))

f <- apply(y, 1, function(x){min(which(x == 1))})
l <- apply(y, 1, function(x){max(which(x == 1))})
z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
for(i in 1:nrow(y)){
  for(j in f[i]:l[i]){
    z[i,j] <- 1
  }
  if(l[i] < ncol(y)){
    for(j in (l[i]+1):ncol(y)){
      z[i,j] <- NA
    }
  }
}
newN <- rep(0, ncol(y))
for(i in 1:ncol(y)){
  newN[i] <- sum(f == i)
}


fake_js_data <- list(
  y = y,
  # w = w,
  z = z,
  z.constrain = matrix(1, nrow = nrow(z), ncol = ncol(z)),
  newN = newN
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
  z = replace(z, is.na(z), 1)
)

source("models/survival/cjs_model_null.R")

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

# save(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
apply(alive, 2, sum)
summary(rslt[,c(2:5, 12:15)])
