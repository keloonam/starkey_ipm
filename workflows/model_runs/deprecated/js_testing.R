require(nimble)

n    <- 100
s.f  <- 0.9
s.m  <- 0.7
s.c  <- 0.5
p    <- 0.8
nocc <- 5
naug <- 300

# Sampler variables
ni <- 60000
nt <- 1
nb <- 30000
nc <- 3
na <- 1000

params <- c(
  "sf",
  "sm",
  "sc",
  "p.f",
  "p.m",
  "ep",
  "Nf",
  "Nm",
  "Nt",
  "Ns"
)

alive <- matrix( 0, nrow = n, ncol = nocc)
y     <- matrix( 0, nrow = n, ncol = nocc)
c     <- matrix( 0, nrow = n, ncol = nocc)
s     <- matrix(NA, nrow = n, ncol = nocc)
f     <- sample(1:(nocc-1), replace = T, size = n)
m     <- rbinom(n, 1, 0.5)
drop  <- c()
for(i in 1:nrow(alive)){
  alive[i, f[i]] <- 1
  c[i, f[i]]     <- 1
  
  for(j in 2:ncol(alive)){
    if(c[i,j-1] == 1){
      s[i,j] <- s.c
    }else{
      if(m[i] == 1){
        s[i,j] <- s.m
      }else{
        s[i,j] <- s.f
      }
    }
    if(alive[i,j-1] == 1){
      alive[i,j] <- rbinom(1, 1, s[i,j])
    }
  }
  for(j in 1:ncol(alive)){
    y[i,j] <- rbinom(1, 1, alive[i,j] * p)
  }
  if(sum(y[i,]) == 0){
    drop <- c(drop, i)
  }
}

n.total   <- apply(alive, 2, sum)
n.females <- apply(alive * (1-m), 2, sum)
n.males   <- apply(alive * m, 2, sum)

y <- y[-drop,]
c <- c[-drop,]
m <- m[-drop ]

f <- apply(y, 1, function(x){min(which(x == 1))})
l <- apply(y, 1, function(x){max(which(x == 1))})
w <- matrix(0, nrow = nrow(y), ncol = ncol(y)) # known state alive

for(i in 1:nrow(y)){
  for(j in f[i]:l[i]){
    w[i,j] <- 1
  }
  if(l[i] < ncol(y)){
    for(j in (l[i]+1):ncol(y)){
      w[i,j] <- NA
    }
  }
}
y <- cbind(rep(0, nrow(y)), y)
ca <- cbind(c, rep(0, nrow(c)))
w <- cbind(rep(0, nrow(w)), w)
ca <- ca * y

# z:    1 -- not yet alive;    2 -- adult;    2 -- dead
z <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
for(i in 1:nrow(z)){
  for(t in 2:ncol(z)){
    if(!is.na(w[i,t])){
      if(w[i,t] == 1){
        z[i,t] <- 2
      }
    }
  }
}

y <- rbind(y, matrix(0,  nrow = naug, ncol = ncol(y)))
z <- rbind(z, matrix(NA, nrow = naug, ncol = ncol(z)))
ca <- rbind(ca, matrix(0,  nrow = naug, ncol = ncol(ca)))
m <- c(m, rep(NA, naug))

for(i in 1:nrow(z)){
  if(any(!is.na(z[i,]))){
    if(any(c[i,] == 1)){
      tmp <- which(c[i,] == 1) - 1
      z[i,1:tmp] <- 1
    }
  }
}

z_init <- z
z_init[,1] <- 1
for(i in 1:nrow(z)){
  for(t in 1:ncol(z)){
    if(is.na(z_init[i,t])){
      z_init[i,t] <- 2
    }
  }
}
for(i in 1:nrow(y)){
  for(t in 1:ncol(y)){
    if(y[i,t] == 0){
      y[i,t] <- 2
    }
  }
}
h <- matrix(0, nrow = nrow(y), ncol = ncol(y))
z[,1] <- NA
avail <- matrix(1, nrow = nrow(ca), ncol = ncol(ca))
for(i in 1:nrow(ca)){
  if(any(ca[i,] == 1)){
    calf_sess <- which(ca[i,] == 1)
    if(calf_sess > 1){
      avail[i,1:calf_sess] <- 0
    }
  }
}
data <- list(
  y = y,
  m = m,
  z = z,
  c = ca,
  h = h,
  avail = avail
)
nocc <- ncol(y)
constants <- list(
  nocc = ncol(y),
  nind = nrow(y)
  # l    = rep(ncol(y), nrow(y))
)
jags.data <- list(
  y = y,
  m = m,
  h = h,
  z = z,
  c = ca,
  avail = avail,
  # z_con = matrix(1, nrow = nrow(z), ncol = ncol(z)),
  nocc = ncol(y),
  nind = nrow(y)
)
z_init[,1] <- NA

inits <- list(
  sf = runif(nocc-1, 0, 1),
  sm = runif(nocc-1, 0, 1),
  sc = runif(nocc-1, 0, 1),
  pf = runif(nocc-1, 0, 1),
  pm = runif(nocc-1, 0, 1),
  # pc = runif(nocc-1, 0, 1),
  # ec = runif(nocc-1, 0, 1),
  # ea[1] = runif(1, 0, 1),
  sh = runif(1, 0, 1),
  ph = runif(1, 0, 1)
  # z  = z_init
)

source("models/survival/js_multistate_elk_nocalves_irandp.R")

rslt <- nimbleMCMC(
  code = code,
  constants = constants,
  data = data,
  inits = inits,
  monitors = params,
  niter = ni,
  nburnin = nb,
  nchains = nc,
  progressBar = T,
  samplesAsCodaMCMC = T,
  summary = F,
  check = F
)

# require(rjags)
# model <- jags.model(file = "models//survival//js_multistate_elk_jags.txt",
#            data = jags.data,
#            inits = inits,
#            n.chains = nc,
#            n.adapt = na)
# update(model, n.iter = nb)
# rslt <- coda.samples(model, variable.names = params, n.iter = ni)

# save(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
# apply(alive, 2, sum)
# summary(rslt[,c(2:5, 12:15)])
