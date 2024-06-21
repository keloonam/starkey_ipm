# Scratch work for testing if deterministic harvest pushes N upwards
#Environment====================================================================
require(tidyverse); require(mcmcplots); require(rjags)

#Variables======================================================================

ny <- 5
N1 <- c(50, 100)
Rm <- 0.25
Sm <- 0.9
Hm <- 0.1
Pm <- 0.5

#Simulate population without harvest============================================

# Simulate N
N <- matrix(NA, nrow = 2, ncol = ny)
N[,1] <- N1
for(t in 2:ny){
  N[1,t] <- rbinom(1, N[2,t-1], Rm)
  N[2,t] <- rbinom(1, sum(N[1:2,t-1]), Sm)
}

# Simulate cjs data - only "capturing" individuals alive at N1 for simplicity
y <- matrix(0, nrow = sum(N1), ncol = ny)
z <- y
z[,1] <- 1
for(i in 1:nrow(z)){
  for(t in 2:ncol(z)){
    z[i,t] <- rbinom(1, 1, Sm*z[i,t-1])
  }
}
y <- matrix(rbinom(length(z), 1, Pm*z), nrow = nrow(z), ncol = ncol(z))
unseen <- which(apply(y, 1, sum) == 0)
y <- y[-unseen,]
f <- apply(y, 1, which.max)
late_obs <- which(f == ny)
f <- f[-late_obs]; y <- y[-late_obs,]
z_init <- matrix(NA, nrow(y), ncol(y))
for(i in 1:nrow(y)){
  definitely_alive <- min(which(y[i,] == 1)):max(which(y[i,] == 1))
  z_init[i,definitely_alive] <- 1
}
# Simulate ratio data
rdt <- matrix(NA, nrow = ny, ncol = 2)
for(t in 2:nrow(rdt)){
  rdt[t,1] <- rbinom(1, N[1,t], Pm)
  rdt[t,2] <- rbinom(1, N[2,t-1], Pm)
}

# Simulate count data
c <- matrix(rpois(length(N), N), nrow = nrow(N), ncol = ncol(N))

#Build JAGS Objects=============================================================

jd <- list(
  rdt = rdt,
  y = y,
  c = c,
  ny = ny,
  ni = nrow(y),
  f = f,
  z = z_init
)
ji <- list(
  R = c(NA, rep(0.5, ny-1)),
  S = c(NA, rep(0.9, ny-1))
)
jp <- c("N", "S", "R", "P")

#Run the model==================================================================
jm <- jags.model(
  file = "models//misc//toy_ipm.txt",
  data = jd,
  inits = ji,
  n.chains = 3,
  n.adapt = 1000
)
update(jm, n.iter = 5000)

rslt_nh <- coda.samples(
  jm,
  variable.names = jp,
  n.iter = 25000,
  thin = 5
)

mcmcplot(rslt_nh)

trueN_nh <- N

#Variables======================================================================

ny <- 5
N1 <- c(50, 100)
Rm <- 0.5
Sm <- 0.9
Hm <- 0.25
Pm <- 0.5

#Simulate population without harvest============================================

# Simulate N
N <- matrix(NA, nrow = 2, ncol = ny)
N[,1] <- N1
nh <- rep(0, ny)
for(t in 2:ny){
  n2t <- rbinom(1, sum(N[1:2,t-1]), Hm)
  N[1,t] <- rbinom(1, N[2,t-1], Rm)
  N[2,t] <- rbinom(1, n2t, Sm)
  nh[t] <- n2t - N[2,t]
}

# Simulate cjs data - only "capturing" individuals alive at N1 for simplicity
y <- matrix(0, nrow = sum(N1), ncol = ny)
z <- y
z[,1] <- 1
for(i in 1:nrow(z)){
  for(t in 2:ncol(z)){
    z[i,t] <- rbinom(1, 1, Sm*z[i,t-1])
  }
}
y <- matrix(rbinom(length(z), 1, Pm*z), nrow = nrow(z), ncol = ncol(z))
unseen <- which(apply(y, 1, sum) == 0)
y <- y[-unseen,]
f <- apply(y, 1, which.max)
late_obs <- which(f == ny)
f <- f[-late_obs]; y <- y[-late_obs,]
z_init <- matrix(NA, nrow(y), ncol(y))
for(i in 1:nrow(y)){
  definitely_alive <- min(which(y[i,] == 1)):max(which(y[i,] == 1))
  z_init[i,definitely_alive] <- 1
}
# Simulate ratio data
rdt <- matrix(NA, nrow = ny, ncol = 2)
for(t in 2:nrow(rdt)){
  rdt[t,1] <- rbinom(1, N[1,t], Pm)
  rdt[t,2] <- rbinom(1, N[2,t-1], Pm)
}

# Simulate count data
c <- matrix(rpois(length(N), N), nrow = nrow(N), ncol = ncol(N))

#Build JAGS Objects=============================================================

jd <- list(
  rdt = rdt,
  y = y,
  c = c,
  ny = ny,
  ni = nrow(y),
  f = f,
  z = z_init,
  h = nh
)
ji <- list(
  R = c(NA, rep(0.5, ny-1)),
  S = c(NA, rep(0.9, ny-1))
)
jp <- c("N", "S", "R", "P")

#Run the model==================================================================
jm <- jags.model(
  file = "models//misc//toy_harvest_ipm.txt",
  data = jd,
  inits = ji,
  n.chains = 3,
  n.adapt = 1000
)
update(jm, n.iter = 5000)

rslt_wh <- coda.samples(
  jm,
  variable.names = jp,
  n.iter = 25000,
  thin = 5
)

mcmcplot(rslt_wh)

trueN_wh <- N


