require(nimble)

code <- nimbleCode({
  # Model
  for (i in 1:BATCHES) {
    for (j in 1:SAMPLES) {
      y[i,j] ~ dnorm(mu[i], sd = sigma.within);
    }
    mu[i] ~ dnorm(theta, sd = sigma.between);
  }
  
  # Priors
  theta ~ dnorm(0.0, 1.0E-10);
  sigma.within ~ dunif(0, 100)
  sigma.between ~ dunif(0, 100)
})

model <- nimbleModel(code, constants = list(BATCHES = 6, SAMPLES = 5))
data <- matrix(c(1545, 1540, 1595, 1445, 1595, 1520, 1440, 1555, 1550, 
                 1440, 1630, 1455, 1440, 1490, 1605, 1595, 1515, 1450, 1520, 1560, 
                 1510, 1465, 1635, 1480, 1580, 1495, 1560, 1545, 1625, 1445), nrow = 6)

mcmc_conf <- configureMCMC(model, monitors = c("theta"), print = T, enableWAIC = F)
mcmc <- buildMCMC(model)
comp_model <- compileNimble(model)
c_mcmc <- compileNimble(mcmc, project = model)

samples <- runMCMC(c_mcmc, niter = 1000, nburnin = 500, thin = 1, nchains = 3)
dplyr::bind_rows(lapply(samples, dplyr::as_tibble))
lapply(lapply(samples, as.mcmc), summary)
