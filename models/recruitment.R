code <- nimbleCode({
  # Priors
  b0 ~ dlogis(0, 1)
  bt ~ dlogis(0, 1)
  bm ~ dlogis(0, 1)
  sdr ~ dunif(0, 5)
  for(t in 1:n_years){
    er[t] ~ dnorm(0, sd = sdr)
    logit(R[t]) <- b0 + bt*xt[t] + bm*xm[t] + er[t]
  } # t
  
  # Model
  for(t in 1:n_years){
    n_calf[t] ~ dbinom(R[t], n_cow[t])
  } # t
  
  # Goodness-of-fit
  for(t in 1:n_years){
    n_cnew[t] ~ dbinom(R[t], n_cow[t])
    c_diff[t] <- n_calf[t] - n_cnew[t]
  }
  mn_c_diff <- mean(c_diff[1:n_years])
  sd_c_diff <- sd(n_cnew[1:n_years]) - sd(n_calf[1:n_years])
})