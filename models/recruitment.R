code <- nimbleCode({
  # Priors
  for(t in 1:n_years){
    R[t] ~ dunif(0,1)
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