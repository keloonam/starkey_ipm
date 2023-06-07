nimble_code <- nimble::nimbleCode({
  
  # Priors
  b0_mean ~ dlogis(0, 1)
  b0_sd   ~ T(dnorm(0, sd = 10), 0, 10)
  b_lact_mean ~ dlogis(0, 1)
  b_lact_sd   ~ T(dnorm(0, sd = 10), 0, 10)
  for(t in 1:ny){
    b0[t]     ~ dnorm(b0_mean, sd = b0_sd)
    b_lact[t] ~ dnorm(b_lact_mean, sd = b_lact_sd)
  }
  b_dens_open ~ dlogis(0, 1)
  b_pdsi_open ~ dlogis(0, 1)
  b_dens_lact ~ dlogis(0, 1)
  b_pdsi_lact ~ dlogis(0, 1)
  
  # GLM
  for(i in 1:nr){
    logit(p[i]) <- b0[yr[i]] + 
      b_lact[yr[i]] * lac[i] +
      b_dens_open   * elk[i] * (1 - lac[i]) +
      b_pdsi_open   * pdi[i] * (1 - lac[i]) +
      b_dens_lact   * elk[i] *      lac[i]  +
      b_pdsi_lact   * pdi[i] *      lac[i]

  # Model  
    preg[i] ~ dbin(p[i], nobs[i])
  }
}
)
