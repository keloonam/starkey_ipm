nimble_code <- nimble::nimbleCode({
  
  # Priors
  b0_mean        ~ dlogis(0, 1)
  b_yng_mean     ~ dlogis(0, 1)
  b_old_mean     ~ dlogis(0, 1)
  b_lac_prm_mean ~ dlogis(0, 1)
  b_lac_yng_mean ~ dlogis(0, 1)
  b_lac_old_mean ~ dlogis(0, 1)
  b0_sd          ~ T(dnorm(0, sd = 10), 0, 10)
  b_yng_sd       ~ T(dnorm(0, sd = 10), 0, 10)
  b_old_sd       ~ T(dnorm(0, sd = 10), 0, 10)
  b_lac_prm_sd   ~ T(dnorm(0, sd = 10), 0, 10)
  b_lac_yng_sd   ~ T(dnorm(0, sd = 10), 0, 10)
  b_lac_old_sd   ~ T(dnorm(0, sd = 10), 0, 10)
  logit(preg_yng_mean) <- b0_mean + b_yng_mean
  for(t in 1:ny){
    b0[t]        ~ dnorm(b0_mean,        sd = b0_sd)
    b_yng[t]     ~ dnorm(b_yng_mean,     sd = b_yng_sd)
    b_old[t]     ~ dnorm(b_old_mean,     sd = b_old_sd)
    b_lac_prm[t] ~ dnorm(b_lac_prm_mean, sd = b_lac_prm_sd)
    b_lac_yng[t] ~ dnorm(b_lac_yng_mean, sd = b_lac_yng_sd)
    b_lac_old[t] ~ dnorm(b_lac_old_mean, sd = b_lac_old_sd)
    logit(preg_young[t]) <- b0[t] + b_yng[t]
    logit(a_preg_lac_prime[t]) <- b0[t] + b_lac_prm[t]
  }
  bden_dry_prm ~ dlogis(0, 1)
  bden_dry_yng ~ dlogis(0, 1)
  bden_dry_old ~ dlogis(0, 1)
  bpdi_dry_prm ~ dlogis(0, 1)
  bpdi_dry_yng ~ dlogis(0, 1)
  bpdi_dry_old ~ dlogis(0, 1)
  bden_lac_prm ~ dlogis(0, 1)
  bden_lac_yng ~ dlogis(0, 1)
  bden_lac_old ~ dlogis(0, 1)
  bpdi_lac_prm ~ dlogis(0, 1)
  bpdi_lac_yng ~ dlogis(0, 1)
  bpdi_lac_old ~ dlogis(0, 1)
  
  # GLM
  for(i in 1:nr){
    logit(p[i]) <- b0[yr[i]] + 
      b_yng[yr[i]]                             * yng[i] +
      b_old[yr[i]]                             * old[i] +
      b_lac_prm[yr[i]]          *      lac[i]  * prm[i] +
      b_lac_yng[yr[i]]          *      lac[i]  * yng[i] +
      b_lac_old[yr[i]]          *      lac[i]  * old[i] +
      bden_dry_prm * elk[yr[i]] * (1 - lac[i]) * prm[i] +          
      bden_dry_yng * elk[yr[i]] * (1 - lac[i]) * yng[i] +
      bden_dry_old * elk[yr[i]] * (1 - lac[i]) * old[i] +
      bpdi_dry_prm * pdi[yr[i]] * (1 - lac[i]) * prm[i] +          
      bpdi_dry_yng * pdi[yr[i]] * (1 - lac[i]) * yng[i] +
      bpdi_dry_old * pdi[yr[i]] * (1 - lac[i]) * old[i] + 
      bden_lac_prm * elk[yr[i]] *      lac[i]  * prm[i] +          
      bden_lac_yng * elk[yr[i]] *      lac[i]  * yng[i] +
      bden_lac_old * elk[yr[i]] *      lac[i]  * old[i] +
      bpdi_lac_prm * pdi[yr[i]] *      lac[i]  * prm[i] +          
      bpdi_lac_yng * pdi[yr[i]] *      lac[i]  * yng[i] +
      bpdi_lac_old * pdi[yr[i]] *      lac[i]  * old[i]

  # Model  
    preg[i] ~ dbin(p[i], nobs[i])
    rsdl[i] <- preg[i] / nobs[i] - p[i]
  }
}
)
