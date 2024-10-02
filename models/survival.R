code <- nimbleCode({
  # cjs code for Starkey elk with js abundance estimate added
  # 18-August-2022
  
  #Priors=======================================================================
  for(t in 1:(nocc)){
    # time varying effects on
    # detection probability (p), and survival probability (s)
    # including intercepts (0), males (m), and calves (c)
    bp0[t] ~ dlogis(0, 1)
    bpm[t] ~ dlogis(0, 1)
    bs0[t] ~ dlogis(0, 1)
    bsm[t] ~ dlogis(0, 1)
    bsc[t] ~ dlogis(0, 1)
    bsh[t] ~ dlogis(0, 1)
    bph[t] ~ dlogis(0, 1)
  }
  
  for(i in 1:nind){
    for(t in f[i]:l[i]){
      # probabilities to actually use
      logit(p[i,t]) <- bp0[t] + bpm[t]*m[i]*(1-c[i,t]) + bph[t]*h[i,t]
      logit(s[i,t]) <- bs0[t] + bsm[t]*m[i]*(1-c[i,t]) + bsh[t]*h[i,t] + 
        bsc[t]*c[i,t]
    } #i
  } #t
  
  #Process======================================================================
  for(i in 1:nind){
    for(t in (f[i]+1):l[i]){
      # State Process -------- alive >>> z = 1
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- s[i,t-1] * z[i,t-1]
      
      # Observation Process -- seen >>> y = 1
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
  } #i
  
  #Values=======================================================================
  for(t in 1:nocc){
    logit(survival_af[t]) <- bs0[t]
    logit(survival_am[t]) <- bs0[t] + bsm[t]
    logit(survival_ca[t]) <- bs0[t] + bsc[t]
    logit(prob_af[t])     <- bp0[t]
    logit(prob_am[t])     <- bp0[t] + bpm[t]
  }
  
  #Goodness of Fit==============================================================
  for(i in 1:nind){
    for(t in (f[i]+1):l[i]){
      w[i,t] ~ dbern(mu2[i,t])
    }
    n_cap_obs[i] <- sum(y[i,(f[i]+1):l[i]])
    n_cap_sim[i] <- sum(w[i,(f[i]+1):l[i]])
  }
  sd_cap_obs <- sd(n_cap_obs[1:nind])
  sd_cap_sim <- sd(n_cap_sim[1:nind])
  mean_cap_obs <- mean(n_cap_obs[1:nind])
  mean_cap_sim <- mean(n_cap_obs[1:nind])
  pp_mean_difference <- mean_cap_obs - mean_cap_sim
  pp_sd_difference <- sd_cap_obs - sd_cap_sim
})