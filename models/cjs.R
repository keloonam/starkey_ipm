code <- nimbleCode({
  # cjs code for Starkey elk with js abundance estimate added
  # 18-August-2022
  
  #Priors=======================================================================
  for(t in 1:(nocc)){
    # time varying effects on
    # detection probability (p), and survival probability (s)
    # including intercepts (0), males (m), and calves (c)
    # Detection Probability
    bp0[t] ~ dlogis(0,1) # Intercept
    bpm[t] ~ dlogis(0,1) # Male
    bph[t] ~ dlogis(0,1) # Herd
    # Survival
    bs0[t] ~ dlogis(0,1) # Intercept
    bsm[t] ~ dlogis(0,1) # Male
    bsc[t] ~ dlogis(0,1) # Calf
    bsh[t] ~ dlogis(0,1) # Herd
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
  
  #Goodness of Fit==============================================================
  for(i in 1:nind){
    for(t in (f[i]+1):l[i]){
      y_new[i,t] ~ dbern(mu2[i,t])
    }
    obs_i_sum[i] <- sum(y[i,(f[i]+1):l[i]])
    sim_i_sum[i] <- sum(y_new[i,(f[i]+1):l[i]])
    dif_i_sum[i] <- sim_i_sum[i] - obs_i_sum[i]
  }
  
  post_diff_mn <- mean(dif_i_sum[1:nind])
  post_diff_sd <- sd(dif_i_sum[1:nind])
  
  #Values=======================================================================
  for(t in 1:nocc){
    logit(survival_af[t]) <- bs0[t]
    logit(survival_am[t]) <- bs0[t] + bsm[t]
    logit(survival_ca[t]) <- bs0[t] + bsc[t]
    logit(prob_af[t])     <- bp0[t]
    logit(prob_am[t])     <- bp0[t] + bpm[t]
  }
})