nimble_code <- nimbleCode({
  #Observation Models===========================================================
  ##### Abundance #####
  for(i in 1:nNC){ 
    NC_est[i,1] ~ T(dnorm(NC[NC_est_t[i]], sd = NC_est[i,2]), 0, )
  }
  for(i in 1:nNF){
    NF_est[i,1] ~ T(dnorm(NF[NF_est_t[i]], sd = NF_est[i,2]), 0, )
  }
  for(i in 1:nNM){
    NM_est[i,1] ~ T(dnorm(NM[NM_est_t[i]], sd = NM_est[i,2]), 0, )
  }
  
  ##### Counts #####
  sd_afcount ~ dnorm(0, 0.0001)
  tau_afcount <- 1/sd_afcount^2
  for(i in 1:nn_fc){
    af_count[i] ~ T(dnorm(NFaug[af_count_t[i]], tau_afcount), 0, )
  }
  
  sd_amcount ~ dnorm(0, 0.0001)
  tau_amcount <- 1/sd_amcount^2
  for(i in 1:nn_mc){
    am_count[i] ~ T(dnorm(NMaug[am_count_t[i]], tau_amcount), 0, )
  }
})