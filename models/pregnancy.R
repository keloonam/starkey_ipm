nimble_code <- nimble::nimbleCode({
  # GLM
  for(i in 1:nr){
    # Prior
    p[i] ~ dunif(0, 1)
    # Model  
    preg[i] ~ dbin(p[i], nobs[i])
  }
}
)