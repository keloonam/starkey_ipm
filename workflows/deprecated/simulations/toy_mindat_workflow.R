require(R2jags)

cat(
  "model{
    # toy model for testing negative binomial minimum count
    #r ~ dunif(0,5000)
    
    p2 ~ dunif(0,1)
    for(i in 1:nyear){
      p1[i] ~ dunif(0,1)
      N[i] ~ dpois(500)
      M[i] ~ dpois(500)
    }
    for(i in 1:nyear){
      mindat1[i] ~ dbin(p1[i], N[i])
      mindat2[i] ~ dnegbin(p2, M[i])
    }
  }",
  file = here::here("minimum_count_test_model.txt")
)

jags_data <- list(
  nyear = 10,
  mindat1 = c(50,100,200,300,400,500,600,700,800,900),
  mindat2 = c(50,100,200,300,400,500,600,700,800,900)
)

jags(
  data = jags_data,
  parameters.to.save = c("N", "p1", "p2", "M"),
  model.file = here::here("minimum_count_test_model.txt"),
  n.chains = 3,
  n.iter = 5000,
  n.burnin = 1000,
  n.thin = 1
)
