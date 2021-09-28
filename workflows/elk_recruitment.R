# Recruitment Model -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

#Envirnoment====================================================================

require(tidyverse); require(R2jags)
load("data//elk_cjs_data.Rdata")

#Data_prep======================================================================

n_calf <- apply(
  elk_cjs_data$calf * elk_cjs_data$y, 
  2, 
  sum, 
  na.rm = T
)[-c(1,32,33)]

n_cow <- apply(
  elk_cjs_data$y * (1-elk_cjs_data$male) * (1-elk_cjs_data$calf),
  2,
  sum,
  na.rm = T
)[-c(1,32,33)]

jags_data <- list(
  n_cow = n_cow,
  n_calf = n_calf,
  n_years = length(n_cow),
  pr_p = 0.5
)

parms = c("R", "b.0")

rslt <- jags(
  data = jags_data,
  parameters.to.save = parms,
  model.file = "models//recruitment.txt",
  n.chains = 3,
  n.iter = 10000,
  n.thin = 1
)
save(rslt, file = "recruitment_years_2to31.Rdata")
