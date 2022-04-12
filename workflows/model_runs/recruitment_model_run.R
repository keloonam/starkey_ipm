# Recruitment Model -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

#Envirnoment====================================================================

require(tidyverse); require(R2jags)
counts <- read_csv("data//min_n_handle_summaries.csv")
# load("data//elk_cjs_data.Rdata")

#Data_prep======================================================================

# n_calf <- apply(
#   elk_cjs_data$calf * elk_cjs_data$y, 
#   2, 
#   sum, 
#   na.rm = T
# )[-c(1,32,33)]
# 
# n_cow <- apply(
#   elk_cjs_data$y * (1-elk_cjs_data$male) * (1-elk_cjs_data$calf),
#   2,
#   sum,
#   na.rm = T
# )[-c(1,32,33)]

n_calf <- counts$ca
n_cow <- counts$fe_ad

jags_data <- list(
  n_cow = n_cow,
  n_calf = n_calf,
  n_years = length(n_cow)
)

params = c("R")

rslt <- jags(
  data = jags_data,
  parameters.to.save = params,
  model.file = "models//recruitment//recruitment.txt",
  n.chains = 3,
  n.iter = 10000,
  n.thin = 1
)

# mcmcplots::mcmcplot(rslt)
save(rslt, file = "results//recruitment_years_2to31.Rdata")
