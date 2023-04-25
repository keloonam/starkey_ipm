# Kenneth Loonam
# February 2023
# Post hoc analyses of covariate effects on lambda and calf survival

#Environment====================================================================
require(dplyr); require(ggplot2); require(ggsci); require(rjags)
require(cowplot); require(tidyr); require(purrr)
load("results//ipm_result_21apr2023_R&S_pdi.Rdata")
load("data//elk_ipm_data_21apr2023.Rdata")

#Data Prep======================================================================

r <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("_r\\[", names(.))) %>%
  as.matrix()
r_mean <- apply(r, 1, mean)
r_sd   <- apply(r, 1, sd)

NT <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("N_tot", names(.))) %>%
  as.matrix()

NL <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("N_lam", names(.))) %>%
  as.matrix()

lam <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("lambda", names(.))) %>%
  as.matrix()

NS <- matrix(data = NA, nrow = nrow(NT), ncol = ncol(NT))
NS[,1] <- NT[,1]

#Simulate_NS_under_exponential_growth===========================================

for(t in 2:ncol(NS)){
  NS[,t] <- NT[,t-1] * exp(rnorm(1, r_mean, r_sd))
}

#Fit_DD_to_simulated_population=================================================
alpha_S <- beta_S <- sigma_S <- rep(NA, nrow(NS))
for(i in 1:nrow(NS)){
  lm1 <- lm(log(NS[i,-1]) ~ log(NT[i,-ncol(NS)]))
  alpha_S[i] <- coef(lm1)[1]
  beta_S[i]  <- coef(lm1)[2] - 1
  sigma_S[i] <- summary(lm1)$sigma
}

#Fit_DD_to_"real"_population====================================================
# NL is the expected population based on survival and recruitment
# it removes the effects of harvest and management removals/additions
alpha_T <- beta_T <- sigma_T <- rep(NA, nrow(NT))
for(i in 1:nrow(NL)){
  lm2 <- lm(log(NL[i,]) ~ log(NT[i,-ncol(NS)]))
  alpha_T[i] <- coef(lm2)[1]
  beta_T[i]  <- coef(lm2)[2] - 1
  sigma_T[i] <- summary(lm2)$sigma
}

#Visualize_the_results==========================================================

dd_res <- tibble(
  beta = c(beta_S, beta_T),
  source = c(
    rep("simulated exponential growth", length(beta_S)),
    rep("observed", length(beta_T)))
)

ggplot(data = dd_res, aes(x = beta, color = source)) +
  geom_density() +
  theme_classic()
