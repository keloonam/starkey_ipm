require(purrr); require(dplyr); require(tidyr); require(rjags)

fd <- readRDS("data//the_ipm_data.rds")
load("data//elk_ipm_data_21apr2023.Rdata")
load("results//ipm_result_21apr2023_R&S_pdi.Rdata")
full_data <- readRDS("data//the_ipm_data.rds")
data <- summary(rslt)
q    <- data$quantiles
dq <- as_tibble(q) %>%
  mutate(var = dimnames(q)[[1]])
NC_m <- dq %>%
  filter(grepl("N_c", .$var)) %>%
  pull(`50%`)
NF_m <- dq %>%
  filter(grepl("N_f", .$var)) %>%
  pull(`50%`)
NM_m <- dq %>%
  filter(grepl("N_m", .$var)) %>%
  pull(`50%`)
NL_m <- dq %>%
  filter(grepl("N_lam", .$var)) %>%
  pull(`50%`)
NL_m <- c(mean(NL_m), NL_m)
R_m <- dq %>%
  filter(grepl("R\\[", .$var)) %>%
  pull(`50%`) 
R_m <- c(mean(R_m), R_m)
lambda <- dq %>%
  filter(grepl("lambda", .$var)) %>%
  pull(`50%`) 
lambda <- c(mean(lambda), lambda)
R_Bpu_m <- dq %>%
  filter(grepl("R_cg", .$var)) %>%
  pull(`50%`) 
R_Bvy_m <- dq %>%
  filter(grepl("R_wt", .$var)) %>%
  pull(`50%`)
R_Bvm_m <- dq %>%
  filter(grepl("R_wm", .$var)) %>%
  pull(`50%`)
R_Bdd_m <- dq %>%
  filter(grepl("R_dd", .$var)) %>%
  pull(`50%`)
SC_Bpu_m <- dq %>%
  filter(grepl("S_cg", .$var)) %>%
  pull(`50%`) 
SC_Bvy_m <- dq %>%
  filter(grepl("S_wt", .$var)) %>%
  pull(`50%`)
SC_Bvm_m <- dq %>%
  filter(grepl("S_wm", .$var)) %>%
  pull(`50%`)
SC_Bdd_m <- dq %>%
  filter(grepl("S_dd", .$var)) %>%
  pull(`50%`)
SF_m <- dq %>%
  filter(grepl("survival_af", .$var)) %>%
  pull(`50%`)
SM_m <- dq %>%
  filter(grepl("survival_am", .$var)) %>%
  pull(`50%`)
SC_m <- dq %>%
  filter(grepl("survival_ca", .$var)) %>%
  pull(`50%`)
SF_m <- c(mean(SF_m), SF_m)
SM_m <- c(mean(SM_m), SM_m)
SC_m <- c(mean(SC_m), SC_m)
NFaug_m <- NF_m - fd$NFhar + (NC_m/2 * SC_m)
NMaug_m <- NM_m - fd$NMhar + (NC_m/2 * SC_m)
R_B0_m <- dq %>%
  filter(grepl("R_B0", .$var)) %>%
  pull(`50%`)
logit <- function(x){
  log(x/(1-x))
}
R_B0_mean <- dq %>%
  filter(grepl("R_B0", .$var)) %>%
  pull(`50%`) 
R_B0_sd <- sd(logit(R_m))
SC_B0_mean <- dq %>%
  filter(grepl("S_C_mean", .$var)) %>%
  pull(`50%`) %>%
  logit()
NCe_sd <- sqrt(R_m * NF_m * (1 - R_m))
NCe_mean <- R_m * NF_m
NCSe_mean <- (c(0, NC_m[-length(NC_m)]) + full_data$NCman) * SC_m
NCSe_sd <- sqrt(NCSe_mean * (1 - SC_m))
NCaug_m <- NCSe_mean
for(i in 1:length(NCaug_m)){
  if(full_data$NCmin[i] > NCaug_m[i]){
    NCaug_m[i] <- full_data$NCmin[i]
  }
}
NFSe_mean <- (c(0, NF_m[-length(NF_m)]) + full_data$NFman) * SF_m
NFSe_sd <- sqrt(NFSe_mean * (1 - SF_m))
NMSe_mean <- (c(0, NM_m[-length(NM_m)]) + full_data$NMman) * SM_m
NMSe_sd <- sqrt(NMSe_mean * (1 - SM_m))
NCSe_lambda <- c(0, NC_m[-length(NC_m)] * SC_m[-1])
NFe_lambda <- c(0, NF_m[-length(NF_m)] * SF_m[-1] + NCSe_lambda[-1] / 2)
NMe_lambda <- c(0, NM_m[-length(NM_m)] * SM_m[-1] + NCSe_lambda[-1] / 2)
NCe_lambda <- (NFe_lambda + NCSe_lambda) * R_m
N_lambda <- NFe_lambda + NMe_lambda + NCe_lambda
lambda <- dq %>%
  filter(grepl("lambda", .$var)) %>%
  pull(`50%`)
lambda <- c(0, lambda)
Ntot <- NM_m + NF_m + NC_m
sd_afcount <- 80
tau_afcount <- 1/(sd_afcount^2)
sd_amcount <- 30
tau_amcount <- 1/(sd_amcount^2)
ilogit <- function(x){
  return(1/(1+exp(-x)))
}
R0 <- ilogit(
  logit(R_m) - 
    R_Bvy_m * ipm_data$palmer_index - 
    R_Bvm_m * c(0, ipm_data$palmer_index[-length(ipm_data$palmer_index)]) - 
    R_Bdd_m * ipm_data$elk_density -
    R_Bpu_m * ipm_data$cougar_density
  )
SC0 <- ilogit(
  logit(SC_m) - 
    SC_Bvy_m * ipm_data$palmer_index - 
    SC_Bvm_m * c(0, ipm_data$palmer_index[-length(ipm_data$palmer_index)]) - 
    SC_Bdd_m * ipm_data$elk_density -
    SC_Bpu_m * ipm_data$cougar_density
)

inits <- list(
  R = R_m,
  NFaug = NFaug_m,
  NMaug = NMaug_m,
  NCaug = NCaug_m,
  NC = NC_m,
  NF = NF_m,
  NM = NM_m,
  R_B0_mean = mean(logit(R0)), 
  R_B0_sd = sd(logit(R0)), 
  R_Bvy = R_Bvy_m, 
  R_Bvm = R_Bvm_m, 
  R_Bpu = R_Bpu_m, 
  R_Bdd = R_Bdd_m, 
  R_B0 = logit(R0), 
  SC_B0_mean = mean(logit(SC0)), 
  # SF_mean = mean(SF_m), 
  # SM_mean = mean(SM_m), 
  SF_B0_mean = logit(mean(SF_m)),
  SM_B0_mean = logit(mean(SM_m)),
  SC_B0_sd = sd(logit(SC0)), 
  SF_B0_sd = 3*sd(logit(SF_m)), 
  SM_B0_sd = 3*sd(logit(SM_m)), 
  SC_Bvy = SC_Bvy_m, 
  SC_Bvm = SC_Bvm_m, 
  SC_Bpu = SC_Bpu_m, 
  SC_Bdd = SC_Bdd_m, 
  SF = SF_m, 
  SM = SM_m, 
  SC = SC_m, 
  SF_B0 = logit(SF_m), 
  SM_B0 = logit(SM_m), 
  SC_B0 = logit(SC0), 
  SC_Byr = logit(SC_m), 
  NCe_mean = NCe_mean, 
  NCe_sd = NCe_sd, 
  NCSe_mean = NCSe_mean, 
  NCSe_sd = NCSe_sd, 
  NFSe_mean = NFSe_mean, 
  NFSe_sd = NFSe_sd, 
  NMSe_mean = NMSe_mean, 
  NMSe_sd = NMSe_sd, 
  NCSe_lambda = NCSe_lambda, 
  NFe_lambda = NFe_lambda, 
  NMe_lambda = NMe_lambda, 
  NCe_lambda = NCe_lambda, 
  N_lambda = NL_m, 
  LAMBDA = lambda, 
  Ntot = Ntot, 
  sd_afcount = sd_afcount, 
  tau_afcount = tau_afcount, 
  sd_amcount = sd_amcount, 
  tau_amcount = tau_amcount 
)
inits$NCaug[1]  <- NA
inits$NFaug[1]  <- NA
inits$NMaug[1]  <- NA
inits$R_B0[1]   <- NA
inits$SC_B0[1]  <- NA
inits$SF_B0[1]  <- NA
inits$SM_B0[1]  <- NA

saveRDS(inits, file = "data//the_ipm_inits.rds")
rm(list = ls())
