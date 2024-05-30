require(purrr); require(dplyr); require(tidyr); require(rjags)

fd <- readRDS("data//the_ipm_data.rds")
load("data//elk_ipm_data_21apr2023.Rdata")
load("results//ipm_result_21apr2023_R&S_pdi.Rdata")
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
SC_B0_m <- dq %>%
  filter(grepl("S_C_mean", .$var)) %>%
  pull(`50%`) %>%
  logit()

inits <- list(
  R = R_m,
  NFaug = NFaug_m,
  NMaug = NMaug_m,
  NC = NC_ma,
  R_B0_mean = 0,
  R_B0_sd = 1,
  R_Bvy = 0,
  R_Bvm = 0,
  R_Bpu = 0,
  R_Bdd = 0,
  R_B0 = rep(0, cnst$n_year),
  SC_B0_mean = 1,
  SF_mean = .9,
  SM_mean = .8,
  SC_B0_sd = 1,
  SF_B0_sd = 1,
  SM_B0_sd = 1,
  PF_B0_mean = .7,
  PM_B0_mean = .4,
  SC_Bvy = 0,
  SC_Bvm = 0,
  SC_Bpu = 0,
  SC_Bdd = 0,
  S__Bhe = 0,
  P__Bhe = 0,
  SF = c(0.884770525,
         0.889001221,
         0.901899499,
         0.782582815,
         0.969685794,
         0.880314217,
         0.962460604,
         0.918201538,
         0.980822869,
         0.958146548,
         0.935717436,
         0.857003911,
         0.971228656,
         0.823726637,
         0.892128052,
         0.94147829,
         0.94393519,
         0.962483867,
         0.874637305,
         0.8769206,
         0.947995336,
         0.948126502,
         0.905557352,
         0.916445677,
         0.964549789,
         0.927572767,
         0.902715829,
         0.87945314,
         0.964015907,
         0.586101054,
         0.934172145,
         0.945457986,
         0.405536903,
         0.884650128
  ),
  SM = c(0.760445645,
         0.827459361,
         0.608605794,
         0.774641003,
         0.456817942,
         0.853863678,
         0.932848557,
         0.647822222,
         0.786981871,
         0.891836691,
         0.895242527,
         0.83777281,
         0.877889374,
         0.793611417,
         0.532948225,
         0.847230942,
         0.815926834,
         0.607894006,
         0.746633978,
         0.788936987,
         0.866870234,
         0.891211443,
         0.926472189,
         0.81703107,
         0.784293764,
         0.798438117,
         0.860758045,
         0.575583445,
         0.73502042,
         0.809603998,
         0.735831226,
         0.92917555,
         0.325154839,
         0.762322668
  ),
  SC = c(0.737036848,
         0.777284558,
         0.306886163,
         0.91429004,
         0.899786264,
         0.886486318,
         0.947097171,
         0.841966395,
         0.94125792,
         0.91744994,
         0.930535903,
         0.88232726,
         0.934271591,
         0.854962102,
         0.817068509,
         0.809010935,
         0.814618543,
         0.897819545,
         0.668726539,
         0.660063051,
         0.671526593,
         0.727844492,
         0.529458412,
         0.562617538,
         0.635847343,
         0.55438308,
         0.592073328,
         0.338506906,
         0.644386021,
         0.707520546,
         0.629802288,
         0.848658224,
         0.695773043,
         0.734352731
  ),
  S = array(0.9, dim = c(full_data$n_year, 3, 2)),
  PM_B0 = rep(3, cnst$n_year),
  PF_B0 = rep(1, cnst$n_year)
)