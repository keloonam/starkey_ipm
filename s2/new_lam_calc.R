# Kenneth Loonam
# February 2023
# Post hoc analyses of covariate effects on lambda and calf survival

#Environment====================================================================
require(dplyr); require(ggplot2); require(ggsci); require(rjags)
require(cowplot); require(tidyr); require(purrr)

#Data===========================================================================
rs2 <- readRDS("s2//results//ipmrs_27sep2025_spei12.rds")
ipm_data <- readRDS("s2//ipm_data_25sep2024.rds")
rs1 <- readRDS("s2//results//ipmrs_27sep2025_null.rds")

#Functions======================================================================
expit <- function(x){
  (1/(1+exp(-x))) %>% return()
}
logit <- function(x){
  (log(x / (1-x))) %>% return()
}

#Build_lambda_1=================================================================
rsf <- rs2 %>%
  map(as_tibble) %>%
  bind_rows()
rsn <- rs1 %>%
  map(as_tibble) %>%
  bind_rows()

sim_lam <- function(fdt, cov, length_cov, r_cov_name, s_cov_name){
  cov_rng <- seq(min(cov), max(cov), length.out = length_cov)
  lam <- matrix(NA, nrow = nrow(fdt), ncol = length_cov)
  r <- lam
  sf <- lam
  sc <- lam
  rb1 <- fdt %>% pull(r_cov_name)
  sb1 <- fdt %>% pull(s_cov_name)
  rb0 <- fdt %>% pull(R_B0)
  sb0 <- fdt %>% pull(S_C_B0_ps) %>% logit()
  psf <- fdt %>% pull(S_Y_F_B0_ps)
  
  for(i in 1:nrow(lam)){
    r[i,]  <- (rb0[i] + rb1[i] * cov_rng) %>% expit()
    sc[i,] <- (sb0[i] + sb1[i] * cov_rng) %>% expit()
    sf[i,] <- psf[i]
  }
  for(i in 1:nrow(lam)){
    for(j in 1:ncol(lam)){
      tm <- matrix(0, nrow = 2, ncol = 2)
      tm[1,2] <- (r[i,j] / 2) * sf[i,j]
      tm[2,1] <- sc[i,j]
      tm[2,2] <- sf[i,j]
      lam[i,j] <- eigen(tm)$values[1]
    }
  }
  lam_df <- lam %>%
    as_tibble() %>%
    map(quantile, probs = c(0.025, 0.5, 0.975)) %>%
    bind_rows() %>%
    set_names("lci", "mci", "uci") %>%
    mutate(cov = cov_rng)
  cov_dif <- cov_rng[1] - cov_rng[length(cov_rng)]
  lam_dif <- lam[,1] - lam[,ncol(lam)]
  est_slope <- lam_dif / cov_dif
  out <- list(
    lam_df = lam_df,
    est_slope = est_slope
  )
}

cg_sim_dt <- sim_lam(
  fdt = rsf, 
  cov = ipm_data$cdens, 
  length_cov = 10,
  r_cov_name = "R_cg",
  s_cov_name = "S_cg"
)
wm_sim_dt <- sim_lam(
  fdt = rsf, 
  cov = ipm_data$spei12, 
  length_cov = 10,
  r_cov_name = "R_wm",
  s_cov_name = "S_wm"
)
wt_sim_dt <- sim_lam(
  fdt = rsf, 
  cov = ipm_data$spei12, 
  length_cov = 10,
  r_cov_name = "R_wt",
  s_cov_name = "S_wt"
)
dd_sim_dt <- sim_lam(
  fdt = rsf, 
  cov = ipm_data$nelk, 
  length_cov = 10,
  r_cov_name = "R_dd",
  s_cov_name = "S_dd"
)
c1 <- "#0000ff"
c2 <- "#db720f"
dd_sim_dt$lam_df %>%
  ggplot(aes(x = cov, y = mci)) +
  geom_ribbon(aes(ymax = uci, ymin = lci), fill = c1, alpha = 0.25) +
  geom_line(color = c1) +
  theme_classic()
dd_sim_dt$est_slope %>% quantile(c(0.025, 0.5, 0.975))
(dd_sim_dt$est_slope > 0) %>% mean()

#Lambda_2=======================================================================
calc_fe_only_lambda <- function(fs, r, cs, nf, nc){
  n1 <- nf + nc/2
  n2 <- nf*fs + nc*cs/2 + nf*r*fs/2
  out <- n2 / n1
  return(out)
}
lambda_from_n <- matrix(NA, nrow = nrow(rsn), ncol = 35)
s_af <- rsn %>%
  select(grep("survival_af", names(.))) %>%
  as.matrix()
s_ca <- rsn %>%
  select(grep("survival_ca", names(.))) %>%
  as.matrix()
r_es <- rsn %>%
  select(grep("R\\[", names(.))) %>%
  as.matrix()
n_ca <- rsn %>%
  select(grep("N_c\\[", names(.))) %>%
  select(-ncol(.)) %>%
  as.matrix()
n_af <- rsn %>%
  select(grep("N_f", names(.))) %>%
  select(-ncol(.)) %>%
  as.matrix()
for(i in 1:nrow(rsn)){
  for(j in 1:ncol(lambda_from_n))
    lambda_from_n[i,j] <- calc_fe_only_lambda(
      fs = s_af[i,j],
      r = r_es[i,j],
      cs = s_ca[i,j],
      nf = n_af[i,j], 
      nc = n_ca[i,j]
    )
}
lfn_df <- lambda_from_n %>%
  as_tibble() %>%
  map(quantile, probs = c(0.025, 0.5, 0.975)) %>%
  bind_rows() %>%
  set_names("lci", "mci", "uci") %>%
  mutate(yr = 1989:2023)
lfn_df %>%
  ggplot(aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = c2, alpha = 0.25) +
  geom_line(color = c2) +
  geom_hline(yintercept = 1, linetype = 2)
glm(
  lfn_df$mci ~ ipm_data$cdens[-1] + ipm_data$spei12[-1] + ipm_data$spei12[-36] + ipm_data$nelk[-36],
  family = gaussian(link = "log")
) %>% summary()
