pull_demo <- function(dtf, pat, yrs){
  dtf %>%
    select(grep(pat, names(.))) %>%
    map(quantile, probs = c(0.025, 0.5, 0.975)) %>%
    bind_rows() %>%
    set_names("lci", "mci", "uci") %>%
    mutate(yr = yrs) %>%
    return()
}

calc_fe_only_lambda <- function(fs, r, cs, nf, nc){
  n1 <- nf + nc/2
  n2 <- nf*fs + nc*cs/2 + nf*r*fs/2
  out <- n2 / n1
  return(out)
}

build_lfn_tib <- function(rsn){
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
  return(lfn_df)
}

expit <- function(x){
  1/(1+exp(-x))
}
logit <- function(x){
  log(x / (1-x))
}

calc_marg_line <- function(cov, b0, b1){
  dt_line <- matrix(data = NA, nrow = length(cov), ncol = 4)
  dt_line[,1] <- cov
  for(i in 1:length(cov)){
    y <- expit(b0 + cov[i]*b1) 
    dt_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
  }
  dt_line <- as_tibble(dt_line)
  names(dt_line) <- c("val", "lci", "mci", "uci")
  return(dt_line)
}



calc_eigen_lambda <- function(fs, r, cs){
  tm <- matrix(0, nrow = 2, ncol = 2)
  tm[1,2] <- r * fs / 2
  tm[2,1] <- cs
  tm[2,2] <- fs
  out <- eigen(tm)$values[1]
  return(out)
}

build_el_tib <- function(rsn){
  lambda_from_e <- matrix(NA, nrow = nrow(rsn), ncol = 35)
  s_af <- rsn %>%
    select(grep("survival_af", names(.))) %>%
    as.matrix()
  s_ca <- rsn %>%
    select(grep("survival_ca", names(.))) %>%
    as.matrix()
  r_es <- rsn %>%
    select(grep("R\\[", names(.))) %>%
    as.matrix()
  for(i in 1:nrow(rsn)){
    for(j in 1:ncol(lambda_from_e))
      lambda_from_e[i,j] <- calc_eigen_lambda(
        fs = s_af[i,j],
        r = r_es[i,j],
        cs = s_ca[i,j]
      )
  }
  lfn_df <- lambda_from_e %>%
    as_tibble() %>%
    map(quantile, probs = c(0.025, 0.5, 0.975)) %>%
    bind_rows() %>%
    set_names("lci", "mci", "uci") %>%
    mutate(yr = 1989:2023)
  return(lfn_df)
}


sim_lam <- function(fdt, cov, length_cov, r_cov_name, s_cov_name){
  cov_rng <- seq(min(cov) - 1, max(cov) + 1, length.out = length_cov)
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
    mutate(val = cov_rng)
  cov_dif <- cov_rng[1] - cov_rng[length(cov_rng)]
  lam_dif <- lam[,1] - lam[,ncol(lam)]
  est_slope <- lam_dif / cov_dif
  out <- list(
    lam_df = lam_df,
    est_slope = est_slope
  )
  return(out)
}
