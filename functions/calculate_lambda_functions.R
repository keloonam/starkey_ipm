require(dplyr); require(purrr); require(tidyr); require(furrr)
ilogit <- function(x){
  return(1/(1+exp(-x)))
}
calculate_observed_lambda <- function(tm, Nvec){
  N0 <- sum(Nvec, na.rm = T)
  N1 <- sum(tm %*% Nvec)
  out <- c(observed_lambda = N1/N0)
  return(out)
}
calculate_final_lambda <- function(tm){
  out <- c(final_lambda = as.numeric(eigen(tm)$values[1]))
  return(out)
}
calculate_lambdas <- function(dt){
  out <- tibble(
    Lf = calculate_final_lambda(dt$tm),
    Lo = calculate_observed_lambda(dt$tm, dt$nv)
  )
  return(out)
}

calculate_pp <- function(
  b0, bcg, bdd, bwm, byr, cg, dd, wm,
  use_cg = T,
  use_dd = T,
  use_wm = T
){
  real_scale <- b0 + bcg*cg*use_cg + bdd*dd*use_dd + bwm*wm*use_wm + byr
  out <- ilogit(real_scale)
  return(out)
}

calculate_sn <- function(
    b0, bcg, bdd, bwt, byr, cg, dd, wt,
    use_cg = T,
    use_dd = T,
    use_wt = T
){
  real_scale <- b0 + bcg*cg*use_cg + bdd*dd*use_dd + bwt*wt*use_wt + byr
  out <- ilogit(real_scale)
  return(out)
}
calculate_sc <- function(
    b0, bcg, bdd, bwm, bwt, byr, cg, dd, wm, wt,
    use_cg = T,
    use_dd = T,
    use_wt = T,
    use_wm = T
){
  real_scale <- b0 + byr + 
    bcg*cg*use_cg + 
    bdd*dd*use_dd + 
    bwt*wt*use_wt + 
    bwm*wm*use_wm 
  out <- ilogit(real_scale)
  return(out)
}
calculate_sf <- function(
    b0, bcg, bdd, bwm, bwt, byr, cg, dd, wm, wt,
    use_cg = T,
    use_dd = T,
    use_wt = T,
    use_wm = T
){
  real_scale <- b0 + byr + 
    bcg*cg*use_cg + 
    bdd*dd*use_dd + 
    bwt*wt*use_wt + 
    bwm*wm*use_wm 
  out <- ilogit(real_scale)
  return(out)
}
build_elk_tm <- function(PP, SN, SC, SF, H1, H2){
  tm <- matrix(0, 2, 2)
  tm[1,2] <- PP * SN
  tm[2,1] <- SC * H1
  tm[2,2] <- SF * H2
  return(tm)
}
build_tm <- function(stpn, rsdt, cvdt, yrx){
  cv <- cvdt %>% filter(yr == yrx)
  cg <- cv %>% pull(cg)
  dd <- cv %>% pull(dd)
  wm <- cv %>% pull(wm)
  wt <- cv %>% pull(wt)
  
  rs <- rsdt %>% filter(mcmc_step == stpn)
  ppyr <- rs %>% filter(yr == yrx) %>% pull(ppyr)
  ppcg <- rs %>% filter(!is.na(ppcg)) %>% pull(ppcg)
  ppdd <- rs %>% filter(!is.na(ppdd)) %>% pull(ppdd)
  ppwm <- rs %>% filter(!is.na(ppwm)) %>% pull(ppwm)
  ppb0 <- rs %>% filter(!is.na(ppb0)) %>% pull(ppb0)
  pp <- calculate_pp(
    b0 = ppb0, 
    bcg = ppcg,
    bdd = ppdd,
    bwm = ppwm,
    byr = ppyr,
    wm = wm,
    cg = cg,
    dd = dd)
  
  snyr <- rs %>% filter(yr == yrx) %>% pull(snyr)
  sncg <- rs %>% filter(!is.na(sncg)) %>% pull(sncg)
  sndd <- rs %>% filter(!is.na(sndd)) %>% pull(sndd)
  snwt <- rs %>% filter(!is.na(snwt)) %>% pull(snwt)
  snb0 <- rs %>% filter(!is.na(snb0)) %>% pull(snb0)
  sn <- calculate_sn(
    b0 = snb0, 
    bcg = sncg,
    bdd = sndd,
    bwt = snwt,
    byr = snyr,
    wt = wt,
    cg = cg,
    dd = dd)
  
  scyr <- rs %>% filter(yr == yrx) %>% pull(scyr)
  sccg <- rs %>% filter(!is.na(sccg)) %>% pull(sccg)
  scdd <- rs %>% filter(!is.na(scdd)) %>% pull(scdd)
  scwm <- rs %>% filter(!is.na(scwm)) %>% pull(scwm)
  scwt <- rs %>% filter(!is.na(scwt)) %>% pull(scwt)
  scb0 <- rs %>% filter(!is.na(scb0)) %>% pull(scb0)
  sc <- calculate_sc(
    b0 = scb0, 
    bcg = sccg,
    bdd = scdd,
    bwt = scwt,
    bwm = scwm,
    byr = scyr,
    wt = wt,
    wm = wm,
    cg = cg,
    dd = dd)
  
  sfyr <- rs %>% filter(yr == yrx) %>% pull(sfyr)
  sfcg <- rs %>% filter(!is.na(sfcg)) %>% pull(sfcg)
  sfdd <- rs %>% filter(!is.na(sfdd)) %>% pull(sfdd)
  sfwm <- rs %>% filter(!is.na(sfwm)) %>% pull(sfwm)
  sfwt <- rs %>% filter(!is.na(sfwt)) %>% pull(sfwt)
  sfb0 <- rs %>% filter(!is.na(sfb0)) %>% pull(sfb0)
  sf <- calculate_sf(
    b0 = sfb0, 
    bcg = sfcg,
    bdd = sfdd,
    bwt = sfwt,
    bwm = sfwm,
    byr = sfyr,
    wt = wt,
    wm = wm,
    cg = cg,
    dd = dd)
  
  sh1 <- rs %>% filter(yr == yrx) %>% mutate(sh1 = 1 - h1) %>% pull(sh1)
  sh2 <- rs %>% filter(yr == yrx) %>% mutate(sh2 = 1 - h2) %>% pull(sh2)
  
  out <- build_elk_tm(PP = pp, SN = sn, SC = sc, SF = sf, H1 = sh1, H2 = sh2)
  return(out)
}

add_abundance_to_list <- function(tm, nv){
  nvn <- matrix(c(nv$nc, nv$nf), nrow = 2, ncol = 1)
  out <- list(tm = tm, nv = as.vector(as.matrix(nv)))
}
