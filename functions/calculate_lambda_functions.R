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
    Lo = calculate_observed_lambda(dt$tm, dt$nv[1:2])
  )
  return(out)
}

calculate_pp <- function(
  b0o, b0p, b0y, bdd, bwm, byr, nf, nof, npf, nyf, dd, wm
){
  pp_old <- ilogit(b0o + bdd*dd + bwm*wm + byr)
  pp_prm <- ilogit(b0p + bdd*dd + bwm*wm + byr)
  pp_yng <- ilogit(b0y + bdd*dd + bwm*wm + byr)
  out <- (pp_old*nof + pp_prm*npf + pp_yng*nyf) / nf
  return(out)
}

calculate_sn <- function(
    b0, bcg, bdd, byr, cg, dd
){
  real_scale <- b0 + bcg*cg + bdd*dd + byr
  out <- ilogit(real_scale)
  return(out)
}
calculate_sc <- function(
    b0, bcg, bwm, bwt, byr, cg, wm, wt
){
  real_scale <- b0 + byr + bcg*cg + bwt*wt + bwm*wm 
  out <- ilogit(real_scale)
  return(out)
}
calculate_sf <- function(
    b0, bwm, bwt, byr, wm, wt
){
  real_scale <- b0 + byr + bwt*wt + bwm*wm 
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
  rs <- rsdt %>% filter(mcmc_step == stpn)
  cv <- cvdt %>% filter(yr == yrx)
  if("cg" %in% names(cv)){
    cg <- cv %>% pull(cg)
  }else{
    cg <- rs %>% filter(yr == yrx) %>% pull(cges)
  }
  dd <- cv %>% pull(dd)
  wm <- cv %>% pull(wm)
  wt <- cv %>% pull(wt)
  
  
  ppyr <- rs %>% filter(yr == yrx) %>% pull(ppyr)
  nf <- rs %>% filter(yr == yrx) %>% pull(nf)
  nof <- rs %>% filter(yr == yrx) %>% pull(nof)
  npf <- rs %>% filter(yr == yrx) %>% pull(npf)
  nyf <- rs %>% filter(yr == yrx) %>% pull(nyf)
  ppdd <- rs %>% filter(!is.na(ppdd)) %>% pull(ppdd)
  ppwm <- rs %>% filter(!is.na(ppwm)) %>% pull(ppwm)
  ppb0o <- rs %>% filter(!is.na(ppb0o)) %>% pull(ppb0o)
  ppb0p <- rs %>% filter(!is.na(ppb0p)) %>% pull(ppb0p)
  ppb0y <- rs %>% filter(!is.na(ppb0y)) %>% pull(ppb0y)
  
  
  pp <- calculate_pp(
    b0o = ppb0o, 
    b0p = ppb0p, 
    b0y = ppb0y, 
    bdd = ppdd,
    bwm = ppwm,
    byr = ppyr,
    nf  = nf,
    nof = nof,
    npf = npf,
    nyf = nyf,
    wm  = wm,
    dd  = dd)
  
  snyr <- rs %>% filter(yr == yrx) %>% pull(snyr)
  sncg <- rs %>% filter(!is.na(sncg)) %>% pull(sncg)
  sndd <- rs %>% filter(!is.na(sndd)) %>% pull(sndd)
  snb0 <- rs %>% filter(!is.na(snb0)) %>% pull(snb0)
  sn <- calculate_sn(
    b0 = snb0, 
    bcg = sncg,
    bdd = sndd,
    byr = snyr,
    cg = cg,
    dd = dd)
  
  scyr <- rs %>% filter(yr == yrx) %>% pull(scyr)
  sccg <- rs %>% filter(!is.na(sccg)) %>% pull(sccg)
  scwm <- rs %>% filter(!is.na(scwm)) %>% pull(scwm)
  scwt <- rs %>% filter(!is.na(scwt)) %>% pull(scwt)
  scb0 <- rs %>% filter(!is.na(scb0)) %>% pull(scb0)
  sc <- calculate_sc(
    b0 = scb0, 
    bcg = sccg,
    bwt = scwt,
    bwm = scwm,
    byr = scyr,
    wt = wt,
    wm = wm,
    cg = cg)
  
  sfyr <- rs %>% filter(yr == yrx) %>% pull(sfyr)
  sfwm <- rs %>% filter(!is.na(sfwm)) %>% pull(sfwm)
  sfwt <- rs %>% filter(!is.na(sfwt)) %>% pull(sfwt)
  sfb0 <- rs %>% filter(!is.na(sfb0)) %>% pull(sfb0)
  sf <- calculate_sf(
    b0 = sfb0, 
    bwt = sfwt,
    bwm = sfwm,
    byr = sfyr,
    wt = wt,
    wm = wm)
  
  sh1 <- rs %>% filter(yr == yrx) %>% mutate(sh1 = 1 - h1) %>% pull(sh1)
  sh2 <- rs %>% filter(yr == yrx) %>% mutate(sh2 = 1 - h2) %>% pull(sh2)
  
  out <- build_elk_tm(PP = pp, SN = sn, SC = sc, SF = sf, H1 = sh1, H2 = sh2)
  return(out)
}

add_abundance_to_list <- function(tm, nv){
  nvn <- matrix(c(nv$nc, nv$nf), nrow = 2, ncol = 1)
  out <- list(tm = tm, nv = as.vector(as.matrix(nv)))
}
