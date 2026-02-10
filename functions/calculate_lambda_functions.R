require(dplyr); require(purrr); require(tidyr); require(furrr); require(rlang)
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
  b0o, b0p, b0y, bdd, bwm, bwt, bcg, byr, nf, nof, npf, nyf, dd, wm, wt, cg
){
  pp_old <- ilogit(b0o + bdd*dd + bwm*wm + bwt*wt + bcg*cg + byr)
  pp_prm <- ilogit(b0p + bdd*dd + bwm*wm + bwt*wt + bcg*cg + byr)
  pp_yng <- ilogit(b0y + bdd*dd + bwm*wm + bwt*wt + bcg*cg + byr)
  out <- (pp_old*nof + pp_prm*npf + pp_yng*nyf) / nf
  return(out)
}
calculate_sx <- function(
    b0, bcg, bdd, bwt, bwm, byr, cg, dd, wt, wm
){
  real_scale <- b0 + bcg*cg + bdd*dd + bwt*wt + bwm*wm + byr
  ilogit(real_scale) %>% return()
}
pull_vrx <- function(dt, vrx){
  vrx <- enquo(vrx)
  if(as_name(vrx) %in% names(dt)){
    out <- dt %>%
      filter(!is.na(!!vrx)) %>% pull(!!vrx)
  }else{
    out <- 0
  }
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
  ppdd <- rs %>% pull_vrx(ppdd)
  ppwm <- rs %>% pull_vrx(ppwm)
  ppcg <- rs %>% pull_vrx(ppcg)
  ppwt <- rs %>% pull_vrx(ppwt)
  ppb0o <- rs %>% pull_vrx(ppb0o)
  ppb0p <- rs %>% pull_vrx(ppb0p)
  ppb0y <- rs %>% pull_vrx(ppb0y)
  
  pp <- calculate_pp(
    b0o = ppb0o, 
    b0p = ppb0p, 
    b0y = ppb0y, 
    bdd = ppdd,
    bwm = ppwm,
    bwt = ppwt,
    bcg = ppcg,
    byr = ppyr,
    nf  = nf,
    nof = nof,
    npf = npf,
    nyf = nyf,
    wm  = wm,
    dd  = dd,
    wt = wt,
    cg = cg)
  
  snyr <- rs %>% filter(yr == yrx) %>% pull(snyr)
  sncg <- rs %>% pull_vrx(sncg)
  sndd <- rs %>% pull_vrx(sndd)
  snwm <- rs %>% pull_vrx(snwm)
  snwt <- rs %>% pull_vrx(snwt)
  snb0 <- rs %>% pull_vrx(snb0)
  sn <- calculate_sx(
    b0 = snb0, 
    bcg = sncg,
    bdd = sndd,
    bwt = snwt,
    bwm = snwm,
    byr = snyr,
    cg = cg,
    dd = dd,
    wt = wt,
    wm = wm)
  
  scyr <- rs %>% filter(yr == yrx) %>% pull(scyr)
  sccg <- rs %>% pull_vrx(sccg)
  scdd <- rs %>% pull_vrx(scdd)
  scwm <- rs %>% pull_vrx(scwm)
  scwt <- rs %>% pull_vrx(scwt)
  scb0 <- rs %>% pull_vrx(scb0)
  sc <- calculate_sx(
    b0 = scb0, 
    bcg = sccg,
    bdd = scdd,
    bwt = scwt,
    bwm = scwm,
    byr = scyr,
    cg = cg,
    dd = dd,
    wt = wt,
    wm = wm)
  
  sfyr <- rs %>% filter(yr == yrx) %>% pull(sfyr)
  sfcg <- rs %>% pull_vrx(sfcg)
  sfdd <- rs %>% pull_vrx(sfdd)
  sfwm <- rs %>% pull_vrx(sfwm)
  sfwt <- rs %>% pull_vrx(sfwt)
  sfb0 <- rs %>% pull_vrx(sfb0)
  sf <- calculate_sx(
    b0 = sfb0, 
    bcg = sfcg,
    bdd = sfdd,
    bwt = sfwt,
    bwm = sfwm,
    byr = sfyr,
    cg = cg,
    dd = dd,
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
