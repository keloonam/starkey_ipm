pull_param_vec <- function(param, fdt, stu){
  fdt %>% 
    mutate(target_variable = .[[param]]) %>%
    group_by(stpn) %>% 
    summarise(out = mean(target_variable, na.rm = T)) %>%
    pull(out) %>%
    return()
}
calculate_derivative <- function(param, eq, mu){
  tmp <- eval(deriv(Leq, param), envir = mu)
  attributes(tmp)$gradient %>% 
    as.vector() %>%
    tibble(val = .) %>%
    mutate(param = param) %>%
    mutate(stpn = 1:nrow(.)) %>%
    return()
}
calculate_elasticity_with_covariances <- function(
    target_step, fdt, sensitivities
  ){
  vcv <- fdt %>%
    filter(stpn == target_step) %>%
    select(nfp, ncp, sh1, sh2, cg, dd, wt, wm, cm, dm, ppyr, snyr, scyr, sfyr) %>%
    as.matrix() %>%
    var()
  sns <- sensitivities %>%
    filter(stpn == target_step) %>%
    pivot_wider(names_from = param, values_from = val) %>%
    select(nfp, ncp, sh1, sh2, cg, dd, wt, wm, cm, dm, ppyr, snyr, scyr, sfyr) %>%
    as.matrix() %>% as.vector()
  cnt <- vcv * outer(sns, sns)
  out <- rowSums(cnt) %>%
    as_tibble() %>%
    mutate(stpn = target_step) %>%
    mutate(param = c(
      "nfp", "ncp", "sh1", "sh2", "cg", "dd", "wt", "wm", "cm", "dm", "ppyr", "snyr", 
      "scyr", "sfyr"
    ))
  return(out)
}

calculate_elasticity_without_covariances <- function(
    target_step, fdt, sensitivities
){
  vcv <- fdt %>%
    filter(stpn == target_step) %>%
    select(nfp, ncp, sh1, sh2, cg, dd, wt, wm, cm, dm, ppyr, snyr, scyr, sfyr) %>%
    as.matrix() %>%
    var()
  sns <- sensitivities %>%
    filter(stpn == target_step) %>%
    pivot_wider(names_from = param, values_from = val) %>%
    select(nfp, ncp, sh1, sh2, cg, dd, wt, wm, cm, dm, ppyr, snyr, scyr, sfyr) %>%
    as.matrix() %>% as.vector()
  cnt <- diag(vcv) * diag(outer(sns, sns))
  out <- cnt %>%
    as_tibble() %>%
    mutate(stpn = target_step) %>%
    mutate(param = c(
      "nfp", "ncp", "sh1", "sh2", "cg", "dd", "wt", "wm", "cm", "dm", "ppyr", "snyr", 
      "scyr", "sfyr"
    ))
  return(out)
}

get_dfs_mns <- function(target_step, fdt){
  fdt %>% 
    filter(stpn == target_step) %>%
    mutate(
      dpp = ppyr - lag(ppyr),
      dsn = snyr - lag(snyr),
      dsc = scyr - lag(scyr),
      dsf = sfyr - lag(sfyr),
      dh1 = sh1 - lag(sh1),
      dh2 = sh2 - lag(sh2),
      dnf = nfp - lag(nfp),
      dnc = ncp - lag(ncp),
      dcg = cg - lag(cg),
      dcm = cm - lag(cm),
      ddd = dd - lag(dd),
      ddm = dm - lag(dm),
      dwt = wt - lag(wt),
      dwm = wm - lag(wm)
    ) %>%
    mutate(
      ppyr = (ppyr + lag(ppyr))/2,
      snyr = (snyr + lag(snyr))/2,
      scyr = (scyr + lag(scyr))/2,
      sfyr = (sfyr + lag(sfyr))/2,
      sh1 = (sh1 + lag(sh1))/2,
      sh2 = (sh2 + lag(sh2))/2,
      nfp = (nfp + lag(nfp))/2,
      ncp = (ncp + lag(ncp))/2,
      cg = (cg + lag(cg))/2,
      cm = (cm + lag(cm))/2,
      dd = (dd + lag(dd))/2,
      dm = (dm + lag(dm))/2,
      wt = (wt + lag(wt))/2,
      wm = (wm + lag(wm))/2
    ) %>%
    mutate(yrx = lag(yrx)) %>%
    filter(!is.na(yrx)) %>%
    select(
      stpn, yrx, 
      dpp, dsn, dsc, dsf, dh1, dh2, dnf, dnc, dcg, dcm, ddd, ddm, dwt, dwm,
      cg, cm, dd, dm, wt, wm, nfp, ncp, sh1, sh2,
      ppb0, ppyr, ppcg,       ppwm, ppdd, 
      snb0, snyr, sncg, snwt, sndd,
      scb0, scyr, sccg, scwt, scwm, scdd,
      sfb0, sfyr, sfcg, sfwt, sfwm, sfdd
    ) %>%
    return()
}

pull_diff_mat <- function(xdt, dpx){
  xdt %>%
    mutate(dp = !!sym(dpx)) %>%
    select(stpn, yrx, dp) %>%
    pivot_wider(names_from = yrx, values_from = dp) %>%
    select(-stpn) %>%
    as.matrix()
}

pull_mn_mat <- function(dpx, dtx){
  dtx %>%
    mutate(param = !!sym(dpx)) %>%
    select(stpn, yrx, param) %>%
    pivot_wider(names_from = yrx, values_from = param) %>%
    select(-stpn) %>%
    as.matrix() %>%
    return()
}

calc_deriv_ann <- function(param, eq, mu){
  tmp <- eval(deriv(eq, param), envir = mu)
  nstp <- nrow(mu[[1]])
  nyrs <- ncol(mu[[1]])
  attributes(tmp)$gradient %>% 
    as.vector() %>%
    tibble(val = .) %>%
    mutate(param = param) %>%
    mutate(stpn = rep(1:nstp, nyrs)) %>%
    mutate(yrx = sort(rep(2:(nyrs+1), nstp))) %>%
    pivot_wider(names_from = yrx, values_from = val) %>%
    select(-c(stpn, param)) %>%
    as.matrix() %>%
    return()
}
