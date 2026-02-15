require(tidyverse); require(rjags)
ipmdt <- readRDS("data//ipm_data_03feb2026.rds")
lgstrs <- readRDS("results//fbipm_rslt_06feb2026_full_lgst.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
meanrs <- readRDS("results//fbipm_rslt_06feb2026_full_mean.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
normrs <- readRDS("results//fbipm_rslt_06feb2026_full_norm.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
odfwrs <- readRDS("results//fbipm_rslt_06feb2026_full_odfw.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
rcnsrs <- readRDS("results//fbipm_rslt_06feb2026_full_rcns.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()

rslt_list <- list(
  lgst = lgstrs,
  mean = meanrs,
  norm = normrs,
  odfw = odfwrs,
  rcns = rcnsrs
)

to_odds <- function(x){
  return((exp(x) - 1) * 100)
}

param_summ <- function(fdt, param){
  fdt %>% 
    pull({{param}}) %>%
    to_odds() %>%
    quantile(c(0.025, 0.5, 0.975)) %>%
    return()
}
expit <- function(x){
  1/(1+exp(-x))
}
ann_mean <- function(fdt, param_pattern){
  fdt %>%
    select(grep(param_pattern, names(.))) %>%
    as.matrix() %>% colMeans() %>% mean() %>%
    return()
}

lgstrs %>% param_summ(SNB_dd); lgstrs %>% param_summ(SNB_cg)
meanrs %>% param_summ(SNB_dd); meanrs %>% param_summ(SNB_cg)
normrs %>% param_summ(SNB_dd); normrs %>% param_summ(SNB_cg)
odfwrs %>% param_summ(SNB_dd); odfwrs %>% param_summ(SNB_cg)
rcnsrs %>% param_summ(SNB_dd); rcnsrs %>% param_summ(SNB_cg)

full_param_summ <- function(fdt_list, param){
  rs1 <- fdt_list$lgst %>% param_summ({{param}})
  rs2 <- fdt_list$mean %>% param_summ({{param}})
  rs3 <- fdt_list$norm %>% param_summ({{param}})
  rs4 <- fdt_list$odfw %>% param_summ({{param}})
  rs5 <- fdt_list$rcns %>% param_summ({{param}})
  
  rbind(rs1, rs2, rs3, rs4, rs5) %>% colMeans() %>% round() %>% return()
}
full_param_summ(rslt_list, SFB_wm)
ann_mean(meanrs, "N_f\\[") * ann_mean(meanrs, "Pbar")
sncg <- meanrs %>% pull(SNB_cg)
snb0 <- meanrs %>% pull(SN_B0)
expit(snb0 + sncg*ipmdt$cg_mean_estimate[1]) %>% mean()
expit(snb0 + sncg*ipmdt$cg_mean_estimate[20]) %>% mean()
188*0.73 - 188*0.46
