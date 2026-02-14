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
rs1 <- lgstrs %>% param_summ(PB_dd); lgstrs %>% param_summ(PB_wm)
rs2 <- meanrs %>% param_summ(PB_dd); meanrs %>% param_summ(PB_wm)
rs3 <- normrs %>% param_summ(PB_dd); normrs %>% param_summ(PB_wm)
rs4 <- odfwrs %>% param_summ(PB_dd); odfwrs %>% param_summ(PB_wm)
rs5 <- rcnsrs %>% param_summ(PB_dd); rcnsrs %>% param_summ(PB_wm)

full_param_summ <- function(fdt_list, param){
  rs1 <- fdt_list$lgst %>% param_summ({{param}})
  rs2 <- fdt_list$mean %>% param_summ({{param}})
  rs3 <- fdt_list$norm %>% param_summ({{param}})
  rs4 <- fdt_list$odfw %>% param_summ({{param}})
  rs5 <- fdt_list$rcns %>% param_summ({{param}})
  
  rbind(rs1, rs2, rs3, rs4, rs5) %>% colMeans() %>% round() %>% return()
}
full_param_summ(rslt_list, PB_dd)
