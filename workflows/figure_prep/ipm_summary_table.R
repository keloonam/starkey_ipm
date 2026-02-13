require(tidyverse)
rs_raw <- readRDS("results//fbipm_rslt_06feb2026_full_mean.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
expit <- function(x){
  out <- 1/(1+exp(-x))
  return(out)
}
build_tibble_row <- function(fdt, param){
  if(grepl("sd", param)){
    x <- fdt %>% pull(param)
  }else{
    x <- fdt %>% pull(param) %>% expit()
  }
  tibble(
    Parameter = param,
    LCRI = round(quantile(x, 0.025), 2),
    Median = round(quantile(x, 0.5), 2),
    MCRI = round(quantile(x, 0.975), 2)
  ) %>%
    return()
}

params_oi <- c(
  "P_B0[1]", "P_B0[2]", "P_B0[3]", "P_sd", "SN_B0", "SN_sd", "SC_B0", "SC_sd", 
  "SF_B0", "SF_sd", "SM_B0", "SM_sd", 
  "H_B0[1,1]", "H_B0[2,1]", "H_B0[1,2]", "H_B0[2,2]",
  "sd_afcount", "sd_amcount"
)

map(params_oi, build_tibble_row, fdt = rs_raw) %>%
  bind_rows() %>%
  write_csv(file = "figures//ipm_demographic_summary_table.csv")
