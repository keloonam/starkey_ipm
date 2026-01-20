require(tidyverse)
use_years <- 1989:2023

ilogit <- function(x){
  1/(1+exp(-x))
}

logit <- function(x){
  log(x/(1-x))
}

pull_parameter_series <- function(demo_rate, use_years){
  frs <- readRDS("results//mscr_rslt.rds") %>%
    map(as_tibble) %>%
    bind_rows()  %>%
    select(grep(demo_rate, names(.)))
  cntr <- names(frs)
  names(cntr) <- as.character(use_years)
  frs %>%
    rename(all_of(cntr)) %>%
    pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
    mutate(param = demo_rate) %>%
    return()
}

long_results <- map(
  .x = c("Pf", "Pm", "Snaf", "Snam", "Snca", "Sham", "Shaf", "Shjm", "Shjf"),
  .f = pull_parameter_series, 
  use_years = use_years) %>%
  bind_rows()

rslt_summ <- long_results %>%
  mutate(rsp = logit(val)) %>%
  group_by(param, yr) %>%
  summarise(mn = mean(rsp), tau = 1/sd(rsp)^2)

saveRDS(rslt_summ, file = "results//mscr_rslt_summ.rds")

# long_results %>%
#   filter(param == "Shaf") %>%
#   mutate(val = logit(val)) %>%
#   group_by(yr, param) %>%
#   summarise(
#     lci = quantile(val, 0.025), 
#     mci = quantile(val, 0.5), 
#     uci = quantile(val, 0.975)) %>%
#   ggplot(aes(x = yr, y = mci, col = param)) +
#   geom_pointrange(aes(ymin = lci, ymax = uci)) +
#   theme_classic()
