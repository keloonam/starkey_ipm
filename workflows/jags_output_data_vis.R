# IPM data visualization
# Kenneth Loonam
# April 2021

#Variables======================================================================

data_file <- "results//ipm_result.Rdata"

#Environment====================================================================

require(tidyverse); require(ggplot2); require(mcmcplots); require(R2jags)

#Data===========================================================================

load(data_file)

data_raw <- rslt$BUGSoutput$summary

s_af <- data_raw[grep("survival_af", dimnames(data_raw)[[1]]), ] %>%
  as_tibble(.) %>%
  mutate(year = (1:nrow(.)) + 1987)

s_am <- data_raw[grep("survival_am", dimnames(data_raw)[[1]]), ] %>%
  as_tibble(.) %>%
  mutate(year = (1:nrow(.)) + 1987)

s_ca <- data_raw[grep("survival_ca", dimnames(data_raw)[[1]]), ] %>%
  as_tibble(.) %>%
  mutate(year = (1:nrow(.)) + 1987)

recruit <- data_raw[grep("recruitment", dimnames(data_raw)[[1]]), ] %>%
  as_tibble(.) %>%
  mutate(year = (1:nrow(.)) + 1987)



N_af <- data_raw[grep("N_af", dimnames(data_raw)[[1]]), ] %>%
  as_tibble(.) %>%
  mutate(year = (1:nrow(.)) + 1987)

N_am <- data_raw[grep("N_am", dimnames(data_raw)[[1]]), ] %>%
  as_tibble(.) %>%
  mutate(year = (1:nrow(.)) + 1987)

N_ca <- data_raw[grep("N_ca", dimnames(data_raw)[[1]]), ] %>%
  as_tibble(.) %>%
  mutate(year = (1:nrow(.)) + 1987)

#Functions======================================================================

plot_fn <- function(x, sz = 2){
  plot_data <- data_raw[grep(x, dimnames(data_raw)[[1]]), ] %>%
    as_tibble(.) %>%
    mutate(year = (1:nrow(.)) + 1987)
  
  ggplot(plot_data) +
    geom_point(aes(x = year, y = mean), size = sz) + 
    geom_line(aes(x = year, y = mean), size = sz - .75) +
    geom_line(aes(x = year, y = `2.5%`), color = "cyan3", size = sz - .75) +
    geom_line(aes(x = year, y = `97.5%`), color = "cyan3", size = sz - .75) + 
    theme_light() +
    ylab(x) +
    xlab("Year") + 
    theme(text = element_text(size = 18)) +
    guides(color = "none", size = "none")
}

#Examples=======================================================================

plot_fn("N_tot")
plot_fn("N_af")
plot_fn("survival_af")
plot_fn("recruitment")
plot_fn("survival_af")

