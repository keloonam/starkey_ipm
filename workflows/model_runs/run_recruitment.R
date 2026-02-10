# Run the cjs model
# Kenneth Loonam
# August 2024

parameter <- "spei_12m" #  

# Load packages
require(nimble); require(tidyverse)
scl <- function(x){
  out <- (x-mean(x))/sd(x)
  return(out)
}

cvdt <- readRDS("data//unscaled_covariates.rds")

# Specify results location
result_file <- paste0("results//recruitment_test//", parameter, ".rds")

# Load the data
cd <- readRDS("data//recruitment_data.rds")
cd$data$xt <- cvdt %>%
  filter(yr %in% cd$constants$years) %>%
  pull(parameter) %>%
  scl()
cd$data$xm <- cvdt %>%
  mutate(prm = lag(get(parameter))) %>%
  filter(yr %in% cd$constants$years) %>%
  pull(prm) %>%
  scl()

# Load the model
source("models//recruitment.R")

# Set trace monitors
params <- c(
  "b0", "bt", "bm", "R", "er", "mn_c_diff", "sd_c_diff"
)

# Run it
rslt <- nimbleMCMC(
  code        = code,
  constants   = cd$constants,
  data        = cd$data,
  monitors    = params,
  inits       = cd$initial_values,
  niter       = 25000,
  nburnin     = 10000,
  nchains     = 3,
  progressBar = T,
  summary     = F,
  check       = F
)  

saveRDS(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)

tmp_summ <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select("bt", "bm", "b0", "mn_c_diff", "sd_c_diff") %>%
  pivot_longer(1:ncol(.), names_to = "param") %>%
  group_by(param) %>%
  summarise(
    lci = quantile(value, 0.025),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.975)
  ) %>%
  ungroup() %>%
  mutate(covariate = parameter)

full_summ <- readRDS("results//recruitment_test//cov_rslt_summ.rds") %>%
  bind_rows(tmp_summ)
saveRDS(full_summ, "results//recruitment_test//cov_rslt_summ.rds")

rm(list = ls())

readRDS("results//recruitment_test//cov_rslt_summ.rds") %>%
  filter(covariate %in% c(
    "mm_precip", "ndvi", "pdi_september", "spei_12m", "temp_mean"
  )) %>%
  filter(param %in% c("bm", "bt")) %>%
  ggplot(aes(x = param, y = mci, color = covariate)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge(width = 0.5)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "11") +
  xlab("Timing") + ylab("Estimate") +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = c("Prior Year", "Current Year")) +
  scale_color_discrete(labels = c(
    "Precipitation",
    "NDVI",
    "PDI", 
    "SPEI",
    "Temperature"
  ))
ggsave(
  "figures//recruitment_covariate_test_fig.png",
  width = 15,
  height = 8,
  units = "cm",
  dpi = 600
)

bucvdt <- cvdt %>%
  filter(yr > 1987) %>%
  select(yr, pdi_september, mm_precip, temp_mean, ndvi, spei_12m)
bucvdt %>%
  ggplot(aes(x = yr, y = spei_12m)) +
  geom_point() +
  geom_smooth(method = "lm")
plot(lm(bucvdt$yr ~ bucvdt$spei_12m))

plot_corr <- function(dt, x, y){
  xl <- case_when(
    x == "yr" ~ "Year",
    x == "pdi_september" ~ "PDI",
    x == "mm_precip" ~ "Precipitation (mm)",
    x == "temp_mean" ~ "Temperature (C)",
    x == "ndvi" ~ "NDVI",
    x == "spei_12m" ~ "SPEI"
  )
  yl <- case_when(
    y == "yr" ~ "Year",
    y == "pdi_september" ~ "PDI",
    y == "mm_precip" ~ "Precipitation (mm)",
    y == "temp_mean" ~ "Temperature (C)",
    y == "ndvi" ~ "NDVI",
    y == "spei_12m" ~ "SPEI"
  )
  xd <- dt %>% pull(x)
  yd <- dt %>% pull(y)
  r2 <- round(summary(lm(yd ~ xd))$r.squared, 3)
  out <- dt %>%
    ggplot(aes(x = get(x), y = get(y))) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() + xlab(xl) + ylab(yl) +
    labs(title = bquote({R^2} == .(r2)))
  return(out)
}
yearspei <- bucvdt %>% plot_corr(x = "yr", y = "spei_12m") + 
  xlab(element_blank())
yearndvi <- bucvdt %>% plot_corr(x = "yr", y = "ndvi") + 
  xlab(element_blank()) +
  ylab(element_blank()) + theme(axis.text.y = element_blank())
yearprcp <- bucvdt %>% plot_corr(x = "yr", y = "mm_precip") + 
  xlab(element_blank()) +
  ylab(element_blank()) + theme(axis.text.y = element_blank())
yeartemp <- bucvdt %>% plot_corr(x = "yr", y = "temp_mean") + 
  xlab(element_blank()) +
  ylab(element_blank()) + theme(axis.text.y = element_blank())
yearpdsi <- bucvdt %>% plot_corr(x = "yr", y = "pdi_september") + 
  xlab(element_blank()) +
  ylab(element_blank()) + theme(axis.text.y = element_blank())

speindvi <- bucvdt %>% plot_corr(x = "spei_12m", y = "ndvi") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank())
speiprcp <- bucvdt %>% plot_corr(x = "spei_12m", y = "mm_precip") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank())
speitemp <- bucvdt %>% plot_corr(x = "spei_12m", y = "temp_mean") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank())
speipdsi <- bucvdt %>% plot_corr(x = "spei_12m", y = "pdi_september")

ndviprcp <- bucvdt %>% plot_corr(x = "ndvi", y = "mm_precip") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylab(element_blank()) + theme(axis.text.y = element_blank())
ndvitemp <- bucvdt %>% plot_corr(x = "ndvi", y = "temp_mean") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylab(element_blank()) + theme(axis.text.y = element_blank())
ndvipdsi <- bucvdt %>% plot_corr(x = "ndvi", y = "pdi_september") +
  ylab(element_blank()) + theme(axis.text.y = element_blank())

prcptemp <- bucvdt %>% plot_corr(x = "mm_precip", y = "temp_mean") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylab(element_blank()) + theme(axis.text.y = element_blank())
prcppdsi <- bucvdt %>% plot_corr(x = "mm_precip", y = "pdi_september") +
  ylab(element_blank()) + theme(axis.text.y = element_blank())

temppdsi <- bucvdt %>% plot_corr(x = "temp_mean", y = "pdi_september") +
  ylab(element_blank()) + theme(axis.text.y = element_blank())

emptyplt <- ggplot(data = NULL, aes(x = 1:10, y = 1:10)) + theme_classic() +
  ylab(element_blank()) + xlab(element_blank()) +
  theme(
    axis.line = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank())
require(cowplot)
plot_grid(
  speindvi, emptyplt, emptyplt, emptyplt,
  speiprcp, ndviprcp, emptyplt, emptyplt,
  speitemp, ndvitemp, prcptemp, emptyplt,
  speipdsi, ndvipdsi, prcppdsi, temppdsi,
  nrow = 4, ncol = 4
)
ggsave(
  "figures//climate_covariate_correlation.png",
  width = 18,
  height = 18,
  units = "cm", 
  dpi = 600,
  scale = 1.25
)

