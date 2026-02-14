require(nimble); require(tidyverse)

# Environment===================================================================
prism <- read_csv("data//prism_data_starkey_centroid.csv") %>%
  set_names("dt", "mm_precip", "t_min", "t_mean", "t_max") %>%
  mutate(dt = mdy(dt)) %>%
  mutate(yr = year(dt), mn = month(dt)) %>%
  filter(mn %in% c(5,6,7,8,9)) %>%
  group_by(yr) %>%
  summarise(
    mm_precip = sum(mm_precip),
    temp_mean = mean(t_mean)
  ) %>%
  ungroup() %>%
  # mutate(precip = scl(mm_precip)) %>%
  # mutate(temp = scl(t_mean)) %>%
  select(yr, mm_precip, temp_mean)


scl <- function(x){
  ((x - mean(x)) / sd(x)) %>% return()
}

spei1 <- read_csv("data//spei.csv") %>%
  mutate(dt = my(DATA)) %>%
  mutate(yr = year(dt)) %>%
  mutate(mn = month(dt)) %>%
  mutate(yr = case_when(
    yr > 2024 ~ yr - 100,
    T ~ yr
  )) %>%
  mutate(spei_1m = SPEI_1) %>%
  mutate(spei_3m = SPEI_3) %>%
  mutate(spei_6m = SPEI_6) %>%
  mutate(spei_12m = SPEI_12) %>%
  mutate(spei_24m = SPEI_24) %>%
  mutate(spei_48m = SPEI_48) %>%
  filter(mn == 9) %>%
  filter(yr < 2024) %>%
  select(yr, spei_1m, spei_3m, spei_6m, spei_12m, spei_24m, spei_48m) %>%
  arrange(yr) %>%
  mutate(rng = "a")

spei2 <- spei1 %>% filter(yr > 1987) %>%
  mutate(rng = "b")
spei <- bind_rows(spei1, spei2)

spei_a <- spei1 %>%
  mutate(yr = yr - mean(yr)) %>%
  select(yr, spei_12m) %>% as.matrix()
spei_b <- spei2 %>%
  mutate(yr = yr - mean(yr)) %>%
  select(yr, spei_12m) %>% as.matrix()


nim_code <- nimbleCode({
  b0_a ~ dnorm(0, sd = 10)
  b1_a ~ dnorm(0, sd = 10)
  sd_a ~ dunif(0, 10)
  b0_b ~ dnorm(0, sd = 10)
  b1_b ~ dnorm(0, sd = 10)
  sd_b ~ dunif(0, 10)
  
  for(i in 1:n_a){
    spei_a[i,2] ~ dnorm(mu_a[i], sd_a)
    mu_a[i] <- b0_a + b1_a * spei_a[i,1]
  }
  for(i in 1:n_b){
    spei_b[i,2] ~ dnorm(mu_b[i], sd_b)
    mu_b[i] <- b0_b + b1_b * spei_b[i,1]
  }
})

nim_data = list(
  spei_a = spei_a,
  spei_b = spei_b
)
nim_cnst = list(
  n_a = nrow(spei_a),
  n_b = nrow(spei_b)
)

nim_rslt <- nimbleMCMC(
  code = nim_code,
  constants = nim_cnst,
  data = nim_data,
  monitors = c("b0_a", "b0_b", "b1_a", "b1_b", "sd_a", "sd_b"),
  thin = 1, 
  niter = 11000,
  nburnin = 1000,
  nchains = 3
)
mcmcplots::mcmcplot(nim_rslt)

rs <- nim_rslt %>% map(as_tibble) %>% bind_rows()

ln_a <- tibble(yr = spei1$yr, lci = NA, mci = NA, uci = NA, stg = "a")
for(i in 1:nrow(spei_a)){
  mu_a <- spei_a[i,1] * rs$b1_a + rs$b0_a
  ln_a$lci[i] <- quantile(mu_a, 0.025)
  ln_a$mci[i] <- quantile(mu_a, 0.5)
  ln_a$uci[i] <- quantile(mu_a, 0.975)
}


ln_b <- tibble(yr = spei2$yr, lci = NA, mci = NA, uci = NA, stg = "b")
for(i in 1:nrow(spei_b)){
  mu_b <- spei_b[i,1] * rs$b1_b + rs$b0_b
  ln_b$lci[i] <- quantile(mu_b, 0.025)
  ln_b$mci[i] <- quantile(mu_b, 0.5)
  ln_b$uci[i] <- quantile(mu_b, 0.975)
}


ldt <- bind_rows(ln_a, ln_b)
climate_plot <- ggplot(data = ldt, aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci, color = stg, fill = stg), alpha = 0.25) +
  geom_line(data = ldt, aes(x = yr, y = mci, , color = stg)) +
  scale_color_manual(values = c("#b06205", "#000077")) +
  scale_fill_manual(values = c("#b06205", "#000077")) +
  geom_point(dat = spei, aes(x = yr, y = spei_12m)) +
  theme_classic() + theme(legend.position = "none") +
  xlab("Year") + ylab("SPEI") + labs(title = "b - Water balance")

#Now_Pumas======================================================================

cvdt <- readRDS("data//unscaled_covariates.rds") %>%
  select(yr, pd_reconstruction, pd_odfw_est, pd_logistic, pd_full_mean, pd_mortalities)
norm_rslt <- readRDS("results//fbipm_rslt_06feb2026_best_norm.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>% 
  select(grep("cdrng", names(.)))
names(norm_rslt) <- as.character(1988:2023)
get_mnsd <- function(dt, param){
  x <- dt %>%
    filter(yr > 1993 & yr < 2013) %>%
    pull(param)
  out <- c(mn = mean(x), sd = sd(x))
  return(out)
}
rescale_column <- function(dt, param){
  mnsd <- get_mnsd(dt, param)
  x <- dt %>% pull(param)
  out <- (x - mnsd[1]) / mnsd[2]
  return(out)
}

puma <- cvdt %>%
  mutate(
    rcns = rescale_column(., "pd_reconstruction"),
    lgst = rescale_column(., "pd_logistic"),
    odfw = rescale_column(., "pd_odfw_est")
  ) %>%
  select(yr, rcns, lgst, odfw) %>%
  mutate(mnes = rowMeans(select(., 2:4), na.rm = T)) %>%
  filter(yr > 1987) %>%
  pivot_longer(2:ncol(.), names_to = "prm", values_to = "val") %>%
  mutate(parameter = case_when(
    prm == "rcns" ~ "Reconstruction",
    prm == "lgst" ~ "Logistic Growth",
    prm == "odfw" ~ "ODFW Estimate",
    prm == "mnes" ~ "Index Mean"
  ))

lngdt <- norm_rslt %>%
  mutate(stpn = 1:nrow(.)) %>%
  pivot_longer(1:(ncol(.)-1), names_to = "yr", values_to = "val")
rngmn <- lngdt %>% filter(yr > 1993 & yr < 2013) %>% pull(val) %>% mean()
rngsd <- lngdt %>% filter(yr > 1993 & yr < 2013) %>% pull(val) %>% sd()

cgnorm <- lngdt %>%
  mutate(y = (val - rngmn) / rngsd) %>%
  group_by(yr) %>%
  summarise(
    lci = quantile(y, 0.025),
    val = quantile(y, 0.5),
    uci = quantile(y, 0.975)
  ) %>%
  mutate(yr = as.numeric(yr))
cgnormadd <- cgnorm %>%
  mutate(prm = "norm", parameter = "Normal Draw") %>%
  select(yr, prm, val, parameter)
puma_plot <- puma %>%
  bind_rows(cgnormadd) %>%
  filter(prm %in% c("rcns", "odfw", "lgst", "mnes", "norm")) %>%
  ggplot(aes(x = yr, y = val, color = parameter, linetype = parameter)) +
  # geom_point() +
  geom_line() +
  geom_ribbon(
    data = cgnorm, 
    aes(ymin = lci, ymax = uci), 
    alpha = 0.2,
    color = "#ccc",
    fill = "#E69F00",
    linetype = 1,
    linewidth = 0.1
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "Reconstruction" = "#800080",
    "Random" = "#E69F00",
    "ODFW" = "#009E73",
    "Mean" = "#D55E00",
    "Logistic" = "#0072B2"
  )) +
  theme(legend.title = element_blank()) +
  xlab("Year") + ylab("Scaled Value") + labs(title = "a - Puma indices") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(
    legend.position = c(.35, .85))
  

#Now_combine_and_save_plots=====================================================

require(cowplot)
plot_grid(
  puma_plot, climate_plot,
  nrow = 2, ncol = 1,
  rel_heights = c(1, 1)
)
ggsave(
  "figures//main_text_covariate_fig.tiff",
  width = 8,
  height = 14,
  units = "cm", 
  dpi = 600,
  scale = 1.2
)

