require(rjags); require(tidyverse)
cvdt <- readRDS("data//unscaled_covariates.rds") %>%
  filter(yr > 1987)
cvdt$pd_odfw_est[is.na(cvdt$pd_odfw_est)] <- min(cvdt$pd_odfw_est, na.rm = T)
ipm_rslt <- readRDS("results/ipm_rslt_04feb2026_null_fb.rds")

nelk <- ipm_rslt %>% map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>%
  select(grep("Ntot", names(.))) %>% as.matrix() %>% colMeans()
ncow <- ipm_rslt %>% map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>%
  select(grep("N_f\\[", names(.))) %>% as.matrix() %>% colMeans()
cvdt$nelk <- nelk
cgnorm <- readRDS("results//fbipm_rslt_06feb2026_full_norm.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>% 
  select(grep("cdrng", names(.))) %>%
  as.matrix() %>% colMeans()
cvdt$pd_norm <- cgnorm
scl <- function(x){
  out <- (x - mean(x)) / sd(x)
  return(out)
}

plot_corr <- function(dt, x, y){
  xl <- case_when(
    x == "spei_12m" ~ "SPEI",
    x == "nelk" ~ "Elk Abundance"
  )
  yl <- case_when(
    y == "pd_reconstruction" ~ "Reconstruction",
    y == "pd_odfw_est" ~ "ODFW Estimate",
    y == "pd_full_mean" ~ "Index Mean",
    y == "pd_logistic" ~ "Logistic Growth",
    y == "spei_12m" ~ "SPEI",
    y == "nelk" ~ "Elk Density",
    y == "pd_norm" ~ "Normal Draw"
  )
  xd <- dt %>% pull(x) %>% scl()
  yd <- dt %>% pull(y) %>% scl()
  r2 <- round(summary(lm(yd ~ xd))$r.squared, 3)
  if(xl == "SPEI"){
    out <- dt %>%
      ggplot(aes(x = get(x), y = get(y))) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() + xlab(xl) + ylab(yl) +
      labs(title = bquote({R^2} == .(r2)))
  }else{
    out <- dt %>%
      ggplot(aes(x = get(x), y = scl(get(y)))) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() + xlab(xl) + ylab(yl) +
      labs(title = bquote({R^2} == .(r2)))
  }
  return(out)
}
nelkspei <- cvdt %>% plot_corr(x = "nelk", y = "spei_12m") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank())
nelklgst <- cvdt %>% plot_corr(x = "nelk", y = "pd_logistic") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank())
nelkmean <- cvdt %>% plot_corr(x = "nelk", y = "pd_full_mean") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank())
nelknorm <- cvdt %>% plot_corr(x = "nelk", y = "pd_norm") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank())
nelkodfw <- cvdt %>% plot_corr(x = "nelk", y = "pd_odfw_est") + 
  xlab(element_blank()) + theme(axis.text.x = element_blank())
nelkrcns <- cvdt %>% plot_corr(x = "nelk", y = "pd_reconstruction")

speilgst <- cvdt %>% plot_corr(x = "spei_12m", y = "pd_logistic") + 
  xlab(element_blank()) + theme(axis.text = element_blank()) + 
  ylab(element_blank())
speimean <- cvdt %>% plot_corr(x = "spei_12m", y = "pd_full_mean") + 
  xlab(element_blank()) + theme(axis.text = element_blank()) + 
  ylab(element_blank())
speinorm <- cvdt %>% plot_corr(x = "spei_12m", y = "pd_norm") + 
  xlab(element_blank()) + theme(axis.text = element_blank()) + 
  ylab(element_blank())
speiodfw <- cvdt %>% plot_corr(x = "spei_12m", y = "pd_odfw_est") + 
  xlab(element_blank()) + theme(axis.text = element_blank()) + 
  ylab(element_blank())
speircns <- cvdt %>% plot_corr(x = "spei_12m", y = "pd_reconstruction") +
  ylab(element_blank()) + theme(axis.text.y = element_blank())

emptyplt <- ggplot(data = NULL, aes(x = 1:10, y = 1:10)) + theme_classic() +
  ylab(element_blank()) + xlab(element_blank()) +
  theme(
    axis.line = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank())
require(cowplot)
plot_grid(
  nelkspei, emptyplt,
  nelklgst, speilgst,
  nelkmean, speimean,
  nelknorm, speinorm,
  nelkodfw, speiodfw,
  nelkrcns, speircns,
  nrow = 6, ncol = 2
)
ggsave(
  "figures//ipm_covariate_correlation.png",
  width = 18,
  height = 25,
  units = "cm", 
  dpi = 600,
  scale = 1.25
)
rm(list = ls())

