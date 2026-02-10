require(tidyverse); require(rjags); require(cowplot)
tag <- "full_rcns"
rslt_file <- paste0("results//fbipm_rslt_06feb2026_", tag, ".rds")
save_file <- paste0("figures//ipm_residual_plots_", tag, ".png")

frs <- readRDS(rslt_file) %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
fdt <- readRDS("data//ipm_data_03feb2026.rds")

cgcov <- case_when(
  tag == "full_lgst" ~ fdt$cg_logistic_growth,
  tag == "full_mean" ~ fdt$cg_mean_estimate,
  tag == "full_norm" ~ NA,
  tag == "full_odfw" ~ fdt$cg_odfw_estimate,
  tag == "full_rcns" ~ fdt$cg_reconstruction
)
if(tag == "full_norm"){
  cgcov <- readRDS(rslt_file) %>%
    map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>%
    select(grep("cdrng", names(.))) %>%
    as.matrix() %>% colMeans()
}

cgrange <- c(min(cgcov), max(cgcov))
wtrange <- c(min(fdt$spei12), max(fdt$spei12))
ddrange <- c(min(fdt$nelk), max(fdt$nelk))

cgxl <- c(cgrange[1] - diff(cgrange)*0.01, cgrange[2] + diff(cgrange)*0.01)
ddxl <- c(ddrange[1] - diff(ddrange)*0.01, ddrange[2] + diff(ddrange)*0.01)
wtxl <- c(wtrange[1] - diff(wtrange)*0.01, wtrange[2] + diff(wtrange)*0.01)

ppcolor <- "#0072B2"
sncolor <- "#D55E00"
sccolor <- "#800080"
sfcolor <- "#009E73"

#c(t, r, b, l) - travels clockwise from top
xapm <- c(0,   0.1, 0.2, -.1)
yapm <- c(0.2, 0,   -0.2,   0.2)
appm <- c(0.2, 0.1, -.2, -0.1)

rename_cols <- function(dtx, new_names){
  stopifnot(ncol(dtx) == length(new_names))
  names(dtx) <- new_names
  return(dtx)
}

pull_vr <- function(dtx, str_pattern, new_names = as.character(1989:2023)){
  dtx %>%
    select(grep(str_pattern, names(.))) %>%
    rename_cols(new_names) %>%
    mutate(stpn = 1:nrow(.)) %>%
    pivot_longer(1:(ncol(.)-1), names_to = "yr", values_to = "val")
}

build_plot_data <- function(dtx, str_pattern, ipm_dt = fdt, drop_year = 1){
  dtx %>% pull_vr(str_pattern) %>%
    mutate(yr = as.numeric(yr)) %>%
    summarise(
      lci = quantile(val, 0.025),
      mci = quantile(val, 0.5),
      uci = quantile(val, 0.975),
      .by = yr
    ) %>%
    mutate(dd = ipm_dt$nelk[-drop_year]) %>%
    mutate(cg = cgcov[-drop_year]) %>%
    mutate(wt = ipm_dt$spei12[-drop_year]) %>%
    mutate(wm = ipm_dt$spei12[-length(ipm_dt$spei12)]) %>%
    return()
}

plot_residual <- function(plot_data, covariate, clr, xl, appm = appm){
  plot_data %>%
    mutate(xvar := {{covariate}}) %>%
    ggplot(aes(x = xvar, y = mci, ymin = lci, ymax = uci)) +
    geom_pointrange(color = clr) +
    theme_classic() +
    xlim(xl) + #ylim(yl) +
    xlab(element_blank()) + ylab(element_blank()) +
    theme(axis.text = element_blank()) +
    geom_hline(yintercept = 0, linetype = "11", color = "#bbb") %>%
    return()
}
ppdt <- frs %>% build_plot_data("P_Byr")
sndt <- frs %>% build_plot_data("SN_Byr")
scdt <- frs %>% build_plot_data("SC_Byr")
sfdt <- frs %>% build_plot_data("SF_Byr")

ppcg <- plot_residual(ppdt, cg, ppcolor, cgxl) +
  theme(plot.margin = unit(appm, "cm"))
ppwm <- plot_residual(ppdt, wm, ppcolor, wtxl) +
  theme(plot.margin = unit(appm, "cm"))
ppdd <- plot_residual(ppdt, dd, ppcolor, ddxl) +
  theme(plot.margin = unit(appm, "cm"))

sncg <- plot_residual(sndt, cg, sncolor, cgxl) +
  theme(plot.margin = unit(appm, "cm"))
snwt <- plot_residual(sndt, wt, sncolor, wtxl) +
  theme(plot.margin = unit(appm, "cm"))
snwm <- plot_residual(sndt, wm, sncolor, wtxl) +
  theme(plot.margin = unit(appm, "cm"))
sndd <- plot_residual(sndt, dd, sncolor, ddxl) +
  theme(plot.margin = unit(appm, "cm"))

sccg <- plot_residual(scdt, cg, sccolor, cgxl) +
  theme(plot.margin = unit(appm, "cm"))
scwm <- plot_residual(scdt, wm, sccolor, wtxl) +
  theme(plot.margin = unit(appm, "cm"))
scwt <- plot_residual(scdt, wt, sccolor, wtxl) +
  theme(plot.margin = unit(appm, "cm"))
scdd <- plot_residual(scdt, dd, sccolor, ddxl) +
  theme(plot.margin = unit(appm, "cm"))

sfcg <- plot_residual(sfdt, cg, sfcolor, cgxl) +
  theme(plot.margin = unit(appm, "cm"))
sfwm <- plot_residual(sfdt, wm, sfcolor, wtxl) +
  theme(plot.margin = unit(appm, "cm"))
sfwt <- plot_residual(sfdt, wt, sfcolor, wtxl) +
  theme(plot.margin = unit(appm, "cm"))
sfdd <- plot_residual(sfdt, dd, sfcolor, ddxl) +
  theme(plot.margin = unit(appm, "cm"))

ppyp <- ppdt %>% ggplot(aes(x = cg, y = mci, ymin = lci, ymax = uci)) +
  theme_classic() +
  xlab(element_blank()) +
  theme(
    axis.text.x = element_blank(), 
    axis.line = element_blank(), 
    axis.ticks = element_blank()
  ) + ylab("Pregnancy") +
  theme(plot.margin = unit(yapm, "cm"))
snyp <- sndt %>% ggplot(aes(x = cg, y = mci, ymin = lci, ymax = uci)) +
  theme_classic() +
  xlab(element_blank()) +
  theme(
    axis.text.x = element_blank(), 
    axis.line = element_blank(), 
    axis.ticks = element_blank()
  ) + ylab("Neonate Survival") +
  theme(plot.margin = unit(yapm, "cm"))
scyp <- scdt %>% ggplot(aes(x = cg, y = mci, ymin = lci, ymax = uci)) +
  theme_classic() +
  xlab(element_blank()) +
  theme(
    axis.text.x = element_blank(), 
    axis.line = element_blank(), 
    axis.ticks = element_blank()
  ) + ylab("Calf Survival") +
  theme(plot.margin = unit(yapm, "cm"))
sfyp <- sfdt %>% ggplot(aes(x = cg, y = mci, ymin = lci, ymax = uci)) +
  theme_classic() +
  xlab(element_blank()) +
  theme(
    axis.text.x = element_blank(), 
    axis.line = element_blank(), 
    axis.ticks = element_blank()
  ) + ylab("Adult Survival") +
  theme(plot.margin = unit(yapm, "cm"))

cgxp <- ppdt %>% ggplot(aes(x = cg, y = mci, ymin = lci, ymax = uci)) +
  theme_classic() +
  ylab(element_blank()) +
  theme(
    axis.text.y = element_blank(), 
    axis.line = element_blank(), 
    axis.ticks = element_blank()
  ) + xlab("Puma Density") + xlim(cgxl) +
  theme(plot.margin = unit(xapm, "cm"))
ddxp <- ppdt %>% ggplot(aes(x = dd, y = mci, ymin = lci, ymax = uci)) +
  theme_classic() +
  ylab(element_blank()) +
  theme(
    axis.text.y = element_blank(), 
    axis.line = element_blank(), 
    axis.ticks = element_blank()
  ) + xlab("Elk Density") + xlim(ddxl) +
  theme(plot.margin = unit(xapm, "cm"))
wtxp <- ppdt %>% ggplot(aes(x = dd, y = mci, ymin = lci, ymax = uci)) +
  theme_classic() +
  ylab(element_blank()) +
  theme(
    axis.text.y = element_blank(), 
    axis.line = element_blank(), 
    axis.ticks = element_blank()
  ) + xlab("SPEI (t)") + xlim(wtxl) +
  theme(plot.margin = unit(xapm, "cm"))
wmxp <- ppdt %>% ggplot(aes(x = dd, y = mci, ymin = lci, ymax = uci)) +
  theme_classic() +
  ylab(element_blank()) +
  theme(
    axis.text.y = element_blank(), 
    axis.line = element_blank(), 
    axis.ticks = element_blank()
  ) + xlab("SPEI (t-1)") + xlim(wtxl) +
  theme(plot.margin = unit(xapm, "cm"))
blxy <- ggplot() + 
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

plot_grid(
  ppyp, ppdd, blxy, ppwm, ppcg,
  snyp, sndd, snwt, snwm, sncg,
  scyp, scdd, scwt, scwm, sccg,
  sfyp, sfdd, sfwt, sfwm, sfcg,
  blxy, ddxp, wtxp, wmxp, cgxp,
  rel_widths = c(0.25, 1, 1, 1, 1),
  rel_heights = c(1, 1, 1, 1, 0.25),
  align = "none",
  ncol = 5,
  label_size = 2)

ggsave(
  save_file,
  width = 18,
  height = 18,
  scale = 1.2,
  units = "cm",
  dpi = 600)

rm(list = ls())
