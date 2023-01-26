#Header=========================================================================
# Kenneth Loonam
# January 2023
# Figures to support IPM publication

#Environment====================================================================
# Packages
require(dplyr); require(ggplot2); require(ggsci); require(rjags)

# Results to load
load("results//ipm_result_05jan2023_R_pdis.Rdata")

# Derived result summaries
data <- summary(rslt)
q    <- data$quantiles

#To Do List=====================================================================

# IPM explanation figure
# Demographic rates through time, multi-panel figure
# Posterior Distributions
# Lambda regression confidence intervals
# Convergence plots

#