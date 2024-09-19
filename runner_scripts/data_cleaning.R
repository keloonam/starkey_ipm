# Run scripts to clean all of the data necessary for the ipm
# Kenneth Loonam
# August 2024

#Packages=======================================================================
require(dplyr); require(tidyr); require(purrr); require(readr)
require(lubridate); require(nimble)

# Clean the capture and handling data
ah_fp <- "data//animal_handling.csv" # file path for animal handling csv
start_year <- 1987 # Keep these expansive - they are for data prep only
end_year <- 2024 # Only shorten if you want to exclude years we know things from
source("functions//handling_data_prep_functions.R")
source("workflows//data_prep//handling.R")

# Prepare cjs data
yr_range <- c(1988, 2023)
capture_handling_data <- "data//capture_handling_data.rds"
source("workflows//data_prep//cjs_data_prep.R")

# Clean CJS results
# This requires running the cjs model runner script at least once
# "runner_scripts//run_cjs.R"
results_file <- "results//cjs_rslt_06sep2024.rds"
# This option toggles removal of years with questionable estimates
# The years flagged for removal have known causes for bad estimates, such as:
  # e.g. individuals removed without ids recorded
remove_bad_years <- T
source("workflows//data_prep//cjs_summarize_results.R")

# Prepare recruitment data
yr_range <- c(1988, 2023)
capture_handling_data <- "data//capture_handling_data.rds"
source("workflows//data_prep//recruitment_data_prep.R")

# Clean Recruitment results
# This requires running the recruitment model runner script at least once
# runner_scripts//run_recruitment.R
results_file <- "results//recruitment_rslt_10sep2024.rds"
recruitment_data <- "data//recruitment_data.rds" # needed for calendar years
source("workflows//data_prep//recruitment_summarize_results.R")

# Prepare misc. data
# This script hard codes data objects and combines them with data from cap/hand
# The hard coded data is from capture and handling summaries and includes:
  # number of harvested animals
  # management movement numbers
  # minimum known alive (feed-ground attendance)
# Where appropriate, those results are compared to cap/hand results
# If they both contain information on a value, the larger number is used
  # I decided that counting extra elk is harder than missing some in a count
capture_handling_data <- "data//capture_handling_data.rds"
cjs_data <- "data//cjs_data.rds"
cjs_results_summary <- "results//cjs_summary.rds"
yr_range <- c(1988, 2023)
source("workflows//data_prep//misc_data_prep.R")

# Prepare the covariate data
# This script contains hard-coded values to limit file numbers
# It also references multiple .csvs and will fail without those
# This script is also entirely optional if you are not using covariates
source("workflows//data_prep//covariate_data_prep.R")

# Prepare the full IPM data
year_range <- 1988:2023
max_age <- 25
climate_cov <- "pdi_growing" # pdi_growing, pdi_september, temp, precip, ndvi
# PDI = palmer drought index
# Temp, precip, and ndvi are means over may-september
puma_cov <- "pd_logis" # pd_recon, pd_morts, pd_odfwe, pd_logis
# Recon = reconstruction
# Morts = mortalities
# ODFWE = regression estimate from ODWF
# logis = logistic growth model build from reconstruction
# only logis covers the fully 1987 to 2023 timeframe currently
source("workflows//data_prep//ipm_data_bundling.R")

# Prepare the IPM initial values
source("workflows//data_prep//ipm_build_inits.R")
