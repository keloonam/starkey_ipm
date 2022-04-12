# Project runner script
# Kenneth Loonam
# March 2022

#Dependencies===================================================================

# tidyverse
# csv of animal and animal handling tables saved as "data//handling.csv"

#Elk_data_cleaning==============================================================

# vector of years from t = 1 to n -- real world year values
years <- 1988:2021

# code for the species in the capture and handling database
species <- "E"

# save file, all other data prep scripts point to this file, hard coded
elk_data_location <- "data//elk_data.Rdata"

# Run script to build elk data tibbles and save as a list in data folder
source("workflows//data_prep//elk_data_cleaning.R")

#Harvest_data===================================================================

target_herd <- "main"
years <- 1988:2021
source("workflows//data_prep//elk_harvest.R")
rm(list = ls())

#Elk_min_n======================================================================

target_herd <- "main"
years <- 1988:2021
elk_data_location <- "data//elk_data.Rdata"
feedground_data_location <- "data//min_n_handle_summaries.csv"
data_destination <- "data//elk_minimum_count_data.Rdata"
source("workflows//data_prep//elk_min_n.R")
rm(list = ls())

#CJS_model_run==================================================================

# do this in CJS script!!
# "source(workflows//model_runs//cjs_model_run_nimble_cip.R)"
# rm(list = ls())

#Ratio_recruitment_data=========================================================

source("workflows//model_runs//recruitment_model_run.R")
rm(list = ls())

#IPM_data_prep==================================================================

source("workflows//data_prep//elk_ipm_data_prep.R")
rm(list = ls())
# You are on the ratio data. You need to filter the years with bad estimates.
# After that, you are moving on through the ipm data prep.