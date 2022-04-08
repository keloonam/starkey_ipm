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

#combine feedground and capture and handling into one min n for elk and mule deer