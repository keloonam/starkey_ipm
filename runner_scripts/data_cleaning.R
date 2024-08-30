# Run scripts to clean all of the data necessary for the ipm
# Kenneth Loonam
# August 2024

#Packages=======================================================================
require(dplyr); require(tidyr); require(purrr); require(readr)
require(lubridate)

# Clean the capture and handling data
ah_fp <- "data//animal_handling.csv" # file path for animal handling csv
start_year <- 1987 # Keep these expansive - they are for data prep only
end_year <- 2024 # Only shorten if you want to exclude years we know things from
source("functions//handling_data_prep_functions.R")
source("workflows//data_prep//handling.R")

# Prepare cjs data
yr_range <- c(1988, 2023)
