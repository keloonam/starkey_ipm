# Build the initial values for the ipm based on the data
dtlist <- readRDS("data//ipm_data.rds")
inits <- list(
  sd_R = 0.1,
  R_YR = c(NA, rep(0, dtlist$nyr-1)),
  sd_SC = 0.1,
  sd_SF = 0.1,
  sd_SM = 0.1,
  SC_B0 = 3,
  SM_B0 = 3,
  SF_B0 = 3,
  SC_YR = c(NA, rep(0, dtlist$nyr-1)),
  SF_YR = c(NA, rep(0, dtlist$nyr-1)),
  SM_YR = c(NA, rep(0, dtlist$nyr-1))
)

saveRDS(inits, "data//ipm_inits.rds")
rm(list = ls()[-which(ls() == "yr_range")])