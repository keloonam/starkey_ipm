# Build the initial values for the ipm based on the data
dtlist <- readRDS("data//ipm_data.rds")



inits <- list(
  sd_R = 0.1,
  R_B0 = 0,
  R_WT = 0,
  R_WM = 0,
  R_CG = 0,
  R_DD = 0,
  R_YR = c(NA, rep(0, dtlist$nyr-1)),
  sd_SC = 0.1,
  sd_SF = 0.1,
  sd_SM = 0.1,
  SC_B0 = 3,
  SM_B0 = 3,
  SF_B0 = 3,
  SC_WT = 0,
  SC_WM = 0,
  SC_CG = 0,
  SC_DD = 0,
  SC_YR = c(NA, rep(0, dtlist$nyr-1)),
  SF_YR = c(NA, rep(0, dtlist$nyr-1)),
  SM_YR = c(NA, rep(0, dtlist$nyr-1))
)

saveRDS(inits, "data//ipm_inits.rds")
rm(list = ls())