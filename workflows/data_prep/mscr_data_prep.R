# Prepare data for the cjs run
# Kenneth Loonam
# August 2024

#Environment====================================================================

source("functions//cjs_data_prep_functions.R")
dtlist <- readRDS(capture_handling_data)

ids_to_remove <- find_removals(
  data_list = dtlist, 
  yr_range = yr_range)

ya <- to_matrix(
  x = dtlist$capture_history, 
  yr_range = yr_range, 
  removals = ids_to_remove
)

yh <- to_matrix(
  x = dtlist$harvest_year, 
  yr_range = yr_range, 
  removals = ids_to_remove
)

yt <- ((ya + yh) > 0)*1

# f <- find_first_value_position(y_old)

# z <- build_z(y_old, f)

male <- as.numeric(to_vector(
  dtlist$sex, 
  removals = ids_to_remove, 
  target_column = "sex"
) == "M")

# herd <- (to_matrix(
#   dtlist$herd_assignment,
#   yr_range = yr_range,
#   removals = ids_to_remove
# ) != "main") * 1

yc <- (to_matrix(
  dtlist$annual_age,
  yr_range = yr_range,
  removals = ids_to_remove
) == 1) * 1

# l <- dtlist$harvest_year %>%
#   to_matrix(yr_range = yr_range, removals = ids_to_remove) %>%
#   find_first_value_position()

ma.af <- matrix(
  0, 
  nrow = ncol(yt) - 1, 
  ncol = ncol(yt) * 2 - 1)
dimnames(ma.af) <- list(
  as.character(1988:2022), 
  c(paste0(sort(rep(as.character(1989:2023), 2)), c("a","h")), "unseen")
  )
ma.am <- ma.jf <- ma.jm <- ma.af
for(i in 1:nrow(yt)){
  # Start here Kenneth
  # You need to write custom m-arrays for each age/sex class
  # The marray function from the ipm book might be useful
  sn <- which(yt[i,] != 0)
  suppressWarnings(hr <- min(which(yh[i,] != 0)))
  hr[hr==Inf] <- 0
  for(j in 1:(length(sn)-1)){
    rid <- sn[j]
    cid <- 2 * (sn[j+1] - 1) - 1 + (sn[j+1] == hr)
    
    afi <- (1-male[i]) * (1-yc[i,rid])
    ami <-    male[i]  * (1-yc[i,rid])
    jfi <- (1-male[i]) *    yc[i,rid]
    jmi <-    male[i]  *    yc[i,rid]
    
    ma.af[rid,cid] <- ma.af[rid,cid] + afi
    ma.am[rid,cid] <- ma.am[rid,cid] + ami
    ma.jf[rid,cid] <- ma.jf[rid,cid] + jfi
    ma.jm[rid,cid] <- ma.jm[rid,cid] + jmi
    if(sum(c(afi, ami, jfi, jmi)) > 1){
      print("Warning: double counting!!!")
    }
  }
}
r.af <- rep(0, nrow(ma.af))
names(r.af) <- as.character(1988:2022)
r.am <- r.jf <- r.jm <- r.af
for(t in 1:nrow(ma.af)){
  r.af[t] <- sum(ya[,t] * (1-male) * (1-yc[,t]))
  r.am[t] <- sum(ya[,t] *    male  * (1-yc[,t]))
  r.jf[t] <- sum(ya[,t] * (1-male) *    yc[,t])
  r.jm[t] <- sum(ya[,t] *    male  *    yc[,t])
}
ma.af[,ncol(ma.af)] <- r.af - rowSums(ma.af[,-ncol(ma.af)])
ma.am[,ncol(ma.am)] <- r.am - rowSums(ma.am[,-ncol(ma.am)])
ma.jf[,ncol(ma.jf)] <- r.jf - rowSums(ma.jf[,-ncol(ma.jf)])
ma.jm[,ncol(ma.jm)] <- r.jm - rowSums(ma.jm[,-ncol(ma.jm)])

nt <- nrow(ma.af) + 1
full_cjs_data <- list(
  constants = list(
    nt = nt
  ),
  data = list(
    ma.af = ma.af,
    ma.am = ma.am,
    ma.jf = ma.jf,
    ma.jm = ma.jm,
    r.af = r.af,
    r.am = r.am,
    r.jf = r.jf,
    r.jm = r.jm
  ),
  initial_values = list(
    Pf = rep(0.5, nt-1),
    Pm = rep(0.5, nt-1),
    Snaf = rep(0.5, nt-1),
    Snam = rep(0.5, nt-1),
    Snca = rep(0.5, nt-1),
    Shaf = rep(0.5, nt-1),
    Sham = rep(0.5, nt-1),
    Shjf = rep(0.5, nt-1),
    Shjm = rep(0.5, nt-1)
  )
)

saveRDS(full_cjs_data, "data//cjs_data.rds")
rm(list = ls()[-which(ls() == "yr_range")])
