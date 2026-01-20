require(dplyr); require(mcmcplots); require(purrr); require(rjags)
fd <- list(
  r1 = readRDS("results//ipm//tests//08july_rsltA1.rds"),
  r2 = readRDS("results//ipm//tests//08july_rsltA2.rds"),
  r3 = readRDS("results//ipm//tests//08july_rsltA3.rds")
) 
# mcmcplots::mcmcplot(fd)

add_run_id <- function(x, id){
  x %>%
    mutate(run_id = id) %>%
    return()
}

rslt <- fd %>%
  map(as.matrix) %>% map(as_tibble) %>%
  map2(.x = ., .y = c(1,2,3), .f = add_run_id) %>%
  bind_rows()
