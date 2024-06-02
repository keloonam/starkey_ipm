tau <- 0.001
rnorm(10000, mean = 0, sd = 100) %>% abs() %>% hist()

require(purrr)

rs <- rslt %>% 
  map(as_tibble) %>%
  bind_rows() 
rs %>%
  pull(`SF_B0_mean`) %>% unique()
