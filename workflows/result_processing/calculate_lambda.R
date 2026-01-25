source("functions//calculate_lambda_functions.R")
tags <- paste0("bst_", c("lgst", "mean", "norm", "odfw", "rcns"))
for(tag in tags){
  save_file_name <- paste0("results//full_lambda_tibble_", tag, ".rds")
  tmls_file_name <- paste0("results//transition_matrix_list_", tag, ".rds")
  mtdt_file_name <- paste0("data//tm_list_metadata_", tag, ".rds")
  
  tmna_list <- readRDS(tmls_file_name)
  meta_data <- readRDS(mtdt_file_name)
  
  plan(multisession, workers = 16)
  lambda_results <- future_map(
    .x = tmna_list,
    .f = calculate_lambdas
  ) %>%
    bind_rows()
  cbind(lambda_results, meta_data) %>%
    saveRDS(save_file_name)
}