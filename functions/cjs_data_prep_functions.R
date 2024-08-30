find_removals <- function(data_list, yr_range){
  
  never_captured <- data_list$capture_history %>%
    select(id, as.character(yr_range[1]:yr_range[2])) %>%
    mutate(n_caps = rowSums(across(where(is.numeric)))) %>%
    filter(n_caps == 0) %>%
    pull(id)
  
  never_in_mains <- data_list$herd_assignment %>%
    select(id, as.character(yr_range[1]:yr_range[2])) %>%
    mutate(
      across(
        as.character(yr_range[1]:yr_range[2]), 
        function(x)as.numeric(x == "main"))) %>%
    mutate(n_main = rowSums(across(where(is.numeric)))) %>%
    filter(n_main == 0) %>%
    pull(id)
}

to_matrix <- function(x){
  x %>%
    arrange(id) %>%
    select(as.character(yr_range[1]:yr_range[2])) %>%
    as.matrix() %>%
    return()
}

