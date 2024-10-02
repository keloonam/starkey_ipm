'%!in%' <- function(x,y)!('%in%'(x,y))


find_removals <- function(data_list, yr_range){
  
  never_captured <- data_list$capture_history %>%
    select(id, as.character(yr_range)) %>%
    mutate(n_caps = rowSums(across(where(is.numeric)), na.rm = T)) %>%
    filter(n_caps == 0) %>%
    pull(id)
  
  never_in_mains <- data_list$herd_assignment %>%
    select(id, as.character(yr_range)) %>%
    mutate(
      across(
        as.character(yr_range), 
        function(x)as.numeric(x == "main"))) %>%
    mutate(n_main = rowSums(across(where(is.numeric)), na.rm = T)) %>%
    filter(n_main == 0) %>%
    pull(id)
  
  f <- to_matrix_wo_filter(data_list$capture_history, yr_range = yr_range) %>%
    find_first_value_position()
  l <- to_matrix_wo_filter(data_list$harvest_year, yr_range = yr_range) %>%
    find_first_value_position()
  
  first_is_terminal <- data_list$capture_history$id[f==l]
  
  out <- unique(c(
    first_is_terminal,
    never_in_mains,
    never_captured
  ))
  return(out)
}

find_first_value_position <- function(history_matrix){
  history_matrix[,ncol(history_matrix)] <- 1
  out <- apply(history_matrix, 1, function(x) min(which(x == 1)))
  return(out)
}
  

to_matrix_wo_filter <- function(x, yr_range){
  x %>%
    arrange(id) %>%
    select(as.character(yr_range)) %>%
    as.matrix() %>%
    return()
}

to_matrix <- function(x, yr_range, removals){
  x %>%
    arrange(id) %>%
    filter(id %!in% removals) %>%
    select(as.character(yr_range)) %>%
    as.matrix() %>%
    return()
}

build_z <- function(y, f){
  z <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  l_k <- rep(NA, nrow(y))
  for(i in 1:nrow(y)){
    l_k[i] <- max(which(y[i,] == 1))
    z[i, (f[i]):l_k[i]] <- 1
  }
  return(z)
}

to_vector <- function(x, removals, target_column){
  x %>%
    arrange(id) %>%
    filter(id %!in% removals) %>%
    pull(target_column) %>%
    return()
}
