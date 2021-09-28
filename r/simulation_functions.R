# Simulation Functions
# Kenneth Loonam
# February 2020

#Simulate_Population============================================================

sim_pop_fn <- function(
  n_af, n_am, n_jf, n_jm, af_s, am_s, jf_s, jm_s, recruit, n_years
  ){
  
  n_tot <- sum(n_af, n_am, n_jf, n_jm)
  
  population <- tibble(
    id = 1:n_tot,
    male = c(
      rep(0, n_af), 
      rep(1, n_am), 
      rep(0, n_jf), 
      rep(1, n_jm)),
    stat_0 = c(
      rep("adult", sum(n_af, n_am)), 
      c(rep("juvenile", sum(n_jf, n_jm))))
  ) %>%
    bind_cols(as_tibble(matrix(
      as.character(NA), 
      nrow = n_tot, 
      ncol = n_years, 
      dimnames = list(NULL, paste0("stat_", as.character(1:n_years)))
    )))
  
  for(i in 1:n_years){
    
    population <- population %>%
      mutate(s_p = case_when(
        (.[i+2] == "dead") ~ 0,
        (.[i+2] ==    "adult" & .$male == 0) ~ af_s,
        (.[i+2] ==    "adult" & .$male == 1) ~ am_s,
        (.[i+2] == "juvenile" & .$male == 0) ~ jf_s,
        (.[i+2] == "juvenile" & .$male == 1) ~ jm_s
      )) %>%
      mutate(survived = rbinom(nrow(.), 1, .$s_p))
    
    for(j in 1:nrow(population)){
      if(population$survived[j] == 1){
        population[j,i+3] <- "adult"
      }else{
        population[j,i+3] <- "dead"
      }
    }
    
    n_born <- population %>%
      filter(male == 0) %>%
      filter(.[i+2] == "adult") %>%
      nrow() %>%
      rbinom(n = 1, size = ., prob = recruit)
    
    new_ids <- (max(population$id) + 1):(max(population$id) + n_born)
    
    if(n_born > 0){
      
      new_histories <- tibble(
        id = new_ids,
        male = rbinom(length(new_ids), 1, 0.5),
        stat_0 = rep(NA, n_born)
      ) %>%
        bind_cols(as_tibble(matrix(
          NA, 
          nrow = n_born, 
          ncol = n_years, 
          dimnames = list(NULL, paste0("stat_", as.character(1:n_years)))
        )))
      
      new_histories[,i + 3] <- "juvenile"
      
      population <- population %>%
        select(-s_p, -survived) %>%
        bind_rows(new_histories)
      
    }else{
      population <- population %>%
        select(-s_p, -survived)
    }
  }
  return(population)
}

#Simulate_feed_ground_captures==================================================

is_available <- function(x){
  out <- (!is.na(x) & x != "dead")
  return(out)
}

sim_cap_fn <- function(population, cap_p){
  
  availability <- population %>%
    select(-male, -id) %>%
    as.matrix() %>%
    is_available()
  
  dimnames(availability) <- list(population$id, dimnames(availability[[2]]))
  
  out <- availability * matrix(
    rbinom(
      n = length(availability),
      size = 1,
      p = cap_p
    ),
    nrow = nrow(availability),
    ncol = ncol(availability)
    )
  
  keep <- (apply(out, 1, sum) > 0)
  
  return(out[keep,])
}




