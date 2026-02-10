require(tidyverse)

scl <- function(x){
  ((x-mean(x, na.rm = T))/sd(x, na.rm = T)) %>% return()
}
build_age_at_harvest_data <- function(chdt, use_years, n_age){
  n_har <- array(data = 0, dim = c(n_age, 2, length(use_years)))
  for(i in 1:nrow(chdt$sex)){
    for(j in 1:length(use_years)){
      k <- which(names(chdt$annual_age) == use_years[j])
      if(chdt$harvest_year[i,k] == 1){
        if(chdt$annual_age[i,k] > n_age){
          a <- n_age
        }else(
          a <- pull(chdt$annual_age[i,k])
        )
        # print(c(i,j))
        x <- 1+(chdt$sex$sex[i]=="M")
        n_har[a,x,j] <- n_har[a,x,j] + 1
      }
    }
  }
  n_har[1,,] <- 0
  return(n_har)
}
logit <- function(x){
  out <- log(x/(1-x))
  return(out)
}
expit <- function(x){
  out <- 1/(1+exp(-x))
  return(out)
}