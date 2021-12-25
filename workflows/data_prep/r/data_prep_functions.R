# model data prep functions
# Kenneth Loonam
# February 2020

ch_init_fn <- function(ch, f){
  for(i in 1:nrow(ch)){
    ch[i,1:f[i]] <- NA
  }
  out <- ifelse(!is.na(ch), 1, ch)
  
  return(out)
}
