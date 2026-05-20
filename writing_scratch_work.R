nft <- rsraw %>%
  select(grep("N_f\\[", names(.))) %>%
  as.matrix()
nfo <- rsraw %>%
  select(grep("N_fo\\[", names(.))) %>%
  as.matrix()


prop_old <- nfo/nft
p_diff <- matrix(NA, nrow = ncol(prop_old), ncol = ncol(prop_old))
# for(i in 1:(ncol(prop_old)-1)){
#   for(j in (i+1):ncol(prop_old)){
#     p_diff[i,j] <- mean((prop_old[,i] - prop_old[,j]) > 0, na.rm = T)
#   }
# }


old_elk <- tibble(
  yr = 1988:2023,
  i025 = apply(prop_old, 2, quantile, probs = 0.025),
  i500 = apply(prop_old, 2, quantile, probs = 0.5),
  i975 = apply(prop_old, 2, quantile, probs = 0.975)
)
ggplot(old_elk, aes(x = yr, y = i500)) + 
  geom_ribbon(aes(ymin = i025, ymax = i975), alpha = 0.25) +
  geom_line()

age <- tmp$annual_age[,3:38] %>% as.matrix()
sex <- tmp$sex %>% pull(sex) %>% (function(x)(x == "F")*1)
kal <- tmp$known_alive[,3:38] %>% as.matrix()
inm <- tmp$herd_assignment[,3:38] %>% as.matrix() %>% (function(x)(x == "main")*1)

old_elk <- (age*sex*kal*inm) %>% 
  as_tibble() %>% 
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(!is.na(val) & val > 0) %>%
  mutate(old = val > 14) %>%
  group_by(yr) %>%
  summarise(
    mean_age = mean(val),
    sd_age = sd(val),
    prop_old = mean(old)
  ) %>%
  print(n = 35) %>%
  pull(prop_old)


x <- runif(100000, 0, 100) %>% abs() %>% (function(x)1/x^2)
hist(x)
tau <-quantile(x, 0.1)
rnorm(100000, 0, 1/sqrt(tau)) %>% hist()
