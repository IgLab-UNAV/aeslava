
library(tidyverse)
library(broom)

divalfa <- diversidad_alfa

# Simpson
simpson <- dplyr::select(divalfa, c(Group, Simpson))
simpson %>% 
  aov(Simpson~Group, data=.) %>% 
  tidy() %>%
  slice(1) %>% 
  pull(statistic) -> t_original
replicate(n = 10000,
          expr = {
            simpson %>% 
              mutate(Group = sample(Group)) %>% 
              aov(Simpson~Group, data=.) %>%
              tidy() %>% 
              slice(1) %>% 
              pull(statistic)
          }) -> t_permu
(sum(t_permu > t_original)+1)/10001

# Shannon
shannon <- dplyr::select(divalfa, c(Group, Shannon))
shannon %>% 
  aov(Shannon~Group, data=.) %>% 
  tidy() %>% 
  slice(1) %>% 
  pull(statistic) -> t_original
replicate(n = 10000,
          expr = {
            shannon %>% 
              mutate(Group = sample(Group)) %>% 
              aov(Shannon~Group, data=.) %>%
              tidy() %>% 
              slice(1) %>% 
              pull(statistic)
          }) -> t_permu
(sum(t_permu > t_original)+1)/10001
