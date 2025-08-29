
library(tidyverse)
library(broom)

divalfa <- diversidad_alfa

# Simpson
simpson <- dplyr::select(divalfa, c(Genotype, Simpson))
simpson %>% 
  t.test(Simpson~Genotype, alternative = 'greater', data=.) %>% 
  tidy() %>% 
  pull(statistic) -> t_original
replicate(n = 10000,
          expr = {
            simpson %>% 
              mutate(Genotype = sample(Genotype)) %>% 
              t.test(Simpson~Genotype, alternative = 'greater', data=.) %>%
              tidy() %>% 
              pull(statistic)
          }) -> t_permu
(sum(t_permu > t_original)+1)/10001

# Shannon
shannon <- dplyr::select(divalfa, c(Genotype, Shannon))
shannon %>% 
  t.test(Shannon~Genotype, alternative = 'greater', data=.) %>% 
  tidy() %>% 
  pull(statistic) -> t_original
replicate(n = 10000,
          expr = {
            shannon %>% 
              mutate(Genotype = sample(Genotype)) %>% 
              t.test(Shannon~Genotype, alternative = 'greater', data=.) %>%
              tidy() %>% 
              pull(statistic)
          }) -> t_permu
(sum(t_permu > t_original)+1)/10001
