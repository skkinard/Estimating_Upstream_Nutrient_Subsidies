# SC_mixing_model_calc_region
# written by: sean Kinard
# last edit: 2023-09-17
#----------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------
source('03_public/toolkit.R') # load packages and helper-functions
library(simmr)

d <- read_csv('03_public/output/isotope_CNS_2020_01_clean.csv')

d <- d %>% 
  filter(! species %in% bad_sources) %>%
  filter(species !='Periphyton') %>%
  add_sitevars() %>%
  mutate(distcat = ifelse(baydist_km<25, '<25 km', '>25 km')) %>%
  add_is_dam()

source("03_public/calculation/20_SC_mixing_model_function_region.R")
source("03_public/calculation/20_SC_mixing_model_function_calallen.R")

#----------------------------------------------------------------------------
# CS mixing model: site_code
#----------------------------------------------------------------------------
# blank table
mmix_site <- tibble(
  deviance = numeric(),
  Fresh = numeric(),
  Estuarine = numeric(),
  m_group = character(),
  site = character(),
  statistic = character() )

mmix_site <- run_simmr_region(my_data = d,
                              my_group = 'site_code') %>%
  full_join(mmix_site) %>%
  mutate(dataset = 'site')

#----------------------------------------------------------------------------
# CS mixing model: estuary distance categories
#----------------------------------------------------------------------------
# blank table
mmix_distcat <- tibble(
  deviance = numeric(),
  Fresh = numeric(),
  Estuarine = numeric(),
  m_group = character(),
  site = character(),
  statistic = character() )

mmix_distcat <- run_simmr_region(my_data = d,
                          my_group = 'distcat') %>%
    full_join(mmix_distcat) %>%
  mutate(dataset = 'distcat') 

#----------------------------------------------------------------------------
# CS mixing model: transient_type
#----------------------------------------------------------------------------
# blank table
mmix_transient <- tibble(
  deviance = numeric(),
  Aquatic = numeric(),
  Terrestrial = numeric(),
  m_group = character(),
  site = character(),
  statistic = character() )

mmix_transient <- run_simmr_region(my_data = filter(d, is_dam == 'No Dam'),
                            my_group = 'is_diadromous') %>%
    full_join(mmix_transient) %>%
  mutate(dataset = 'is_diadromous')

#----------------------------------------------------------------------------
# CS mixing model: Callallen
#----------------------------------------------------------------------------
# blank table
mmix_calallen <- tibble(
  deviance = numeric(),
  Aquatic = numeric(),
  Terrestrial = numeric(),
  m_group = character(),
  site = character(),
  statistic = character() )

mmix_calallen <- run_simmr_calallen(my_data = d,
                               my_group = 'is_dam') %>%
  mutate(dataset = 'calallen')

#----------------------------------------------------------------------------
# CS mixing model: distant Dams
#----------------------------------------------------------------------------
# blank table
mmix_ddam <- tibble(
  deviance = numeric(),
  Aquatic = numeric(),
  Terrestrial = numeric(),
  m_group = character(),
  site = character(),
  statistic = character() )

mmix_ddam <- run_simmr_region(my_data = filter(d, ! site_code %in% c('UN', 'LN')),
                                   my_group = 'is_dam') %>%
  full_join(mmix_ddam) %>%
  mutate(dataset = 'calallen',
         m_group = ifelse(m_group == 'Dam', 'Above Other', 'None'))

#----------------------------------------------------------------------------
# Combine all mmix
#----------------------------------------------------------------------------
mmix_all <- full_join(mmix_site, mmix_distcat) %>%
  full_join(mmix_transient) %>%
  full_join(mmix_calallen) %>%
  full_join(mmix_ddam) %>%
  select(Estuarine, statistic, m_group, dataset) %>%
  mutate(Estuarine = 100*Estuarine) %>%
  pivot_wider(names_from = statistic, values_from = Estuarine)

write_csv(mmix_all, '03_public/output/CS_mix_out_region.csv')
#----------------------------------------------------------------------------
# End SC_mixing_model_calc_universal