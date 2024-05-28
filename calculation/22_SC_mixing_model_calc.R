# CH_mixing_model_calc
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
  filter(species !='Periphyton')

source("03_public/calculation/20_SC_mixing_model_function.R")
#----------------------------------------------------------------------------
# CS mixing model: Guilds
#----------------------------------------------------------------------------
# blank table
mmix_guild <- tibble(
  deviance = numeric(),
  Fresh = numeric(),
  Estuarine = numeric(),
  m_group = character(),
  site = character(),
  statistic = character() )

# sites
msites <- d %>% add_sitevars() %>% filter(site_type != 'Estuary') %>%
  pull(site_code) %>% unique()

#  test <- run_simmr(my_data=d, my_location='TR', my_group='order')

# loop mixing models across sites and aggregate statistics and quantile info
for(i in msites) {
  mmix_guild <- run_simmr(my_data = d,
                          my_location = i,
                          my_group = 'order') %>%
    full_join(mmix_guild) }

mmix_order <- mmix_guild %>%
  mutate(dataset = 'order')

#----------------------------------------------------------------------------
# CS mixing model: transient_type
#----------------------------------------------------------------------------
# blank table
mmix_trophic <- tibble(
  deviance = numeric(),
  Aquatic = numeric(),
  Terrestrial = numeric(),
  m_group = character(),
  site = character(),
  statistic = character() )

# loop mixing models across sites and aggregate statistics and quantile info
for(i in msites) {
  mmix_trophic <- run_simmr(my_data = d,
                            my_location = i,
                            my_group = 'transient_type') %>%
    full_join(mmix_trophic) }

mmix_transient <- mmix_trophic %>%
  mutate(dataset = 'transient_type')

#----------------------------------------------------------------------------
# CS mixing model: species
#----------------------------------------------------------------------------

d %>% filter(guild %in% c('Fish', 'Invertebrate')) %>%
  group_by(site_code, species) %>% 
  summarise(count=length(species)) %>%
  pivot_wider(names_from=site_code, values_from=count) %>%
  print(n=35)

common_species <- c('clarkii', 'megalotis', 'macrochirus', 'pugio', 'affinis')

# blank table
mmix_species <- tibble(
  deviance = numeric(),
  Aquatic = numeric(),
  Terrestrial = numeric(),
  m_group = character(),
  site = character(),
  statistic = character() )

# loop mixing models across sites and aggregate statistics and quantile info
for(i in msites) {
  mmix_species <- run_simmr(my_data = d,
                            my_location = i,
                            my_group = 'species') %>%
    full_join(mmix_species) }

mmix_species <- mmix_species %>%
  mutate(dataset = 'species')

#----------------------------------------------------------------------------
# Combine all mmix
#----------------------------------------------------------------------------
mmix_all <- full_join(mmix_order, mmix_transient) %>%
  full_join(mmix_species) %>%
  select(Estuarine, site, statistic, m_group, dataset) %>%
  mutate(Estuarine = 100*Estuarine) %>%
  pivot_wider(names_from = statistic, values_from = Estuarine)

write_csv(mmix_all, '03_public/output/CS_mix_out.csv')
#----------------------------------------------------------------------------
# End CS_mixing_model_calc