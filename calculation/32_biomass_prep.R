# 15_biomass_prep
# Sean Kinard
# 2024-05-21

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source(here::here('03_public', 'toolkit.R')) # load packages and helper-functions

library(bio.infer)

# load data
fb <- read_csv(here('03_local_files', 'data', 'community', 
                    'fish_biomass_density_alldates.csv'))
ib <- read_csv(here('03_local_files', 'data', 'community', 
                    'invert_biomass_density_alldates.csv'))
fspec <- read_csv(here('03_local_files', 'data', 'community', 
                       'fish_species_with_transient.csv'))

#------------------------------------------------------------------------------
# Data Prep
#------------------------------------------------------------------------------
# merge check (by = join_by('lowest_taxon'))
colnames(fb)
colnames(fspec)

## match format
fb <- fb %>%
  mutate(lowest_taxon = str_replace_all(lowest_taxon, '_', '. ')) %>%
  mutate(lowest_taxon = str_replace_all(lowest_taxon, 'P. crocro', 'R. crocro'))

## fish lowest taxon match?
fb_species <- fb %>% pull(lowest_taxon) %>% unique()
keep(fb_species, ! fb_species %in% fspec$lowest_taxon)

## invertebrate lowest taxon match?
ib_species <- ib %>% pull(family) %>% unique()
keep(ib_species, ! ib_species %in% fspec$family)

i_unclass <-  tibble(family = keep(ib_species, ! ib_species %in% fspec$family))

i_unclass <- i_unclass %>%
  mutate(id = row_number()) %>%
  rename(lowest_taxon = family) %>%
  select(id, lowest_taxon) %>%
  as.data.frame() %>%
  get.taxonomic() %>%
  mutate(order = str_to_title(ORDER), 
         family = str_to_title(FAMILY)) %>%
  select(order, family, lowest_taxon) %>%
  as_tibble()
  
# some unclassified invertebrates are euryhaline
i_euryhaline <- c('Ceratopogonidae', 'Tabanidae', 'Hydrochidae', 'Hydracarina', 
                  'Hydrometridae', 'Belostomatidae', 'Nematoda')

# add classifications
i_class <- i_unclass %>%
  mutate(transient_type = ifelse(family %in% i_euryhaline, 'Euryhaline', 'Freshwater'),
         transient = ifelse(transient_type != 'Freshwater', 'Transient', 'Resident'),
         is_diadromous = transient_type,
         guild = ifelse(family %in% c('Bufonidae', 'Ranidae'), 'Vertebrate', 
                        'Invertebrate'))

# write_csv(i_unclass, '03_public/output/unclassified_invertebrates.csv')


# Merge biomass and species classification data
fspec2 <- full_join(fspec, i_class) %>% 
  mutate(genus = ifelse(family == 'Chironomidae', NA, family) ) %>%
  unique()

# check invertebrate dates match fish sample event dates
f_dates_missing_in_i <- anti_join(
  fb %>%
    select(site_code, collection_period) %>%
    unique(),
  ib %>%
    select(site_code, collection_period) %>%
    unique()) %>%
  arrange(site_code, collection_period)
# mostly 2019 samples missing in invertebrate data

# group dates by annual quarter
fb <- fb %>%
  filter(collection_period > ymd('2017-06-01')) %>%
  mutate(year = year(collection_period),
         quarter = quarter(collection_period),
         collection_period = paste(year, quarter, sep='-Q')) %>%
  select(-year, -quarter) %>% 
  group_by(site_code, collection_period, lowest_taxon) %>%
  summarize(biomass_mu = mean(biomass, na.rm=T),
            biomass_sd = sd(biomass, na.rm = T),
            biomass_n = length(biomass)) %>%
  arrange(desc(biomass_mu))

ib <- ib %>%
  filter(collection_period > ymd('2017-06-01')) %>%
  mutate(year = year(collection_period),
         quarter = quarter(collection_period),
         collection_period = paste(year, quarter, sep='-Q')) %>%
  select(-year, -quarter) %>% 
  group_by(site_code, collection_period, family) %>%
  summarize(biomass_mu = mean(biomass, na.rm=T,
            biomass_sd = sd(biomass, na.rm = T),
            biomass_n = length(biomass)))

# combine fish and invertebrate biomass data
cb <- full_join(fb, ib%>%rename(lowest_taxon=family))

cb %>%
  filter(! is.na(biomass_mu)) %>%
  arrange(site_code, collection_period, desc(biomass_mu))

cb_taxons <- cb %>% pull(lowest_taxon) %>% unique()

missing <- keep(cb_taxons,! cb_taxons  %in% fspec2$lowest_taxon)

# sites where P. pugio were not distinguished from P. ohione
cb %>%
  filter(lowest_taxon %in% missing) %>%
  pull(site_code) %>%
  unique()

# possible taxonomic info
fspec2 %>%
  filter(str_detect(family, 'Palae'))

# sites where no macros were found (anecdotally remembered)
no_macros_found <- c('EM', 'GC', 'PD', 'SF', 'TR', 'WM')

# fill in resolution for Palaemonidae in cb
cb <- cb %>%
  mutate(lowest_taxon = case_when(
    lowest_taxon == 'Palaemonidae' & site_code %in% no_macros_found ~ 'P. pugio',
    lowest_taxon == 'Palaemonidae' ~ 'M. ohione',
    T ~ lowest_taxon
  ))

d <- left_join(cb, fspec2)

# correct biomass percent from guild-specific to entire community
d <- d %>%
  filter(!is.na(biomass_mu)) %>%
  ungroup() %>%
  group_by(site_code, collection_period) %>%
  mutate(biomass_etot = sum(biomass_mu, na.rm=T)) %>%
  ungroup() %>%
  mutate(biomass_percent = biomass_mu/biomass_etot*100) %>%
  select(contains('biomass'), everything()) %>%
  arrange(site_code, collection_period, lowest_taxon) %>%
  filter(! is.na(biomass_percent))

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(d, here('03_public', 'output', 
                  'fish_and_invert_biomass_with_taxonomic_and_transient_info.csv'))

#------------------------------------------------------------------------------
# End 15_biomass_prep
