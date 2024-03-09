# 16_clean_isotope_CNS
# Sean Kinard
# 2023-07-24

# summarize environmental variables
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('toolkit.R')

d <- read_csv('data/isotope_raw/TERRG_CNS_combined.csv') 

d_species <- read_csv('data/community/fish_species_with_transient.csv') %>%
  mutate(lowest_taxon=str_replace_all(lowest_taxon, " ", ""))


#------------------------------------------------------------------------------
# Formatting
#------------------------------------------------------------------------------
d <- d %>% r_friendly_colnames() %>%
  rename(site_code = site, 
         collection_date=collection, 
         guild=type,
         lowest_taxon=species) %>%
  select(-details, -tray, -uid, -project) %>%
  select(site_code, collection_date, guild, lowest_taxon, max_forklength, 
         my_rep, carbon, nitrogen, sulfur, c_percent, n_percent, s_percent) 

#------------------------------------------------------------------------------
# match lowest_taxon
#------------------------------------------------------------------------------
d <- d %>%
  mutate(lowest_taxon = case_when(
    lowest_taxon %in% c('cambarid', 'p.clarkii') ~ 'P.clarkii',
    lowest_taxon == 'csapidus' ~ 'C.sapidus',
    lowest_taxon == 'macrobrachium' ~ 'M.ohione',
    lowest_taxon %in% c('mysid', 'palaemonid') ~ 'P.pugio',
    lowest_taxon == 'thiarid' ~ 'melanoides',
    lowest_taxon == 'unknown' ~ 'coenagrionid',
    lowest_taxon == 'L.hybrid' ~ 'L.macrochirus',
    lowest_taxon == 'M.berylina' ~ 'M.beryllina',
    lowest_taxon == 'M.cephalus' ~ 'L.macrochirus',
    TRUE ~ lowest_taxon),
    lowest_taxon = str_to_title(lowest_taxon))

# merge to taxonomic key
d <- left_join(d%>%select(-guild), d_species) %>%
  select(lowest_taxon, transient, genus, everything())

#----------------------------------------------------------------------------
# remove suspicious data
#----------------------------------------------------------------------------
d <- d %>% filter(c_percent > 1)
d <- d %>% mutate(
  is_outlier = ifelse(sulfur < -2 & site_code == 'NB', 'yes', 'no')) %>%
  filter(is_outlier != 'yes') %>%
  select(-is_outlier)
#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(d %>% filter(! is.na(sulfur)), 
          'analysis/output/isotope_CNS_2020_01_clean.csv')

write_csv(d, 'analysis/output/isotope_CN_2020_01_05_09_clean.csv')

