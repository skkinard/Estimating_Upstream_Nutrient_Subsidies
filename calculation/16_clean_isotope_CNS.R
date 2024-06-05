# 16_clean_isotope_CNS
# Sean Kinard
# 2023-07-24

# summarize environmental variables
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('03_public/toolkit.R')

d <- read_csv('03_local_files/data/isotope_raw/TERRG_CNS_combined.csv') 

d_species <- read_csv('03_local_files/data/community/fish_species_with_transient.csv') %>%
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

bad <- d %>%
  filter(is.na(guild)) 

bad <- bad %>%
  mutate(lowest_taxon = case_when(
    lowest_taxon %in% c("Melanoides", "C.sapidus") ~ lowest_taxon,
    ! str_detect(lowest_taxon, 'id') ~ lowest_taxon,
    str_detect(lowest_taxon, 'id') ~ paste(lowest_taxon, 'ae', sep=''),
    T ~ NA)) %>%
  mutate(id = row_number())

#----------------------------------------------------------------------------
# Fix missing taxonomic information for mispelled invertebrates
#----------------------------------------------------------------------------
library(bio.infer)

add_tax <- bad %>%
  select(id, lowest_taxon) %>%
  as.data.frame() %>%
  get.taxonomic() %>%
  mutate(order = str_to_title(ORDER), 
         family = str_to_title(FAMILY)) %>%
  select(order, family, lowest_taxon) %>%
  as_tibble()

fixed <- add_tax %>% 
  mutate(
    family = case_when(
      lowest_taxon == 'Belastomatidae' ~ 'Belastomatidae',
      lowest_taxon == 'Cyrenidae' ~ 'Corbiculidae',
      lowest_taxon == 'Dysticidae' ~ 'Dytiscidae',
      lowest_taxon == 'Hirvdinidae' ~ 'Hirudinidae',
      lowest_taxon == 'Rharrisii' ~ 'Panopeidae',
      lowest_taxon == 'Simaliidae=' ~ 'Simuliidae',
      lowest_taxon == 'Tanyponidae' ~ 'Chironomidae',
      lowest_taxon == 'Velliidae' ~ 'Veliidae',
      lowest_taxon == 'C.sapidus' ~ 'Portunidae', 
      T ~ family),
    order = case_when(
      lowest_taxon == 'Belastomatidae' ~ 'Hemiptera',
      lowest_taxon == 'Cyrenidae' ~ 'Veneroida',
      lowest_taxon == 'Dysticidae' ~ 'Coleoptera',
      lowest_taxon == 'Hirvdinidae' ~ 'Hirudinea',
      lowest_taxon == 'Rharrisii' ~ 'Decapoda',
      lowest_taxon == 'Simaliidae=' ~ 'Diptera',
      lowest_taxon == 'Tanyponidae' ~ 'Diptera',
      lowest_taxon == 'Velliidae' ~ 'Hemiptera',
      lowest_taxon == 'C.sapidus' ~ 'Decapoda', 
      T ~ order ),
    transient_type = ifelse(lowest_taxon == 'Rharrisii', 'Euryhaline', 
                            'Freshwater'),
    guild = 'Invertebrate') %>% unique()

bad_fixed <- bad %>%
  select(-c(family, order, guild, transient_type)) %>%
  left_join(fixed)

d <-d %>% 
  filter(!is.na(guild)) %>%
  full_join(bad_fixed)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(d %>% filter(! is.na(sulfur)), 
          '03_public/output/isotope_CNS_2020_01_clean.csv')

write_csv(d, '03_public/output/isotope_CN_2020_01_05_09_clean.csv')
