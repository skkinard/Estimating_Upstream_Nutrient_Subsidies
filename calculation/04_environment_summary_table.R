# 04_site_table
# Sean Kinard
# 2023-07-17

# summarize environmental variables
#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
source('toolkit.R')

dl <- read_csv("data/environment/environment_longterm.csv")

ddist <- read_csv('data/environment/site_distances.csv')

dw <- read_csv("data/environment/MDN_site_lat_lon.csv") %>%
  r_friendly_colnames() %>%
  rename(lat=site_lat, lon=site_lon)

ds <- read_csv("data/environment/environment_shortterm.csv") %>%
  filter(collection_period > ymd('2019-12-31') &
           collection_period < ymd('2020-12-31')) %>%
  select(-c('nitrate', 'phosphate', 'ammonia', 'doc'))

dn <- read_csv('data/environment/nutrient_TERRG.csv') %>%
  rename(nitrate=no3n, phosphate=ortho_p, ammonia=nh4_n) %>%
  create_combo_period() %>%
  select(-site_period, -collection_date)

#------------------------------------------------------------------------------
# Long-term Environment
#------------------------------------------------------------------------------
lte_watershed <- dl %>%
  select(site_code, any_of(location), any_of(climate), any_of(landuse)) %>%
  right_join(dw %>% select(site_code, staid) ) %>%
  filter(site_code %in% my_sites) %>%
  pretty_titles()
  
lte_flow <- dl %>% 
  filter(site_code %in% my_sites) %>%
  select(site_code, contains('q_')) %>%
  pretty_titles()
colnames(lte_flow) <- c('Site Code', 
                        colnames(select(lte_flow, contains('Q'))) %>%
  str_to_upper())

lte_water_quality <- dl %>% 
  filter(site_code %in% my_sites) %>%
  select(site_code, any_of(water_quality)) %>%
  pretty_titles()

lte_geomorph <- dl %>% 
  filter(site_code %in% my_sites) %>%
  select(site_code, any_of(geomorph)) %>% pretty_titles()

lte_algae <- dl %>% 
  filter(site_code %in% my_sites) %>%
  select(site_code, any_of(algae)) %>% pretty_titles()

lte_distances <- ddist  %>% pretty_titles()

lte_table_long <- left_join(lte_watershed, lte_flow) %>% 
  left_join(lte_water_quality) %>%
  left_join(lte_geomorph) %>%
  left_join(lte_algae) %>%
  left_join(lte_distances) %>%
  filter(`Site Code` %in% my_sites) %>%
  select("USGS STAID", everything()) %>%
  column_to_rownames(var = "Site Code") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var='Variable') %>%
  as_tibble()

caption_lte_summary <- 'Table of long-term environmental variables for sampling locations. Values represent 20 year averages. Flow is average annual discharge, HFPP is the proportion of the annual discharge that is 3x higher than the Flow, LFPP is the proportion of discharge below the 25th percentile, Flashiness is the cumulative changes in day to day discharge divided by cumulative annual discharge, Season approximates the degree to which the flow varies during the course of a single year.'

# site location
dw <- dw %>% 
  arrange(desc(site_type), lon) %>%
  select(site_type, everything()) %>%
  pretty_titles() %>%
  rename(Latitude=Lat, Longitude=Lon)

#------------------------------------------------------------------------------
# short term data
#------------------------------------------------------------------------------

# merge nutrient to habitat
d_ste <- full_join(ds, dn) %>% 
  pivot_longer(cols=-c('site_code', 'collection_period'), 
               names_to='xname', values_to='xvalue') %>%
  group_by(site_code, collection_period, xname) %>%
  dplyr::summarize(xvalue=mean(xvalue, na.rm=T)) %>%
  pivot_wider(names_from=xname, values_from=xvalue)

# impute missing with interpolation within site
d_ste <- d_ste %>% impute_interpolation()

# aggregate data near stable isotope collections (iso_periods)
d_ste <- d_ste %>%
  create_iso_period() %>%
  select(-collection_period) %>%
  pivot_longer(cols=-c('site_code', 'iso_period'), 
               names_to='xname', values_to='xvalue') %>%
  group_by(site_code, iso_period, xname) %>%
  dplyr::summarize(xvalue=mean(xvalue, na.rm=T)) %>%
  pivot_wider(names_from=xname, values_from=xvalue)


ste_qtr_flow <- d_ste %>% select(site_code, contains('q2wk')) %>%
  pretty_titles()
colnames(lte_flow) <- c('Site Code', 
                        colnames(select(lte_flow, contains('Q'))) %>%
                          str_to_upper())

ste_qtr_water_quality <- d_ste %>% select(site_code, any_of(water_quality)) %>%
  pretty_titles()

ste_qtr_geomorph <- d_ste %>% select(site_code, any_of(geomorph)) %>% pretty_titles()

ste_qtr_algae <- d_ste %>% select(site_code, any_of(algae)) %>% pretty_titles()

#------------------------------------------------------------------------------
# Export lte Tables
#------------------------------------------------------------------------------
lte_summary_tables <- list(lte_table_long,
                           lte_watershed,
                           lte_flow,
                           lte_water_quality,
                           lte_geomorph,
                           lte_algae,
                           lte_distances)

names(lte_summary_tables) <- c('lte_table_long',
                               'lte_watershed',
                               'lte_flow',
                               'lte_water_quality',
                               'lte_geomorph',
                               'lte_algae',
                               'lte_distances')

for (i in 1:length(lte_summary_tables)) {
  my_place <- paste('analysis/output/summary_table_',
                    names(lte_summary_tables[i]), 
                    ".csv", sep='')
  my_object <- lte_summary_tables[[i]]
  write_csv(my_object, my_place) }

write_csv(dw, 'analysis/output/summary_table_site_location.csv')

#------------------------------------------------------------------------------
# Export ste Tables
#------------------------------------------------------------------------------
ste_qtr_summary_tables <- list(d_ste,
                           ste_qtr_flow,
                           ste_qtr_water_quality,
                           ste_qtr_geomorph,
                           ste_qtr_algae)

names(ste_qtr_summary_tables) <- c('ste_qtr_long',
                               'ste_qtr_flow',
                               'ste_qtr_water_quality',
                               'ste_qtr_geomorph',
                               'ste_qtr_algae')

for (i in 1:length(ste_qtr_summary_tables)) {
  my_place <- paste('analysis/output/summary_table_', 
                    names(ste_qtr_summary_tables[i]), 
                    ".csv", sep='')
  my_object <- ste_qtr_summary_tables[[i]]
  write_csv(my_object, my_place) }

#------------------------------------------------------------------------------
# End 04_site_table