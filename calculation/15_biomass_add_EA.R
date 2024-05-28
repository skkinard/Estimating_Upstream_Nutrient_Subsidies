# 15_biomass_add_EA
# Sean Kinard
# 2024-05-21
#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
# Background: This research project quantifies the contributions of estuarine-derived nutrients to coastal streams spanning a natural precipitation gradient in South-Central Texas. It reveals estuarine assimilation in inconspicuous migrants as well as freshwater taxa and tests relationships to climate, geographic, and anthropogenic features.

# In this script, I combine 2 biological data sets (estuarine assimilation (EA) and biomass) to calcuate a species-biomass-weighted estimate for the overall community estuarine assimilation for 9 streams. Then, I test relationships between EA and annual rainfall using regression at various scales of comparison.

# These results corroborate the increased assimilation of estuarine-derived nutrients in more arid climate. Second, they show that this relationship is conserved within both freshwater and euryhaline species. Thus, we intriguingly reveal that freshwater taxa consume more estuarine materials in arid environment, and that euryhaline species consume less estuarine materials in humid environments. These results do not imply the mechanisms of consumption. For example, freshwater fish in an arid environment may consume euryhaline wanderers or they may sojourn to the nearby estuary to consume estuarine materials directly before returning to the freswhater environment where they were caught.

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------

# load core project packages and custom-made helper-functions
source(here::here('03_public', 'toolkit.R'))

# load some script-specific packages
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(here)

# load survey data containing biomass and density for taxa for each sample event
db <- read_csv(here('03_public', 'output', 
                    'fish_and_invert_biomass_with_taxonomic_and_transient_info.csv'))

species_key <- read_csv(here('03_local_files', 'data', 'community',
                             'fish_species_with_transient.csv'))

# load Bayesian mixing model outputs from stable isotope data
# load Estuarine Assimilation Estimates for: order, transient_type,  species   
iso <- read_csv(here('03_public', 'output', 'CS_mix_out.csv'))

# Match isotope data column names to merge survey
iso_order <- iso %>% filter(dataset == 'order') %>%
  rename(order = m_group, EA_order_mu = mean, EA_order_sd = sd, 
         site_code = site) %>%
  select(site_code, order, EA_order_mu, EA_order_sd)

# Extract taxonomic information for merging
iso_species <- iso %>% filter(dataset == 'species') %>%
  rename(species = m_group, EA_species_mu = mean, EA_species_sd= sd, 
         site_code = site) %>%
  select(site_code, species, EA_species_mu, EA_species_sd) %>%
  left_join(species_key) %>%
  filter(!(site_code == 'PL' & species == 'maculatus' & order == 'Gobiiformes'))

# extract EA for transient type: fresh, euryhaline, amphidromous, catadromous
iso_transient <- iso %>% filter(dataset == 'transient_type') %>%
  rename(transient_type = m_group, EA_transient_mu = mean, EA_transient_sd=sd, 
         site_code = site) %>%
  select(site_code, transient_type, EA_transient_mu, EA_transient_sd) 

# merge: survey, EA by order, taxonomic info, EA by transient type
dc <- db %>%
  left_join(iso_order) %>%
  left_join(iso_species) %>%
  left_join(iso_transient) %>%
  arrange(desc(biomass_mu))

#------------------------------------------------------------------------------
# Visualize: widespread taxa
#------------------------------------------------------------------------------

# data prep: find widespread taxonomic orders for regression
order_widespread <- iso_order %>%
  group_by(order) %>%
  summarize(n_samples = length(EA_order_mu)) %>%
  arrange(desc(n_samples)) %>%
  filter(n_samples>8) %>%
  pull(order)

# data prep: find widespread invertebrate species for regression
species_widespread <- iso_species %>%
  group_by(species) %>%
  summarize(n_samples = length(EA_species_mu)) %>%
  arrange(desc(n_samples)) %>%
  filter(n_samples>6) %>%
  pull(species)

# base plot function
vis_widespread <- function(x) {
  x %>%
    ggplot(aes(x=annualrain, y=EA_order_mu)) +
    stat_poly_eq(label.x=.5, label.y=.95, formula=y~x,
                 color='black', use_label(c("adj.R2","p")), size=6) +
    geom_point(size=4, fill='red', shape=21) +
    labs(x='Rainfall (cm/yr)', y='% Estuarine') +
    theme_bw(base_size=20)
}

# visualize EA versus rainfall within widespread orders
p_EA_vs_rain_by_order <- iso_order %>%
  filter(order %in% order_widespread) %>%
  add_rain() %>%
  vis_widespread +
  facet_wrap(~order)
# statistically significant: Centrarchiformes EA decreases with Rainfall
# NS: Cyprinidontiformes, Decapoda
# Interpretation: Within widespread taxanomic orders, centrarchiformes display a linear decline in EA with annual rainfall. Cyprinidontiformes exhibit a similar relationship at a higher significance threshold (p = 0.106), and Decapoda also visually resembles other orders but was crucially missing data from the most arid site.

# visualize EA versus rainfall within widespread species
p_EA_vs_rain_by_species <- iso_species %>%
  filter(species %in% species_widespread) %>%
  add_rain() %>%
  vis_widespread() +
  facet_wrap(~species)
# Statistically significant: Pugio and macrochirus* EA decreases with Rainfall
# NS clarkii

#------------------------------------------------------------------------------
# Visualize: Transient Type
#------------------------------------------------------------------------------
# EA versus rainfall within Transient Type
p_EA_vs_rain_by_transient <- iso_transient %>%
  filter(transient_type %in% c('Amphidromous', 'Freshwater')) %>%
  add_rain() %>%
  ggplot(aes(x=annualrain, y=EA_transient_mu)) +
  facet_wrap(~transient_type) +
  stat_poly_eq(label.x=.5, label.y=.95, formula=y~x,
               color='black', use_label(c("adj.R2","p")), size=6) +
  geom_point(size=4, fill='red', shape=21) +
  geom_smooth(lty=1, lwd=.4, se=F, method='lm') +
  labs(x='Rainfall (cm/yr)', y='% Estuarine') +
  theme_bw(base_size=20)
# Statistically significant: Amphidromous* and Freshwater* EA decreases with Rainfall
# Interpretation: Amphidromous and Freshwater organisms consume more estuarine-derived materials in more arid regions. Transient and catadromous diets do not relate linearly to annual rainfall, which implies stable consumption of estuarine-derived materials regardless of rainfall regime.

# resident versus % biomass transient
dc_vis_tran <- dc %>%
  filter(collection_period %in% c('2019-Q4', '2020-Q1')) %>%
  group_by(site_code, collection_period, transient_type) %>%
  summarize(biomass_percent = sum(biomass_percent, na.rm=T)) %>%
  ungroup() %>%
  left_join(iso_transient %>% 
              filter(transient_type == 'Freshwater') %>%
              rename(EA_fresh = EA_transient_mu) %>%
              select(site_code, EA_fresh) )

vis_tran <-function(X) {
p_ea_fresh_vs_pcnt_biomass_all <-  X %>%
  ggplot(aes(x=biomass_percent, y=EA_fresh)) +
  stat_poly_eq(label.x=.5, label.y=.95, formula=y~x,
               color='black', use_label(c("adj.R2","p")), size=6) +
  geom_point(size=4, fill='red', shape=21) +
  geom_smooth(lty=1, lwd=.4, se=F, method='lm') +
  labs(x='% Community Biomass', y='% EA in Freshwater Taxa') +
  theme_bw(base_size=20) }

# all taxa groups
p_ea_fresh_vs_pcnt_biomass_alldc_vis_tran %>%
  vis_tran() +
facet_wrap(~transient_type)
# Statistically Significant: Euryhaline
# NS: Amphidromous, Catadromous, Freshwater

# only euryhaline
p_ea_fresh_vs_pcnt_biomass_eury <- dc_vis_tran %>%
  filter(transient_type == 'Euryhaline') %>%
  vis_tran()
# Statistically Significant: Euryhaline

#------------------------------------------------------------------------------
# export
#------------------------------------------------------------------------------
write_csv(  dc,   here('03_public', 'output', 'biomass_add_EA.csv'))
#------------------------------------------------------------------------------
# End 15_biomass_add_EA
