# 15_biomass_add_EA
# Sean Kinard
# 2024-05-21
#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
# Background:
## This research project aims to quantify the contributions of estuarine-derived nutrients to coastal streams across a natural precipitation gradient in South-Central Texas. It investigates estuarine assimilation in both inconspicuous migrants and freshwater taxa while examining its relationships with climate, geography, and anthropogenic factors.

# Project Overview:
## In this script, I merge two biological datasets—estuarine assimilation (EA) and biomass—to calculate a species-biomass-weighted estimate for the overall community estuarine assimilation across nine streams. Subsequently, I examine the relationships between EA and annual rainfall using regression analysis at various scales of comparison.

# Key Findings:
## These results confirm the increased assimilation of estuarine-derived nutrients in regions with arid climates. Furthermore, they demonstrate that this relationship holds true for both freshwater and euryhaline species. Notably, the findings reveal a counterintuitive pattern: freshwater taxa exhibit greater consumption of estuarine materials in arid environments, while euryhaline species show reduced assimilation in humid environments. It's important to note that these results do not imply the mechanisms of consumption. For instance, freshwater fish in arid environments may either consume euryhaline wanderers or make periodic trips to nearby estuaries to directly consume estuarine materials before returning to their freshwater habitat.

# Relevance: 
# Our research sheds light on the intricate dynamics of estuarine assimilation in coastal stream ecosystems, spanning a natural precipitation gradient. By quantifying the contributions of estuarine-derived nutrients and their relationships with climatic factors, our findings challenge conventional assumptions. We reveal a nuanced pattern where freshwater taxa exhibit unexpected reliance on estuarine materials in arid environments, while euryhaline species display reduced assimilation in humid regions. These insights not only deepen our understanding of ecosystem dynamics but also have practical implications for managing and conserving coastal ecosystems, highlighting the importance of considering species-specific responses to environmental variability.

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
                    'fish_and_invert_biomass_with_taxonomic_and_transient_info.csv')) %>%
  select(-transient) %>%
  rename(transient = transient_type)

species_key <- read_csv(here('03_local_files', 'data', 'community',
                             'fish_species_with_transient.csv'))

# load Bayesian mixing model outputs from stable isotope data
# load Estuarine Assimilation Estimates for: order, transient_type,  species   
iso <- read_csv(here('03_public', 'output', 'CS_mix_out.csv')) %>%
  mutate(dataset=str_remove_all(dataset, '_type'))

# format isotope data by scales of comparison
extract_iso_group <- function(xgroup) {
  output <- iso %>% filter(dataset == xgroup) %>%
    rename(ZZZ = m_group, EA_ZZZ_mu = mean, EA_ZZZ_sd = sd, 
           site_code = site) %>%
    select(site_code, ZZZ, contains('EA')) %>%
    add_rain()
  
  colnames(output) <- str_replace_all(colnames(output), 'ZZZ', xgroup)
  
  return(output)
}

# Match isotope data column names to merge survey
iso_guild <- extract_iso_group('guild')

# Match isotope data column names to merge survey
iso_order <- extract_iso_group('order')

# Extract taxonomic information for merging
iso_species <- extract_iso_group('species') %>%
  left_join(species_key) %>%
  filter(!(site_code == 'PL' & species == 'maculatus' & order == 'Gobiiformes')) %>%
  select(site_code, species, EA_species_mu, EA_species_sd) %>%
  add_rain()

# extract EA for transient type: fresh, euryhaline, amphidromous, catadromous
iso_transient <- extract_iso_group('transient')

# merge: survey, EA by order, taxonomic info, EA by transient type
dc <- db %>%
  left_join(iso_guild) %>%
  left_join(iso_order) %>%
  left_join(iso_species) %>%
  left_join(iso_transient) %>%
  arrange(desc(biomass_mu))

dc %>% filter(is.na(EA_order_mu)) %>%
  select(site_code, lowest_taxon, order, contains('EA'))

#------------------------------------------------------------------------------
# Setup: EA Versus Rainfall (Taxonomic Groups)
#------------------------------------------------------------------------------

# table function: EA versus Rainfall
table_lm_stats <- function(x) {
  
  lm(formula = EA_XX_mu ~ annualrain, 
     data = x %>% add_rain()) %>%
    summary()  %>%
    broom::tidy()
}

# base plot: EA versus Rainfall
plot_ea_v_rain <- function(x) {

  x %>%
    ggplot(aes(x=annualrain, y=EA_XX_mu)) +
    stat_poly_eq(label.x=.5, label.y=.95, formula=y~x,
                 color='black', use_label(c("adj.R2","p")), size=6) +
    geom_point(size=4, color='blue', fill='skyblue', shape=21, alpha=.5) +
    geom_point(size=4, color='blue', fill=NA, shape=21) +
    labs(x='Rainfall (cm/yr)', y='% Estuarine') +
    theme_bw(base_size=20) +
    geom_smooth(data = . %>% filter(p.value < 0.1 & p.value >= 0.05), 
                method = "lm", se = FALSE, 
                color = "blue", lwd=.5, lty=2)  +
    geom_smooth(data = . %>% filter(p.value < 0.05), 
                method = "lm", se = FALSE, 
                color = "blue", lwd=.5, lty=1)
}

# Generate table and plot for x_group of comparison
ea_vs_rain <-function(x_data, x_group, n_sites=1) {
  temp_data <- x_data
  colnames(temp_data) <- str_replace_all(colnames(temp_data), x_group, 'XX')
  
  # list widespread groups
  widespread <- temp_data %>%
    group_by(XX) %>%
    summarize(n_samples = length(EA_XX_mu)) %>%
    filter(n_samples>n_sites) %>%
    pull(XX)
  
  # Table: Regression Statistics
  t_EA_vs_rain <- temp_data %>%
    filter(XX %in% widespread) %>%
    group_by(XX) %>%
    nest() %>%
    mutate(lm = map(data, table_lm_stats)) %>%
    unnest(lm) %>%
    filter(term == 'annualrain') %>%
    select(-term, -data) 
  
  colnames(t_EA_vs_rain) <- str_replace_all(colnames(t_EA_vs_rain), x_group, 'XX')
  
  # visualize 
  p_EA_vs_rain <- temp_data %>%
    filter(XX %in% widespread) %>%
    left_join(t_EA_vs_rain) %>%
    plot_ea_v_rain() +
    facet_wrap(~XX)
  
  colnames(t_EA_vs_rain) <- str_replace_all(colnames(t_EA_vs_rain), 'XX', x_group)
  
  output <- list(figure = p_EA_vs_rain,
       table = t_EA_vs_rain)
  
  return(output)
}

#------------------------------------------------------------------------------
# EA Versus Rainfall: Within Widespread Orders
#------------------------------------------------------------------------------

widespread_orders <- ea_vs_rain(x_data = iso_order, 
                                x_group= 'order', 
                                n_sites = 8)

widespread_orders$table
widespread_orders$figure

# Notes
## Statistically significant: EA decreases with Rainfall for Centrarchiformes.
## NS: Cyprinidontiformes, Decapoda

# Interpretation
## The analysis indicates a significant negative correlation between estuarine assimilation (EA) and annual rainfall within the taxonomic order Centrarchiformes. Additionally, while Cyprinidontiformes exhibit a similar relationship, it falls just short of statistical significance (p = 0.106). Notably, Decapoda also demonstrates a pattern resembling other orders, but data from the most arid site were missing, potentially affecting the analysis. These findings suggest that, within widespread taxonomic orders, centrarchiformes display a clear linear decline in EA with increasing annual rainfall, indicating species-specific responses to environmental variables. However, further investigation is warranted to account for missing data and confirm the observed trends across all taxa.


#------------------------------------------------------------------------------
# EA versus rainfall within widespread species
#------------------------------------------------------------------------------

widespread_species <- ea_vs_rain(x_data = iso_species,
                                 x_group= 'species',
                                 n_sites = 6)

widespread_species$table
widespread_species$figure

# Notes
## Statistically significant: EA decreases with Rainfall for Palaemonetes pugio and Lepomis macrochirus.
## NS: Procambarus clarkii

# Interpretation
## The analysis reveals a significant negative correlation between estuarine assimilation (EA) and annual rainfall for two widespread species: Lepomis macrochirus (bluegill sunfish) with a p-value below 0.1 threshold and Palaemonetes pugio (daggerblade grass shrimp) with a p-value below 0.05 threshold. In contrast, EA in other species like Procambarus clarkii did not exhibit a linear correlation with annual rainfall. These findings suggest that estuarine exploitation varies among species, with Pugio and macrochirus showing a decreased reliance on estuarine-derived materials in regions with higher rainfall. This indicates that estuarine foraging behavior is species-specific and not universally inherent within the community.


#------------------------------------------------------------------------------
# EA versus rainfall within Transient Type
#------------------------------------------------------------------------------

widespread_transient <- ea_vs_rain(x_data = iso_transient,
                                   x_group= 'transient', 
                                   n_sites = 5)

widespread_transient$table
widespread_transient$figure

# Statistically significant: EA decreases with Rainfall for Amphidromous and Freshwater organisms.

# Interpretation:
## The analysis reveals a significant negative correlation between estuarine assimilation (EA) and rainfall for both Amphidromous and Freshwater organisms. This suggests that these organisms tend to consume more estuarine-derived materials in regions with lower rainfall, possibly due to increased concentration of nutrients and salinity during dry periods. In contrast, transient and catadromous organisms show no linear relationship with annual rainfall, indicating a consistent consumption of estuarine-derived materials irrespective of rainfall patterns.

  
#------------------------------------------------------------------------------
# resident EA versus % biomass transient
#------------------------------------------------------------------------------

# table function: EA versus Rainfall
table_lm_stats2 <- function(x) {
  
  lm(formula = EA_fresh ~ biomass_percent, data = x) %>%
    summary()  %>%
    broom::tidy()
}

plot_ea_v_tb <-function(X) {
  X %>%
    ggplot(aes(x=biomass_percent, y=EA_fresh)) +
    stat_poly_eq(label.x=.5, label.y=.95, formula=y~x,
                 color='black', use_label(c("adj.R2","p")), size=6) +
    geom_point(size=4, color='blue', fill='skyblue', shape=21, alpha=.5) +
    geom_point(size=4, color='blue', fill=NA, shape=21) +
    labs(x='Transient Biomass (%)', y="Resident EA (%)") +
    theme_bw(base_size=20) +
    geom_smooth(data = . %>% filter(p.value < 0.1 & p.value >= 0.05), 
                method = "lm", se = FALSE, 
                color = "blue", lwd=.5, lty=2)  +
    geom_smooth(data = . %>% filter(p.value < 0.05), 
                method = "lm", se = FALSE, 
                color = "blue", lwd=.5, lty=1) }

# data prep:
d_ea_v_tb <- dc %>%
  filter(collection_period %in% c('2019-Q4', '2020-Q1')) %>%
  group_by(site_code, collection_period, transient) %>%
  summarize(biomass_percent = sum(biomass_percent, na.rm=T)) %>%
  ungroup() %>%
  mutate(present = ifelse(near(biomass_percent, 0), 0, 1)) %>%
  group_by(transient, collection_period) %>%
  mutate(n_sites = sum(present)) %>%
  ungroup() %>%
  select(-present) %>%
  add_rain() %>%
  left_join(iso_transient %>% 
              filter(transient == 'Freshwater') %>%
              rename(EA_fresh = EA_transient_mu) %>%
              select(site_code, EA_fresh) )

# Table: Regression Statistics
t_EA_vs_tb <- d_ea_v_tb %>%
  group_by(transient) %>%
  nest() %>%
  mutate(lm = map(data, table_lm_stats2)) %>%
  unnest(lm) %>%
  filter(term == 'biomass_percent') %>%
  select(-term, -data) 

# figure: Resident EA Versus % Transient Biomass
p_ea_v_tb <- d_ea_v_tb %>%
  left_join(t_EA_vs_tb) %>%
  plot_ea_v_tb() +
  facet_wrap(~transient)

# Notes
## Statistically Significant: Euryhaline
## NS: Amphidromous, Catadromous, Freshwater

# Interpretation:
## The analysis reveals a statistically significant positive relationship between resident estuarine assimilation and the proportion of euryhaline biomass within the community, comprising both fish and invertebrates. 
## Two potential causal mechanisms can be inferred: 
### (1) Euryhaline species might contribute to the deposition of estuarine-derived nutrients into freshwater ecosystems.
### (2) Alternatively, the presence of euryhaline species in a system may correlate with an increased likelihood of freshwater species traveling to nearby estuaries to exploit estuarine-derived nutrients.

#------------------------------------------------------------------------------
# export
#------------------------------------------------------------------------------
write_csv(  dc,   here('03_public', 'output', 'biomass_add_EA.csv'))
#------------------------------------------------------------------------------
# End 15_biomass_add_EA