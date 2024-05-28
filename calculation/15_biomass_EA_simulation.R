# 15_biomass_correction
# Sean Kinard
# 2024-05-21

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source(here::here('03_public', 'toolkit.R')) # load packages and helper-functions
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(here)

# load data
dc <- read_csv(here('03_public', 'output', 'biomass_add_EA.csv'))

# prepare for simulations
d_empty <- dc %>%
  mutate(mu = case_when(
    is.na(EA_species_mu) & is.na(EA_order_mu) ~ EA_transient_mu,
    is.na(EA_species_mu) ~ EA_order_mu,
    T ~ EA_species_mu ),
    sd = case_when(
      is.na(EA_species_sd) & is.na(EA_order_sd) ~ EA_transient_sd,
      is.na(EA_species_sd) ~ EA_order_sd,
      T ~ EA_species_sd )) %>%
  select(-starts_with('EA')) %>%
  rename(EA_mu = mu, EA_sd = sd) %>%
  select(site_code, collection_period, lowest_taxon, EA_mu, EA_sd, biomass_mu,
         biomass_percent) %>%
  group_by(site_code, collection_period, lowest_taxon) %>%
  nest()

#----------------------------------------------------------------------------
# Simulate Communities
#----------------------------------------------------------------------------
sim_ea <- function(xdata) {
  # number of pulls = the integer of relative biomass (as proportion) multiplied by 200. In other words n=grams of lowest taxon per 200 grams of biomass for that sample event. In this way there arre 200 total draws for each event, and the draws per species are based on their biomass relative to the biomass of fish and invertebrates for that event.
  r_pulls <- xdata$biomass_percent*2%>%ceiling()
  # take random samples from a normal distribution
  # mean and sd are taken from simmr aa estimates for species > guild
  EA_sim <- rnorm(n = r_pulls, mean = xdata$EA_mu, sd = xdata$EA_sd)
  output <- tibble(EA = EA_sim)
  return(output) }

d_fill <- d_empty %>%
  mutate(sim_tib = map(data, sim_ea)) %>%
  unnest(sim_tib) %>%
  unnest(data)

# Add taxonomic and transient categories
d_fill <- left_join(
  d_fill, 
  dc %>% 
    select(lowest_taxon, order, family, genus, species, guild, 
           transient, transient_type, is_diadromous)%>% 
    unique())

#----------------------------------------------------------------------------
# Visualize Results
#----------------------------------------------------------------------------
# Estuarine Assimilation versus rainfall
p_EA_v_Rainfall <- d_fill %>%
  filter(collection_period %in% c('2020-Q1', '2019-Q4')) %>%
  group_by(site_code, transient_type) %>%
  summarize(EA_mu_corrected = mean(EA, na.rm=T)) %>%
  add_rain() %>%
  ggplot(aes(x = annualrain, y = EA_mu_corrected, fill = transient_type, 
             color = transient_type, shape=transient_type)) +
  geom_point(size = 5, alpha=.4) +
  geom_point(size = 5, fill=NA) +
  scale_shape_manual(values = c(21,22,23,24)) +
  scale_fill_manual(values = c('purple', 'blue', 'green', 'orange')) +
  scale_color_manual(values = c('purple4', 'blue4', 'green4', 'orange4')) +
  theme_bw() +
  labs(y = 'Mean Estuarine Assimilation',
       x = 'Rainfall (cm/yr)')

d_fill

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------

write_csv(dc, here::here('03_public', 'output', 
                         '15_combined_biomass_density.csv'))

#------------------------------------------------------------------------------
# End 15_biomass_correction