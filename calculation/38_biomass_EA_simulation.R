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

# TR missing EA data for invertebrates, calculate site_EA_mu from fish for imputation
dc <- dc %>%
  group_by(site_code) %>%
  mutate(EA_site_mu = mean(EA_guild_mu, na.rm=T),
         EA_site_sd = mean(EA_guild_sd, na.rm=T)) %>%
  group_by(site_code) %>% # regional average sd for imputing missing values (n=1)
  mutate(site_sd = mean(EA_guild_sd, na.rm=T)) %>%
  ungroup()

# prepare for simulations
d_empty <- dc %>%
  mutate(mu = case_when(
    is.na(EA_guild_mu) ~ EA_site_mu,
    is.na(EA_transient_mu) ~ EA_guild_mu,
    is.na(EA_species_mu) ~ EA_transient_mu,
    T ~ EA_species_mu),
    sd = case_when(
      is.na(EA_guild_sd) ~ EA_site_sd,
      is.na(EA_transient_sd) ~ EA_guild_sd,
      is.na(EA_species_sd) ~ EA_transient_sd,
      T ~ EA_species_sd)) %>%
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
  # number of pulls = the integer of relative biomass (as proportion) multiplied by 200. In other words n=grams of lowest taxon per 200 grams of biomass for that sample event. In this way there are 200 total draws for each event, and the draws per species are based on their biomass relative to the biomass of fish and invertebrates for that event.
  r_pulls <- xdata$biomass_percent*200%>%ceiling()
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
# Calculate Biomass-Adjusted Community Mean EA
#----------------------------------------------------------------------------
d_site_EA_mu <- d_fill %>%
  group_by(site_code) %>%
  summarize(site_EA_mu = mean(EA, na.rm=T),
            site_EA_sd = sd(EA, na.rm=T))

d_transient_EA_mu <- d_fill %>%
  group_by(site_code, transient_type) %>%
  summarize(transient_EA_mu = mean(EA, na.rm=T),
            trasient_EA_sd = sd(EA, na.rm=T))

d_transient_EA_mu %>%
  add_rain() %>%
  add_sitevars() %>%
  ggplot(aes(x=baydist_km, y=transient_EA_mu)) +
  facet_wrap(~transient_type) +
  stat_poly_eq(label.x=.5, label.y=.95, formula=y~x,
               color='black', use_label(c("adj.R2","p")), size=4) +
  geom_point(size=3, color='blue', fill='skyblue', shape=21, alpha=.5) +
  geom_point(size=3, color='blue', fill=NA, shape=21) +
  theme_bw(base_size=12) +
  geom_smooth(method = "lm", se = FALSE, 
              color = "blue", lwd=.5, lty=2)
# relationships between EA and annual rain are similar to those observed with the unaltered stable isotope data. However, when using simulated biomass-weighted communities, the statistical significance slightly diminishes. This discrepancy may stem from inadequate coverage in stable isotope sampling during this phase, leading to a bias toward the guild or site community mean EA. Consequently, this dilutes the distinctions observed between categories of transient type.

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



#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(d_fill, here::here('03_public', 'output', 
                             '15_simulated_biomass_weighted_communities_with_EA.csv'))


write_csv(dc, here::here('03_public', 'output', 
                         '15_combined_biomass_density.csv'))

#------------------------------------------------------------------------------
# End 15_biomass_correction