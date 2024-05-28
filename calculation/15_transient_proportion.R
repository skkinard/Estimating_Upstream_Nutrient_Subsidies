# 15_transient_proportion
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
dc <- read_csv(here::here('03_public', 'output', 
                               '15_combined_biomass_density.csv'))

#------------------------------------------------------------------------------
# transient proportion
#------------------------------------------------------------------------------
# total biomass
total_biomass <- dc %>%
  group_by(site_code, collection_period) %>%
  mutate(total_biomass = sum(biomass, na.rm=T)) %>%
  ungroup() %>%
  select(total_biomass, biomass, everything()) %>%
  arrange(desc(biomass))

# transient_type_biomass
transient_type_biomass <- total_biomass %>%
  group_by(site_code, collection_period, transient_type) %>%
  mutate(transient_type_biomass = sum(biomass, na.rm=T)) %>%
  ungroup() %>%
  select(transient_type_biomass, everything())
  
# transient_type_biomass_percent
tb_stats <- transient_type_biomass %>%
  mutate(transient_type_biomass_percent = 
           transient_type_biomass/total_biomass*100) %>%
  select(transient_type_biomass_percent, everything()) %>%
  ungroup()

#------------------------------------------------------------------------------
# visualize
#------------------------------------------------------------------------------
# site averages
p_biomass_v_rain_site <- tb_stats %>%
  group_by(site_code, transient_type) %>%
  summarize(mu_percent_euryhaline = 
              mean(transient_type_biomass_percent, na.rm=T)) %>%
  add_rain() %>%
  filter(transient_type %in% c('Euryhaline', 'Amphidromous')) %>%
  ggplot(aes(x=annualrain, y = mu_percent_euryhaline)) +
  facet_wrap(~transient_type, ncol=2) +
  geom_smooth(method='lm', se=F, lwd=.6, lty=2, show.legend=F, color='red3') +
  geom_point(size=4, color='black') +
  stat_poly_eq(label.x=.5, label.y=.95,
               color='black', use_label(c("adj.R2","p"))) +
  scale_fill_manual(values=my_colors) +
  theme_bw(base_size=14) +
  ylim(c(0,50)) +
  labs(x='Rainfall (cm/yr)',
       y='% Biomass')

# all sample events
p_biomass_v_rain_event <- tb_stats %>%
  group_by(site_code, collection_period, transient_type) %>%
  summarize(mu_percent_euryhaline = 
              mean(transient_type_biomass_percent, na.rm=T)) %>%
  add_rain() %>%
  filter(transient_type %in% c('Euryhaline', 'Amphidromous')) %>%
  ggplot(aes(x=annualrain, y = mu_percent_euryhaline)) +
  facet_wrap(~transient_type, ncol=2) +
  geom_smooth(method='lm', se=F, lwd=.6, lty=2, show.legend=F, color='red3') +
  geom_point(size=4, fill='black', alpha=.2) +
  stat_poly_eq(label.x=.5, label.y=.95,
               color='black', use_label(c("adj.R2","p"))) +
  scale_fill_manual(values=my_colors) +
  theme_bw(base_size=14) +
  ylim(c(0,50)) +
  labs(x='Rainfall (cm/yr)',
       y='% Biomass')

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(tb_stats, here::here('03_public', 'output', 
                               '15_transient_biomass_proportion.csv'))

#------------------------------------------------------------------------------
# End 15_transient_proportion