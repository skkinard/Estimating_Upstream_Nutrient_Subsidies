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
  mutate(total_biomass = sum(biomass_mu, na.rm=T)) %>%
  ungroup() %>%
  select(total_biomass, biomass_mu, everything()) %>%
  arrange(desc(biomass_mu))

# transient_type_biomass
transient_type_biomass <- total_biomass %>%
  group_by(site_code, collection_period, transient_type) %>%
  mutate(transient_type_biomass = sum(biomass_mu, na.rm=T)) %>%
  ungroup() %>%
  select(transient_type_biomass, everything())
  
# transient_type_biomass_percent
tb_stats <- transient_type_biomass %>%
  mutate(transient_type_biomass_percent = 
           transient_type_biomass/total_biomass*100) %>%
  select(transient_type_biomass_percent, everything()) %>%
  ungroup()

#------------------------------------------------------------------------------
# Community Biomass Proportion Vs Rainfall: Setup Functions
#------------------------------------------------------------------------------

# table function: Biomass versus Rainfall
table_lm_stats <- function(x) {
  
  lm(formula = mu_percent_biomass ~ annualrain, 
     data = x %>% add_rain()) %>%
    summary()  %>%
    broom::tidy()
}

# Figure Function: Biomass versus Rainfall
visualize_biomass_vs_rain <- function(x) {
  x %>%
    ungroup() %>%
    ggplot(aes(x=annualrain, y = mu_percent_biomass)) +
    stat_poly_eq(label.x=.5, label.y=.95, formula=y~x,
                 color='black', use_label(c("adj.R2","p")), size=4) +
    geom_point(size=4, color='blue', fill='skyblue', shape=21, alpha=.5) +
    geom_point(size=4, color='blue', fill=NA, shape=21) +
    labs(x='Rainfall (cm/yr)', y='% of Community Biomass') +
    theme_bw(base_size=14) +
    geom_smooth(data = . %>% filter(p.value < 0.1), 
                method = "lm", se = FALSE, 
                color = "blue", lwd=.5, lty=2) +
    geom_smooth(data = . %>% filter(p.value < 0.05), 
                method = "lm", se = FALSE, 
                color = "blue", lwd=.5, lty=1)
} 

#------------------------------------------------------------------------------
# Community Biomass Proportion Vs Rainfall: Site Average: All survey data
#------------------------------------------------------------------------------
d_site_mu <- tb_stats %>%
  group_by(site_code, transient_type) %>%
  summarize(mu_percent_biomass = 
              mean(transient_type_biomass_percent, na.rm=T)) %>%
  add_rain() 

t_site_mu <- d_site_mu %>%
  group_by(transient_type) %>%
  nest() %>%
  mutate(lm = map(data, table_lm_stats)) %>%
  unnest(lm) %>%
  filter(term == 'annualrain') %>%
  select(-term, -data) 

p_site_mu <- d_site_mu %>%
  left_join(t_site_mu) %>%
  visualize_biomass_vs_rain() +
  facet_wrap(~transient_type)

#------------------------------------------------------------------------------
# Community Biomass Proportion Vs Rainfall: (by Year)
#------------------------------------------------------------------------------
d_by_year <- tb_stats %>%
  mutate(year = year(ymd(collection_period))) %>%
  group_by(year, site_code, transient_type) %>%
  summarize(mu_percent_biomass = 
              mean(transient_type_biomass_percent, na.rm=T)) %>%
  add_rain() 

t_by_year <- d_by_year %>%
  group_by(year, transient_type) %>%
  nest() %>%
  mutate(lm = map(data, table_lm_stats)) %>%
  unnest(lm) %>%
  filter(term == 'annualrain') %>%
  select(-term, -data) 

p_by_year <- d_by_year %>%
  left_join(t_by_year) %>%
  visualize_biomass_vs_rain() +
  facet_grid(year~transient_type)

#------------------------------------------------------------------------------
# Community Biomass Proportion Vs Rainfall: (by Quarter)
#------------------------------------------------------------------------------
d_by_quarter <- tb_stats %>%
  mutate(quarter = quarter(ymd(collection_period))) %>%
  group_by(quarter, site_code, transient_type) %>%
  summarize(mu_percent_biomass = 
              mean(transient_type_biomass_percent, na.rm=T)) %>%
  add_rain() %>%
  ungroup()

t_by_quarter <- d_by_quarter %>%
  group_by(quarter, transient_type) %>%
  nest() %>%
  mutate(lm = map(data, table_lm_stats)) %>%
  unnest(lm) %>%
  filter(term == 'annualrain') %>%
  select(-term, -data) %>%
  ungroup()

p_by_quarter <- d_by_quarter %>%
  left_join(t_by_quarter) %>%
  visualize_biomass_vs_rain() +
  facet_grid(quarter~transient_type)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(tb_stats, here::here('03_public', 'output', 
                               '15_transient_biomass_proportion.csv'))

#------------------------------------------------------------------------------
# End 15_transient_proportion