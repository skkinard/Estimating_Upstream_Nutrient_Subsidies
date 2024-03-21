# 24_CS_mix_vis
# written by: Sean Kinard
# last edit: 2023-09-17
#----------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------
source('03_public/toolkit.R') # load packages and helper-function

d <- read_csv('03_public/output/CS_mix_out.csv')

#----------------------------------------------------------------------------
# Plots
#----------------------------------------------------------------------------
d %>%
  filter(dataset == 'transient_type') %>%
  rename(site_code=site) %>%
  add_sitevars() %>%
  filter(site_type == 'Stream') %>%
  add_rain() %>%
  ggplot(aes(x=annualrain, y=mean)) +
  facet_wrap(~m_group) +
  geom_point()

d %>%
  filter(dataset == 'order') %>%
  rename(site_code=site) %>%
  add_sitevars() %>%
  filter(site_type == 'Stream') %>%
  add_rain() %>%
  ggplot(aes(x=annualrain, y=mean)) +
  facet_wrap(~m_group) +
  geom_point()

d %>%
  filter(dataset == 'species') %>%
  rename(site_code=site) %>%
  add_sitevars() %>%
  filter(site_type == 'Stream') %>%
  add_rain() %>%
  ggplot(aes(x=annualrain, y=mean)) +
  facet_wrap(~m_group) +
  geom_point()

#----------------------------------------------------------------------------
# Table
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Export
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# End 24_CS_mix_vis