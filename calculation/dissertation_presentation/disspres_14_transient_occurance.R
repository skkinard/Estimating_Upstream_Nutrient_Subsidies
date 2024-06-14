# 14_diadromy_occurance
# Sean Kinard
# 2023-07-20

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('03_public/toolkit.R') # load packages and helper-functions
library(patchwork)
library(ggpubr)
library(ggpmisc)

d <- read_csv('03_public/output/community_efish_2017_2020.csv')

d <- d %>%
  rowwise() %>%
  mutate(pass_total = sum(pass1, pass2, pass3, na.rm=T)) %>%
  ungroup() %>%
  right_join(d)

d <- d %>%
  group_by(site_code, collection_date) %>%
  dplyr::summarize(ef_event_total = sum(pass_total, na.rm=T)) %>%
  right_join(d)

probability_of_occurance <- d %>%
  ungroup() %>%
  group_by(guild, site_code, collection_date, transient_type) %>%
  dplyr::summarize(pass_total = ifelse(sum(pass_total, na.rm=T) > 0, 1 , 0)) %>%
  ungroup() %>%
  pivot_wider(names_from=transient_type, values_from=pass_total, 
              values_fill = 0) %>%
  pivot_longer(cols=c('Catadromous', 'Amphidromous', 'Amphidromous', 
                      'Anadromous'),
               names_to='transient_type', values_to='detected') %>%
  group_by(site_code, transient_type) %>%
  dplyr::summarize(probability_of_occurance=sum(detected) / length(detected))

extract_lm_stats <- function(xdata) {
  lm(probability_of_occurance~annualrain, data = xdata) %>%
    glance() %>%
    select(r.squared, p.value) }

d_transient_prob <- probability_of_occurance %>%
  add_rain() %>%
  group_by(transient_type) %>%
  nest() %>%
  mutate(lm_stats=map(data, extract_lm_stats)) %>%
  unnest(lm_stats) %>%
  unnest(data) %>%
  fix_site_order() %>%
  mutate(probability_of_occurance=probability_of_occurance*100)

plot_transient_prob <-  d_transient_prob %>%
  filter(transient_type =='Amphidromous') %>%
  ggplot(aes(annualrain, probability_of_occurance)) +
  geom_hline(yintercept=50, lty=2, lwd=.3, color='grey') +
  geom_smooth(data = d_transient_prob %>% filter(p.value < 0.05),
              method='lm', se=F, lwd=.6, lty=2, show.legend=F, color='red3') +
  geom_point(size=4, shape = 21, fill='cyan', alpha=.5) +
  geom_point(size=4, shape = 21, color='cyan', fill=NA) +
  stat_poly_eq(label.x=.5, label.y=.95,
               color='white', use_label(c("adj.R2","p")), size = 6) +
  scale_fill_manual(values=my_colors) +
  dark_theme_grey(base_size=20) +
  ylim(c(0,100)) + 
  labs(x='Annual Rain (cm)', y='Amphidromous Presence (%)')

d_transient_prob_prep <- d_transient_prob %>%
  add_sitevars() %>%
  pivot_longer(cols=c(annualrain, baydist_km), 
               names_to = 'XNAME', values_to = 'XVALUE') %>%
  mutate(XNAME = ifelse(XNAME=='annualrain', 'Rain (cm/yr)', 
                        'Estuary Distance (km)')) %>%
  filter(transient_type=='Amphidromous')

plot_transient_prob2 <- d_transient_prob_prep %>%
  ggplot(aes(XVALUE, probability_of_occurance)) +
  facet_wrap(~XNAME, scales='free') +
  geom_smooth(data = d_transient_prob_prep %>% filter(XNAME == 'Rain (cm/yr)'),
              method='lm', se=F, lwd=.6, lty=2, show.legend=F, color='red3') +
  geom_point(size=4, color='black') +
  stat_poly_eq(label.x=.5, label.y=.95,
               color='black', use_label(c("adj.R2","p"))) +
  scale_fill_manual(values=my_colors) +
  theme_bw(base_size=12) +
  ylim(c(0,100)) + 
  labs(x=element_blank(), y='Detection (%)')

#------------------------------------------------------------------------------
# Proportion of Transient Fish of Total Fish
#------------------------------------------------------------------------------
# Average proportion of transient within caught individuals
d_fraction <- d %>%
  ungroup() %>%
  group_by(guild, site_code, collection_date, transient_type) %>%
  dplyr::summarize(pass_total = sum(pass_total, na.rm=T)) %>%
  left_join(select(d, site_code, collection_date, ef_event_total)) %>%
  mutate(catch_proportion=pass_total/ef_event_total*100) %>%
  filter(transient_type != 'Freshwater')

plot_transient_fraction <- d_fraction %>%
  add_rain() %>%
  fix_site_order() %>%
  ggplot(aes(x=annualrain, y=catch_proportion)) +
  facet_wrap(~transient_type, ncol=1) +
  geom_smooth(method='lm', se=F, lty=2, color='red3', lwd=.5, span=.8) +
  geom_boxplot(aes(fill = site_code), notch=T, width=2) +
  stat_poly_eq(label.x=.99, label.y=.98, 
               color='black', use_label(c("adj.R2", "p"))) +
  scale_fill_manual(values=my_colors) +
  theme_bw(base_size=12) +
  labs(x='Annual Rain (cm)', y='Transient Proportion of Fish (%) [When Present]')

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
# write_csv(d_transient_prob, '03_public/output/d_transient_prob.csv')
# 
# ggsave('03_public/visualization/14_transient_occurance_probability.png',
#        plot=plot_transient_prob, width=9, height=12, units='in')
# 
# ggsave('03_public/visualization/14_transient_occurance_fraction.png',
#        plot=plot_transient_fraction, width=9, height=12, units='in')

#------------------------------------------------------------------------------
# End 14_diadromy_occurance