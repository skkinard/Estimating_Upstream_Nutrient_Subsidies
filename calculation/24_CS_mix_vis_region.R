# 24_CS_mix_vis_region
# written by: Sean Kinard
# last edit: 2023-09-17
#----------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------
source('03_public/toolkit.R') # load packages and helper-function

d <- read_csv('03_public/output/CS_mix_out_region.csv')

d_iso_raw <- read_csv('03_public/output/isotope_CNS_2020_01_clean.csv')

table_UN_prediction <- read_csv('03_public/output/table_UN_prediction.csv')

d_iso <- d_iso_raw %>% filter(! species %in% bad_sources) %>%
  filter(species !='Periphyton')

site_vars <- read_csv('03_local_files/data/isotope_raw/MDN_clean_final.csv') %>%
  r_friendly_colnames() %>%
  rename(site_code=site) %>%
  select(site_code, site_type, pptavg_site, latitude, longitude, nearest_bay, 
         baydist_km, elev_site_m) %>%
  mutate(site_type = str_replace_all(site_type, 'Fresh', 'Stream')) %>%
  unique() %>%
  pivot_longer(cols=-any_of(c('site_code', 'site_type', 'nearest_bay', 
                              'pptavg_site', 'latitude', 'longitude')),
               names_to='xname', values_to='xvalue') %>%
  mutate(xvalue = ifelse(xvalue==-999, NA,xvalue)) %>%
  pivot_wider(names_from=xname, values_from=xvalue) %>%
  arrange(desc(site_type), longitude) %>%
  left_join(read_csv("03_local_files/data/environment/MDN_site_lat_lon.csv") %>%
              r_friendly_colnames() %>%
              select(site_name, site_code, staid)) %>%
  select(site_code, site_name, staid, latitude, longitude, everything())
#----------------------------------------------------------------------------
# Iso Tables
#----------------------------------------------------------------------------
table_iso_sources <- d_iso_raw %>%
  filter(! site_code %in% c('LN', 'UN')) %>%
  add_sitevars() %>%
  filter(guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  mutate(lowest_taxon = case_when(
    lowest_taxon == 'Filamentous_algae'~ 'Algae',
    lowest_taxon == 'Coarse_detritus' ~ 'Detritus',
    lowest_taxon == 'Mix_tv' ~ 'Riparian Leaves',
    T ~ lowest_taxon)) %>%
  group_by(site_type, site_code, lowest_taxon) %>%
  dplyr::summarize(C13 = extract_iso(carbon),
                   S34 = extract_iso(sulfur)) %>%
  ungroup()

table_iso_stream <- d_iso_raw %>%
  add_sitevars() %>%
  filter(! site_code %in% c('LN', 'UN')) %>%
  filter(! guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  group_by(site_type, site_code, guild, family) %>%
  dplyr::summarize(C13 = extract_iso(carbon),
                   S34 = extract_iso(sulfur)) %>%
  ungroup()

table_iso_dam <- d_iso_raw %>%
  add_sitevars() %>%
  filter(site_code %in% c('LN', 'UN')) %>%
  filter(! guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  group_by(site_type, site_code, guild, family) %>%
  dplyr::summarize(C13 = extract_iso(carbon),
                   S34 = extract_iso(sulfur)) %>%
  ungroup()

#----------------------------------------------------------------------------
# Source Stats
#----------------------------------------------------------------------------
d_sources <- d_iso %>%
  filter(guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  pivot_longer(cols=c(carbon, sulfur), 
               names_to = 'x_name', 
               values_to = 'x_val') %>%
  fix_site_order() %>%
  fix_source_order() %>%
  create_site_type() %>%
  filter(site_type %in% c('Stream', 'Estuary')) %>%
  select(guild, species, site_code, site_type, collection_date, x_name, x_val)

# boostrap means and CI
# Bootstrap estimate means and 95% CI
source_boot <- d_sources %>%
  rename(isotope=x_name) %>%
  group_by(site_type, isotope) %>%
  filter(!is.na(x_val)) %>%
  dplyr::summarize(
    mu = estimate_mu(x_val),
    lower = estimate_mu_lower(x_val),
    upper = estimate_mu_upper(x_val) )

site_list <- unique(source_boot$site_type)

plot_source_boot <- source_boot %>%
  ggplot(aes(x=site_type, y=mu)) +
  facet_wrap(~isotope, labeller = iso_facet_lab) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, color='black') +
  geom_point(size=3) +
  geom_text(data=source_boot%>%filter(site_type=='Estuary'),
            aes(y=upper+1.5), color = 'red', label='A', size = 5) +
  geom_text(data=source_boot%>%filter(site_type=='Stream'),
            aes(y=upper+1.5), color = 'red', label='B', size = 5) +
  labs(x=element_blank(), y=expression(delta)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())

table_source_boot <- source_boot

table_summary_stats_sources <- d_sources %>%
  rename(isotope=x_name) %>%
  group_by(site_type, guild, isotope) %>%
  filter(!is.na(x_val)) %>%
  dplyr::summarize(x_mu = mean(x_val, na.rm=T),
                   x_sd = sd(x_val, na.rm=T),
                   x_n = length(x_val),
                   x_min = min(x_val, na.rm=T),
                   x_max = max(x_val, na.rm=T)) %>%
  ungroup() %>%
  mutate(iso_abr = substr(isotope, 1,1)%>%str_to_upper(),
         site_label = ifelse(site_type=='Estuary', 'Estuary', 'Stream'),
         plot_label = paste(site_label, iso_abr, sep='-'))

#----------------------------------------------------------------------------
# scatterplot transient
#----------------------------------------------------------------------------
d_source_stats <- d_sources %>%
  rename(isotope=x_name) %>%
  group_by(site_type, isotope) %>%
  filter(!is.na(x_val)) %>%
  dplyr::summarize(x_mu = mean(x_val, na.rm=T),
                   x_sd = sd(x_val, na.rm=T)) %>%
  ungroup() %>%
  pivot_wider(names_from=isotope, values_from = x_mu:x_sd)

colnames(d_source_stats) <- str_replace_all(colnames(d_source_stats), 'x_', '')

plot_scatter_transient <- d_iso %>%
  filter(guild %in% c('Fish', 'Invertebrate')) %>%
  ggplot() +
  geom_errorbar(data=d_source_stats, 
                aes(x=mu_carbon, 
                    ymin=mu_sulfur-sd_sulfur, 
                    ymax=mu_sulfur+sd_sulfur), width=0) +
  geom_errorbarh(data=d_source_stats, 
                aes(y=mu_sulfur, 
                    xmin=mu_carbon-sd_carbon, 
                    xmax=mu_carbon+sd_carbon), height=0) +
  geom_point(aes(x=carbon, y=sulfur, 
                 fill=is_diadromous, shape = is_diadromous, color=is_diadromous),
             size=3, alpha=.5) +
  geom_point(aes(x=carbon, y=sulfur, shape = is_diadromous, color=is_diadromous),
             size=3, fill=NA) +
  geom_label(data=d_source_stats, 
             aes(y=mu_sulfur,x=mu_carbon, label=site_type)) +
  scale_shape_manual('', values=21:24) +
  scale_fill_manual('', values=c('cyan3', 'blue', 'red')) +
  scale_color_manual('', values=c('cyan3', 'blue', 'red')) +
  labs(x=expression(paste(delta)^13*C),
       y=expression(paste(delta)^34*S)) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.24, 0.85)) +
  theme(legend.background = element_rect(colour = 'grey50', fill = 'white', 
                                         linetype='solid')) +
  theme(legend.title = element_blank(),
        legend.margin = margin(3, 3, 3, 3),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm")) +
  coord_equal()


#----------------------------------------------------------------------------
# scatterplot Dam
#----------------------------------------------------------------------------
plot_scatter_dam <- d_iso %>%
  filter(site_code %in% c("UN", "LN")) %>%
  add_is_dam() %>%
  mutate(is_dam = ifelse(is_dam == 'Dam', 'Above Dam', 'Below Dam')) %>%
  filter(guild %in% c('Fish', 'Invertebrate')) %>%
  ggplot() +
  geom_errorbar(data=d_source_stats, 
                aes(x=mu_carbon, 
                    ymin=mu_sulfur-sd_sulfur, 
                    ymax=mu_sulfur+sd_sulfur), width=0) +
  geom_errorbarh(data=d_source_stats, 
                 aes(y=mu_sulfur, 
                     xmin=mu_carbon-sd_carbon, 
                     xmax=mu_carbon+sd_carbon), height=0) +
  geom_text_repel(aes(x=carbon, y=sulfur, label=lowest_taxon)) +
  geom_point(aes(x=carbon, y=sulfur),
             size=3, shape=21, fill='white') +
  geom_point(aes(x=carbon, y=sulfur, fill=is_dam),
             size=3, shape=21, alpha=.8) +
  geom_label(data=d_source_stats, 
             aes(y=mu_sulfur,x=mu_carbon, label=site_type)) +
  scale_shape_manual('', values=21:24) +
  scale_fill_manual('', values=c('red', 'blue')) +
  labs(x=expression(paste(delta)^13*C),
       y=expression(paste(delta)^34*S)) +
  theme_bw(base_size=20) +
  theme(legend.position = c(0.2, 0.85)) +
  theme(legend.background = element_rect(colour = 'grey50', fill = 'white', 
                                         linetype='solid')) +
  theme(legend.title = element_blank())
  

#----------------------------------------------------------------------------
# Transient x Estuarine
#----------------------------------------------------------------------------
d_prep <- d %>%
  filter(dataset == 'is_diadromous') %>%
  mutate(m_group = fct_relevel(m_group, c('Freshwater', 'Catadromous',
                                          'Amphidromous', 'Esturaine')))

plot_transient_estuarine <- d_prep %>%
  ggplot(aes(x=m_group, y=mean)) +
  geom_hline(yintercept=50, lty=2, lwd=.5, color='grey') +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width = 0.2) +
    geom_point(size=4) +
  geom_text(data = filter(d_prep, m_group == 'Freshwater'),
            aes(y=`97.5%`+5), label='A', color='red', size= 5) +
  geom_text(data = filter(d_prep, m_group == 'Diadromous'),
            aes(y=`97.5%`+5), label='A', color='red', size= 5) +
  geom_text(data = filter(d_prep, m_group == 'Euryhaline'),
            aes(y=`97.5%`+5), label='B', color='red', size= 5) +
  labs(x=element_blank(), y='% Estuarine') +
  ylim(c(0,100)) +
  theme_bw(base_size=20)

vis_mix_CI <- function(xdata) {
  xdata %>%
    ggplot(aes(x=m_group, y=mean)) +
    geom_hline(yintercept=50, lty=2, lwd=.5, color='grey') +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width = 0.2) +
    geom_point(size=4) +
    geom_text(data = filter(xdata, mean > 50),
              aes(y=`97.5%`+5), label='*', color='red', size= 10) +
    labs(x=element_blank(), y='% Estuarine') +
    ylim(c(0,100)) +
    theme_bw(base_size=15)
}

d_vis_mix_CI <- d %>%
  mutate(m_group = case_when(m_group == 'No Dam' ~ 'Below C.', 
                             m_group == 'Dam' ~ 'Above C.',
                             T ~ m_group),
    m_group = fct_relevel(
    m_group, 
    c('Freshwater', 'Diadromous', 'Euryhaline', 'Below C.'))) %>%
  group_by(dataset) %>%
  nest() %>%
  mutate(mix_CI = map(data, vis_mix_CI))

plot_rain_dist <- d_vis_mix_CI$data[[1]] %>%
  rename(site_code=m_group) %>%
  add_sitevars() %>%
  add_rain() %>%
  add_is_dam() %>%
  filter(is_dam == 'No Dam') %>%
  pivot_longer(cols=c(baydist_km, annualrain), 
               names_to='X_NAME', values_to = 'X_VALUE' ) %>%
  mutate(X_NAME = str_replace_all(X_NAME, 'annualrain', 'Rainfall (cm/yr)'),
         X_NAME = str_replace_all(X_NAME, 'baydist_km', 'Distance (km)')) %>%
  ggplot(aes(x=X_VALUE, y=mean)) +
  facet_wrap(~ X_NAME, scales='free_x') +
  stat_poly_eq(label.x=.5, label.y=.95, formula=y~log(x),
               color='black', use_label(c("adj.R2","p")), size=6) +
  geom_point(size=4, fill='red', shape=21) +
  scale_fill_viridis_c('Distance (km)', direction=-1) +
  labs(x=element_blank(), y='% Estuarine') +
  ylim(c(0,100)) +
  theme_bw(base_size=20)
  
plot_mix_multi <- d_vis_mix_CI$mix_CI[[2]] +
    ggtitle('Estuary Distance') +
  d_vis_mix_CI$mix_CI[[3]] +
    ggtitle('Transient Type') +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y=element_blank()) +
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1)) +
  d_vis_mix_CI$mix_CI[[4]] +
    ggtitle('Dams') +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y=element_blank()) +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))

table_transient_estuarine <- d %>%
  mutate(m_group = case_when(m_group == 'No Dam' ~ 'Below C.D.', 
                             m_group == 'Dam' ~ 'Above C.D.',
                             T ~ m_group),
         m_group = fct_relevel(
           m_group, 
           c('Amphidromous', 'Catadromous', 'Freshwater', 'Below C.D.'))) %>%
  mutate(dataset = case_when(
    dataset == 'site' ~ 'Site',
    dataset == 'distcat' ~ 'Distance',
    dataset == 'transient_type' ~ 'Transient',
    dataset == 'calallen' ~ 'Dam')) %>%
  rename(transient_type = m_group, 
         estuarine=mean, lower=`2.5%`, upper=`97.5%`,
         Comparison = dataset) %>%
  select(Comparison, transient_type, estuarine, lower, upper)

#----------------------------------------------------------------------------
# Dam x Estuarine
#----------------------------------------------------------------------------
d_prep_cal <- d %>%
  filter(dataset == 'calallen') %>%
  filter(m_group %in% c('Dam', 'No Dam')) %>%
  mutate(m_group = ifelse(m_group == 'Dam', 'Above Dam', 'Below Dam')) %>%
  left_join(table_UN_prediction %>%
              mutate(m_group = ifelse(site_code == 'UN', 'Above Dam', 'Below Dam')) %>%
              select(m_group, pred)) %>%
  rename(Observed=mean, Predicted=pred) %>%
  pivot_longer(cols=any_of(c('Observed', 'Predicted')), 
               names_to='type', values_to='x_val') %>%
  mutate(`2.5%` = ifelse(type=='Predicted', NA, `2.5%`),
         `97.5%` = ifelse(type=='Predicted', NA, `97.5%`))

plot_dam_estuarine <- d_prep_cal %>%
  ggplot(aes(x=m_group)) +
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`, group=type), 
                width = 0.2, position=position_dodge(width=.1)) +
  geom_point(aes(y=x_val, shape=type, fill=type),
             size=4, position=position_dodge(width=.1)) +
  scale_fill_manual(values=c('black', 'red')) +
  scale_shape_manual(values=c(21, 22)) +
  labs(x=element_blank(), y='% Estuarine') +
  ylim(c(0,100)) +
  theme_bw(base_size=15) +
  theme(legend.position = c(0.25, 0.85)) +
  theme(legend.background = element_rect(colour = 'grey50', fill = 'white', 
                                         linetype='solid')) +
  theme(legend.title = element_blank(),
        legend.margin = margin(3, 3, 3, 3),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))

#----------------------------------------------------------------------------
# End 24_CS_mix_vis