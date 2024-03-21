# 08_environment_short_term_preview
# Sean Kinard
# 2023-07-18
#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('03_public/toolkit.R')

ds <- read_csv("03_local_files/data/environment/environment_shortterm.csv") %>%
  filter(collection_period > ymd('2019-12-31') &
           collection_period < ymd('2020-12-31')) %>%
  select(-c('nitrate', 'phosphate', 'ammonia', 'doc'))

dn <- read_csv('03_local_files/data/environment/nutrient_TERRG.csv') %>%
  rename(nitrate=no3n, phosphate=ortho_p, ammonia=nh4_n) %>%
  create_combo_period() %>%
  select(-site_period, -collection_date)

# merge nutrient to habitat
d_ste <- full_join(ds, dn) %>% 
  pivot_longer(cols=-c('site_code', 'collection_period'), 
               names_to='xname', values_to='xvalue') %>%
  group_by(site_code, collection_period, xname) %>%
  dplyr::summarize(xvalue=mean(xvalue, na.rm=T)) %>%
  pivot_wider(names_from=xname, values_from=xvalue)

# impute missing with interpolation within site
d_ste <- d_ste %>% impute_interpolation()

index_vars <- c('site_code', 'collection_period')
remove_index <-function(df) {
  select(df, -any_of(index_vars)) }

#------------------------------------------------------------------------------
# natural log transformation
#------------------------------------------------------------------------------
ln_vars <- d_ste %>%
  remove_index() %>%
  pivot_longer(cols=everything(), names_to='xvar', values_to='xval') %>%
  group_by(xvar) %>%
  dplyr::summarize(range_over_mean = range(xval, na.rm=T)/median(xval, na.rm=T)) %>%
  filter(range_over_mean > 50) %>%
  pull(xvar)

# visual verfication
plot_hist_lnvars_raw <- d_ste %>% 
  select(index_vars, any_of(ln_vars)) %>%
  pivot_longer(cols=-any_of(index_vars), 
               names_to='xvar', values_to='xval') %>%
  ggplot(aes(xval)) +
  facet_wrap(~xvar, scales='free') +
  geom_histogram()

# log transformation
d_ste_ln <- d_ste %>% 
  select(index_vars, any_of(ln_vars)) %>%
  pivot_longer(cols=-any_of(index_vars), 
               names_to='xvar', values_to='xval') %>%
  mutate(xval = log(xval),
         xval = ifelse(is.infinite(xval), NA, xval)) %>%
  pivot_wider(names_from='xvar', values_from='xval')
colnames(d_ste_ln) <- paste('ln_', colnames(d_ste_ln), sep='')
d_ste_ln <- d_ste_ln %>% rename(site_code = ln_site_code,
                                collection_period = ln_collection_period)
d_ste <- d_ste %>% select(-any_of(ln_vars)) %>% left_join(d_ste_ln)

ln_names <- paste('ln_', ln_vars, sep='')

# visual verfication
plot_hist_lnvars_transformed <- d_ste %>% 
  select(index_vars, contains('ln_')) %>%
  pivot_longer(cols=-any_of(index_vars), 
               names_to='xvar', values_to='xval') %>%
  ggplot(aes(xval)) +
  facet_wrap(~xvar, scales='free') +
  geom_histogram()

#------------------------------------------------------------------------------
# Root-transformation
#------------------------------------------------------------------------------
sqrt_vars <- d_ste %>%
  remove_index() %>%
  select(-contains('ln_')) %>%
  pivot_longer(cols=everything(), names_to='xvar', values_to='xval') %>%
  group_by(xvar) %>%
  dplyr::summarize(range_over_mean = range(xval, na.rm=T)/median(xval, na.rm=T)) %>%
  filter(range_over_mean > 10) %>%
  pull(xvar)

# visual verfication
plot_hist_sqrtvars_raw <- d_ste %>% 
  select(index_vars, any_of(sqrt_vars)) %>%
  pivot_longer(cols=-any_of(index_vars), 
               names_to='xvar', values_to='xval') %>%
  ggplot(aes(xval)) +
  facet_wrap(~xvar, scales='free') +
  geom_histogram()

# log transformation
d_ste_sqrt <- d_ste %>% 
  select(index_vars, any_of(sqrt_vars)) %>%
  pivot_longer(cols=-any_of(index_vars), 
               names_to='xvar', values_to='xval') %>%
  mutate(xval = log(xval),
         xval = ifelse(is.infinite(xval), NA, xval)) %>%
  pivot_wider(names_from='xvar', values_from='xval')
colnames(d_ste_sqrt) <- paste('sqrt_', colnames(d_ste_sqrt), sep='')
d_ste_sqrt <- d_ste_sqrt %>% rename(site_code = sqrt_site_code,
                                collection_period = sqrt_collection_period)
d_ste <- d_ste %>% select(-any_of(sqrt_vars)) %>% left_join(d_ste_sqrt)

ln_names <- paste('sqrt_', sqrt_vars, sep='')

# visual verfication
plot_hist_sqrtvars_transformed <- d_ste %>% 
  select(index_vars, contains('sqrt_')) %>%
  pivot_longer(cols=-any_of(index_vars), 
               names_to='xvar', values_to='xval') %>%
  ggplot(aes(xval)) +
  facet_wrap(~xvar, scales='free') +
  geom_histogram()

# restore ln_vars to categories
climate <- restore_category(climate)
landuse <- restore_category(landuse)
flow <- restore_category(flow)
water_quality <- restore_category(water_quality)
geomorph <- restore_category(geomorph)
algae <- restore_category(algae)

#------------------------------------------------------------------------------
# autocorrelation
#------------------------------------------------------------------------------
table_ste_correlation <- d_ste %>% 
  correlate() %>%
  shave() %>% 
  fashion()

plot_ste_correlation <- d_ste %>% 
  correlate() %>%
  shave() %>%
  rplot(colours = c("#F24405FF", "black", "#00FFFFFF")) +
  dark_theme_grey() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

#------------------------------------------------------------------------------
# densiogram
#------------------------------------------------------------------------------
ex_hist <- function(my_data) {
  my_data %>% 
    pivot_longer(cols=-site_code, names_to='xvar', values_to='xval') %>%
    ggplot(aes(xval, fill=site_code, color =site_code)) +
    geom_density(fill = NA) +
    geom_density(color = NA, alpha=.1)+
    facet_wrap(~xvar, scales='free') +
    dark_theme_grey(base_size = 14) +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    xlab(element_blank()) +
    ylab(element_blank())}

plot_ste_densiogram_flow <- d_ste %>% 
  select(site_code, any_of(flow)) %>%
  ex_hist() # need to log transform many q variables

plot_ste_densiogram_waterquality <- d_ste %>% 
  select(site_code, any_of(water_quality)) %>%
  ex_hist() # need to log transform conductivity

plot_ste_densiogram_geomorph <- d_ste %>% 
  select(site_code, any_of(geomorph)) %>%
  ex_hist()

plot_ste_densiogram_algae <- d_ste %>% 
  select(site_code, any_of(algae)) %>%
  ex_hist()

# -----------------------------------------------------------------------------
# violin
# -----------------------------------------------------------------------------
ex_violin <- function(my_data) {
  my_data %>%
    fix_site_order() %>%
    pivot_longer(cols= -any_of(index_vars),
                 names_to = "xvar",
                 values_to = "xval") %>%
    ggplot(aes(x=site_code, y=xval, fill=site_code)) +
    facet_wrap(~xvar, scale = 'free') +
    geom_violin(show.legend = F) +
    ylab(element_blank())+
    scale_fill_manual(values = my_colors) +
    labs(fill='Site') +
    xlab(element_blank()) +
    ylab(element_blank()) +
    dark_theme_grey(base_size=14) }

plot_ste_violin_flow <- select(d_ste, site_code, any_of(flow)) %>% ex_violin()
plot_ste_violin_geomorph <- d_ste %>% 
  select(site_code, any_of(geomorph), any_of(ln_vars)) %>% ex_violin()
plot_ste_violin_algae <-  select(d_ste, site_code, any_of(algae)) %>% ex_violin()
plot_ste_violin_water_quality <- select(d_ste, site_code, any_of(water_quality)) %>%
  ex_violin()

# -----------------------------------------------------------------------------
# time series
# -----------------------------------------------------------------------------
time_series <- function(my_data, x_span) {
  my_data %>%
    add_rain() %>%
    pivot_longer(cols = -any_of(c('annualrain', index_vars)),
                 names_to = "xvar",
                 values_to = "xval") %>%
    fix_site_order() %>%
    ggplot(aes(x=collection_period, y=xval, fill= site_code)) +
    facet_wrap(~xvar, scale = 'free') +
    geom_point(size=3, color = 'black', shape =21, alpha=.5,
               show.legend = F) +
    geom_smooth(aes(color = site_code),
                method = "loess", size = 1.3, se=FALSE, lty=1, lwd=.4, 
                span=x_span,
                show.legend = T) +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    ylab(element_blank()) +
    xlab(element_blank()) +
    dark_theme_grey(base_size=14) +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) }

plot_ste_timeseries_flow <- d_ste %>% 
    select(any_of(index_vars), any_of(flow)) %>%
    time_series(x_span=.4)

plot_ste_timeseries_water_quality <- d_ste %>% 
  select(any_of(index_vars), any_of(water_quality)) %>%
    time_series(x_span=.4)

plot_ste_timeseries_geomorph <- d_ste %>% 
    select(any_of(index_vars), any_of(geomorph)) %>%
    time_series(x_span=.4)

plot_ste_timeseries_algae <- d_ste %>% 
    select(any_of(index_vars), any_of(algae)) %>%
    time_series(x_span=.4)

# -----------------------------------------------------------------------------
# time series (SCALED)
# -----------------------------------------------------------------------------
scaling <- function(my_data) {
  x <- my_data
  
  old_cols <- colnames(x)
  
  for (i in 3:length(colnames(my_data))) {
    
    new_col <- paste(colnames(x[i]), 
                     'scaled',
                     sep='_')
    
    x <- x %>%
      select(1,2,i) %>%
      pivot_wider(names_from = site_code, 
                  values_from = 3) %>%
      mutate(AR = scale(AR),
             EM = scale(EM),
             GC = scale(GC),
             MR = scale(MR),
             PD = scale(PD),
             PL = scale(PL),
             SF = scale(SF),
             TR = scale(TR),
             WM = scale(WM) ) %>%
      pivot_longer(cols= AR:WM,
                   names_to = "site_code",
                   values_to = new_col,
                   values_drop_na = T) %>%
      right_join(x) %>%
      select(all_of(old_cols), new_col)
    
    old_cols <- colnames(x)  }
  
  x <- select(x, site_code, collection_period, contains('scaled')) %>%
    add_rain()
  
  return(x) } 

sitescaled_timeseries <- function(my_data, x_span) {
  my_data %>%
    scaling %>%
    add_rain() %>%
    pivot_longer(cols = -any_of(c('annualrain', index_vars)),
                 names_to = "xvar",
                 values_to = "xval") %>%
    fix_site_order() %>%
    ggplot(aes(x=collection_period, y=xval, fill= site_code)) +
    facet_wrap(~xvar, scale = 'free') +
    geom_point(size=3, color = 'black', shape =21, alpha=.5,
               show.legend = F) +
    geom_smooth(aes(group=NA),
                method = "loess", size = 1.3, se=FALSE, 
                linetype=1, span=x_span, color = 'chartreuse3',
                show.legend = F) +
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    ylab(element_blank()) +
    xlab(element_blank()) +
    dark_theme_grey(base_size=14) +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) }

plot_ste_sitescaled_timeseries_flow <- d_ste %>% 
  select(any_of(index_vars), any_of(flow)) %>%
  sitescaled_timeseries(x_span=.1)

plot_ste_sitescaled_timeseries_water_quality <- d_ste %>% 
    select(any_of(index_vars), any_of(water_quality)) %>%
  sitescaled_timeseries(x_span=.1)

plot_ste_sitescaled_timeseries_geomorph <- d_ste %>% 
    select(any_of(index_vars), any_of(geomorph)) %>%
  sitescaled_timeseries(x_span=.1)

plot_ste_sitescaled_timeseries_algae <- d_ste %>% 
    select(any_of(index_vars), any_of(algae)) %>%
  sitescaled_timeseries(x_span=.1)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(d_ste, 
          '03_public/output/environment_short_term_post.csv')
write_csv(table_ste_correlation,
          '03_public/output/environment_short_term_correlation.csv')


ggsave(filename = '03_public/visualization/08_plot_ste_correlation.png',
       plot_ste_correlation,
       width=9, height=9, units="in")

# density figures
ggsave(filename = '03_public/visualization/08_plot_ste_densiogram_flow.png',
       plot_ste_densiogram_flow,
       width=9, height=9, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_densiogram_waterquality.png',
       plot_ste_densiogram_waterquality,
       width=9, height=9, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_densiogram_geomorph.png',
       plot_ste_densiogram_geomorph,
       width=9, height=6, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_densiogram_algae.png',
       plot_ste_densiogram_algae,
       width=9, height=3, units="in")

# violin plots
ggsave(filename = '03_public/visualization/08_plot_ste_violin_flow.png',
       plot_ste_violin_flow,
       width=9, height=9, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_violin_water_quality.png',
       plot_ste_violin_water_quality,
       width=9, height=9, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_violin_geomorph.png',
       plot_ste_violin_geomorph,
       width=9, height=6, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_violin_algae.png',
       plot_ste_violin_algae,
       width=9, height=3, units="in")

# timeseries plots
ggsave(filename = '03_public/visualization/08_plot_ste_timeseries_flow.png',
       plot_ste_timeseries_flow,
       width=15, height=12, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_timeseries_water_quality.png',
       plot_ste_timeseries_water_quality,
       width=9, height=9, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_timeseries_geomorph.png',
       plot_ste_timeseries_geomorph,
       width=9, height=6, units="in")
ggsave(filename = '03_public/visualization/08_plot_ste_timeseries_algae.png',
       plot_ste_timeseries_algae,
       width=9, height=3, units="in")

# sitescaled_timeseries plots
ggsave(
  filename = '03_public/visualization/08_plot_ste_sitescaled_timeseries_flow.png',
  plot_ste_sitescaled_timeseries_flow,
  width=9, height=9, units="in")
ggsave(
  filename =
    '03_public/visualization/08_plot_ste_sitescaled_timeseries_water_quality.png',
  plot_ste_sitescaled_timeseries_water_quality,
  width=9, height=9, units="in")
ggsave(
  filename =
    '03_public/visualization/08_plot_ste_sitescaled_timeseries_geomorph.png',
  plot_ste_sitescaled_timeseries_geomorph,
  width=9, height=6, units="in")
ggsave(
  filename =
    '03_public/visualization/08_plot_ste_sitescaled_timeseries_algae.png',
  plot_ste_sitescaled_timeseries_algae,
  width=9, height=3, units="in")

#------------------------------------------------------------------------------
# End 08_environment_long_term_preview