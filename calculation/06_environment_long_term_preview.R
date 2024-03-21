# 06_environment_long_term_preview
# Sean Kinard
# 2023-06-20

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('03_public/toolkit.R')

lterm <- read_csv('03_public/output/summary_table_lte_table_long.csv') %>% 
  column_to_rownames('Variable') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('site_code') %>%
  as_tibble() %>%
  r_friendly_colnames() %>%
  select(-usgs_staid, -lat, -lon)

#------------------------------------------------------------------------------
# autocorrelation
#------------------------------------------------------------------------------
table_lte_correlation <- lterm %>% 
  correlate() %>%
  shave() %>% 
  fashion()
  
plot_lte_correlation_plot <- lterm %>% 
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
  ggplot(aes(xval)) +
  geom_density(color = 'chartreuse3', fill = NA) +
  geom_density(fill = 'chartreuse3', color = NA, alpha=.2)+
  facet_wrap(~xvar, scales='free') +
  dark_theme_grey(base_size = 14) +
  xlab(element_blank()) +
  ylab(element_blank()) }

plot_lte_histogram_landuse <- lterm %>% 
  select(site_code, contains('land')) %>%
  ex_hist() # no log transforms necessary

plot_lte_histogram_flow <- lterm %>% 
  select(site_code, any_of(flow)) %>%
  ex_hist() # need to log transform many q variables

plot_lte_histogram_waterquality <- lterm %>% 
  select(site_code, any_of(water_quality)) %>%
  ex_hist() # need to log transform conductivity

plot_lte_histogram_geomorph <- lterm %>% 
  select(site_code, any_of(geomorph)) %>%
  ex_hist()

plot_lte_histogram_distance <- lterm %>% 
  select(site_code, any_of(distance)) %>%
  ex_hist()

#------------------------------------------------------------------------------
# lte scatterplots
#------------------------------------------------------------------------------
scatter_plot <- function(my_data) {
  my_data %>%
    pivot_longer(cols = -site_code,
                 names_to = "variable",
                 values_to = "value") %>%
    fix_site_order() %>%
    add_rain() %>%
    ggplot(aes(x=annualrain, y=value)) +
    facet_wrap(~variable, scale = 'free', ncol=2) +
    geom_smooth(method = "loess", se=FALSE, color="grey30", linetype=2,
                size=1.2, span=.8) +
    geom_point(aes(fill=site_code), color = 'black', size = 4, shape=21) +
    scale_fill_manual(values=my_colors) +
    labs(fill='Site') +
    xlab('Rain (cm/yr)') +
    ylab(element_blank()) +
    dark_theme_grey(base_size=14) }

plot_lte_scatter_landuse <- lterm %>%
  select(site_code, contains('land')) %>%
  scatter_plot()

plot_lte_scatter_flow <- lterm %>% 
  select(site_code, any_of(flow)) %>%
  scatter_plot() +
  facet_wrap(~variable, scale = 'free', ncol=3)

plot_lte_scatter_waterquality <- lterm %>% 
  select(site_code, any_of(water_quality)) %>%
  scatter_plot() +
  facet_wrap(~variable, scale = 'free', ncol=3)

plot_lte_scatter_geomorph <- lterm %>% 
  select(site_code, any_of(geomorph)) %>%
  scatter_plot()

plot_plot_lte_scatter_distance <- lterm %>% 
  select(site_code, any_of(distance)) %>%
  scatter_plot()

#------------------------------------------------------------------------------
# Export tables
#------------------------------------------------------------------------------
write_csv(lterm, '03_public/output/environment_long_term_post.csv')
write_csv(table_lte_correlation, '03_public/output/table_lte_correlation.csv')

# #------------------------------------------------------------------------------
# # Export figures
# #------------------------------------------------------------------------------
my_objects <- ls()
my_figure_names <- my_objects[str_detect(my_objects, 'plot_')]

my_figures <- list(
  plot_lte_correlation_plot, plot_lte_histogram_distance, 
  plot_lte_histogram_flow, plot_lte_histogram_geomorph, 
  plot_lte_histogram_landuse, plot_lte_histogram_waterquality, 
  plot_lte_scatter_flow, plot_lte_scatter_geomorph, plot_lte_scatter_landuse,
  plot_lte_scatter_waterquality, plot_plot_lte_scatter_distance)

names(my_figures) <- my_figure_names

for (i in 1:length(my_figures)) {
  my_place <- paste('03_public/visualization/06_', names(my_figures[i]), ".png", sep='')
  my_object <- my_figures[[i]]
  ggsave(my_place,
         plot = my_object,
         width = 10,
         height = 10,
         units = c("in")) }

#------------------------------------------------------------------------------
# End 06_environment_long_term_preview