# 18_marine_end_members
# Project: Tracking Marine Nutrients
# stable isotopes: carbon, nitrogen, sulfur
# written by: sean Kinard
# last edit: 2023-07-28
#----------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------
source('toolkit.R') # load packages and helper-functions

d <- read_csv('analysis/output/isotope_CNS_2020_01_clean.csv')

d <- d %>% filter(! species %in% bad_sources)
#----------------------------------------------------------------------------
# Sources
#----------------------------------------------------------------------------
d_sources <- d %>%
  filter(guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  pivot_longer(cols=c(carbon, sulfur), 
               names_to = 'x_name', 
               values_to = 'x_val') %>%
  fix_site_order() %>%
  fix_source_order() %>%
  create_site_type() %>%
  filter(site_type %in% c('Stream', 'Estuary')) %>%
  select(guild, species, site_code, site_type, collection_date, x_name, x_val)

#----------------------------------------------------------------------------
# Boxplot Source Comparisons
#----------------------------------------------------------------------------
boxplot_sources <- function(df) {
  df %>%
    ggplot(aes(site_code, x_val, fill = species) ) +
    geom_boxplot() +
    scale_fill_manual('Type', values=my_colors[c(1,3,5,7,9)])   +
    guides(fill = guide_legend(override.aes=list(shape=22))) +
    labs(shape = NULL) +
    labs(fill = NULL) +
    xlab(element_blank()) +
    ylab(expression(paste(delta)~('\u2030'))) +
    theme_bw(base_size=14) +
    facet_wrap(~x_name, labeller = iso_facet_lab, scales='free', 
               ncol=1) }

plot_stream_source_v_type <- d_sources %>%
  filter(site_type=='Stream') %>%
  boxplot_sources() +
  ggtitle('Streams')

plot_estuary_source_v_type <- d_sources %>%
  filter(site_type=='Estuary') %>%
  boxplot_sources() +
  ggtitle('Estuary')

plot_dam_source_v_type <- d %>%
  filter(guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  pivot_longer(cols=c(carbon, sulfur), 
               names_to = 'x_name', 
               values_to = 'x_val') %>%
  fix_site_order() %>%
  fix_source_order() %>%
  create_site_type() %>%
  filter(site_type=='Dam') %>%
  boxplot_sources() +
  ggtitle('Dam')
  
plot_CS_marine_vs_fresh_sources <- d %>%
  filter(guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  pivot_longer(cols=c(carbon, sulfur, nitrogen), 
               names_to = 'x_name', 
               values_to = 'x_val') %>%
  fix_site_order() %>%
  fix_source_order() %>%
  create_site_type() %>%
  filter(site_type %in% c('Stream', 'Estuary')) %>%
  #filter(! species %in% bad_sources) %>%
  select(guild, species, site_code, site_type, collection_date, x_name, x_val) %>%
  ggplot(aes(site_type, x_val) ) +
  geom_boxplot(color='black', fill='grey80',alpha=.75) +
  stat_compare_means(label.y.npc = "bottom", label.x.npc = "left") +
  xlab(element_blank()) +
  ylab(expression(paste(delta)~('\u2030'))) +
  theme_bw(base_size=14) +
  facet_wrap(~x_name, labeller = iso_facet_lab, scales='free', ncol=3)

plot_CS_marine_vs_fresh_source_types <- d_sources %>%
  filter(x_name != 'nitrogen') %>%
  ggplot(aes(site_type, x_val, fill=species) ) +
  geom_boxplot(color='black', alpha=.75) +
  stat_compare_means(label.y.npc = "bottom", label.x.npc = "left") +
  scale_fill_manual('Type', values=my_colors[c(1,3,5,7,9)])   +
  xlab(element_blank()) +
  ylab(expression(paste(delta)~('\u2030'))) +
  theme_bw(base_size=14) +
  facet_wrap(~x_name, labeller = iso_facet_lab, scales='free', ncol=1)

#----------------------------------------------------------------------------
# Table Source Comparisons
#----------------------------------------------------------------------------
stream_aov <- function(df) {
  aov(permil ~ site_code, data = df) }

stream_tidy <- function(aov_object) {
  tidy(aov_object) %>%
    filter(term=='site_code') }

stream_THSD <- function(aov_object) {
  TukeyHSD(aov_object, conf.level=0.95) %>%
    tidy() %>% 
    filter(adj.p.value<.05) }

stream_boxplot <- function(df) {
  df %>%
    ggplot(aes(x=site_code, y=permil)) +
    geom_boxplot(fill='grey70', color='black') +
    stat_summary(fun.y=mean, geom="point", 
                 shape=23, size=4, fill="white", color="black") +
    stat_compare_means(label.y.npc = "bottom", label.x.npc = "left") +
    xlab(element_blank()) +
    ylab(expression(paste(delta)~('\u2030'))) +
    theme_bw(base_size=14) }

source_v_site <- d_sources %>%
  filter(site_type=='Stream') %>%
  rename(isotope=x_name, permil=x_val) %>%
  group_by(isotope, guild) %>%
  nest() %>%
  mutate(aov_results = map(data, stream_aov)) %>%
  mutate(aov_summary = map(aov_results, stream_tidy)) %>%
  mutate(THSD = map(aov_results, stream_THSD)) %>%
  mutate(plot_box = map(data, stream_boxplot))

plot_source_v_site_aq_c <- source_v_site$plot_box[[1]] + 
  ylab(expression(paste(delta)~C^13~('\u2030'))) +
  geom_hline(yintercept=-19, lty=2, color='blue') +
  geom_hline(yintercept=-24, lty=2, color='blue') +
  geom_hline(yintercept=-20, lty=2, color='red') +
  geom_hline(yintercept=-35, lty=2, color='red') 

plot_source_v_site_aq_s <- source_v_site$plot_box[[2]] +
  ylab(expression(paste(delta)~S^34~('\u2030'))) +
  geom_hline(yintercept=17, lty=2, color='blue') +
  geom_hline(yintercept=21, lty=2, color='blue') +
  geom_hline(yintercept=2, lty=2, color='red') +
  geom_hline(yintercept=6, lty=2, color='red') 

table_thsd_source_v_site <- source_v_site %>% select(guild, isotope, THSD) %>%
  unnest(THSD) %>% select(-null.value, -contains('conf'))

#----------------------------------------------------------------------------
# literature 
#----------------------------------------------------------------------------
# Fry: Sulfur: "continental vegetation seems to average near +2 to +6 permille over large areas and is quite distinct from the +17 to +21 permille values of marine plankton and seaweeds."

# Fry: Carbon: oceanic plankton = -19:-24 permille , freshwater = -20:-45 permille, terrestrial veg-=-28 permlile

### MacAvoy, S. E., S. A. Macko, S. P. McIninch, and G. C. Garman. “Marine Nutrient Contributions to Freshwater Apex Predators.” Oecologia 122, no. 4 (March 1, 2000): 568–73. https://doi.org/10.1007/s004420050980.
## Resident freshwater C,S,N -25.6, 4.8, 12.2
## Anadromos Alosa     C,S,N -18.6, 18.4, 12.9
## Ictalurus furcatus  C,S,N -21.7, 10.6, 16.3

## Hicks, Brendan J., Mark S. Wipfli, Dirk W. Lang, and Maria E. Lang. “Marine-Derived Nitrogen and Carbon in Freshwater-Riparian Food Webs of the Copper River Delta, Southcentral Alaska.” Oecologia 144, no. 4 (August 1, 2005): 558–69. https://doi.org/10.1007/s00442-005-0035-2.
##             No spawners (C, N)        Spawners present
### Fish       -35:-29 , 4:8             -35:-20 , 6:15
### Aq inv     -40:-30 , 1:7             -34:-29,  2:9
### Aq veg     -31:-28 , 2:4             -30:-28   1:2
### Te Veg     -28:-27 , -2:-1           -28:-26   1:2

#----------------------------------------------------------------------------
# Literature comparison
#----------------------------------------------------------------------------
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

# boostrap means and CI
# Bootstrap estimate means and 95% CI
source_boot <- d_sources %>%
  rename(isotope=x_name) %>%
  group_by(site_type, guild, isotope) %>%
  filter(!is.na(x_val)) %>%
  dplyr::summarize(
    mu = estimate_mu(x_val),
    lower = estimate_mu_lower(x_val),
    upper = estimate_mu_upper(x_val) )

site_list <- unique(source_boot$site_type)

plot_sulfur_v_lit <- source_boot %>%
  filter(isotope=='sulfur') %>%
  ggplot(aes(x=site_type, y=mu, fill=guild)) +
  geom_hline(yintercept=17, lty=2, color='blue') +
  geom_hline(yintercept=21, lty=2, color='blue') +
  geom_hline(yintercept=2, lty=2, color='red') +
  geom_hline(yintercept=6, lty=2, color='red') +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, color='black',
                position=position_dodge(width=.5)) +
  geom_point(shape=23, size=5, color='black',
             position=position_dodge(width=.5)) +
  scale_fill_manual(values=c('green', 'brown')) +
  labs(x=element_blank(), y=expression(paste(delta)~S^34~('\u2030'))) +
  theme_bw(base_size=12) +
  theme(legend.title=element_blank()) +
  coord_flip()

plot_carbon_v_lit <- source_boot %>%
  filter(isotope=='carbon') %>%
  ggplot(aes(x=site_type, y=mu, color=guild, fill=guild)) +
  geom_hline(yintercept=-19, lty=2, color='blue') +
  geom_hline(yintercept=-24, lty=2, color='blue') +
  geom_hline(yintercept=-20, lty=2, color='red') +
  geom_hline(yintercept=-45, lty=2, color='red') +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, color='black',
                position=position_dodge(width=.5)) +
  geom_point(shape=23, size=5, color='black',
             position=position_dodge(width=.5)) +
  scale_fill_manual(values=c('green', 'brown')) +
  labs(x=element_blank(), y=expression(paste(delta)~C^13~('\u2030'))) +
  theme_bw(base_size=12) +
  theme(legend.title=element_blank()) +
  coord_flip()

plot_sources_v_lit <- plot_carbon_v_lit + plot_sulfur_v_lit + 
  plot_layout(guides = 'collect')

# Sulfur: 
  # Estuary algae (+18.3 +- 1.3) and overlap expected marine values (17-21).
  # Stream algae (+7.5 +- 4.3) overlap freshwater sulfur signatures (+2 to +6).

# Carbon:
  # Estuary algae (-19.5) is within expected marine carbon13 values (-19 to -24)
  # Stream Terrestrial sources (-27.5+-5.6) are 8 permille depleted in comparison to Estuary aquatic sources (-19.5+-3.8)

#----------------------------------------------------------------------------
# Bootstrap Source Stats
#----------------------------------------------------------------------------
# Source Extant based on bootstrap 95% confidence intervals (Liberal Mean)
s_marine <- table_summary_stats_sources %>% 
  filter(isotope=='sulfur') %>%
  filter(guild == 'Aquatic_Producer') %>%
  filter(site_type == 'Estuary') %>%
  pull(x_max)

s_terrestrial <- table_summary_stats_sources %>% 
  filter(isotope=='sulfur') %>%
  filter(guild == 'Aquatic_Producer') %>%
  filter(site_type == 'Stream') %>%
  pull(x_min)

c_marine <- table_summary_stats_sources %>% 
  filter(isotope=='carbon') %>%
  filter(guild == 'Aquatic_Producer') %>%
  filter(site_type == 'Estuary') %>%
  pull(x_max)

c_terrestrial <- table_summary_stats_sources %>% 
  filter(isotope=='carbon') %>%
  filter(guild == 'Aquatic_Producer') %>%
  filter(site_type == 'Stream') %>%
  pull(x_min)

table_end_members <- tibble(
  type = c('Estuary', 'Terrestrial'),
  sulfur = c(s_marine, s_terrestrial),
  carbon = c(c_marine, c_terrestrial))

aq_source_table <- d %>%
  filter(guild %in% c('Aquatic_Producer')) %>%
  fix_site_order() %>%
  fix_source_order() %>%
  create_site_type() %>%
  filter(site_type %in% c('Stream', 'Estuary')) %>%
  mutate(site_type = str_replace_all(site_type, 'Estuary', 'Estuary')) %>%
  group_by(site_type, guild) %>%
  dplyr::summarize(c_mu = mean(carbon, na.rm=T),
                   c_sd = sd(carbon, na.rm=T),
                   s_mu = mean(sulfur, na.rm=T),
                   s_sd = sd(sulfur, na.rm=T),) %>%
  ungroup()


#----------------------------------------------------------------------------
# Universal Baseline
#----------------------------------------------------------------------------
d_prep <- d %>%
  fix_site_order() %>%
  create_site_type() %>%
  mutate(mixture_type=ifelse(
    guild %in% c('Fish', 'Invertebrate'), 'Animal', 'Source'))

plot_scatter_all <-  d_prep %>%
  ggplot(aes(carbon, sulfur, 
             shape = site_type, fill=site_code) ) +
  facet_wrap(~mixture_type, ncol=1) +
  geom_hline(yintercept=s_marine, lty=2, color='cyan3') +
  geom_hline(yintercept=s_terrestrial, lty=2, color='red3') +
  geom_vline(xintercept=c_marine, lty=2, color='cyan3') +
  geom_vline(xintercept=c_terrestrial, lty=2, color='red3') +
  geom_rect(data=aq_source_table%>% filter(site_type=='Estuary'),
                aes(x = c_mu, y = s_mu, shape = NA,
                  xmin= c_mu-c_sd,
                    xmax=c_mu+c_sd,
                    ymin=s_mu-s_sd,
                    ymax=s_mu+s_sd), fill='cyan', alpha=.3) +
  geom_rect(data=aq_source_table%>% filter(site_type=='Stream'),
            aes(x = c_mu, y = s_mu, shape = NA,
                xmin= c_mu-c_sd,
                xmax=c_mu+c_sd,
                ymin=s_mu-s_sd,
                ymax=s_mu+s_sd), fill='red', alpha=.3) +
  geom_point(size=3, alpha=.75) +
  scale_shape_manual('Site Type', values=c(22, 21, 24)) +
  scale_fill_manual('Site', values=my_sites_colors) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab(expression(paste(delta)~C^13~('\u2030'))) +
  ylab(expression(paste(delta)~S^34~('\u2030'))) +
  theme_bw(base_size=14) +
  geom_text(aes(x=-15,y=s_terrestrial+.7, label='Fresh Min'), size = 4, color='red3') +
  geom_text(aes(x=-15,y=s_marine+.7, label='Estuary Max'), size = 4, color='cyan3')

c_aov <- aov(carbon~site_code, 
    data = d_prep %>% filter(site_type == 'Stream' &
                             mixture_type=='Source'))
c_aov %>% tidy()

s_aov <- aov(sulfur~site_code, 
    data = d_prep %>% filter(site_type == 'Stream' &
                               mixture_type=='Source'))
s_aov %>% tidy()

s_thsd <- TukeyHSD(s_aov, conf.level=.95)
plot(s_thsd, las = 2)

# sulfure signatures within sources differ from regional averages at TR and GC
#----------------------------------------------------------------------------
# Site-Specific Baseline
#----------------------------------------------------------------------------

scatterplot_site <- function(df) {
  
  output <- df %>%
    ggplot(aes(carbon, sulfur, 
               shape = mixture_type) ) +
    geom_hline(yintercept=s_marine, lty=2, lwd=.5,alpha=.5, color='cyan3') +
    geom_hline(yintercept=s_terrestrial, lty=2,lwd=.5,alpha=.5, color='red3') +
    geom_vline(xintercept=c_marine, lty=2,lwd=.5,alpha=.5, color='cyan3') +
    geom_vline(xintercept=c_terrestrial, lty=2,lwd=.5,alpha=.5, color='red3') +
    geom_rect(data = df%>% slice(1),
              aes(x = estuary_source_c_mu, 
                  y = estuary_source_s_mu, 
                  shape = NA,
                  xmin= estuary_source_c_mu-estuary_source_c_sd,
                  xmax= estuary_source_c_mu+estuary_source_c_sd,
                  ymin= estuary_source_s_mu-estuary_source_s_sd,
                  ymax= estuary_source_s_mu+estuary_source_s_sd), 
              fill='cyan', alpha=.3) +
    geom_rect(data = df%>% slice(1),
              aes(x = stream_source_c_mu, 
                  y = stream_source_s_mu, 
                  shape = NA,
                  xmin= stream_source_c_mu-stream_source_c_sd,
                  xmax= stream_source_c_mu+stream_source_c_sd,
                  ymin= stream_source_s_mu-stream_source_s_sd,
                  ymax= stream_source_s_mu+stream_source_s_sd), 
              fill='red', alpha=.3) +
    geom_text(aes(x=-30,y=s_terrestrial-.8, label='Fresh'), size = 4, color='red3') +
    geom_text(aes(x=-15,y=s_marine+.8, label='Estuary'), size = 4, color='cyan3') +
    geom_point(shape=21, size=2) +
    xlim(c(-35,-10)) +
    ylim(c(-5,25)) +
    xlab(expression(paste(delta)~C^13~('\u2030'))) +
    ylab(expression(paste(delta)~S^34~('\u2030'))) +
    theme_bw(base_size=14) +
    geom_text(aes(x=-22,y=24, label=site_code2[1]), size = 5)
  
  return(output) }

estuary_source_stats <- d %>%
  filter(guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  add_sitevars() %>%
  filter(site_type %in% c('Estuary')) %>%
  group_by(site_code) %>%
  dplyr::summarize(estuary_source_c_mu=mean(carbon,na.rm=T),
                   estuary_source_c_sd=sd(carbon, na.rm=T),
                   estuary_source_s_mu=mean(sulfur,na.rm=T),
                   estuary_source_s_sd=sd(sulfur, na.rm=T)) %>%
  ungroup()

stream_source_stats <- d %>%
  filter(guild %in% c('Aquatic_Producer', 'Terrestrial_Producer')) %>%
  add_sitevars() %>%
  filter(site_type %in% c('Stream')) %>%
  group_by(site_code) %>%
  dplyr::summarize(stream_source_c_mu=mean(carbon,na.rm=T),
                   stream_source_c_sd=sd(carbon, na.rm=T),
                   stream_source_s_mu=mean(sulfur,na.rm=T),
                   stream_source_s_sd=sd(sulfur, na.rm=T)) %>%
  ungroup()

my_scatterplots <- d %>%
  add_sitevars() %>%
  filter(site_type %in% c('Stream')) %>%
  left_join(estuary_source_stats %>%
              rename(nearest_bay = site_code)) %>%
  left_join(stream_source_stats) %>%
  fix_site_order() %>%
  arrange(site_code) %>%
  mutate(site_code2 = site_code) %>%
  group_by(site_code) %>%
  mutate(mixture_type=ifelse(
    guild %in% c('Fish', 'Invertebrate'), 'Animal', 'Source')) %>%
  filter(mixture_type=='Animal') %>%
  nest() %>%
  mutate(my_scatter = map(data, scatterplot_site))

remove_x <- function(my_plot) {
  my_plot + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank()) }

remove_y <- function(my_plot) {
  my_plot + theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank()) }

remove_xy <- function(my_plot) {
  my_plot + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank()) }

plot_sc_TR <- my_scatterplots$my_scatter[[1]] %>% remove_x() + ylab(element_blank())
plot_sc_SF <- my_scatterplots$my_scatter[[2]] %>% remove_xy()
plot_sc_AR <- my_scatterplots$my_scatter[[3]] %>% remove_xy()
plot_sc_MR <- my_scatterplots$my_scatter[[4]] %>% remove_x() 
plot_sc_PD <- my_scatterplots$my_scatter[[5]] %>% remove_xy()
plot_sc_PL <- my_scatterplots$my_scatter[[6]] %>% remove_xy()
plot_sc_GC <- my_scatterplots$my_scatter[[7]] + xlab(element_blank()) + ylab(element_blank())
plot_sc_WM <- my_scatterplots$my_scatter[[8]] %>% remove_y()
plot_sc_EM <- my_scatterplots$my_scatter[[9]] %>% remove_y() + xlab(element_blank())

plot_scatter_sites <- (plot_sc_TR + plot_sc_SF + plot_sc_AR) /
  (plot_sc_MR + plot_sc_PD +plot_sc_PL) /
  (plot_sc_GC + plot_sc_WM +plot_sc_EM) +
  plot_layout(guides='collect')

plot_scatter_dam <- my_scatterplots$my_scatter[[10]] +
my_scatterplots$my_scatter[[11]]%>% remove_y()

# Local sulfur values are enriched at TR and PD, sites that qualify as potential marsh or wetland habitat, which may be attributed to local redox by soil bacteria.

# Otherwise, with few exceptions, local mixtures are highly clustered near local stable isotope signatures. Estuarine influence or atmospheric deposition may influence local source signatures.

#----------------------------------------------------------------------------
# End 18_marine_end_members