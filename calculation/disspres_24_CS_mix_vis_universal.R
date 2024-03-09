# 24_CS_mix_vis_universal
# written by: Sean Kinard
# last edit: 2023-09-17
#----------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------
source('toolkit.R') # load packages and helper-function

d <- read_csv('analysis/output/CS_mix_out_universal.csv') %>%
  rename(site_code = site) %>%
  add_sitevars() %>%
  add_is_dam() %>%
  add_rain() 

d_spe <- read_csv('analysis/output/isotope_CNS_2020_01_clean.csv') %>%
  select(order, family, species, guild, transient_type) %>%
  unique()

d_transient_prob <- read_csv('analysis/output/d_transient_prob.csv')

#----------------------------------------------------------------------------
# Bay Distance
#----------------------------------------------------------------------------
d_baydist <- d %>%
  add_sitevars() %>%
  filter(site_type == 'Stream') %>%
  filter(is_dam == 'No Dam') %>%
  add_rain() %>%
  unique() %>%
  filter(dataset == 'transient_type')

plot_bay_dist <- d_baydist  %>%
  ggplot(aes(x=baydist_km, y=mean)) +
  geom_point(aes(shape=m_group, fill=m_group), alpha=.5, size=3, shape=21) +
  geom_point(aes(shape=m_group, color=m_group), 
             size=3, fill=NA, shape=21) +
  stat_poly_line(se=F, method='lm', formula=y~log(x),
                 color='grey85', lwd=.75, lty=2,) +
  stat_poly_eq(label.x=.5, label.y=.95, formula=y~log(x),
               color='white', use_label(c("eq", "adj.R2","p")),
               eq.x.rhs="ln(x)", size = 6) +
  scale_shape_manual('', values=21:24) +
  scale_fill_manual('', values=c('red', 'yellow',  'cyan')) +
  scale_color_manual('', values=c('red', 'yellow', 'cyan')) +
  labs(fill=element_blank()) +
  dark_theme_grey(base_size=20) +
    theme(legend.position = c(0.75, 0.7)) +
    theme(legend.background = element_rect(colour = 'grey50', fill = 'grey5', 
                                           linetype='solid')) +
    theme(legend.title = element_blank(),
          legend.margin = margin(3, 3, 3, 3),
          legend.spacing.x = unit(0, "mm"),
          legend.spacing.y = unit(0, "mm")) +
  labs(y="% Estuarine", x="Estuary Distance (km)")

plot_bay_dist2 <- d_baydist  %>%
  ggplot(aes(x=annualrain, y=mean)) +
  geom_point(aes(shape=m_group, fill=m_group), alpha=.5, size=3) +
  geom_point(aes(shape=m_group, fill=m_group), size=3, fill=NA) +
  stat_poly_line(se=F, method='lm', formula=y~log(x),
                 color='grey25', lwd=.75, lty=2,) +
  stat_poly_eq(label.x=.5, label.y=.95, formula=y~x,
               color='black', use_label(c("eq", "adj.R2","p"))) +
  scale_shape_manual('', values=21:24) +
  scale_fill_manual('', values=c('red', 'yellow', 'blue', 'cyan')) +
  labs(fill=element_blank()) +
  theme_bw(base_size=14) +
  theme(legend.position = c(0.75, 0.7)) +
  theme(legend.background = element_rect(colour = 'grey50', fill = 'white', 
                                         linetype='solid')) +
  theme(legend.title = element_blank(),
        legend.margin = margin(3, 3, 3, 3),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm")) +
  labs(y="% Estuarine", x="Rain (cm/yr)")

#----------------------------------------------------------------------------
# EA vs amphidromous Detection Probability
#----------------------------------------------------------------------------
d_eadp <- d %>%
  add_sitevars() %>%
  filter(site_type == 'Stream') %>% 
  filter(m_group=='Freshwater') %>%
  select(site_code, mean) %>%
  rename(EA=mean) %>%
  left_join(d_transient_prob%>%filter(transient_type=='Amphidromous')%>%
              select(site_code, probability_of_occurance) %>%
              rename(Amp_Det_Pro = probability_of_occurance)) %>%
  na.omit()

lm_eadp <- lm(EA~Amp_Det_Pro, data = d_eadp)
summary(lm_eadp)

lm_stats_eadp = glance(lm_eadp)

plot_eadp <- d_eadp %>%
ggplot(aes(EA, Amp_Det_Pro)) +
  geom_point(fill='cyan', alpha=.5, shape = 21, size = 3, show.legend = T) +
  geom_point(color='cyan', fill=NA, shape = 21, size = 3, show.legend = T) +
  geom_smooth(method = 'lm',
              color='red', se=F, lwd=.8, lty=2, show.legend=F) +
  labs(y= "Freshwater EA (%)",
       x= "Amphidromous Presence (%)") +
  stat_poly_eq(label.x=.5, label.y=.95,
               color='white', use_label(c("adj.R2","p")),
               size=6) +
  dark_theme_grey(base_size=20)

# -----------------------------------------------------------------------------
# Multiple regression
# -----------------------------------------------------------------------------
my_full_formula <- formula(mean~ annualrain * elev_site_m * baydist_km)
full_astep <- lm(my_full_formula, data =  d_baydist)
astep <- stepAIC(full_astep, direction='both')

d_lm <- tibble(
  x_string = c(
    "mean~ annualrain",
    "mean~ annualrain^2",
    "mean~ baydist_km",
    "mean~ annualrain^2",
    "mean~ elev_site_m",
    "mean~ elev_site_m",
    "mean~ annualrain + baydist_km",
    "mean~ annualrain * baydist_km",
    "mean~ annualrain + elev_site_m",
    "mean~ annualrain * elev_site_m",
    "mean~ annualrain + elev_site_m + baydist_km"))

extract_lm_smry <- function(a_string) {
  a_string %>% 
  formula() %>% 
  lm(data=d_baydist) %>% 
  glance() %>%
  select(r.squared, adj.r.squared, p.value, AIC, nobs) }

table_lm_basdist_AIC <- d_lm %>%
  mutate(lm_smry = map(x_string, extract_lm_smry)) %>%
  unnest(lm_smry) %>%
  arrange(AIC)

make_lm_stats_AIC <- function(xdata) {
  lm(formula(x_string), data = xdata) %>%
    tidy() }

best_formulas <- table_lm_basdist_AIC %>%
  slice(1:3) %>% 
  select(x_string, adj.r.squared, p.value, AIC) %>%
  nest_by(formula = x_string) %>%
  mutate(lm_stats = list(lm(as.formula(formula), data = d_baydist) %>% tidy())) %>%
  ungroup() %>%
  unnest(data) %>%
  rename(Eq.adj.r2=adj.r.squared, Eq.p=p.value) %>%
  unnest(lm_stats) %>%
  select(-std.error, -statistic,-x_string) %>%
  filter(term %in% c('annualrain', 'baydist_km', 'elev_site_m', 'annualrain:baydist_km')) %>%
  mutate(term = str_replace_all(term, 'annualrain', 'Rain'),
         term = str_replace_all(term, 'baydist_km', 'Distance'),
         term = str_replace_all(term, 'elev_site_m', 'Elevation'),
         formula = str_replace_all(formula, 'annualrain', 'Rain'),
         formula = str_replace_all(formula, 'baydist_km', 'Distance'),
         formula = str_replace_all(formula, 'elev_site_m', 'Elevation'),
         formula = str_replace_all(formula, 'mean~ ', 'EA ~ ')) %>%
  pivot_wider(names_from=term, values_from = estimate) %>%
  rename(coef.p=p.value) %>%
  select(-coef.p, everything(), coef.p)
  
rank1 <- table_lm_basdist_AIC %>% slice(1) %>% pull(x_string) %>% formula()
best_fit <- lm(rank1, data = d_baydist)
table_best_fit_equation <- best_fit %>% tidy()

plot_best_fit_interaction <- d_baydist %>%
  create_site_group() %>%
  ggplot(aes(x=baydist_km, y = mean, fill=site_group, color = site_group)) +
  geom_point(alpha = .5, shape = 21, size = 3, show.legend = T) +
  geom_point(aes(color=site_group), 
             fill=NA, shape = 21, size = 3, show.legend = T) +
  geom_smooth(method = 'lm',
              aes(lty=site_group), se=F, lwd=.8, lty=2, show.legend=F) +
  geom_text_repel(aes(label=site_code), color='white', fill=NA) +
  scale_x_log10() +
  scale_color_manual('', values=c('yellow', 'green3', 'cyan'))+
  scale_linetype_manual('', values = 1:3) +
  scale_fill_manual('', values=c('yellow', 'green3', 'cyan')) +
  labs(y= "% Estuarine",
       x= "Estuary Distance (km)") +
  dark_theme_grey(base_size=20) +
  theme(legend.position = c(0.75, 0.75))

#----------------------------------------------------------------------------
# predict UN with Regression
#----------------------------------------------------------------------------
table_UN_prediction <- d %>% 
  filter(site_code %in% c('UN', 'LN')) %>%
  filter(dataset == 'transient_type') %>%
  group_by(site_code) %>%
  summarize(annualrain = mean(annualrain),
            baydist_km = mean(baydist_km),
            elev_site_m = mean(elev_site_m),
            prcnt_est = mean(mean)) %>%
  add_predictions(model=best_fit, var='pred') %>%
  mutate(diff = prcnt_est - pred)

write_csv(table_UN_prediction, 'analysis/output/table_UN_prediction.csv')

#----------------------------------------------------------------------------
# End 24_CS_mix_vis