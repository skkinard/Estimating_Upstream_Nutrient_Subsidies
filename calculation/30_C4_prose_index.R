# 30_C4_prose_index
# written by: Sean Kinard
# last edit: 2023-09-17
#----------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------
source('03_public/calculation/14_transient_occurance.R')
source('03_public/calculation/24_CS_mix_vis_region.R')
source('03_public/calculation/24_CS_mix_vis_universal.R')

#----------------------------------------------------------------------------
#### Sample Numbers (SN)
d_sn <- d_iso_raw %>% add_sitevars()

# - `r SN_tot` total samples analyzed for carbon-13, sulfur-34, and nitrogen-15.
SN_tot <- length(d_sn$carbon)

# - Nitrogen values were excluded from analyses because source signatures could not be distinguished between estuaries and streams (ANOVA *p* = `r SN_aov_n_p`).
aov_nitrogen <- aov(nitrogen ~ site_type*lowest_taxon, 
                    data = d_sn %>%
                    filter(guild == 'Aquatic_Producer'))

aov_nitrogen %>% tidy() %>% filter(term=='site_type') %>% pull(p.value)
SN_aov_n_p <- aov_nitrogen %>% 
  tidy() %>% 
  filter(term=='site_type:lowest_taxon') %>% 
  pull(p.value) %>%
  round(2)

# - `r SN_s_est` primary producer samples were collected from five estuaries, including filamentous algae, detritus, and aquatic macrophytes.
SN_s_est <- d_sn %>%
  filter(guild=='Aquatic_Producer' & site_type != 'Stream') %>% 
  pull(carbon) %>% length()

# - These basal resources plus riparian leaves and periphyton were collected from twelve stream locations, totalling `r SN_s_str` samples.
SN_s_str <- d_sn %>%
  filter(guild=='Aquatic_Producer' & site_type == 'Stream') %>% 
  pull(carbon) %>% length()

# - From streams, a total of `r SN_f` fish and `r SN_i` invertebrate samples were collected and processed.
SN_f <- d_sn %>%
  filter(guild=='Fish') %>% 
  pull(carbon) %>% length()

SN_i <- d_sn %>%
  filter(guild=='Invertebrate') %>% 
  pull(carbon) %>% length()

# - `r SN_f_amp` amphidromous fish were analyzed, including spp....
SN_f_amp <- d_sn %>%
  filter(guild=='Fish' & transient_type=='Amphidromous') %>% 
  pull(carbon) %>% length()

# - `r SN_i_amp` amphidromous invertebrates were analyzed, including spp....
SN_i_amp <- d_sn %>%
  filter(guild=='Invertebrate' & transient_type=='Amphidromous') %>% 
  pull(carbon) %>% length()

# - `r SN_f_cat` catadromous fish were analyzed, including spp....
SN_f_cat <- d_sn %>%
  filter(guild=='Fish' & transient_type=='Catadromous') %>% 
  pull(carbon) %>% length()

# - `r SN_i_cat` catadromous invertebrates were analyzed, including spp....
SN_i_cat <- d_sn %>%
  filter(guild=='Invertebrate' & transient_type=='Catadromous') %>% 
  pull(carbon) %>% length()

# - `r SN_f_ana` Anadromous fish were analyzed, including spp....
SN_f_ana <- d_sn %>%
  filter(guild=='Fish' & transient_type=='Anadromous') %>% 
  pull(carbon) %>% length()

# - `r SN_i_eur` Anadromous invertebrates were analyzed, including spp....
SN_i_ana <- d_sn %>%
  filter(guild=='Invertebrate' & transient_type=='Anadromous') %>% 
  pull(carbon) %>% length()

# - `r SN_f_eur` Euryhaline fish were analyzed, including spp....
SN_f_eur <- d_sn %>%
  filter(guild=='Fish' & transient_type=='Euryhaline') %>% 
  pull(carbon) %>% length()

# - `r SN_i_ana` Euryhaline invertebrates were analyzed, including spp....
SN_i_eur <- d_sn %>%
  filter(guild=='Invertebrate' & transient_type=='Euryhaline') %>% 
  pull(carbon) %>% length()

# - `r SN_f_fresh` freshwater fish were analyzed, including spp...
SN_f_fresh <- d_sn %>%
  filter(guild=='Fish' & transient_type=='Freshwater') %>% 
  pull(carbon) %>% length()

# - `r SN_i_fresh` freshwater invertebrates were analyzed, including spp...
SN_i_fresh <- d_sn %>%
  filter(guild=='Invertebrate' & transient_type=='Freshwater') %>% 
  pull(carbon) %>% length()

#----------------------------------------------------------------------------
### Source Signatures (SS)
d_SS <- d_iso%>%add_sitevars()

# - Estuarine source $\delta$C^13^ and $\delta$S^34^ values average `r SS_s_est_c13_mu`, and `r SS_s_est_s34_mu` respectively which are within the expected ranges from scientific literature (-19 to -24 $\delta$C^13^ and 17 to 21 $\delta$S^34^) (CITATION Fry 2001, Hicks et al. 2005, MacAvoy et al. 2000).
SS_s_est_c13_mu <- d_SS %>% 
  filter(site_type != 'Stream' &
           lowest_taxon %in% my_sources) %>%
  pull(carbon) %>% mean() %>% round(1)

SS_s_est_s34_mu <- d_SS %>% 
  filter(site_type != 'Stream' &
        lowest_taxon %in% my_sources) %>%
  pull(sulfur) %>% mean() %>% round(1)

# - Stream source $\delta$C^13^ and $\delta$S^34^ values `r SS_s_str_c13_mu`, and `r SS_s_str_s34_mu` also conformed to expected ranges from scientific literature (-19 to -24 $\delta$C^13^, and 2 to 6 $\delta$S^34^) (CITATION Fry 2001, Hicks et al. 2005, MacAvoy et al. 2000).

SS_s_str_c13_mu <- d_SS %>% 
  filter(site_type == 'Stream' &
           lowest_taxon %in% my_sources) %>%
  pull(carbon) %>% mean() %>% round(1)

SS_s_str_s34_mu <- d_SS %>% 
  filter(site_type == 'Stream' &
        lowest_taxon %in% my_sources) %>%
  pull(sulfur) %>% mean() %>% round(1)

# - `r animals_inrange` of `r animals_total` animal signatures fell within the range of source values (Figure: Transient $\delta$S^34^ Vs $\delta$C^13^ Scatterplot).
source_min_c13 <-  d_SS %>% 
  filter(lowest_taxon %in% my_sources) %>% pull(carbon) %>% min()
source_min_s34 <- d_SS %>% 
  filter(lowest_taxon %in% my_sources) %>% pull(sulfur) %>% min()
source_max_c13 <-  d_SS %>% 
  filter(lowest_taxon %in% my_sources) %>% pull(carbon) %>% max()
source_max_s34 <- d_SS %>% 
  filter(lowest_taxon %in% my_sources) %>% pull(sulfur) %>% max()

animals_total <- d_SS %>%
  filter(guild %in% c('Fish', 'Invertebrate')) %>%
  pull(carbon) %>% length()

animals_inrange <- d_SS %>%
  filter(guild %in% c('Fish', 'Invertebrate')) %>%
  filter(carbon > source_min_c13 &
           carbon < source_max_c13 &
           sulfur > source_min_s34 &
           sulfur < source_max_s34) %>%
  pull(carbon) %>% length()

animals_outrange <- animals_total - animals_inrange

# - Despite different community compositions, *P. pugio* below the dam averaged +`r dam_ppugio_c13_diff` $\delta$C^13^ and +`r dam_ppugio_s34_diff` $\delta$S^34^ compared to *P. pugio* above the dam. Similarly, *P. pugio* below the dam averaged +`r LN_minus_stream_ppugio_c13` $\delta$C^13^ and +`r dam_ppugio_s34_diff` $\delta$S^34^ compared to *P. pugio* collected from other stream locations.
d_SS_dam <- d_SS %>%
  filter(site_code %in% c('UN', 'LN') &
           guild %in% c('Invertebrate', 'Fish')) %>%
  select(site_code, lowest_taxon) %>%
  unique()

taxa_UN <- d_SS_dam %>% filter(site_code == 'UN') %>% pull(lowest_taxon)
taxa_LN <- d_SS_dam %>% filter(site_code == 'LN') %>% pull(lowest_taxon)
taxa_stream <- d_SS %>%
  filter(site_code %in% my_streams &
           guild %in% c('Invertebrate', 'Fish')) %>%
  pull(lowest_taxon) %>%
  unique()
taxa_LN_and_UN <- d_SS_dam %>%
  filter(site_code == 'LN') %>%
  filter(lowest_taxon %in% taxa_UN) %>%
  pull(lowest_taxon)
taxa_LN_and_stream <- d_SS_dam %>%
  filter(site_code == 'LN') %>%
  filter(lowest_taxon %in% taxa_stream) %>%
  pull(lowest_taxon)

dam_ppugio <- d_SS %>%
  filter(site_code %in% c('UN', 'LN') &
           lowest_taxon %in% taxa_LN_and_UN) %>%
  group_by(site_code) %>%
  dplyr::summarize(c13_mu = mean(carbon),
                   s34_mu = mean(sulfur))

LN_minus_stream <- function(xdata) {
  diff_c13 <- xdata%>%filter(is_stream=='LN')%>%pull(c13_mu) - 
    xdata%>%filter(is_stream=='Stream')%>%pull(c13_mu)
  diff_s34 <- xdata%>%filter(is_stream=='LN')%>%pull(s34_mu) - 
    xdata%>%filter(is_stream=='Stream')%>%pull(s34_mu)
  x_output <- tibble(LN_minus_stream_c13=diff_c13,
                     LN_minus_stream_s34=diff_s34)
  return(x_output) }

dam_stream_diff <- d_SS %>%
  filter(site_code %in% c(my_streams, 'LN')) %>% 
  filter(lowest_taxon %in% taxa_LN_and_stream) %>%
  mutate(is_stream = ifelse(site_code %in% my_streams, 'Stream', 'LN')) %>%
  group_by(is_stream, lowest_taxon) %>%
  dplyr::summarize(c13_mu = mean(carbon),
                   s34_mu = mean(sulfur)) %>%
  ungroup() %>%
  group_by(lowest_taxon) %>%
  nest() %>%
  mutate(LN_stream_diff = map(data, LN_minus_stream)) %>%
  select(-data) %>%
  unnest(LN_stream_diff)

dam_ppugio_c13_diff <- 
  (dam_ppugio%>%filter(site_code=='LN')%>%pull(c13_mu) - 
     dam_ppugio%>%filter(site_code=='UN')%>%pull(c13_mu)) %>% round(1)
dam_ppugio_s34_diff <- 
  (dam_ppugio%>%filter(site_code=='LN')%>%pull(s34_mu) - 
     dam_ppugio%>%filter(site_code=='UN')%>%pull(s34_mu)) %>% round(1)

LN_minus_stream_ppugio_c13 <- dam_stream_diff%>%
  filter(lowest_taxon=='P.pugio') %>%
  pull(LN_minus_stream_c13) %>% round(1)
LN_minus_stream_ppugio_s34 <- dam_stream_diff%>%
  filter(lowest_taxon=='P.pugio') %>%
  pull(LN_minus_stream_s34) %>% round(1)

#----------------------------------------------------------------------------
### Mixing Models (MM)
# - Mixing models estimate high estuarine assimilation in sites with estuary-distances less than 25 km (`r MM_short_mu`%) which contrast minor estuarine contributions (`r MM_long_mu`%) to more distant stream communities.
MM_short_mu <- d_vis_mix_CI %>%
  filter(dataset=='distcat') %>%
  unnest(data) %>%
  filter(m_group == '<25 km') %>%
  pull(mean) %>%
  round(0)

MM_long_mu <- d_vis_mix_CI %>%
  filter(dataset=='distcat') %>%
  unnest(data) %>%
  filter(m_group == '>25 km') %>%
  pull(mean) %>%
  round(0)

# - Diadromous taxa source primarily estuaries (`r MM_amp_mu`), which contrasts low estuarine assimilation in freshwwater (`r MM_fre_mu`) and catadromous (`r MM_cat_mu`) samples.
MM_dia_mu <- d_vis_mix_CI %>%
  filter(dataset=='is_diadromous') %>%
  unnest(data) %>%
  filter(m_group == 'Diadromous') %>%
  pull(mean) %>%
  round(0)

MM_eur_mu <- d_vis_mix_CI %>%
  filter(dataset=='is_diadromous') %>%
  unnest(data) %>%
  filter(m_group == 'Euryhaline') %>%
  pull(mean) %>%
  round(0)

MM_fre_mu <- d_vis_mix_CI %>%
  filter(dataset=='is_diadromous') %>%
  unnest(data) %>%
  filter(m_group == 'Freshwater') %>%
  pull(mean) %>%
  round(0)

# - The stream community below Calallen Dam assimilated mostly estuarine sources (`r MM_LN_mu`%) which differed from values above Calallen Dam (`r MM_UN_mu`%), above other damn (`r MM_otherdam_mu`%), and sites without dams (`r MM_nodam_mu`%).
MM_LN_mu <- d_vis_mix_CI %>%
  filter(dataset=='calallen') %>%
  unnest(data) %>%
  filter(m_group == 'Below C.') %>%
  pull(mean) %>%
  round(0)

MM_UN_mu <- d_vis_mix_CI %>%
  filter(dataset=='calallen') %>%
  unnest(data) %>%
  filter(m_group == 'Above C.') %>%
  pull(mean) %>%
  round(0)

MM_otherdam_mu <- d_vis_mix_CI %>%
  filter(dataset=='calallen') %>%
  unnest(data) %>%
  filter(m_group == 'Above Other') %>%
  pull(mean) %>%
  round(0)

MM_nodam_mu <- d_vis_mix_CI %>%
  filter(dataset=='calallen') %>%
  unnest(data) %>%
  filter(m_group == 'None') %>%
  pull(mean) %>%
  round(0)

#----------------------------------------------------------------------------
### Transient Distribution (TD)

# - Amphidromous species detection probability is negatively related to annual rainfall (R^2^=`r TD_amp_r2`, *p*=`r TD_amp_p`), and catadromous species are only detected in sites with annual rainfall greater than 72 cm.

TD_amp_lm <- d_transient_prob %>% 
  filter(transient_type=='Amphidromous') %>%
  select(r.squared, p.value) %>%
  unique()

TD_amp_r2 <- TD_amp_lm$r.squared%>%round(2)
TD_amp_p <- TD_amp_lm$p.value%>%round(3)

# - Estuary-distance is not a significant predictor of transient detection.

# EA in Freshwater residents is positively related to amphidromous probability detection
freshEA_vs_Amp_r2 <- lm_stats_eadp$r.squared%>%round(2)
freshEA_vs_Amp_p <- lm_stats_eadp$p.value%>%round(3)

#----------------------------------------------------------------------------
### Linear Models (LM)
# - In the most parsimonious model (adjusted R^2^=`r LM_r2`, *p*=`r LM_p`), estuarine assimilation is negatively related to estuary-distance (`r LM_baydist` %/km) and annual rainfall (`r LM_annualrain` %/cm).
LM_r2 <- table_lm_basdist_AIC %>%
  filter(x_string=='mean~ annualrain + baydist_km') %>%
  pull(adj.r.squared) %>% round(2)

LM_p <- table_lm_basdist_AIC %>%
  filter(x_string=='mean~ annualrain + baydist_km') %>%
  pull(p.value) %>% formatC(format = "e", digits = 0)

LM_baydist <- table_best_fit_equation %>%
  filter(term=='baydist_km') %>%
  pull(estimate) %>% round(2)

LM_annualrain <- table_best_fit_equation %>%
  filter(term=='annualrain') %>%
  pull(estimate) %>% round(2)

# - Using a single predictor, estuarine assimilation relates negatively to log-transformed estuary-distance (R^2^=0.56, *p*<.001).
# - The full model ranked #1 and contained an additional variable (elevation), but the associated *p*-value was not significant.

#----------------------------------------------------------------------------
## Calallen Dam Predictions
# - Estuarine assimilation above Calallen dam is `r predict_diff_UN`% less than the predicted value using annual rainfall and estuary-distance.
predict_diff_UN <- table_UN_prediction %>%
  filter(site_code == 'UN') %>%
  pull(diff) %>% round(0)

# - Predicted values matched observations below the dam (+`r predict_diff_LN`% difference).
predict_diff_LN <- table_UN_prediction %>%
  filter(site_code == 'LN') %>%
  pull(diff) %>% round(0)
