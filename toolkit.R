# diadromy_tools
# Sean Kinard
# 2023-07-17

#------------------------------------------------------------------------------
# load packages
#------------------------------------------------------------------------------
library(tidyverse) # dplyr and ggplot
library(lubridate) # handling dates
library(ggrepel) # plot labels
library(patchwork) # combo plots
library(corrr) # autocorrelation
library(tsibble) # time series
library(imputeTS) # time series imputation
library(broom.mixed) # tidy mixed models
library(ggpubr)
library(ggpmisc) # annotate plot models
library(boot)
library(MASS)
library(interactions)
library(tidymodels)
library(modelr)
library(conflicted)
conflict_prefer_all('dplyr', quiet=T)

#------------------------------------------------------------------------------
# vectors
#------------------------------------------------------------------------------

my_sites <- c("TR", "SF", "AR", "MR", "PD", "PL", "GC", "WM", "EM",
              "BB", "NB", "CB","HB", "LB",
              "LN", "UN")
my_streams <- c("TR", "SF", "AR", "MR", "PD", "PL", "GC", "WM", "EM")
my_bays <- c("BB", "NB", "CB","HB", "LB")
my_dams <- c("LN", "UN")

my_colors <- c('#F5191CFF', '#A54E21FF', '#EC7014FF',
               '#FFC72CFF', '#FFF7BCFF', '#66C2A4FF',
               '#006D2CFF', '#6BAED6FF', '#08519CFF')

my_sites_colors <- c('#FFFFA8FF', '#FFFF57FF', '#FFFF05FF',
                     '#FFCB00FF', '#FFB000FF', '#FF9400FF',
                     '#FF7900FF', '#FF5D00FF', '#FF2600FF',
                     '#00C07AFF', "#009F94FF", "#007696FF", 
                     "#2B4681FF", "#4A075CFF", "grey20", "grey80")

my_stream_colors <- my_sites_colors[1:9]
my_bay_colors <- my_sites_colors[10:14]
my_dam_colors <- my_sites_colors[15:16]

# group predictor colnames for easier indexing
location <- c('lat', 'lon')
climate <- c('annualtemp', 'annualrain')
landuse <- c('developedland', 'forestland', 'cropland', 'otherland')
lt_flow <- colnames(read_csv('03_local_files/data/environment/environment_longterm.csv') %>%
                      dplyr::select(contains('q_')))
st_flow <- colnames(read_csv("03_local_files/data/environment/site_flow_2week_stats.csv") %>%
                      dplyr::select(contains('q2wk')))
flow <- c(lt_flow, st_flow)
water_quality <- c('ammonia', 'conductivity', 'd_oxygen',
                   'doc', 'nitrate', 'ph', 'phosphate',
                   'temperature', 'turbidity')
geomorph <- c('gravel', 'sand', 'silt',
              'canopy', 'depth_mx', 'width')
algae <- c('bluegreen_cyano', 'green_algae', 'diatoms')
distance <- c('elevation_m', 'dam_distance_km', 'lake_density_sqkm', 
              'bay_distance_km')	

my_sources<-c("Filamentous Algae", "Periphyton", "Macrophyte", "Coarse Detritus", 
              "Fine Detritus", "Clover", "Grass", "Mixed TV", "Tree Leaf")

# Isotope facet labels with special characters
iso_facet_lab <- as_labeller(c(
  C="delta~C^13", S="delta~S^34", N="delta~N^15",
  Carbon="delta~C^13", Sulfur="delta~S^34", Nitrogen="delta~N^15", 
  carbon="delta~C^13", sulfur="delta~S^34", nitrogen="delta~N^15",
  "SF"="San~Fernando", "TR"="Tranquitas", "AR"="Aransas", "PD"="Perdido",
  "MR"="Mission", "GC"="Garcitas", "PL"="Placedo", "WM"="W.~Mustang",
  "EM"="E.~Mustang", "UN"="U.~Nueces", "LN"="L.~Nueces", "BB"="Baffin~Bay",
  "NB"="Nueces~Bay", "CB"="Coleto~Bay", "HB"="Hynes~Bay", "LB"="Lavaca~Bay",
  "Stream"="Stream", "Bay"="Bay", "Dam"="Dam"),
  default = label_parsed)

# bad sources (used to filter by variable: "species")
bad_sources <- c("Fine Detritus", "Clover", "Grass")

sites_with_dam <- c('PD', 'WM', 'EM', 'UN')

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
fix_site_order <- function(x) { 
  mutate(x, site_code = fct_relevel(site_code, my_sites)) }

fix_source_order <- function(x) { 
  mutate(x, species = fct_relevel(species, my_sources)) }

create_site_type <- function(x) {
  mutate(x, site_type = case_when(
    site_code %in% my_streams ~ "Stream",
    site_code %in% my_bays ~ "Estuary",
    site_code %in% my_dams ~ "Dam")) }

add_rain <- function(x) { 
  read_csv('03_local_files/data/environment/environment_longterm.csv') %>%
    dplyr::select(site_code, annualrain) %>%
    full_join(tibble(
      site_code = c('UN', 'LN'),
      annualrain = c(62.6,62.6) )) %>%
    right_join(x) }

add_sitevars <- function(df) {
  read_csv('03_local_files/data/isotope_raw/MDN_clean_final.csv') %>%
    r_friendly_colnames() %>%
    select(site, site_type, pptavg_site, nearest_bay, baydist_km, hires_lentic_dens,
           npdes_maj_dens, raw_dis_nearest_maj_npdes, raw_dis_nearest_dam,
           elev_site_m) %>%
    rename(site_code = site) %>%
    mutate(site_type = str_replace_all(site_type, 'Fresh', 'Stream')) %>%
    unique() %>%
    pivot_longer(cols=-any_of(c('site_code', 'site_type', 'nearest_bay')),
                 names_to='xname', values_to='xvalue') %>%
    mutate(xvalue = ifelse(xvalue==-999, NA,xvalue)) %>%
    pivot_wider(names_from=xname, values_from=xvalue) %>%
    right_join(df) }

add_is_dam <- function(x) {x %>%
    mutate(is_dam=ifelse(site_code%in%sites_with_dam, 'Dam', 'No Dam'))}

r_friendly_colnames <- function(d) {
  output <- d
  colnames(output) <- str_to_lower(colnames(output))
  colnames(output) <- str_replace_all(colnames(output), ' ', '_')
  colnames(output) <- str_replace_all(colnames(output), '\\.', '_') 
  return(output) }

pretty_x <- function(x) {
  x <- str_replace(x, '_', ' ')
  x <- str_to_title(x)
  return(x) }

pretty_titles <- function(df) {
  x <- df
  colnames(x) <- str_replace_all(colnames(x), '_', ' ')
  colnames(x) <- str_replace_all(colnames(x), 'land', ' land')
  colnames(x) <- str_replace_all(colnames(x), 'annual', 'annual ')
  colnames(x) <- str_replace_all(colnames(x), 'q q', 'q ')
  colnames(x) <- str_to_title(colnames(x))
  colnames(x) <- str_replace_all(colnames(x), 'Staid', 'USGS STAID')
  colnames(x) <- str_replace_all(colnames(x), 'Doc', 'DOC')
  colnames(x) <- str_replace_all(colnames(x), 'Ph', 'pH')
  colnames(x) <- str_replace_all(colnames(x), 'pHo', 'Pho')
  return(x) }

grid_date <- function(x) {
  x %>% group_by(site_code, collection_period) %>%
    dplyr:: summarize(n=length(site_code)) %>%
    arrange(collection_period) %>%
    pivot_wider(names_from = site_code,
                values_from = n) }

restore_category <- function(x) { 
  case_when(
    x %in% ln_vars ~ paste('ln_', x, sep=''),
    x %in% sqrt_vars ~ paste('sqrt_', x, sep=''),
    TRUE ~ x) }

create_site_group <- function(df) {
  df %>% mutate(site_group = case_when(
    site_code %in% c('TR', 'SF', 'AR', 'LN', 'UN') ~ 'Semi-Arid',
    site_code %in% c('MR', 'PD', 'PL') ~ 'Transition',
    site_code %in% c('GC', 'WM', 'EM') ~ 'Sub-Humid'))%>%
    mutate(site_group = fct_relevel(site_group, c('Semi-Arid',
                                                  'Transition',
                                                  'Sub-Humid')))}
create_iso_period <- function(x) {
  jan_start <- ymd('2019-12-01')
  jan_end <- ymd('2020-02-01')
  june_start <- ymd('2020-05-01')
  june_end <- ymd('2020-07-01')
  sep_start <- ymd('2020-08-01')
  sep_end <- ymd('2020-10-01')
  nov_start <- ymd('2020-10-01')
  nov_end <- ymd('2020-12-01')
  
  x %>%
    mutate(iso_period = case_when(
      collection_period > jan_start & collection_period < jan_end ~ 'January',
      collection_period > june_start & collection_period < june_end ~ 'June',
      collection_period > sep_start & collection_period < sep_end ~ 'September',
      collection_period > nov_start & collection_period < nov_end ~ 'November',
      TRUE ~ NA)) %>%
    filter(! is.na(iso_period)) }

# fill missing environmental predictors with linear interpolation
impute_interpolation <- function(df) {
  df %>%
    rowwise() %>%
    as_tsibble(key=site_code,
               index=collection_period) %>%  
    na_ma(weighting="linear") %>%
    ungroup() %>%
    as_tibble() %>%
    pivot_longer(cols=-any_of(c('site_code', 'collection_period')),
                 names_to='xname',
                 values_to = 'xvalue') %>%
    group_by(site_code, collection_period, xname) %>%
    dplyr::summarize(xvalue=mean(xvalue ,na.rm=T)) %>%
    pivot_wider(names_from=xname, values_from = xvalue) %>%
    ungroup() }


create_site_period <- function(d) {
  my_d <- d %>%
    mutate(my_date = Collection_Date) %>%
    separate(my_date, c("Year", "Month", "Day"), sep = '-') %>%
    mutate(Day = as.numeric(Day),
           Month = as.numeric(Month)) %>%
    mutate(Collection_Month = case_when(
      Site_Code == 'AR' & Collection_Date == ymd("2019-03-21") ~ as.character(Month),
      Day <= 20 ~ as.character(Month),
      TRUE ~ as.character(Month + 1))) %>%
    mutate(Collection_Period = ym(paste(Year, as.character(Collection_Month), 
                                        sep = '-'))) %>%
    mutate(site_period = paste(Site_Code, Collection_Period, sep = '_')) %>%
    dplyr::select( - Year, -Month, -Day, -Collection_Month) 
  
  return(my_d)  }

create_terrg_period <- function(d) {
  my_d <- d %>%
    mutate(my_date = Collection_Date) %>%
    separate(my_date, c("Year", "Month", "Day"), sep = '-') %>%
    mutate(Day = as.numeric(Day),
           Month = as.numeric(Month)) %>%
    mutate(Collection_Month = ifelse(Day < 20,
                                     Month, as.character(Month + 1))) %>%
    mutate(Collection_Period = ym(paste(Year, as.character(Collection_Month), 
                                        sep = '-'))) %>%
    mutate(site_period = paste(Site_Code, Collection_Period, sep = '_')) %>%
    dplyr::select( - Year, -Month, -Day, -Collection_Month) 
  
  return(my_d) }

create_combo_period <- function(df) {
  pre_2020 <- filter(df, collection_date < ymd('2020-01-01')) %>%
    rename(Site_Code = site_code,
           Collection_Date = collection_date) %>%
    create_site_period() %>%
    r_friendly_colnames()
  
  post_2020 <- filter(df, collection_date >= ymd('2020-01-01')) %>%
    rename(Site_Code = site_code,
           Collection_Date = collection_date) %>%
    create_terrg_period() %>%
    r_friendly_colnames()
  
  output <- full_join(pre_2020, post_2020) 
  
  return(output) }

boot_mu <- function(x, resample) {
  mean(x[resample]) }

estimate_mu <- function(my_x) {
  my_output <- boot(my_x, boot_mu, R=2000) %>%
    tidy() %>%
    pull(statistic)
  return(my_output) }

estimate_mu_lower <- function(my_x) {
  my_ci <- boot(my_x, boot_mu, R=2000) %>%
    boot.ci(type="norm")
  my_output <- my_ci[[4]][2]
  return(my_output) }

estimate_mu_upper <- function(my_x) {
  my_ci <- boot(my_x, boot_mu, R=2000) %>%
    boot.ci(type="norm")
  my_output <- my_ci[[4]][3]
  return(my_output) }

extract_iso <- function(x_isotope) {
  paste(mean(x_isotope, na.rm=T)%>%
          format(digits=3, nsmall=1),
        ", ",
        if(is.na(sd(x_isotope, na.rm=T))) {"   "} else {
          sd(x_isotope, na.rm=T)%>%
            format(digits=1, nsmall=1)},
        " (",
        length(x_isotope)%>%
          format(digits=1, nsmall=1),
        ")",
        sep= '') }