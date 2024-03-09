# site_map
# Sean Kinard
# 1-25-2023

#------------------------------------------------------------------------------
# setup
#------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidygeocoder)
library(maps)
library(ggrepel)
library(paletteer)
library(grid)
library(gt)  
library (rgdal)
library (RSQLite)

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
map_limits <- function(x) {
  x %>%
    filter(latitude > 27.25) %>%
    filter(latitude < 29.25) %>%
    filter(longitude > -98.25) %>%
    filter(longitude < -95.75) }

extract_my_stream_crossings <- function(x) {
  x_tr <- x %>% filter(str_detect(stream_name, 'tranquitas'))
  x_sf <- x %>% filter(str_detect(stream_name, 'fernando'))
  x_ar <- x %>% filter(str_detect(stream_name, 'aransas'))
  x_pd <- x %>% filter(str_detect(stream_name, 'perdido'))
  x_mr <- x %>% filter(str_detect(stream_name, 'mission'))
  x_gc <- x %>% filter(str_detect(stream_name, 'garcitas'))
  x_pl <- x %>% filter(str_detect(stream_name, 'placedo'))
  x_nr <- x %>% filter(str_detect(stream_name, 'nueces'))
  x_ms <- x %>% filter(str_detect(stream_name, 'mustang'))
  x_cl <- x %>% filter(str_detect(stream_name, 'coleto'))
  x_gu <- x %>% filter(str_detect(stream_name, 'guadalupe'))
  x_nv <- x %>% filter(str_detect(stream_name, 'navidad'))
  x_lv <- x %>% filter(str_detect(stream_name, 'lavaca'))
  my_output <- bind_rows(x_tr, x_sf, x_ar, x_pd, x_mr, x_gc, x_pl, 
                         x_nr, x_ms, x_cl, x_gu, x_nv, x_lv)
  return(my_output) }
#------------------------------------------------------------------------------
# Shape layers
#------------------------------------------------------------------------------
# state data
state_map_data <- ne_states(country = 'united states of america', 
                            returnclass = "sf")
# rivers
rivers_tx <- read_sf(dsn = "data/map//Surface_Water/Surface_Water.shp")

# lakes
lake1 <- read_sf(dsn = "data/map/NHD_H_1210_HU4_Shape/Shape/NHDWaterbody.shp") %>%
  filter(SHAPE_Area > 6e-05)
lake2 <- read_sf(dsn = "data/map/NHD_H_1211_HU4_Shape/Shape/NHDWaterbody.shp") %>%
  filter(SHAPE_Area > 6e-05)

# Flowline
flowline1 <- read_sf(dsn = "data/map/NHD_H_1210_HU4_Shape/Shape/NHDFlowline.shp") %>%
  mutate(geometry = st_zm(geometry)) %>%
  filter(ftype %in% c(558)) %>%
  filter(SHAPE_Leng > 1e-2)
flowline2 <- read_sf(dsn = "data/map/NHD_H_1211_HU4_Shape/Shape/NHDFlowline.shp") %>%
  mutate(geometry = st_zm(geometry)) %>%
  filter(ftype %in% c(558))%>%
  filter(SHAPE_Leng > 1e-2)

# precip
precip <- read_sf(dsn = "data/map/precip1981_2010_a_tx.shp", 
                  layer = "precip1981_2010_a_tx") %>%
  rename('Rainfall (cm/yr)' = 'Precip_cm')

# oceans
oceans <- read_sf(
  dsn = 'data/map/natural_earth_vector/10m_physical/ne_10m_ocean.shp', 
  layer = "ne_10m_ocean")

# urban
urban <- read_sf(
  dsn = 'data/map/natural_earth_vector/50m_cultural/ne_50m_urban_areas.shp', 
  layer = "ne_50m_urban_areas")

#------------------------------------------------------------------------------
# point features
#------------------------------------------------------------------------------
# cities
cities <- read_sf(
  dsn = 'data/map/natural_earth_vector/10m_cultural/ne_10m_populated_places_simple.shp', 
  layer = "ne_10m_populated_places_simple") %>%
  map_limits()

cities <- cities %>% as_tibble() %>%
  select(latitude, longitude, name) %>%
  mutate(f_type='City') %>%
  rename(f_label=name)

# Site labels
d_site <- read_csv('data/environment/MDN_site_lat_lon.csv')
sample_sites <- d_site %>% 
  rename(latitude = site_lat, longitude = site_lon, f_label=site_code) %>%
  select(latitude, longitude, f_label) %>%
  mutate(f_type=ifelse(f_label=='CD', 'Dam', 'Sample Site'))

# Dams
dam_key <- read_csv('data/map/NHD_area_fcode_key.csv')
dam1 <- read_sf(dsn = "data/map/NHD_H_1210_HU4_Shape/Shape/NHDArea.shp") %>%
  left_join(dam_key) %>% 
  filter(feature_label=='Dam-Weir') %>%
  mutate(gnis_name = ifelse(is.na(gnis_name), 'Unnnamed', gnis_name))
dam2 <- read_sf(dsn = "data/map/NHD_H_1211_HU4_Shape/Shape/NHDArea.shp") %>%
  left_join(dam_key) %>%
  filter(feature_label %in% c('Dam-Weir','Bridge'))

dam1_f <- st_coordinates(st_centroid(dam1$geometry))%>% as_tibble() %>%
  bind_cols(dam1) %>%
  select(X,Y,gnis_name) %>%
  rename(longitude=X, latitude=Y, f_label=gnis_name) %>%
  mutate(f_type='Dam') %>%
  filter(str_detect(f_label, 'Soil', negate=T))
  
dam2_f <- st_coordinates(st_centroid(dam2$geometry))%>% as_tibble() %>%
  bind_cols(dam2) %>%
  select(X,Y,gnis_name) %>%
  rename(longitude=X, latitude=Y, f_label=gnis_name) %>%
  mutate(f_type='Dam')

# crossings
crossing_obs <- read_csv('data/map/DatabaseofStrea/crossing_obs_SC_Texas.csv')

my_obs <- crossing_obs %>% 
  extract_my_stream_crossings() %>% 
  select(latitude, longitude, stream_name, crossing_type) %>%
  rename(f_label = stream_name, f_type = crossing_type)

# merge features
my_point_features <- full_join(cities, sample_sites) %>%
  full_join(dam1_f) %>%
  full_join(dam2_f) %>%
  full_join(my_obs) %>%
  mutate(f_type=str_to_title(f_type),
         f_type = fct_relevel(f_type, c("Sample Site", "Dam", "Culvert", 
                                        "Bridge", "City")))
#------------------------------------------------------------------------------
# Site Map
#------------------------------------------------------------------------------
site_map <- ggplot() +
   geom_sf(data = precip, aes(fill=`Rainfall (cm/yr)`), color = NA) + 
   geom_sf(data = state_map_data, fill = 'grey90', lwd = .4, alpha = .1) + 
   geom_sf(data = urban, fill = 'grey60', alpha=.5) +
   geom_sf(data = lake1, color = NA, fill ='#4D7EAAFF') +
   geom_sf(data = lake2, color = NA, fill ='#4D7EAAFF') +
   geom_sf(data = oceans, fill = 'lightblue1') +
   geom_sf(data = rivers_tx, color = '#4D7EAAFF') +
   geom_sf(data = flowline1, color = '#4D7EAAFF') +
   geom_sf(data = flowline2, color = '#4D7EAAFF') +
  geom_point(data=my_point_features, 
             aes(x=longitude, y=latitude, 
                 shape=f_type, color=f_type, size=f_type)) +
  geom_label_repel(data = cities, 
                   aes(x = longitude, y = latitude, label=f_label),
                   size=3, fill='grey75', color='black', alpha=.7,
                   max.overlaps=1, box.padding = .9) + 
   geom_label_repel(data = filter(d_site, site_code!='CD'), 
                    aes(x = site_lon, y = site_lat, label=site_code),
                    size=4, min.segment.length = 0.1, box.padding = .5) +
   geom_text_repel(data = filter(d_site, site_code=='CD'), 
                    aes(x = site_lon, y = site_lat, label=site_name),
                    size=4, max.overlaps=1, box.padding = 1.5) +
  geom_point(data=my_point_features%>%filter(f_type=='Sample Site'), 
             aes(x=longitude, y=latitude),
             shape=21, color='black', fill='white', size=3) +
  geom_point(data=my_point_features%>%filter(f_type=='Dam'), 
             aes(x=longitude, y=latitude),
             shape=22, color='black', fill='red', size=3) +
  geom_point(data=my_point_features%>%filter(f_type=='Culvert'), 
             aes(x=longitude, y=latitude),
             shape=23, color='black', fill='chartreuse2', size=2) +
   coord_sf(xlim = c(-98.25, -95.75), ylim = c(27.25,29.25)) +
   scale_color_manual(values=c('white', 'red', 'chartreuse2', 'blue', 'grey30'), 
                      name=element_blank()) +
   scale_shape_manual(values=c(16,15,18,18,17), name=element_blank()) +
   scale_size_manual(values=c(3,2,3,2,2), name=element_blank()) +
   guides(color = guide_legend(override.aes = list(size = 4))) +
  paletteer::scale_fill_paletteer_c("ggthemes::Brown", 
                                     direction = -1,
                                     limits = c(50,125),
                                     breaks = c(50,75,100,125),
                                     name = "Rainfall \n(cm/yr)",
                                     guide = guide_colorbar(
                                       direction = "vertical",
                                       title.position = "top"))  +          
   theme_classic(base_size = 18)+                                           
   ggsn::north(location = "topleft", scale = 0.8, symbol = 12,              
               x.min = -95.8, x.max = -95.6, y.min = 27.15, y.max = 27.3) + 
   ggsn::scalebar(location = "bottomright", dist = 50,
                  dist_unit = "km", transform = TRUE, 
                  x.min=-97, x.max=-96, y.min=27.2, y.max=27.5,
                  st.bottom = FALSE, height = 0.025,
                  st.dist = 0.2, st.size = 5) +                            
   ylab(element_blank()) +                                                  
   xlab(element_blank()) +                                                  
   theme(legend.position = c(.9,.35),
         legend.background = element_rect(fill="grey75",
                                          size=0.5, linetype="solid", 
                                          colour ="black"),
         legend.key.size = unit(.75, 'cm'),
         legend.title.align = .5,
         legend.title = element_text(size=12),
         legend.text = element_text(size=10))

#------------------------------------------------------------------------------
# Inset
#------------------------------------------------------------------------------
inset <- ggplot() +
  geom_sf(data = oceans, fill = 'lightblue3') +
  geom_sf(data = state_map_data, fill = 'white', size = 3, linewidth=1) +
  geom_rect(aes(xmin = -98.25, xmax = -95.75, ymin = 27.25, ymax = 29.25), 
            color = "red", fill = NA, linewidth=1.5) +
  geom_text(label='Texas', aes(x=-99, y = 31.5), size = 6) +
  xlim(c(-107, -93)) +
  ylim(c(26,37)) +
  labs(x = NULL, y = NULL) +
  theme_test() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.background = element_rect(fill = "white"))

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------

site_map
print(inset, vp = viewport(0.247, 0.858, width = 0.25, height = 0.25))

ggsave('analysis/visualization/site_map.png',
       plot = site_map,
       width = 11,
       height = 11,
       units = c("in"))

ggsave('analysis/visualization/site_map_inset.png',
       plot = inset,
       width = 2,
       height = 2,
       units = c("in"))

caption_site_map <- 'Stable Isotope collection sites (labeled in white), where primary producers, fish, invertebrates, and environmental data were collected in January 2020. An overlay indicates the average annual precipitation (white-tan-green) from USGS PRISM data (1981-2010). Point features crossing streams of interest include dams (red), culverts (green), and bridges (blue). Cities and urband areas (labeled in grey) were included for geographic reference. This map was made with the National Hydrography Dataset and Natural Earth.'

#------------------------------------------------------------------------------
# End site_map