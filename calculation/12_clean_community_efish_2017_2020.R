# clean_community_efish_2017_2020
# Sean Kinard
# 2023-07-20

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('03_public/toolkit.R') # load packages and helper-functions

fish_species <- read_csv("03_local_files/data/community/fish_species_with_transient.csv")

fish_traits <- read_csv('03_local_files/data/community/functional_trait_binary.csv')

d_1719 <- read_csv('03_local_files/data/community/EFish_abundance_2017-2019.csv')
d_20 <- read_csv('03_local_files/data/community/EFish_abundance_2020.csv')
d_pass_RT <- read_csv('03_local_files/data/community/fish_pass_RAPID_TERRG.csv')

my_vars <- c('site_code', 'collection_date', 'collection_period', 
             'common_name', 'genus', 'species',
             'pass1_abun', 'pass2_abun', 'pass3_abun')

#------------------------------------------------------------------------------
# tidy 2019
#------------------------------------------------------------------------------
tidy_terrg_raw <- function(x) {
  output <- x %>%
  r_friendly_colnames() %>%
    mutate(collection_date=str_replace_all(collection_date, '/', '-'),
           collection_date=mdy(collection_date)) %>%
    rename(Collection_Date = collection_date,
           Site_Code=sitecode) %>%
    create_site_period() %>%
    mutate(Site_Code=substr(Site_Code,1,2)) %>%
    r_friendly_colnames() %>%
    select(any_of(my_vars))
  
  colnames(output) <- str_replace_all(colnames(output), '_abun', '')
  
  output <- output %>% 
    mutate(pass1=as.numeric(pass1),
           pass2=as.numeric(pass2),
           pass3=as.numeric(pass3)) 
    
  return(output) }


pass1719 <- d_1719 %>% tidy_terrg_raw() 
pass20 <- d_20 %>% tidy_terrg_raw()

dc <- full_join(pass1719, pass20) %>% arrange(collection_date, site_code)
#------------------------------------------------------------------------------
# Fix genus species levels
#------------------------------------------------------------------------------

# remove "\xff"
dc <- mutate(dc, genus = str_replace_all(genus, "\xff", ""))

# split "genus species" and import species to blank species cells
dc <- dc %>% 
  mutate(
    genus = ifelse(common_name == 'Largemouth bass', 'Micropterus', genus),
    species = ifelse(common_name == 'Largemouth bass', 'salmoides', species)) %>%
  separate(genus, into=c('genus', 'spe'), sep = " ") %>%
  mutate(species = ifelse(is.na(species) & ! is.na(spe), spe, species)) %>%
  select(-spe)

# convert brookside silverside to inland silverside
dc <- dc %>%
  mutate(genus = ifelse(genus %in% c("Labidesthes", "Labideshes"), 
                        "Menidia", genus))

# convert scarlet shiner to red shiner
dc <- dc %>%
  mutate(genus = ifelse(genus == "Lythrurus", 
                        "Cyprinella", genus))

crayfish_names <- c('Crawdad','Crawfish','Crayfish')

dc <- dc %>%
  mutate(species = case_when(
    species == "pulvereus" ~ "grandis",
    species == "zebrinus" ~ "grandis",
    species == "cyanoguttatus" ~ "cyanoguttatus",
    species == "sicculus" ~ "beryllina",
    species == "osseus" ~ "oculatus",
    species == "auritius" ~ "auritus",
    species == "symmetricus" ~ "macrochirus",
    species == "fumeus" ~ "texanus",
    species == "punctulatus" ~ "punctatus",
    species == "atrocaudalis" ~ "texanus",
    species == "chalybaeus" ~ "texanus",
    species == "shumardi" ~ "chrysotus",
    species == "volucellus" ~ "texanus",
    species == "Sicculus" ~ "beryllina",
    common_name == "Largemouth bass" ~ "salmoides",
    TRUE ~ species)) %>%
  mutate(genus = ifelse(species == 'chrysotus', 'Fundulus', genus))

dc <- dc %>%
  mutate(species = case_when(
    common_name %in% crayfish_names ~ 'clarkii',
    common_name == "Brook Silverside" ~ 'beryllina',
    common_name == "Rio Grande Cichlid" ~ 'cyanoguttatus',
    common_name == "Bluegill/Green hybrid***" ~ 'macrochirus',
    common_name == "Warmouth/Green hybrid***" ~ 'cyanellus',
    common_name == "Grass shrimp" ~ 'pugio',
    genus == "Gambusia" ~ "affinis",
    species == "spp." ~ "megalotis",
    species == "auritis" ~ "auritus",
    species == "auritius" ~ "auritus",
    species == "margi0tus" ~ "marginatus",
    species == "megaloits" ~ "megalotis",
    species == "punctulatus" ~ "punctatus",
    species == "sicculus" ~ "beryllina",
    species == "variegatus" ~ "variegatus",
    genus == 'Macrobrachium' ~ 'ohione',
    substr(species, 1, 4) == "punc" ~ "punctatus",
    TRUE ~ species)) %>%
  mutate(genus = case_when(
    common_name %in% crayfish_names ~ 'Procambarus',
    common_name == "Grass shrimp" ~ 'Palaemonetes',
    species == "macrochirus" ~ 'Lepomis',
    species == "cyanellus" ~ 'Lepomis',
    genus == "Dormitor" ~ "Dormitator",
    genus == "Herichthys" ~ "Herichthys",
    genus == "Labideshes" ~ 'Menidia',
    genus == "Labidesthes" ~ "Menidia",
    genus == "Lepomis" ~ "Lepomis",
    genus == "Minidia" ~ "Menidia",
    substr(species, 1, 4) == "punc" ~ "Ictalurus",
    TRUE ~ genus))

dc <- dc %>% mutate(genus = case_when(
  common_name == "Rio Grande Cichlid" ~ "Herichthys",
  common_name == "Tadpole madtongue" ~ 'Noturus',
  common_name == "Pimephales sp." ~ 'Pimephales',
  common_name == "Fry" ~ 'Gambusia',
  common_name == "mottled killifish" ~ 'Fundulus',
  common_name == "mullet" ~ 'Agonostomus',
  common_name == "anchovy" ~ 'Anchoa',
  common_name == "Menhaden" ~ 'Brevoortia',
  common_name == "killifish" ~ 'Fundulus',
  common_name == "Bass" ~ 'Micropterus',
  common_name == "Topminnow" ~ 'Fundulus',
  common_name == "Uni. Carp" ~ 'Cyprinus',
  common_name == "Unidentified flatfish" ~ 'Trinectes',
  common_name == "Unidentified minnow" ~ 'Pimephales',
  common_name == "Unidentified stripe" ~ 'Fundulus',
  common_name == "Black tailed shiner hybrid" ~ 'Notropis',
  common_name == "blackstripe shiner" ~ 'Notropis',
  common_name == "Gambusia blackspot" ~ 'Gambusia',
  common_name == "Gambusia blackspot?" ~ "Gambusia",
  common_name == "shiner" ~ "Cyprinella",
  common_name == "Mullet" ~ "Mugil",
  common_name == "Tilapia" ~ "Oreochromis",
  common_name == "Uni. Darter" ~ "Etheostoma",
  common_name == "Western Mosquitofish" ~ "Gambusia",
  common_name == "Ribbon shiner" ~ "Notropis",
  common_name == "Tailight shiner" ~ "Cyprinella",
  common_name == "Spotted Bass" ~ "Micropterus",
  common_name == "Warmouth Sunfish" ~ "Lepomis",
  common_name == "Bluegill hybrid" ~ "Lepomis",
  common_name == "Blacktail/red shiner Hybrid" ~ "Notropis",
  common_name == "Unidentified stripe" ~ "Fundulus",
  common_name == "bass" ~ "Micropterus",
  common_name == "Unknown Shiner" ~ "Cyprinella",
  common_name == "Spotted Gar" ~ "Lepisosteus",
  common_name == "Amazon Molly" ~ "Poecilia",
  common_name == "Molly" ~ "Poecilia",
  common_name == "Smallmouth Buffalo" ~ "Ictiobus",
  common_name == "Fiddler crab" ~ "Minuca",
  common_name == "Freshwater drum" ~ "Rhonciscus",
  common_name %in% c("Macrobrachium", 
                     "Giant Freshwater Prawn", 
                     "Giant River Prawn",
                     "Marine shrimp")  ~ "Macrobrachium",
  TRUE ~ genus)) %>%
  mutate(species = case_when(
    common_name == "Tadpole madtongue" ~ 'gyrinus',
    common_name == "Pimephales sp." ~ 'vigilax',
    common_name == "Fry" ~ 'affinis',
    common_name == "mottled killifish" ~ 'grandis',
    common_name == "mullet" ~ 'monticola',
    common_name == "anchovy" ~ 'mitchilli',
    common_name == "Menhaden" ~ 'patronus',
    common_name == "killifish" ~ 'grandis',
    common_name == "Bass" ~ 'salmoides',
    common_name == "Topminnow" ~ 'notatus',
    common_name == "Uni. Carp" ~ 'carpio',
    common_name == "Unidentified flatfish" ~ 'maculatus',
    common_name == "Unidentified minnow" ~ 'vigilax',
    common_name == "Unidentified stripe" ~ 'grandis',
    common_name == "Black tailed shiner hybrid" ~ 'texanus',
    common_name == "blackstripe shiner" ~ 'texanus',
    common_name == "Gambusia blackspot" ~ 'affinis',
    common_name == "shiner" ~ "venusta",
    common_name == "Mullet" ~ "cephalus",
    common_name == "Tilapia" ~ "aureus",
    common_name == "Uni. Darter" ~ "gracile",
    common_name == "Western Mosquitofish" ~ "affinis",
    common_name == "Mosquito fish" ~ "affinis",
    common_name == "Gambusia spp." ~ "affinis",
    common_name == "Gambusia blackspot?" ~ "affinis",
    common_name == "Spotted Bass" ~ "salmoides",
    common_name == "Tailight shiner" ~ "venusta",
    common_name == "Bluntnose Minnow" ~ "vigilax",
    common_name == "Green sunfish" ~ "cyanellus",
    common_name == "Warmouth Sunfish" ~ "gulosus",
    common_name == "Bluegill hybrid" ~ "macrochirus",
    common_name == "Blacktail/red shiner Hybrid" ~ "texanus",
    common_name == "Unidentified stripe" ~ "grandis",
    common_name == "Bluntnose Minnow"~ "vigilax",
    common_name == "bass"~ "salmoides",
    common_name == "Taillight shiner" ~ "texanus",
    common_name == "Tadpole Madtom Catfish" ~ "gyrinus",
    common_name == "Blue Crab" ~ "sapidus",
    common_name == "Unknown Shiner" ~ "lutrensis",
    common_name == "Redbreast Sunfish" ~ "auritus",
    common_name == "Bullhead Catfish" ~ "melas",
    common_name == "Amazon Molly" ~ "latipinna",
    common_name == "Molly" ~ "latipinna",
    common_name == "Smallmouth Buffalo" ~ "bubalus",
    common_name == "Fat Sleeper" ~ "maculatus",
    common_name == "Fiddler crab" ~ "longisignalis",
    common_name == "Freshwater drum" ~ "crocro",
    common_name %in% c("Macrobrachium", 
                       "Giant Freshwater Prawn", 
                       "Giant River Prawn",
                       "Marine shrimp")  ~ "ohione",
    TRUE ~ species))

# change striped mullet to mountain mullet at PL & MR on voucher dates
dc <- dc %>%
  mutate(genus = ifelse(site_code == 'PL' & 
                          collection_date %in% c(
                            '2017-10-31', '2017-11-11', '2017-12-05', '2018-03-02', 
                            '2018-04-01', '2018-05-10', '2018-06-07') &
                          genus == 'Mugil', "Agonostomus", genus),
         species = ifelse(
           site_code == 'PL' & 
             collection_date %in% c(
               '2017-10-31', '2017-11-11', '2017-12-05', '2018-03-02', 
               '2018-04-01', '2018-05-10', '2018-06-07') &
             species == "cephalus", "monticola", species) )

# Change N.texanus at AR to C.lutrensis
dc <- dc %>%
  mutate(genus = case_when(
    site_code == "AR" &
      species == 'texanus' ~ "Cyprinella",
    TRUE ~ genus),
    species = case_when(
      site_code == "AR" &
        species == 'texanus' ~ "lutrensis",
      TRUE ~ species))

dc <- filter(dc, ! common_name %in% c('Unidentified', 'Unidentified?', 'tadpole',
                                      'Tadpoles', 'tiger grub', 'Tadpole'))

#------------------------------------------------------------------------------
# fix bad lepomis
#------------------------------------------------------------------------------
# major_Lepomis (most-caught sunfish per site)
major_Lepomis <- dc %>%
  filter(genus=='Lepomis') %>%
  mutate(sumpass = pass1+pass2+pass3) %>%
  group_by(site_code, species) %>%
  dplyr::summarize(grand_sum = sum(sumpass, na.rm=T)) %>%
  group_by(site_code) %>%
  mutate(max_pass=max(grand_sum)) %>%
  mutate(is_max_pass = ifelse(near(grand_sum, max_pass), 'yes', 'no')) %>%
  filter(is_max_pass == 'yes') %>%
  select(site_code, species)

sites_megalotis <- major_Lepomis %>% 
  filter(species == 'megalotis') %>% pull(site_code)
sites_auritus <- major_Lepomis %>% 
  filter(species == 'auritus') %>% pull(site_code)
sites_macrochirus <- major_Lepomis %>% 
  filter(species == 'macrochirus') %>% pull(site_code)
sites_cyanellus <- major_Lepomis %>% 
  filter(species == 'cyanellus') %>% pull(site_code)

# convert juvi and unkown sf to most-caught sunfish
bad_lepomis <- c("Sunfish juvi", "Juvenile sunfish", "Juvenile Sunfish", 
                 "Juvi. Sunfish", "Unidentified sunfish", 
                 "Longear/Bluegill sunfish", "Bluegill hybrid",
                 "Warmouth/ longear", "Juvenile sunfish")

dc <- dc %>%
  mutate(genus = ifelse(common_name %in% bad_lepomis, "Lepomis", genus),
    species = case_when(
      common_name %in% bad_lepomis &
        site_code %in% sites_megalotis ~ 'megalotis',
      common_name %in% bad_lepomis &
        site_code %in% sites_auritus ~ 'auritus',
      common_name %in% bad_lepomis &
        site_code %in% sites_macrochirus ~ 'macrochirus',
      common_name %in% bad_lepomis &
        site_code %in% sites_cyanellus ~ 'cyanellus',
      TRUE ~ species))

# re-impodc taxanomic based on clean genus and species
dc <- dc %>% 
  select(contains("pass"), site_code, collection_date, genus, species, 
         common_name) %>%
  left_join(fish_species)

#------------------------------------------------------------------------------
# merge with taxanomic
#------------------------------------------------------------------------------
dc2 <- dc %>% left_join(fish_species)

#------------------------------------------------------------------------------
# Export
#------------------------------------------------------------------------------
write_csv(dc2, '03_public/output/community_efish_2017_2020.csv')

#------------------------------------------------------------------------------
# end clean_community_efish_2017_2020