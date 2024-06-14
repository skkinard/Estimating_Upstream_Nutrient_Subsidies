# 000_update_biocalc
# runs fish R scripts in order
# Sean Kinard
# last update: 2024-06-08

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('03_public/toolkit.R') # load packages and helper-functions

my_scripts <- list.files("03_public/calculation")
my_scripts <- my_scripts[str_detect(my_scripts,"update", negate=T)]
my_scripts <- my_scripts[str_detect(my_scripts,"dissertation_presentation", negate=T)]
my_scripts <- my_scripts[str_detect(my_scripts,"environment", negate=T)]
my_scripts <- my_scripts[str_detect(my_scripts,"site_map", negate=T)]
my_scripts <- sort(my_scripts, decreasing=F)

#------------------------------------------------------------------------------
# Source scripts
#------------------------------------------------------------------------------
for(i in 1:length(my_scripts)) {
  source(paste("03_public/calculation/", my_scripts[i], sep=''))
  rm(list=ls())
  my_scripts <- list.files("03_public/calculation")
  my_scripts <- my_scripts[str_detect(my_scripts,"update", negate=T)]
  my_scripts <- my_scripts[str_detect(my_scripts,"dissertation_presentation", negate=T)]
  my_scripts <- my_scripts[str_detect(my_scripts,"environment", negate=T)]
  my_scripts <- my_scripts[str_detect(my_scripts,"site_map", negate=T)]
  my_scripts <- sort(my_scripts, decreasing=F) }

#------------------------------------------------------------------------------
# End 000_update_environment