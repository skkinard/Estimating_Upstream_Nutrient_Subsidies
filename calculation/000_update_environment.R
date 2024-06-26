# 000_update_environment
# runs fish R scripts in order
# Sean Kinard
# last update: 2023-06-23

#------------------------------------------------------------------------------
# Setup
#------------------------------------------------------------------------------
source('03_public/toolkit.R') # load packages and helper-functions

my_scripts <- list.files("analysis/calculation")
my_scripts <- my_scripts[str_detect(my_scripts,"environment", negate=F)]
my_scripts <- my_scripts[str_detect(my_scripts,"update", negate=T)]
my_scripts <- sort(my_scripts, decreasing=F)

#------------------------------------------------------------------------------
# Source scripts
#------------------------------------------------------------------------------
for(i in 1:length(my_scripts)) {
  source(paste("exploration/calculation/", my_scripts[i], sep=''))
  rm(list=ls())
  my_scripts <- list.files("analysis/calculation")
  my_scripts <- my_scripts[str_detect(my_scripts,"environment", negate=F)]
  my_scripts <- my_scripts[str_detect(my_scripts,"update", negate=T)]
  my_scripts <- sort(my_scripts, decreasing=F)  }

#------------------------------------------------------------------------------
# End 000_update_environment