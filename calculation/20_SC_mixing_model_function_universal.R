# Create a function that runs simmr using carbon and sulfur for a given location, and grouping variable.

# Uses universal baseline (average stream vs average estuary source signatures)

run_simmr_uni <- function(my_data, my_location, my_group) {
  
  # fish and invertebrates @ location
  m_tibble <- my_data %>%
    add_sitevars() %>%
    filter(site_code == my_location &
             guild %in% c('Fish', 'Invertebrate'))
  
  # format mixtures
  m_matrix <- m_tibble %>% 
    select(sulfur, carbon) %>%
    as.matrix()
  
  # format mixture grouping variable
  m_groups <- m_tibble %>%
    pull(my_group) %>%
    as.matrix()
  
  # Source stats
  s_stats <-my_data %>%
    add_sitevars() %>%
      filter(guild %in% c('Aquatic_Producer', 'Terrestrial_Producer'),
             site_type %in% c('Stream', 'Estuary')) %>%
    filter( ! site_code %in% c('UN', 'LN')) %>%
    group_by(site_type) %>%
      summarise(sulfur_x = mean(sulfur),
                sulfur_sd = sd(sulfur),
                carbon_x = mean(carbon),
                carbon_sd = sd(carbon)) %>%
    mutate(source_type = case_when(
      site_type =='Stream' ~ 'Fresh',
      site_type =='Estuary' ~ 'Estuarine'))
  
  # source name vector
  s_names <- s_stats %>% pull(source_type) %>% as.matrix()
  
  # source mean sulfur vector
  s_means <- s_stats %>% select(sulfur_x, carbon_x) %>% as.matrix()
  
  # source standard deviation vector
  s_sds <- s_stats %>% select(sulfur_sd, carbon_sd) %>% as.matrix()
  
  # final prep and check for mixing model input        
  simmr_in =  simmr_load(mixtures = m_matrix,
                         group = m_groups,
                         source_names = s_names,
                         source_means = s_means,
                         source_sds = s_sds )
  
  ####################################
  # run mixing model
  ####################################
  
  # Run mixing model with informative priors
  simmr_out <- simmr_mcmc(simmr_in)
  
  ####################################
  # Extract data
  ####################################
  
  # extract statistics and quantile statsitics for each group
  sumstat <- summary(simmr_out, 
                     type = c('statistics', 'quantiles'),
                     group = 1:length(levels(as.factor(m_groups))))
  
  # statistics if more than one group
  if (m_tibble %>% pull(my_group)%>%unique()%>%length() > 1) {
    # Combine stats within list from different groups
    my_stats <- as.data.frame(do.call(cbind, sumstat$statistics)) 
    # create vector of group names
    my_stats_groups <- sort(rep(levels(simmr_out$input$group),2))
    # rename columns to include groups
    colnames(my_stats) <- paste(my_stats_groups, 
                                colnames(my_stats), 
                                sep="_")
    # transpose dataframe
    my_stats <- t(my_stats)
    # add group, site, and statistic columns
    my_stats <- mutate(as_tibble(my_stats), 
                       m_group = my_stats_groups,
                       site = rep(my_location, length(rownames(my_stats))),
                       statistic = rep(c('mean', 'sd'), 
                                       length(levels(simmr_out$input$group))) )
    # quantiles
    my_q <- as.data.frame(do.call(cbind, sumstat$quantiles))
    # create vector of group names
    my_q_groups <- sort(rep(levels(simmr_out$input$group),5))
    # rename columns to include groups
    colnames(my_q) <- paste(my_q_groups, 
                            colnames(my_q), 
                            sep="_")
    # transpose dataframe
    my_q <- t(my_q)
    # add group, site, and statistic columns
    my_q <- mutate(as_tibble(my_q), 
                   m_group = as.factor(my_q_groups),
                   site = rep(my_location, length(rownames(my_q))),
                   statistic = rep(c('2.5%', 
                                     '25%', 
                                     '50%', 
                                     '75%', 
                                     '97.5%'), 
                                   length(levels(simmr_out$input$group))) )
    
    mix_output <- full_join(my_stats, my_q)
    
    
  }
  
  # statistics if only one group
  if (m_tibble %>% pull(my_group)%>%unique()%>%length() == 1) {
    # extract statistics
    my_stats <- as.data.frame(sumstat$statistics)
    # create vector of group names
    my_stats_group <- levels(simmr_out$input$group)
    # rename columns to include groups
    colnames(my_stats) <- paste(my_stats_group, 
                                colnames(my_stats), 
                                sep="_")
    # transpose dataframe
    my_stats <- t(my_stats)
    # add group, site, and statistic columns
    my_stats <- mutate(as_tibble(my_stats), 
                       m_group = my_stats_group,
                       site = my_location,
                       statistic = c('mean', 'sd') )
    # extract quantiles
    my_q <- as.data.frame(sumstat$quantiles)
    # create vector of group names
    my_q_group <- levels(simmr_out$input$group)
    # rename columns to include groups
    colnames(my_q) <- paste(my_q_group, 
                            colnames(my_q), 
                            sep="_")
    # transpose dataframe
    my_q <- t(my_q)
    # add group, site, and statistic columns
    my_q <- mutate(as_tibble(my_q), 
                   m_group = my_q_group,
                   site = my_location,
                   statistic = c('2.5%', 
                                 '25%', 
                                 '50%', 
                                 '75%', 
                                 '97.5%') )
    
    mix_output <- full_join(my_stats, my_q)
  }
  
  ####################################
  # Diagnostics
  ####################################
  
  # Check for algorithm convergence (values should be 1)
  convergence_diagnostic <- summary(simmr_out, type = 'diagnostics')
  
  # Preview input data
  plot(simmr_in) + ggtitle(paste('Preview', my_location, sep='-'))
  
  # Check Priors
  prior_viz(simmr_out)
  
  # check model fit ( values should generally lie within the 50% intervals)
  # posterior_diagnostic <- posterior_predictive(simmr_out) 
  ## Error in `[.data.frame`(y_post_pred_out, , 3) : undefined columns selected
  # posterior_diagnostic_print <- print(posterior_diagnostic)
  
  # density plot
  plot(simmr_out, type = 'density') + ggtitle(paste('Source Density', my_location, sep='-'))
  
  # Source histograms (diagnoal), contour plots of relationships between sources (upper), correlations between sources (lower)
  # large negative correlations indicate that the model cannot discern between the two scources
  # large positive correlations are also possible when mixture data lie in a polygon consisting of multiple competing sources
  plot(simmr_out, type = 'matrix')
  
  ####################################
  # plots and Results
  ####################################
  compare_sources(simmr_out,
                  source_names = c('Fresh',
                                   'Estuarine'))
  
  if (m_tibble %>% pull(my_group)%>%unique()%>%length() > 1) {
    compare_groups(simmr_out,
                   source = 'Estuarine',
                   groups = 1:length(levels(as.factor(m_groups)))) }
  
  return(mix_output)
}
