# remove reactivity bias from titertable
remove_reactivity_bias_logtiter <- function(table) {
  row_names <- rownames(table)
  table <- sapply(as.data.frame(table), function(x) {
    x - (mean(x, na.rm = T) - mean(table, na.rm = T))
  })
  
  rownames(table) <- row_names
  
  return(table)
}

# Function for accounting for a sample's indivdual affect
# calculate GMT per sample, subtract serum group GMT from that to get reactivity bias
# subtract sample's reactivity bias from sample titrations
adjust_individual_effect <- function(data, dil_stepsize = 0){
  
  # calculate gmt per sample
  data %>%
    group_by(sr_name) %>%
    mutate(gmt_sample = meantiter::mean_titers(titer, method ="bayesian", dilution_stepsize = dil_stepsize)$mean) ->titers_variation_adjusted 
  
  # calculate gmt per serum group
  titers_variation_adjusted %>%
    group_by(sr_group) %>%
    mutate(gmt_sr_group = meantiter::mean_titers(titer, method ="bayesian", dilution_stepsize = dil_stepsize)$mean) ->titers_variation_adjusted 
  
  
  # adjust serum reactivity bias
  titers_variation_adjusted %>%
    mutate(reactivity_bias = gmt_sample - gmt_sr_group,
           logtiter = logtiter - reactivity_bias,
          # logtiter = logtiter -gmt_sample,
           titer = 2^logtiter*10) -> titers_variation_adjusted
  
  return(titers_variation_adjusted)
  
}