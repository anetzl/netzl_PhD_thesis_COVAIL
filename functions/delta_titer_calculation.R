# function to get delta data
delta_titer_data <- function(data, log_scale = T, titer_thresh = 40, to_long = T, v1 = "V1", v2 = "V3"){
  
  # do delta v3 - v1 log scale
  # plot titer differences of v3 - v1
  if(to_long){
    diff_data <-  data_to_long(data, titer_thresh)
  } else {
    diff_data <- data
  }
  
  diff_data_v1 <- diff_data %>% filter(visit_code == v1) 
  diff_data_v3 <- diff_data %>% filter(visit_code == v2)
  
  #get common antigens that do not have NA
  common_ags <- intersect(diff_data_v1$ag_name[!is.na(diff_data_v1$titer)], diff_data_v3$ag_name[!is.na(diff_data_v3$titer)])
  
  # get common subjids
  common_subj <- intersect(diff_data_v1$subjid, diff_data_v3$subjid)
  
  diff_data_v3 <- diff_data_v3 %>% filter(subjid %in% common_subj) %>%
    filter(ag_name %in% common_ags)
  diff_data_v1 <- diff_data_v1 %>% filter(subjid %in% common_subj) %>%
    filter(ag_name %in% common_ags)
  
  diff_data_v1 <- diff_data_v1[order(diff_data_v1$subjid, diff_data_v1$ag_name),]
  diff_data_v3 <- diff_data_v3[order(diff_data_v3$subjid, diff_data_v3$ag_name),]
  
  # check if they are identical in ag names and patient ID
  same <- identical(diff_data_v1$ag_name, diff_data_v3$ag_name) & identical(diff_data_v1$subjid, diff_data_v3$subjid)
  
  if(!same){
    warning(paste(v1, "and",v2,"data frames are not the same, simple subtraction does not work"))
    
    return(-1)
  } else {
    
    # now subtract titers
    if(log_scale){
      diff_data_v3$logtiter <- diff_data_v3$logtiter - diff_data_v1$logtiter
      diff_data_v3$titer <- 2^diff_data_v3$logtiter*10
    } else {
      diff_data_v3$titer <- diff_data_v3$titer - diff_data_v1$titer
      diff_data_v3$logtiter <- sapply(diff_data_v3$titer, function(x) {
        ifelse(x < 0, log2(1/(abs(x)/10)), log2(x/10))
        })
    }
    diff_data_v3$sr_group <- gsub(v2, "NA", diff_data_v3$sr_group)
    diff_data_v3$visit_diff <- paste(v2, "-", v1)
    
    return(diff_data_v3)
  }
  
}

delta_titer_gmts <- function(gmt_data){

# testing here
data_long <- data_to_long(sr_group_data,40, antigens = plot_antigens) %>%
  arm_code_wo_vacc()
sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, 40, cols_to_keep = c("arm_code", "visit_code", "inf_code",
                                                                   "age_code", "v_manuf_code"), pre_adj = TRUE) %>%
    filter(!all_below_thresh)

data_wide <- sr_group_gmt_plotdata %>%
ungroup() %>%
  select(!upper:lower) %>%
  select(!titer) %>%
  select(!sr_group) %>%
  pivot_wider(names_from = visit_code, values_from = logtiter)

View(data_wide)

}