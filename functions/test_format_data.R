# function to format csv data, to add sr group and sr name

format_data <- function(data, sr_group_code){
  rownames(sr_group_code) <- sr_group_code$full
 
  # fill NA values in data
  no_na_cohort <- c(which(data$Cohort!=""), nrow(data)+1)
  
  for(p in c(1:(length(no_na_cohort)-1))) {
    
    data$Cohort[c(no_na_cohort[p]:(no_na_cohort[(p+1)]-1))] <- data$Cohort[no_na_cohort[p]]
  }
  
  # fill code names
  data$arm_code <- sr_group_code[data$Cohort, "code"]
  data$age_code <- sr_group_code[data$age, "code"]
  data$inf_code <- sr_group_code[data$Immune.status, "code"]
  data$vacc_code <- sr_group_code[data$Vaccination, "code"]
  data$variant_code <- sr_group_code[data$Variant.of.prior.infection, "code"]
  
  # now add sr name
  data$sr_name <- unlist(lapply(1:nrow(data), function(x){
    paste(data$arm_code[x], data$age_code[x], data$inf_code[x], data$vacc_code[x], data$variant_code[x], x, sep = "-")
  }))
  
  # potentially make numeric in here 
  
  return(data)
}

group_data_sr_group <- function(data, by_arm = T, by_age = F, by_infection = F, by_variant =F, by_vaccine = F) {
  
  # add sr group with the pattern:
  # arm-age-infected-variant-vaccine
  data$sr_group <- ""
  
  if(by_arm) {
    data$sr_group <- paste0(data$sr_group, data$arm_code)
  } else {
    data$sr_group <- paste0(data$sr_group, "NA")
  }
  
  if(by_age) {
    data$sr_group <- paste(data$sr_group, data$age_code, sep = "-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_infection) {
    data$sr_group <- paste(data$sr_group, data$inf_code, sep = "-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_variant){
    data$sr_group <- paste(data$sr_group, data$variant_code, sep ="-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_vaccine){
    data$sr_group <- paste(data$sr_group, data$vacc_code, sep ="-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  return(data)
}


filter_data_sr_group <- function(data, target_arm = NA, target_age = NA, target_infection = NA, target_variant =NA, target_vaccine = NA,
                                 antigens = c("D614G","B.1.1.7", "B.1.1.7+E484K", "B.1.351", "P.1.1", "B.1.617.2", "BA.1", "BA.2")) {
  
  # filter target frame
  filtered <- data
  if(!is.na(target_arm)) {
    filtered %>% filter(arm_code == target_arm) -> filtered
  } 
  if(!(is.na(target_age))) {
    filtered %>% filter(age_code == target_age) -> filtered
  } 
  if(!(is.na(target_infection))) {
    filtered %>% filter(inf_code == target_infection) -> filtered
  }
  if(!(is.na(target_variant))){
    filtered %>% filter(variant_code == target_variant) -> filtered
  }
  if(!(is.na(target_vaccine))){
    filtered %>% filter(vacc_code == target_vaccine) -> filtered
  }
  
  sr_group <- paste(target_arm, target_age, target_infection, target_variant, target_vaccine, sep = "-")
  
  filtered$sr_group <- sr_group
  
  return(filtered)
}

sr_group_log_titers_from_data <- function(data, target_arm = NA, target_age = NA, target_infection = NA, target_variant =NA, target_vaccine = NA,
                                          antigens = c("D614G","B.1.1.7", "B.1.1.7+E484K", "B.1.351", "P.1.1", "B.1.617.2", "BA.1", "BA.2")){
  
  # filter target frame
  filtered <- filter_data_sr_group(data, target_arm, target_age, target_infection, target_variant, target_vaccine, antigens)
  
  filtered_long <- filtered
  filtered %>% 
    select("sr_name", antigens) %>%
    column_to_rownames("sr_name") -> filtered
  
  filtered <- t(filtered)
  filtered <- gsub(",", ".", filtered)
  filtered <- apply(filtered, 1:2, function(x) as.numeric(x))
  filtered <- log2(filtered/10)
  
  return(list(log_titers = filtered, full_sr_group = filtered_long))
}