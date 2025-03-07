#  Do line plots per serum group
rm(list = ls())
library(meantiter)
library(tidyverse)
library(ggplot2)
library(dplyr)

source(file.path("functions", "format_data.R"))
source(file.path("functions","scales.R"))
source(file.path("functions","gmt_calculation.R"))
source(file.path("functions","gmt_fold_change.R"))
source(file.path("functions","titerlineplot.R"))
source(file.path("functions","sr_group_color_functions.R"))
source(file.path("functions","remove_reactivity_bias.R"))
source(file.path("functions","delta_titer_calculation.R"))

source(file.path("code", "00_read_in_data.R"))

BO_D_ajudstment <- c("wo_2dose") 
day_visno <- c("D1","D29", "D91")

remove_all_breakthrough_sub <- FALSE
mark_undetected_bt <- FALSE

# format the data
data <- format_data(data, sr_group_code, lab = lab, combine_b_o_arm = FALSE, remove_all_bfl = FALSE, remove_all_subj_entries = FALSE,
                    arm_code_by = ifelse(as.numeric(day) >= 181 | date == "18APR2023", "trtactual", "trttrue"), stage = stage, keep_bfls = TRUE, only_target_visno = FALSE, 
                    target_visno = day_visno,
                    bo_adjust = NA, remove_all_oos = FALSE, keep_all_oosboost = TRUE, mark_undetected_BT = mark_undetected_bt)


data_b <- data

arm_code_order <- rev(arm_code_order)

for(d_adj in BO_D_ajudstment){
  
  # make figure dir
  if(correct_magnitude){
    figure_dir <- file.path("figures",date, paste0("stage", stage), "titer_lineplots",paste0(d_adj,"_b+o_adj"), "reactivity_adjusted")
    suppressWarnings(dir.create(figure_dir, recursive = TRUE))
    
  } else {
    figure_dir <- file.path("figures",date, paste0("stage", stage), "titer_lineplots", paste0(d_adj,"_b+o_adj"))
    suppressWarnings(dir.create(figure_dir, recursive = TRUE))
    
  }
  
  # now drop B+O, B+O arm day 29 make copy of B+O arm day 
  data <- BO_format_visit_code(data_b, d_adj)
  
  # filter out non responders
  if(lab != "Montefiori"){
    non_responders <- data %>% filter(D614G == 20) %>% filter(B.1.617.2 == 20) %>%
      filter(BA.1 == 20) %>% filter(B.1.351 == 20) %>% filter(`BA.4/5` == 20) %>% pull(subjid)
    
    data <- data %>% filter(!(subjid %in% non_responders))
    warning(paste0(length(unique(non_responders)), " subjects were filtered out because they were non responders"))
    
  }
  
  # remove oosboost and breakthroughs here
  #-------- Mark breakthroughs and remove all out of study boosts
  
  data <- remove_data_after_breakthrough(data, target_visits = day_visno, remove_all_entries = remove_all_breakthrough_sub, only_target_visits = TRUE, bo_adjust = d_adj)
  
  # remove people with an out of study booster
  data <- remove_data_after_oos(data, target_visits = day_visno, remove_all_entries = remove_all_breakthrough_sub, 
                                only_target_visits = TRUE, bo_adjustment = d_adj, only_mark_data = FALSE)
  
  # add info to group by
  # default grouping is by study arm 
  sr_group_data <- group_data_sr_group(data, by_arm = T, by_visit = T, by_strata = F, by_age = F, by_infection = T)
  
  
  # make data wide again
  if(correct_magnitude){
    sr_group_data <- adjust_individual_effect(data_to_long(sr_group_data)) %>%
      select(!logtiter:reactivity_bias) %>%
      pivot_wider(names_from = ag_name, values_from = "titer")
  }
  

  sr_group_data$age_code <- "combined"
  
  sr_group_data$arm_code <- factor(sr_group_data$arm_code, levels = arm_code_order)
  sr_group_data$age_code <- factor(sr_group_data$age_code, levels = age_code_order)
  sr_group_data$inf_code <- factor(sr_group_data$inf_code, levels = inf_code_order)
  
  # subset to idvls for which we have V3 data
  # v3_ids <- sr_group_data %>% filter(visit_code == day_visno[day]) %>% pull(subjid) %>% unique()
  
  # sr_group_data <- sr_group_data %>% filter(subjid %in% v3_ids)
  # create sr_group order for idvl panels
  sr_groups <- unique(sr_group_data$sr_group)
  sr_groups <- c(sr_groups[grepl("-inf", sr_groups)], sr_groups[grepl("non_inf", sr_groups)])
  
  sr_order_df <- sr_group_data[order(sr_group_data$arm_code, sr_group_data$age_code),]
  sr_order_arm_age <- unique(sr_order_df$sr_group)
  

  # ---------------------------- V1, V3 comparison -----------------------------------
  sr_group_order <- sr_order_arm_age
  
  sr_group_data$x_pos_by <- paste0(sr_group_data$visit_code, sr_group_data$inf_code)
  
  sr_group_data <- sr_group_data %>%
    arm_code_wo_vacc() %>%
    filter(visit_code %in% day_visno)

  for(v_manuf in "M"){
          
          sr_group_data_temp <- sr_group_data %>%
            filter(v_manuf_code == v_manuf)
          
          height_targets <- 2 + 1.5*length(unique(sr_group_data_temp$visit_code))
          plot_width <- length(unique(sr_group_data_temp$visit_code))
        
          plots <- titerlineplot_dodge(sr_group_data_temp, sr_group_colors, titer_thresh = 40, antigens = plot_antigens,
                                       facet_n_row = 3, sr_group_order = sr_group_order, gmt_facetter = "visit_code", color_by = c("arm_code"),
                                       line_by = "visit_code", shape_by = "inf_code",
                                       x_position_by = "age_code", cols_to_keep = c("arm_code","age_code", "x_pos_by", "inf_code", "v_manuf_code", "visit_code"), 
                                       show_group_count = T,
                                       gmt_grid_row = "inf_code",
                                       gmt_grid_col = "visit_code",
                                       dodge_group = "inf_code")$gmt + theme(legend.position = "none")
          
        ggsave(file.path(figure_dir, paste0(v_manuf, "_gmts_age_all_both_inf_by_visit.png")), plot = plots, dpi = 300, width = 1 + 2*plot_width, height = height_targets)
        }
  
  
}
