#  Do line plots per serum group
# do B+O, B+O arm v8 as day 29 (= equal to V4 of other arms)
# and day 57 as day 1, as this corresponds to the pre of the 2nd B+O dose
# drop normal day 29 of B+O, B+O arm V4
# and make copy of V8 and put it as V5 (day 91) as it is last recorded time point
rm(list = ls())
library(meantiter)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(dplyr)


source(file.path("functions", "format_data.R"))
source(file.path("functions","scales.R"))
source(file.path("functions","gmt_calculation.R"))
source(file.path("functions","gmt_fold_change.R"))
source(file.path("functions","titerlineplot.R"))
source(file.path("functions","sr_group_color_functions.R"))
source(file.path("functions","remove_reactivity_bias.R"))

# read in data
source(file.path("code", "00_read_in_data.R"))


day_visno <- c("D1","D29", "D91")

mark_undetected_bt <- TRUE

# format the data
data <- format_data(data, sr_group_code, lab = lab, combine_b_o_arm = FALSE, remove_all_bfl = FALSE, remove_all_subj_entries = FALSE,
                    arm_code_by = ifelse(as.numeric(day) >= 181, "trtactual", "trttrue"), stage = stage, keep_bfls = TRUE, only_target_visno = FALSE, 
                    target_visno = day_visno,
                    bo_adjust =NA, remove_all_oos = FALSE, keep_all_oosboost = TRUE, mark_undetected_BT = mark_undetected_bt)

data_b <- data
ag_order <- seq(1, 2*length(plot_antigens), 2)
names(ag_order) <- plot_antigens



# adapt colors, visit_codes
sr_group_colors["bD1", "Color"] <- sr_group_colors["inf", "Color"]
sr_group_colors["D1", "Color"] <- sr_group_colors["non_inf", "Color"]
visit_code_order <- unique(c(visit_code_order, "bD1"))
inf_code_order <- c("non_inf", "non_inf_breakthrough", "inf", "inf_breakthrough", "breakthrough", "inf_oosboost", "non_inf_oosboost", "oosboost")


BO_D_ajudstment <- c("wo_2dose")

for(d_adj in BO_D_ajudstment){
  
  # now drop B+O, B+O arm day 29 make copy of B+O arm day 
  data <- BO_format_visit_code(data_b, d_adj)
  
  
  # filter out non responders
  if(lab != "Montefiori"){
    non_responders <- data %>% filter(D614G == 20) %>% filter(B.1.617.2 == 20) %>%
    filter(BA.1 == 20) %>% filter(B.1.351 == 20) %>% filter(`BA.4/5` == 20) %>% pull(subjid)
  
    data <- data %>% filter(!(subjid %in% non_responders))
    warning(paste0(length(unique(non_responders)), " subjects were filtered out because they were non responders"))
  
  }
  
  #-------- Mark breakthroughs and remove all out of study boosts
  
  data <- mark_breakthrough_data(data, target_visits = day_visno, only_target_visits = TRUE, bo_adjust = d_adj)
  
  # remove people with an out of study booster
  data <- remove_data_after_oos(data, target_visits = day_visno, remove_all_entries = TRUE, 
                                only_target_visits = TRUE, bo_adjustment = d_adj, only_mark_data = TRUE)
  
  # want to filter to visit BEFORE breakthrough, btfl_visno == VISNO gets me the visit AFTER breakhtrough
  # and need to remove breakthroughs of D1 for 2 dose (Day 57)
  data <- data %>%
    filter(btfl_visno != 3.5) %>%
    filter(btfl_visno != 1)
  
  
  # add info to group by
  # default grouping is by study arm 
  # change here visit code briefly for sr group
  visno_order <- sort(unique(data$VISNO))
  
  # count those that got a BT and oosboost as oosboost after oosboost
  data <- data %>%
    mutate(oosboost = ifelse(VISNO >= oosboost_visno, TRUE, FALSE),
          breakthrough = ifelse(VISNO >= btfl_visno & VISNO < oosboost_visno, TRUE, FALSE),
          visit_code = ifelse(VISNO >= btfl_visno, ifelse(VISNO >= oosboost_visno, paste0("o",visit_code), paste0("b",visit_code)), ifelse(VISNO >= oosboost_visno, paste0("o",visit_code), visit_code)))
  
  sr_group_data <- group_data_sr_group(data, by_arm = T, by_visit = T, by_strata = F, by_age = F, by_infection = T)

  # now make visit code normal again
  sr_group_data$visit_code <- gsub("b|o", "", sr_group_data$visit_code)

  sr_group_data$age_code <- "combined"
  
  sr_group_data_b <- sr_group_data
 
  # cycle through infected & non infected
  
    
    if(correct_magnitude){
      figure_dir <- file.path("figures",date, paste0("stage", stage), ifelse(lab == "Montefiori","Montefiori" ,"."), 
       ".", "titer_lineplots",paste0(d_adj,"_b+o_adj"), "reactivity_adjusted", ifelse(mark_undetected_bt,"breakthroughs_w_undetected", "breakthroughs"))
      suppressWarnings(dir.create(figure_dir, recursive = TRUE))
      
    } else {
      figure_dir <- file.path("figures",date, paste0("stage", stage), ifelse(lab == "Montefiori","Montefiori" ,"."), 
        ".", "titer_lineplots", paste0(d_adj,"_b+o_adj"), ifelse(mark_undetected_bt,"breakthroughs_w_undetected", "breakthroughs"))
      suppressWarnings(dir.create(figure_dir, recursive = TRUE))
      
    }
  
    # combine sr groups by age and and filter for inf_stat
    sr_group_data <- sr_group_data_b 
   
    sr_group_data$arm_code <- factor(sr_group_data$arm_code, levels = arm_code_order)
    sr_group_data$age_code <- factor(sr_group_data$age_code, levels = age_code_order)
    sr_group_data$inf_code <- factor(sr_group_data$inf_code, levels = inf_code_order)

    # ---------------------------- V1, V3 comparison -----------------------------------
    
    sr_groups <- unique(sr_group_data$sr_group)

    vacc_manuf_data <- sr_group_data %>%
      arm_code_wo_vacc() %>%
      filter(visit_code %in% day_visno)

    vacc_manuf_data <- data_to_long(vacc_manuf_data, titer_threshold, plot_antigens)

    vacc_manuf_data$arm_code <- factor(vacc_manuf_data$arm_code, levels = arm_code_order)
    vacc_manuf_data$visit_code <- factor(vacc_manuf_data$visit_code, levels = visit_code_order)
   
    vacc_manuf_data_b <- vacc_manuf_data
    # now do it for each vaccine
    for(v_manuf in "M"){
      
      vacc_manuf_data <- vacc_manuf_data_b %>%
        filter(v_manuf_code == v_manuf) %>%
        filter(visit_code != "D1")
      
      
      # change colouring according to breakthrough, by making hack with inf code
      vacc_manuf_data <- vacc_manuf_data %>%
        mutate(inf_code = paste0(inf_code, ifelse(breakthrough, "_breakthrough", ifelse(oosboost, "_oosboost", ""))),
               inf_code = factor(inf_code, levels = inf_code_order))
      
    
      plot_width <- 0.5+2*length(unique(vacc_manuf_data$inf_code))
      height_targets <- 2 + 1.5*length(unique(vacc_manuf_data$visit_code))
      
      # # now look at count of BT vs non-BT
      # sr_gmt <- sr_group_gmt_calc(vacc_manuf_data, 40, cols_to_keep = c("arm_code", "visit_code", "inf_code",
      #                                                                   "age_code", "v_manuf_code", "breakthrough", "oosboost"))
      # 
      # 
      # ## Now I want to plot post breakthrough titers by arm, see if we have arm specific differences here
      # sr_gmt %>%
      #   ungroup() %>%
      #   select(arm_code, visit_code, count, inf_code, breakthrough, oosboost) %>%
      #   unique() %>%
      #   group_by(arm_code, visit_code) %>%
      #   mutate(total_count = sum(count)) %>%
      #   ungroup() %>%
      #   mutate(count_label = paste0(ifelse(breakthrough, "BT n=", ifelse(oosboost, "OOSB n=", "non-BT n=")), count),
      #          y = ifelse(breakthrough, 13, ifelse(oosboost, 11, 12)),
      #          x = ag_order["BA.4/5"],
      #          # inf_code = inf_stat,
      #          age_code = "combined") -> count_table
      
    
      
      titerlineplot_dodge(vacc_manuf_data,sr_group_colors = sr_group_colors, 
                          titer_thresh = titer_threshold, antigens = plot_antigens,
                          color_by = c("arm_code"), shape_by = "age_code",
                          x_position_by = "age_code", cols_to_keep = c("arm_code", "visit_code", "inf_code",
                                                                       "age_code", "v_manuf_code", "breakthrough"), 
                          show_group_count = T,
                          show_mean_line = F, to_long = F, line_by = "visit_code",
                          show_gmt_label = F,
                          gmt_grid_row = "visit_code", gmt_grid_col = "inf_code", dodge_group = "inf_code")$gmt + 
        theme(legend.position = "none")-> plot
      
      
      
      ggsave(file.path(figure_dir, paste0(v_manuf, "_gmts_age_all_both_inf_by_visit_postBT.png")), plot = plot, dpi = 300, width = plot_width, height = height_targets)
      
      # plot$gmt + theme(legend.position = "none") + 
      #   geom_text(data = count_table, 
      #             aes(x = x, 
      #                 y = y,
      #                 label = count_label, 
      #                 color = arm_code),
      #             size = 3,
      #             hjust = 1) -> plot_breakthrough
      # 
      # plot_breakthrough
      
    }
}
    