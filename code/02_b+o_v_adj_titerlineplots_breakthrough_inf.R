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


BO_D_ajudstment <- c("wo_2dose") 
day_visno <- c("D1","D29", "D91", "D181")

remove_all_breakthrough_sub <- FALSE
mark_undetected_bt <- TRUE

# format the data
data <- format_data(data, sr_group_code, lab = lab, combine_b_o_arm = FALSE, remove_all_bfl = FALSE, remove_all_subj_entries = FALSE,
                    arm_code_by = ifelse(as.numeric(day) >= 181 | date == "18APR2023", "trtactual", "trttrue"), stage = stage, keep_bfls = TRUE, only_target_visno = FALSE, 
                    target_visno = day_visno,
                    bo_adjust = NA, remove_all_oos = FALSE, keep_all_oosboost = TRUE, mark_undetected_BT = mark_undetected_bt)


data_b <- data %>%
  filter(visit_code %in% day_visno)
ag_order <- seq(1, 2*length(plot_antigens), 2)
names(ag_order) <- plot_antigens


# adapt colors, visit_codes
sr_group_colors["bD1", "Color"] <- sr_group_colors["inf", "Color"]
sr_group_colors["D1", "Color"] <- sr_group_colors["non_inf", "Color"]
visit_code_order <- unique(c(visit_code_order, "bD1"))

BO_D_ajudstment <- c("wo_2dose") #c("wo_2dose") 

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
  data <- remove_data_after_oos(data, target_visits = day_visno, remove_all_entries = FALSE, 
                                only_target_visits = TRUE, bo_adjustment = d_adj, only_mark_data = TRUE)

  
  # remove data only AFTER an out of study booster, keep time points before
  data %>%
    filter(VISNO < oosboost_visno) -> data
  
  data <- data %>%
    mutate(breakthrough = ifelse(btfl_visno <1000, TRUE, FALSE))
  
  # want to filter to visit BEFORE breakthrough, btfl_visno == VISNO gets me the visit AFTER breakhtrough
  # and need to remove breakthroughs of D1 for 2 dose (Day 57)
  data <- data %>%
    filter(btfl_visno != 3.5) %>%
    filter(btfl_visno != 1)
  


  
  # add info to group by
  # default grouping is by study arm 
  # change here visit code briefly for sr group
  visno_order <- sort(unique(data$VISNO))
 
  data <- data %>%
    mutate(btfl_visno = ifelse(breakthrough, visno_order[match(btfl_visno, visno_order)-1], btfl_visno),
           btfl_visno = ifelse((arm_code == "B+O, B+O:M" & btfl_visno == 3), 1, btfl_visno),
           visit_code = ifelse(btfl_visno == VISNO, paste0("b",visit_code), visit_code))
  
  
  # Remove data for breakthrough subjects after breakthrough
  data <- data %>%
    filter(VISNO <= btfl_visno)
  
  sr_group_data <- group_data_sr_group(data, by_arm = T, by_visit = T, by_strata = F, by_age = F, by_infection = T)

  # now make visit code normal again
  sr_group_data$visit_code <- gsub("b", "", sr_group_data$visit_code)

  sr_group_data$age_code <- "combined"
  
  sr_group_data_b <- sr_group_data
  
  # cycle through infected & non infected
  for(inf_stat in inf_names){
    
    if(correct_magnitude){
      figure_dir <- file.path("figures",date, paste0("stage", stage), ifelse(lab == "Montefiori","Montefiori" ,"."), 
        ifelse(combine_inf_stat,"." ,inf_names[inf_stat]), "titer_lineplots",paste0(d_adj,"_b+o_adj"), "reactivity_adjusted", "breakthroughs")
      suppressWarnings(dir.create(figure_dir, recursive = TRUE))
      
    } else {
      figure_dir <- file.path("figures",date, paste0("stage", stage), ifelse(lab == "Montefiori","Montefiori" ,"."), 
        ifelse(combine_inf_stat,"." ,inf_names[inf_stat]), "titer_lineplots", paste0(d_adj,"_b+o_adj"), "breakthroughs")
      suppressWarnings(dir.create(figure_dir, recursive = TRUE))
      
    }
  
    # combine sr groups by age and and filter for inf_stat
    sr_group_data <- sr_group_data_b %>% 
      filter(inf_code == inf_stat)
    
    sr_group_data$arm_code <- factor(sr_group_data$arm_code, levels = arm_code_order)
    sr_group_data$age_code <- factor(sr_group_data$age_code, levels = age_code_order)
    sr_group_data$inf_code <- factor(sr_group_data$inf_code, levels = inf_code_order)

    # ---------------------------- V1, V3 comparison -----------------------------------
    
    sr_groups <- unique(sr_group_data$sr_group)

    vacc_manuf_data <- sr_group_data %>%
      arm_code_wo_vacc()

    vacc_manuf_data <- data_to_long(vacc_manuf_data, titer_threshold, plot_antigens)

    vacc_manuf_data$arm_code <- factor(vacc_manuf_data$arm_code, levels = arm_code_order)
    
    vacc_manuf_data_b <- vacc_manuf_data
    
    # now do it for each vaccine
    for(v_manuf in "M"){
      
      vacc_manuf_data <- vacc_manuf_data_b %>%
        filter(v_manuf_code == v_manuf)
      
      unique_arms <- length(unique(vacc_manuf_data$arm_code))
      plot_width <- 0.5+2*unique_arms
      
      # calculate non breakthrough gmt
      gmt_data <- vacc_manuf_data[!(grepl("bD", vacc_manuf_data$sr_group)),] %>%
        sr_group_gmt_calc(., titer_threshold, cols_to_keep = c("arm_code", "visit_code", "inf_code","age_code", "breakthrough", "v_manuf_code")) %>%
        filter(!all_below_thresh) %>%
        filter(visit_code != tail(day_visno, n=1)) %>%
        mutate(sr_group = gsub("-D", "-bD", sr_group))
      
      # calculate difference between non breakthrough gmts and gmts
      fc_gmts <- vacc_manuf_data %>%
        sr_group_gmt_calc(., titer_threshold, cols_to_keep = c("arm_code", "visit_code", "inf_code","age_code", "v_manuf_code")) %>%
        filter(!all_below_thresh) %>%
        filter(visit_code != tail(day_visno, n=1)) %>%
        select(!lower:titer) %>%
        select(!count) %>%
        ungroup() %>%
        mutate(breakthrough = grepl("bD", sr_group),
               sr_group = gsub("-D", "-bD", sr_group)) %>%
        pivot_wider(names_from = breakthrough, values_from = logtiter) %>%
        mutate(log_fc = `TRUE` - `FALSE`,
               fc = ifelse(log_fc < 0, paste0("-", round(2^abs(log_fc), )), round(2^log_fc,1)))
      
      # now look at count of old vs young
      
      vacc_manuf_data %>%
        select(sr_group, subjid, stratalab, arm_code, visit_code, inf_code) %>%
        mutate(breakthrough = grepl("bD", sr_group),
               sr_group = gsub("-D", "-bD", sr_group)) %>%
        group_by(sr_group, visit_code, arm_code, stratalab, inf_code) %>% 
        mutate(total_count = length(unique(subjid))) %>%
        filter(breakthrough) %>%
        group_by(sr_group, visit_code, arm_code, stratalab, inf_code, breakthrough, total_count) %>% 
        summarize(count = length(unique(subjid))) %>%
        mutate(count_rel = round(count/total_count,2)) %>%
        select(!total_count) %>%
        mutate(stratalab = ifelse(grepl("Older", stratalab), ">65", "<65")) %>%
        pivot_wider(names_from = stratalab, values_from = c(count, count_rel)) -> count_table
      
      if(TRUE %in% grepl(">65", colnames(count_table))){
        count_table %>% mutate(ag_name = "D614G",
                               y = 0,
                               count_label = paste0("N=",`count_<65`, "; ", `count_>65`, " (", `count_rel_<65`, "; ", `count_rel_>65`,")"),
                               count_label = gsub("NA", "", count_label))-> count_table
      } else {
        count_table %>% mutate(ag_name = "D614G",
                               y = 0,
                               count_label = paste0("N=",`count_<65`, "; (", `count_rel_<65`, "; )"),
                               count_label = gsub("NA", "", count_label))-> count_table
      }
      
      titerlineplot_dodge(vacc_manuf_data[grep("bD", vacc_manuf_data$sr_group),],sr_group_colors = sr_group_colors, 
                          titer_thresh = titer_threshold, antigens = plot_antigens,
                          color_by = c("arm_code"),
                          x_position_by = "age_code", cols_to_keep = c("arm_code", "visit_code", "inf_code",
                                                                       "age_code", "breakthrough", "v_manuf_code"), show_group_count = F,
                          show_mean_line = F, to_long = F, line_by = "visit_code",
                          show_gmt_label = T,
                          gmt_grid_row = "visit_code", gmt_grid_col = "arm_code", dodge_group = "breakthrough", thinner_pre = FALSE)-> plot
      
      gmt_data$visit_code <- factor(gmt_data$visit_code, levels = day_visno)
      fc_gmts$visit_code <- factor(fc_gmts$visit_code, levels = day_visno)
      count_table$visit_code <- factor(count_table$visit_code, levels = day_visno)
    
      plot$all + 
        theme(legend.position = "none") + 
        geom_line(data = gmt_data, 
                  aes(x = ag_name,
                      y = logtiter,
                      group = sr_group),
                  color = "grey50",
                  linetype = "dotted") + 
        geom_pointrange(data = gmt_data, 
                        aes(x = ag_name,
                            y = logtiter,
                            ymin = lower,
                            ymax = upper,
                            group = sr_group),
                        color = "grey50",
                        fill = "white",
                        size = 0.4) + 
        geom_text(data = fc_gmts, 
                  aes(x = ag_name, 
                      y = rep(13.6, nrow(fc_gmts)),
                      label = fc),
                  color = "grey20",
                  size = 2.5) + 
        geom_text(data = count_table, 
                  aes(x = ag_name, 
                      y = y,
                      label = count_label),
                  color = "black",
                  size = 3,
                  hjust = 0) -> plot_breakthrough
      
      height_targets <- 2 + 2*length(unique(fc_gmts$visit_code))
      
      ggsave(file.path(figure_dir, paste0(v_manuf, "_idvl_breakthroughs_by_arm_visit.png")), plot_breakthrough, width = plot_width, height = height_targets, dpi = 300)
      
      
      #     # make plot of pre landscapes
      sr_group_data %>%
        mutate(visit_code = ifelse(breakthrough, gsub("D", "bD", visit_code), visit_code)) %>%
        group_data_sr_group(., by_arm = T, by_visit = T, by_strata = F, by_age = F, by_infection = T) %>%
        # mutate(visit_code = gsub("bD", "D", visit_code)) %>%
        arm_code_wo_vacc() %>%
        filter(v_manuf_code == v_manuf) %>%
        data_to_long(., titer_threshold, plot_antigens) -> test
      
      
      titerlineplot_dodge(test %>% 
                            mutate(inf_code = ifelse(breakthrough, "breakthrough", "non_inf")) %>%
                            filter(visit_code %in% c("bD1", "D1")),sr_group_colors = sr_group_colors, 
                          titer_thresh = titer_threshold, antigens = plot_antigens,
                          color_by = c("inf_code"),
                          x_position_by = "breakthrough", cols_to_keep = c("arm_code", "visit_code", "inf_code",
                                                                           "age_code", "breakthrough", "v_manuf_code"), show_group_count = F,
                          show_mean_line = F, to_long = F, line_by = "age_code",
                          gmt_facetter = "arm_code",
                          gmt_grid_row = "age_code", gmt_grid_col = "arm_code", dodge_group = "breakthrough")-> plot
      plot$all + theme(legend.position = "none",
                       strip.text.y = element_text(size = 8, colour = "white")) -> plot
    
      ggsave(file.path(figure_dir, paste0(v_manuf, "_idvl_breakthroughs_by_arm_D1.png")), plot, width = plot_width, height = 4.5, dpi = 300)
      
    }
    
    

  }
}

