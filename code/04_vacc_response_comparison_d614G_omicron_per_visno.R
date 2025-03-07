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
library(patchwork)
library(gridExtra)
library(grid)
library(Racmacs)

source(file.path("functions", "format_data.R"))
source(file.path("functions","scales.R"))
source(file.path("functions","gmt_calculation.R"))
source(file.path("functions","gmt_fold_change.R"))
source(file.path("functions","titerlineplot.R"))
source(file.path("functions","sr_group_color_functions.R"))
source(file.path("functions","remove_reactivity_bias.R"))
source(file.path("functions","map_longinfo.R"))

do_x_axis_time_plots <- FALSE

# read in map 2x Vax data
map <- read.acmap("data/maps/map_ndsubset_no_outliers_slope_adjusted.ace")
map_long <- long_map_info(map)

map_long %>%
  filter(sr_group %in% c("2x mRNA-1273", "2x mRNA-1273 new")) %>%
  mutate(sr_group = "2x mRNA-1273") %>%
  filter(ag_name == "D614G")-> map_vax

d614g_gmt <- sr_group_gmt_calc(map_vax, 10) %>%
  mutate(standardise_encounters = "2x Vax",
         exact_covail = "D29",
         binned_covail = "D29",
         Study = "Wilks",
         OmicronVariant = "D614G",
         log_fold_change = logtiter,
         vaccine_manufacturer = "mRNA",
         standardised_assay = "PV",
         standardised_pseudo = "Lentiviral")


# read in map 2x Vax data
forest_data <- read.csv("data/titer_data/Netzl_et_al_data_for_imprinting_response.csv")
forest_data %>%
  mutate(mean_month = mean_days/30,
         rounded_month = ifelse(mean_month < 0.75, 0.5, round(mean_month))) %>%
  filter(Comparator.antigen == "D614G") -> forest_data




# do binning by time post. COVAIL has 2w, 1m, 3m, 6m, 9m
month_covail <- c(#"0.5" = "D15",
                  "1" = "D29",
                  "3" = "D91",
                  "6" = "D181",
                  "9" = "D271")

forest_data %>%
  mutate(exact_covail = month_covail[as.character(forest_data$rounded_month)]) -> forest_data

# do binning by 0.75-1.25 of month covail of mean month

covail_bins <- sapply(forest_data$mean_month, function(x){
  
  for(m in names(month_covail)){
    
    m_num <- as.numeric(m)
    if(x > 0.5*m_num & x < 1.5*m_num){
      
      return(month_covail[m])
    
    }
    
  }
  
  return(NA)
  
})

forest_data$binned_covail <- unlist(covail_bins)
forest_data <- plyr::rbind.fill(forest_data, d614g_gmt)

forest_data %>%
  group_by(OmicronVariant, standardise_encounters, exact_covail) %>%
  summarise(gmt = Rmisc::CI(na.omit(log_fold_change))["mean"],
            upper = Rmisc::CI(na.omit(log_fold_change))["upper"],
            lower = Rmisc::CI(na.omit(log_fold_change))["lower"],
            n = length(na.omit(log_fold_change)),
            visit_code = factor(exact_covail, levels = visit_code_order)) %>%
  mutate(logtiter = gmt,
         ag_name = OmicronVariant,
         arm_code = standardise_encounters,
         inf_code = "non_inf",
         x = ag_order[ag_name],
         text = paste0("n(", ag_name, " in ", arm_code, ")=",n)) %>%
  filter(!is.na(visit_code)) %>%
  unique()-> gmt_exact

forest_data %>%
  group_by(OmicronVariant, standardise_encounters, binned_covail) %>%
  summarise(gmt = Rmisc::CI(na.omit(log_fold_change))["mean"],
            upper = Rmisc::CI(na.omit(log_fold_change))["upper"],
            lower = Rmisc::CI(na.omit(log_fold_change))["lower"],
            n = length(na.omit(log_fold_change)),
            visit_code = factor(binned_covail, levels = visit_code_order)) %>%
  mutate(logtiter = gmt,
         ag_name = OmicronVariant,
         arm_code = standardise_encounters,
         inf_code = "non_inf",
         x = ag_order[ag_name],
         text = paste0("n(", ag_name, " in ", arm_code, ")=",n)) %>%
  filter(!is.na(visit_code)) %>%
  unique() -> gmt_binned

# here for only PV
forest_data %>%
  filter(standardised_assay == "PV") %>%
  group_by(OmicronVariant, standardise_encounters, exact_covail) %>%
  summarise(gmt = Rmisc::CI(na.omit(log_fold_change))["mean"],
            upper = Rmisc::CI(na.omit(log_fold_change))["upper"],
            lower = Rmisc::CI(na.omit(log_fold_change))["lower"],
            n = length(na.omit(log_fold_change)),
            visit_code = factor(exact_covail, levels = visit_code_order)) %>%
  mutate(logtiter = gmt,
         ag_name = OmicronVariant,
         arm_code = standardise_encounters,
         inf_code = "non_inf",
         x = ag_order[ag_name],
         text = paste0("n(", ag_name, " in ", arm_code, ")=",n)) %>%
  filter(!is.na(visit_code)) %>%
  unique()-> gmt_exact_pv

forest_data %>%
  filter(standardised_assay == "PV") %>%
  group_by(OmicronVariant, standardise_encounters, binned_covail) %>%
  summarise(gmt = Rmisc::CI(na.omit(log_fold_change))["mean"],
            upper = Rmisc::CI(na.omit(log_fold_change))["upper"],
            lower = Rmisc::CI(na.omit(log_fold_change))["lower"],
            n = length(na.omit(log_fold_change)),
            visit_code = factor(binned_covail, levels = visit_code_order)) %>%
  mutate(logtiter = gmt,
         ag_name = OmicronVariant,
         arm_code = standardise_encounters,
         inf_code = "non_inf",
         x = ag_order[ag_name],
         text = paste0("n(", ag_name, " in ", arm_code, ")=",n)) %>%
  filter(!is.na(visit_code)) %>%
  unique() -> gmt_binned_pv

gmts <- list(#"exact" = gmt_exact,
             "binned" = gmt_binned,
             #"exact_pv" = gmt_exact_pv,
             "binned_pv" = gmt_binned_pv)


# saveRDS(gmts, "data/titer_data/Netzl_et_al_data_for_imprinting_response_binned_gmts.rds")

sr_group_colors_light["D91", "Color"] <- "grey70"
sr_group_colors_light["pre", "Color"] <- "grey70"
    
arm_code_order <- c("pre", arm_code_order)
# Read in the data

source(file.path("code", "00_read_in_data.R"))

BO_D_ajudstment <- c("wo_2dose") #c("wo_2dose", "D85D147")
day_visno <- c("D1", "D29", "D91")

remove_all_breakthrough_sub <- FALSE
mark_undetected_bt <- FALSE

# format the data
data <- format_data(data, sr_group_code, lab = lab, combine_b_o_arm = FALSE, remove_all_bfl = FALSE, remove_all_subj_entries = FALSE,
                    arm_code_by = ifelse(as.numeric(day) >= 181 | date == "18APR2023", "trtactual", "trttrue"), stage = stage, keep_bfls = TRUE, only_target_visno = FALSE, 
                    target_visno = day_visno,
                    bo_adjust = NA, remove_all_oos = FALSE, keep_all_oosboost = TRUE, mark_undetected_BT = mark_undetected_bt)

data_bb <- data

plot_ags <- list("all" = c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5") ,
              "sub" = c("D614G", "BA.1"))


for(name_plot_ags in names(plot_ags)){
  
  
  if(name_plot_ags == "sub"){
    
    # compare only P, P+O & O arm
    data_bb %>%
      filter(arm_code %in% c("O:M", "P+O:M")) -> data_b
    text_angle <- 45
    
    
  } else {
    
    data_b <- data_bb
    text_angle <- 90
    
  }
  
  plot_antigens <- plot_ags[[name_plot_ags]]
  
  ag_order <- seq(1, 2*length(plot_antigens), 2)
  names(ag_order) <- plot_antigens
  
  
  inf_names <- "non_inf"
  
  for(d_adj in BO_D_ajudstment){
    
    # now drop B+O, B+O arm day 29 make copy of B+O arm day 
    data <- BO_format_visit_code(data_b, d_adj)  %>%
      filter(visit_code %in% day_visno)
    
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
    
    
    
    if(remove_all_breakthrough_sub){
      if(date == "18APR2023" & lab == "Montefiori"){
        sanofi_data <- subset_to_visit_subjects(data, c("D15", "D91"))
        
        data <- rbind(sanofi_data, subset_to_visit_subjects(data, day_visno)) %>%
          unique()
      } else {
        data <- subset_to_visit_subjects(data, day_visno)
      }     
    }
    
    if(combine_inf_stat){
      data$inf_code <- "non_inf"
      inf_names <- c(non_inf = "non_inf")
    }
    # add info to group by
    # default grouping is by study arm 
    sr_group_data <- group_data_sr_group(data, by_arm = T, by_visit = T, by_strata = F, by_age = F, by_infection = T)
    sr_group_data$age_code <- "combined"
    
    # make data wide again
    if(correct_magnitude){
      sr_group_data <- adjust_individual_effect(data_to_long(sr_group_data)) %>%
        select(!logtiter:reactivity_bias) %>%
        pivot_wider(names_from = ag_name, values_from = "titer")
    }
    
    sr_group_age <- group_data_sr_group(data, by_arm = T, by_visit = T, by_strata = F, by_age = T, by_infection = T)
    
    # make data wide again
    if(correct_magnitude){
      sr_group_age <- adjust_individual_effect(data_to_long(sr_group_age)) %>%
        select(!logtiter:reactivity_bias) %>%
        pivot_wider(names_from = ag_name, values_from = "titer")
      
    }
    
    sr_group_data <- rbind(sr_group_age, sr_group_data)
    sr_group_data_b <- sr_group_data
    
    # cycle through infected & non infected
    for(inf_stat in inf_names){
      
      figure_dir <- file.path("figures", date, paste0("stage", stage), "additional_figures", ifelse(lab == "Montefiori","Montefiori" ,"."), 
                              "d614g_response_comparison", name_plot_ags)
      suppressWarnings(dir.create(figure_dir, recursive = TRUE))
      
      
      # combine sr groups by age and and filter for inf_stat
      sr_group_data <- sr_group_data_b %>% 
        filter(inf_code == inf_stat)
      
      sr_group_data$arm_code <- factor(sr_group_data$arm_code, levels = arm_code_order)
      sr_group_data$age_code <- factor(sr_group_data$age_code, levels = age_code_order)
      sr_group_data$inf_code <- factor(sr_group_data$inf_code, levels = inf_code_order)
      
      # ---------------------------- V1, V3 comparison -----------------------------------
      
      sr_group_data <- sr_group_data %>% filter(age_code == "combined")
      
      sr_groups <- unique(sr_group_data$sr_group)
      
      #------------------------ Vaccine manufacturer comparison --------------------------------
      vacc_manuf_data <- sr_group_data %>%
        arm_code_wo_vacc() %>%
        filter(v_manuf_code == "M")
      
      
      n_plot_antigens <- length(plot_antigens)/5
      nrow_facet <- length(unique(vacc_manuf_data$v_manuf_code))
      ncol_facet <- length(unique(vacc_manuf_data$arm_code))
      n_days <- length(unique(vacc_manuf_data$visit_code))
      
      plot_width <- 0.5+2*ncol_facet
      unique(vacc_manuf_data$visit_code)
      
      
      forest_data <- forest_data %>%
        mutate(logtiter = log_fold_change,
               ag_name = OmicronVariant,
               arm_code = standardise_encounters,
               inf_code = "non_inf",
               x = ag_order[ag_name],
               visit_code = factor(binned_covail, visit_code_order)) %>%
        filter(!is.na(visit_code))
      
      for(gmt_name in names(gmts)){
        
        if(grepl("pv", gmt_name)){
          
          forest_data_sub_pv <- forest_data %>%
            filter(standardised_assay == "PV")
        } else {
          forest_data_sub_pv <- forest_data
        }
        
        gmt_b <- gmts[[gmt_name]] #%>%
        # mutate(logtiter = gmt,
        #        ag_name = OmicronVariant,
        #        arm_code = standardise_encounters,
        #        inf_code = "non_inf",
        #        x = ag_order[ag_name],
        #        text = paste0("n(", ag_name, " in ", arm_code, ")=",n)) %>%
        # filter(!is.na(visit_code)) %>%
        # unique()
        
        
        ba1 <- "wo"
        
        for(v_target in c("D29", "D91")){
          
          gmt <- gmt_b %>%
            filter(visit_code == v_target) %>%
            mutate(alpha_val = arm_code)
          
          forest_data_sub <- forest_data_sub_pv %>%
            filter(visit_code == v_target) %>%
            mutate(alpha_val = arm_code)
          
          v_data <- vacc_manuf_data %>%
            filter(visit_code %in% c("D1", v_target)) %>%
            mutate(arm_code = ifelse(visit_code == "D1", "pre", arm_code),
                   visit_code = v_target)
          
          p_width <- 2.5
          plot<-titerlineplot_dodge_alpha_pre(v_data,sr_group_colors = sr_group_colors, titer_thresh = titer_threshold, antigens = plot_antigens,
                                    color_by = "arm_code", gmt_facetter = "visit_code",
                                    x_position_by = "age_code", cols_to_keep = c("arm_code", "visit_code", "inf_code",
                                                                                 "age_code", "v_manuf_code"), show_group_count = TRUE,
                                    show_mean_line = F, to_long = T, line_by = "visit_code",
                                    gmt_grid_row = "v_manuf_code", gmt_grid_col = "visit_code", 
                                    dodge_group = "arm_code",
                                    dodge_width = 0.1,
                                    add_group_count_to_legend = FALSE)$gmt  +
            theme(legend.position = "none")
        
          # 2x Vax gmt
          plot2x <- plot +  
            geom_hline(data = gmt %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                         filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = sr_group_colors["P", "Color"]) + 
            geom_text(data = gmt %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                        filter(ag_name == "D614G")%>%
                        unique(), aes(label = text), x =  max(ag_order), y = 1, hjust = 1, size = 2.5, color = sr_group_colors["P", "Color"]) 
          
          plot2x + 
            geom_rect(data = gmt %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                        filter(ag_name == "D614G"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.1, color = NA, fill = sr_group_colors["P", "Color"]) -> plot2x_shade
          
          ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_shade.png")), plot = plot2x_shade, dpi = 300, width = p_width, height = 1.7+2*nrow_facet)
          
          # # 2x Vax Duke      
          # plot <- plot + 
          #   geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("2x Vax")) %>%
          #                filter(ag_name == "D614G") %>%
          #                filter(Study == "Wilks"), aes(yintercept = logtiter), alpha = 0.6, linetype = "solid", color = sr_group_colors["P", "Color"])
          # 
          # ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_duke.png")), plot = plot, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
          # 
          # 3x Vax gmt
          plot3x <- plot2x + 
            geom_hline(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                         filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = "grey40") + 
            geom_text(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                        filter(ag_name == "D614G") %>%
                        unique(), aes(label = text), x =  max(ag_order), y = 0.5, hjust = 1, size = 2.5, color = "grey40")
          
          ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_3xVax.png")), plot = plot3x, dpi = 300, width =2, height = 1.7+2*nrow_facet)
          
          
          plot3x_idvls <- plot3x + 
            geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                         filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = "grey30", alpha = 0.2)
          
          ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_3xVax_idvls.png")), plot = plot3x_idvls, dpi = 300, width = p_width, height = 1.7+2*nrow_facet)
          
          
          plot3x_shade <- plot2x_shade +
            geom_rect(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                        filter(ag_name == "D614G"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.2, color = NA, fill = "grey60") + 
            geom_hline(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                         filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = "grey60") + 
            geom_text(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                        filter(ag_name == "D614G") %>%
                        unique(), aes(label = text), x =  max(ag_order), y = 0.5, hjust = 1, size = 2.5, color = "grey40")
          
          ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_shade_3xVax_shade.png")), plot = plot3x_shade, dpi = 300, width = p_width, height = 1.7+2*nrow_facet)
          
          
          plot3x_shade_only <- plot2x +
            geom_rect(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                        filter(ag_name == "D614G"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.2, color = NA, fill = "grey60") + 
            geom_hline(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                         filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = "grey60") + 
            geom_text(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                        filter(ag_name == "D614G") %>%
                        unique(), aes(label = text), x =  max(ag_order), y = 0.5, hjust = 1, size = 2.5, color = "grey40")
          
          ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_3xVax_shade.png")), plot = plot3x_shade_only, dpi = 300, width = p_width, height = 1.7+2*nrow_facet)
          
          
          # 2x Vax idvls
          plot_2x_idvls <- plot2x + 
            geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                         filter(ag_name == "D614G") %>%
                         filter(!is.na(visit_code)), aes(yintercept = logtiter), alpha = 0.2, linetype = "solid", color = sr_group_colors["P", "Color"])
          ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_idvls.png")), plot = plot_2x_idvls, dpi = 300, width = p_width, height = 1.7+2*nrow_facet)
          
          
          plot_2x_idvls <- plot3x + 
            geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                         filter(ag_name == "D614G") %>%
                         filter(!is.na(visit_code)), aes(yintercept = logtiter), alpha = 0.2, linetype = "solid", color = sr_group_colors["P", "Color"])
          ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_idvls_3x.png")), plot = plot_2x_idvls, dpi = 300, width = p_width, height = 1.7+2*nrow_facet)
          
          plot_2x_idvls <- plot3x + 
            geom_rect(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                        filter(ag_name == "D614G"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.2, color = NA, fill = "grey60") + 
            geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                         filter(ag_name == "D614G") %>%
                         filter(!is.na(visit_code)), aes(yintercept = logtiter), alpha = 0.2, linetype = "solid", color = sr_group_colors["P", "Color"])
          ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_idvls_3xshade.png")), plot = plot_2x_idvls, dpi = 300, width =2, height = 1.7+2*nrow_facet)
          
          if(gmt_name == "binned"){
            
            ba1 <- "w"
            plot_ba1 <- plot3x + 
              geom_hline(data = gmt %>% filter(standardise_encounters %in% c("BA.1 conv")) %>%
                           filter(ag_name == "BA.1"), aes(yintercept = logtiter), linetype = "solid", color = sr_group_colors["O", "Color"]) + 
              geom_text(data = gmt %>% filter(standardise_encounters %in% c("BA.1 conv")) %>%
                          filter(ag_name == "BA.1")%>%
                          unique(), aes(label = text), x =  max(ag_order), y = 0, hjust = 1, size = 2.5, color = sr_group_colors["O", "Color"]) 
            
            plot_ba1shade <- plot_ba1 + 
              geom_rect(data = gmt %>% filter(standardise_encounters %in% c("BA.1 conv")) %>%
                          filter(ag_name == "BA.1"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.1, color = NA, fill = sr_group_colors["O", "Color"])
            
            ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_3xVax_shade.png")), plot = plot_ba1shade, dpi = 300, width = p_width, height = 1.7+2*nrow_facet)
            
            plot_ba1_idvls <- plot_ba1 + 
              geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("BA.1 conv")) %>%
                           filter(ag_name == "BA.1") %>%
                           filter(!is.na(visit_code)), aes(yintercept = logtiter), alpha = 0.2, linetype = "solid", color = sr_group_colors["O", "Color"])
            
            ggsave(file.path(figure_dir, paste0(v_target,"_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_3xVax_idvls.png")), plot = plot_ba1_idvls, dpi = 300, width = p_width, height = 1.7+2*nrow_facet)
            
            
          }
          
          
        }
        
        
        
      }
      
      
    }
    
  }
  
 
  

}

