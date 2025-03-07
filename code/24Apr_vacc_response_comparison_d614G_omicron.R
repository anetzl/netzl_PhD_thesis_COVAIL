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


# saveRDS(gmts, "data/titer_data/Netzl_et_al_data_for_imprinting_response_binned_gmts.csv")

sr_group_colors_light["D91", "Color"] <- "grey70"
# Read in the data

source(file.path("code", "00_read_in_data.R"))

BO_D_ajudstment <- c("wo_2dose") #c("wo_2dose", "D85D147")
day_visno <- c("D1", "D29", "D91", "D181")

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
        arm_code_wo_vacc()
      
      
      n_plot_antigens <- length(plot_antigens)/5
      nrow_facet <- length(unique(vacc_manuf_data$v_manuf_code))
      ncol_facet <- length(unique(vacc_manuf_data$arm_code))
      n_days <- length(unique(vacc_manuf_data$visit_code))
      
      plot_width <- 0.5+2*ncol_facet
      unique(vacc_manuf_data$visit_code)
      
      # x_position_by = "v_manuf_code"
      # titerlineplot_dodge(vacc_manuf_data,sr_group_colors = sr_group_colors, titer_thresh = titer_threshold, antigens = plot_antigens,
      #                     color_by = c("v_manuf_code"),
      #                     x_position_by = "age_code", cols_to_keep = c("arm_code", "visit_code", "inf_code",
      #                                                                  "age_code", "v_manuf_code"), show_group_count = F,
      #                     show_mean_line = F, to_long = T, line_by = "visit_code",
      #                     gmt_grid_row = "visit_code", gmt_grid_col = "arm_code", dodge_group = "v_manuf_code")$gmt + theme(legend.position = "none") -> plot
      # plot
      # stop()
      #   geom_hline(yintercept = d614g_gmt$logtiter, color = "grey40", linetype = "dashed")-> plot
      # 
      # ggsave(file.path(figure_dir, paste0("V_by_manuf_gmts_age_all_by_arm_visit.png")), plot = plot, dpi = 300, width = plot_width, height = 2+2*n_days)
      
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
          
          forest_data_sub <- forest_data %>%
            filter(standardised_assay == "PV")
        } else {
          forest_data_sub <- forest_data
        }
        
        gmt <- gmts[[gmt_name]] #%>%
        # mutate(logtiter = gmt,
        #        ag_name = OmicronVariant,
        #        arm_code = standardise_encounters,
        #        inf_code = "non_inf",
        #        x = ag_order[ag_name],
        #        text = paste0("n(", ag_name, " in ", arm_code, ")=",n)) %>%
        # filter(!is.na(visit_code)) %>%
        # unique()
        
        
        ba1 <- "wo"
        
        plot<-titerlineplot_dodge(vacc_manuf_data,sr_group_colors = sr_group_colors, titer_thresh = titer_threshold, antigens = plot_antigens,
                                  color_by = "arm_code", gmt_facetter = "visit_code",
                                  x_position_by = "age_code", cols_to_keep = c("arm_code", "visit_code", "inf_code",
                                                                               "age_code", "v_manuf_code"), show_group_count = TRUE,
                                  show_mean_line = F, to_long = T, line_by = "visit_code",
                                  gmt_grid_row = "v_manuf_code", gmt_grid_col = "visit_code", 
                                  dodge_group = "arm_code",
                                  dodge_width = 0.1,
                                  add_group_count_to_legend = FALSE)$gmt + theme(legend.position = "none",
                                                                                 axis.text.x = element_text(angle = text_angle))
        
        # 2x Vax gmt
        plot2x <- plot +  
          geom_hline(data = gmt %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                       filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = sr_group_colors["P", "Color"]) + 
          geom_text(data = gmt %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                      filter(ag_name == "D614G")%>%
                      unique(), aes(label = text), x = ag_order["BA.1"], y = 1, hjust = 1, size = 2.5, color = sr_group_colors["P", "Color"]) 
        
        plot2x + 
          geom_rect(data = gmt %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                      filter(ag_name == "D614G"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.1, color = NA, fill = sr_group_colors["P", "Color"]) -> plot2x_shade
        
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_shade.png")), plot = plot2x_shade, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
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
                       filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = "grey60") + 
          geom_text(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                      filter(ag_name == "D614G") %>%
                      unique(), aes(label = text), x = ag_order["BA.1"], y = 0.5, hjust = 1, size = 2.5, color = "grey40")
        
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_3xVax.png")), plot = plot3x, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        
        plot3x_idvls <- plot3x + 
          geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                       filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = "grey60", alpha = 0.2)
        
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_3xVax_idvls.png")), plot = plot3x_idvls, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        
        plot3x_shade <- plot2x_shade +
          geom_rect(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                      filter(ag_name == "D614G"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.2, color = NA, fill = "grey60") + 
          geom_hline(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                       filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = "grey60") + 
          geom_text(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                      filter(ag_name == "D614G") %>%
                      unique(), aes(label = text), x = ag_order["BA.1"], y = 0.5, hjust = 1, size = 2.5, color = "grey40")
        
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_shade_3xVax_shade.png")), plot = plot3x_shade, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        
        plot3x_shade_only <- plot2x +
          geom_rect(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                      filter(ag_name == "D614G"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.2, color = NA, fill = "grey60") + 
          geom_hline(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                       filter(ag_name == "D614G"), aes(yintercept = logtiter), linetype = "solid", color = "grey60") + 
          geom_text(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                      filter(ag_name == "D614G") %>%
                      unique(), aes(label = text), x = ag_order["BA.1"], y = 0.5, hjust = 1, size = 2.5, color = "grey40")
        
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_3xVax_shade.png")), plot = plot3x_shade_only, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        
        # 2x Vax idvls
        plot_2x_idvls <- plot2x + 
          geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                       filter(ag_name == "D614G") %>%
                       filter(!is.na(visit_code)), aes(yintercept = logtiter), alpha = 0.2, linetype = "solid", color = sr_group_colors["P", "Color"])
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_idvls.png")), plot = plot_2x_idvls, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        
        plot_2x_idvls <- plot3x + 
          geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                       filter(ag_name == "D614G") %>%
                       filter(!is.na(visit_code)), aes(yintercept = logtiter), alpha = 0.2, linetype = "solid", color = sr_group_colors["P", "Color"])
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_idvls_3x.png")), plot = plot_2x_idvls, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        plot_2x_idvls <- plot3x + 
          geom_rect(data = gmt %>% filter(standardise_encounters %in% c("3x Vax")) %>%
                      filter(ag_name == "D614G"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.2, color = NA, fill = "grey60") + 
          geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("2x Vax")) %>%
                       filter(ag_name == "D614G") %>%
                       filter(!is.na(visit_code)), aes(yintercept = logtiter), alpha = 0.2, linetype = "solid", color = sr_group_colors["P", "Color"])
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_idvls_3xshade.png")), plot = plot_2x_idvls, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        ba1 <- "w"
        plot_ba1 <- plot3x + 
          geom_hline(data = gmt %>% filter(standardise_encounters %in% c("BA.1 conv")) %>%
                       filter(ag_name == "BA.1"), aes(yintercept = logtiter), linetype = "solid", color = sr_group_colors["O", "Color"]) + 
          geom_text(data = gmt %>% filter(standardise_encounters %in% c("BA.1 conv")) %>%
                      filter(ag_name == "BA.1")%>%
                      unique(), aes(label = text), x = ag_order["BA.1"], y = 0, hjust = 1, size = 2.5, color = sr_group_colors["O", "Color"]) 
        
        plot_ba1shade <- plot_ba1 + 
          geom_rect(data = gmt %>% filter(standardise_encounters %in% c("BA.1 conv")) %>%
                      filter(ag_name == "BA.1"), aes(ymin = lower, ymax = upper, xmin = -Inf, xmax = Inf), alpha = 0.1, color = NA, fill = sr_group_colors["O", "Color"])
        
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_3xVax_shade.png")), plot = plot_ba1shade, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        plot_ba1_idvls <- plot_ba1 + 
          geom_hline(data = forest_data_sub %>% filter(standardise_encounters %in% c("BA.1 conv")) %>%
                       filter(ag_name == "BA.1") %>%
                       filter(!is.na(visit_code)), aes(yintercept = logtiter), alpha = 0.2, linetype = "solid", color = sr_group_colors["O", "Color"])
        
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_",ba1, "_ba1_", gmt_name, "_2xVax_3xVax_idvls.png")), plot = plot_ba1_idvls, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        
        next()
        
        
        ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_visit_", gmt_name, "_ba1_2xVax.png")), plot = plot, dpi = 300, width = 0.5+2*n_days, height = 2+2*nrow_facet)
        
        if(do_x_axis_time_plots){
          
          gmt %>%
            filter(arm_code %in% c("BA.1 conv", "2x Vax")) %>% 
            ggplot(aes(x = visit_code, y = logtiter, color = arm_code, group = arm_code)) + 
            geom_line() + 
            geom_point() +
            facet_wrap(~ag_name)
          
          
          visit_code_plot <- seq(1, 2*length(visit_code_order), 2)
          names(visit_code_plot) <- visit_code_order
          
          # here facet by antigen, x is visit code
          
          gmt$x <- as.numeric(gsub("D", "", gmt$visit_code))
          
          gmt_sub <- rbind(gmt %>% filter(arm_code %in% c("2x Vax", "WT conv")) %>%
                             filter(ag_name %in% c("D614G")),
                           gmt %>% filter(arm_code %in% c("BA.1 conv")) %>%
                             filter(ag_name %in% c("BA.1")))
          
          sr_group_colors["2x Vax", "Color"] <- "skyblue"
            sr_group_colors["WT conv", "Color"] <- "darkblue"
              sr_group_colors["BA.1 conv", "Color"] <- "pink"
                
              temp_colors <- sr_group_colors$Color
              names(temp_colors) <- rownames(sr_group_colors)
              
              color_labels <- c("D614G in 2x Vax","BA.1 in BA.1 infected","BA.1 in BA.1 Vax COVAIL","BA.1 in Prototype COVAIL","BA.1 in BA.1/Prototype COVAIL","D614G in Ancestral infected")
              names(color_labels) <- c("2x Vax", "BA.1 conv", "O", "P", "P+O", "WT conv")
              
              plot<-titerlineplot_dodge(vacc_manuf_data,sr_group_colors = sr_group_colors, titer_thresh = titer_threshold, antigens = "BA.1",
                                        color_by = "arm_code", gmt_facetter = "NA",
                                        x_position_by = "arm_code", cols_to_keep = c("arm_code", "visit_code", "inf_code",
                                                                                     "age_code", "v_manuf_code"), 
                                        show_group_count = FALSE,
                                        show_mean_line = F, to_long = T, line_by = "visit_code",
                                        gmt_grid_row = "v_manuf_code", gmt_grid_col = "ag_name", 
                                        dodge_group = "arm_code",
                                        dodge_width = 0.1,
                                        add_group_count_to_legend = FALSE,
                                        x_axis_var = "visit_code",
                                        plot_grouping_var = "arm_code")$gmt + theme(legend.position = "none")
              
              plot + 
                geom_line(data = gmt_sub, aes(group = arm_code), position = position_dodge(width = 0.1)) + 
                geom_errorbar(data = gmt_sub,
                              aes(ymin = lower,
                                  ymax = upper,
                                  group = "arm_group"), # usually visit code
                              width = 0,
                              position = position_dodge(width = 0.1)
                ) + 
                geom_point(data = gmt_sub, aes(group = arm_code), position = position_dodge(width = 0.1)) + 
                guides(shape = "none", fill = "none")+
                scale_color_manual(values = temp_colors, name = "Data", labels = color_labels)+
                theme(legend.position = "right",
                      strip.text.x = element_blank()) -> plot
              
              ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_ag_", gmt_name, ".png")), plot = plot, dpi = 300, width = 7, height = 2+2*nrow_facet)
              
              
              
              
              plot<-titerlineplot_dodge(vacc_manuf_data,sr_group_colors = sr_group_colors, titer_thresh = titer_threshold, antigens = "BA.1",
                                        color_by = "arm_code", gmt_facetter = "NA",
                                        x_position_by = "arm_code", cols_to_keep = c("arm_code", "visit_code", "inf_code",
                                                                                     "age_code", "v_manuf_code"), 
                                        show_group_count = FALSE,
                                        show_mean_line = F, to_long = T, line_by = "arm_code",
                                        gmt_grid_row = "v_manuf_code", gmt_grid_col = "ag_name", 
                                        dodge_group = "arm_code",
                                        dodge_width = 0.1,
                                        add_group_count_to_legend = FALSE,
                                        x_axis_var = "visit_code",
                                        plot_grouping_var = "arm_code")$gmt + theme(legend.position = "none")
              plot + 
                geom_errorbar(data = gmt_sub,
                              aes(ymin = lower,
                                  ymax = upper,
                                  group = "arm_group"), # usually visit code
                              width = 0,
                              position = position_dodge(width = 0.1)
                ) + 
                geom_point(data = gmt_sub, aes(group = arm_code), position = position_dodge(width = 0.1)) + 
                guides(shape = "none", fill = "none")+
                scale_color_manual(values = temp_colors, name = "Data", labels = color_labels)+
                theme(legend.position = "right",
                      strip.text.x = element_blank()) -> plot
              
              ggsave(file.path(figure_dir, paste0("V_gmts_age_all_by_ag_", gmt_name, "_no_line.png")), plot = plot, dpi = 300, width = 7, height = 2+2*nrow_facet)
              
              
        }
        
        
        
      }
      
      # stop()
      # 
      # 
      # 
      # plot<-titerlineplot_dodge(vacc_manuf_data %>% filter(visit_code %in% day_visno),sr_group_colors = sr_group_colors_light, titer_thresh = titer_threshold, antigens = plot_antigens,
      #                           color_by = c("arm_code", "visit_code"), gmt_facetter = "arm_code",
      #                           x_position_by = "age_code", cols_to_keep = c("arm_code", "visit_code", "inf_code",
      #                                                                        "age_code", "v_manuf_code"), show_group_count = T,
      #                           show_mean_line = F, to_long = T, line_by = "visit_code",
      #                           gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code", dodge_group = "visit_code",
      #                           dodge_width = 0.1,
      #                           add_group_count_to_legend = TRUE)$gmt
      # plot
      # +
      #   geom_hline(yintercept = d614g_gmt$logtiter, color = "grey40", linetype = "dashed") + 
      #   guides(shape = "none") +
      #   theme(legend.position = "top")
      # 
      # ggsave(file.path(figure_dir, paste0("V_gmts_age_all_visit_by_arm.png")), plot = plot, dpi = 300, width = plot_width, height = 2.5+2*nrow_facet)
      # 
      
    }
    
  }
  
 
  

}





# Read in xbb.1.5 data
xbb_data <- read.csv("data/titer_data/CONFIDENTIAL_samT_placeholder_GMT_titers.csv") %>%
  mutate(logtiter = log2(titer/10))

xbb_data %>%
  filter(lab == "moderna") %>%
  filter(ag %in% c("d614g", "xbb.1.5")) %>%
  filter(!grepl("all", sr)) %>%
  mutate(pre_titer = grepl("pre", sr),
         prior_inf = ifelse(grepl("yes", sr), "Prior infection", "No prior infection")) -> xbb_data

xbb_data %>%
  ggplot(aes(x = ag, y = logtiter, color = pre_titer, group = sr)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) + 
  geom_point(color = "white", size = 1) +
  geom_hline(yintercept = d614g_gmt$logtiter, color = "grey40", linetype = "dashed") + 
  scale_y_continuous(labels = function(x) 2^x*10,
                     name = "Titer",
                     breaks = 0:13,
                     limits = c(0, 13)) + 
  scale_x_discrete(labels = c("d614g" = "D614G", "xbb.1.5" = "XBB.1.5"),
                   name = "Antigen variant",
                   expand = c(0, .2)) +
  scale_color_manual(values = c("TRUE" = "#ffa500", "FALSE" = "#ff6500"),
                     labels = c("TRUE" = "Pre", "FALSE" = "2w"),
                     name = "") +
  facet_wrap(~prior_inf) +
  titerplot_theme() +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 6),
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = x_axis_tick_size),
    axis.text.y = element_text(size = 10)
  ) -> p

ggsave(file.path(figure_dir, paste0("Moderna_xbb15_vax.png")), plot = p, dpi = 300, width = 5, height = 2.5+2)



