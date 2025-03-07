# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(r3js)
library(ablandscapes)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

source(file.path("functions", "format_data.R"))
source(file.path("functions", "sams_landscape_functions.R"))
source(file.path("functions", "landscape_functions.R"))



add_plane_to_data3js <- function(data3js, plane_height, plane_color){
  
  
  # Add the titer 50 plane
  x_grid <- seq(from = plot_xlim[1] + 0.85, to = plot_xlim[2], by = 0.5)
  y_grid <- seq(from = plot_ylim[1] + 0.85, to = plot_ylim[2], by = 0.5)
  z_grid <- matrix(plane_height, length(x_grid), length(y_grid))
  
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = plane_color,
    opacity = 0.2
    #  toggle = "Titer 50",
    # wireframe = FALSE
  )
  
  return(data3js)
}


sr_group_colors["B+O:M", "Color"] <- sr_group_colors["B+O, B+O:M", "Color"]
sr_group_colors["B+O", "Color"] <- sr_group_colors["B+O, B+O:M", "Color"]
sr_group_colors["B+O, B+O:M", "Color"] <- "#e9cf4b"
sr_group_colors["B+O, B+O", "Color"] <- "#e9cf4b"
  
set.seed(100)

#------------ Read in GMTs
gmts_test <- readRDS( "data/titer_data/Netzl_et_al_data_for_imprinting_response_binned_gmts.rds")
gmts <- readRDS( "data/titer_data/Netzl_et_al_data_for_imprinting_response_binned_gmts.rds")

#------------------- Map
recalculate_gmt <- FALSE
only_gmt_landscapes <- TRUE
pre_opacity <- 1
# set optimization nr to fit BA.4/5 position; optimization nr 56 gives upper BA.4/5 for presubmission map; 57 for map with new delta sera
opti_nr <- 1

ag_plot_names <- c("D614G" = "D614G", "B.1.617.2" = "Delta", "B.1.351" = "Beta", "BA.1" = "BA.1", "BA.2.12.1" = "BA.2.12.1", "BA.4/BA.5" = "BA.4/BA.5")
# this map is the one without the new delta sera, October 20202
map <- read.acmap(file.path("data", "maps", "map_ndsubset_no_outliers_slope_adjusted.ace"))
agShown(map)[agNames(map) %in% c("BA.1.1", "BA.1+A484K","BA.3", "BA.2")] <- FALSE
#align pre submission map
map <- realignMap(map, map)


# srOutlineWidth(map) <- 0.5
# lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
# 
# agSize(map)[agNames(map) %in% c("BA.4/BA.5", "BA.2")] <- 12
# agFill(map)[agNames(map) %in% c("BA.4/BA.5", "BA.2")] <- agFill(map)[agNames(map) =="BA.1.1"]
# 
# pdf("figures/06OCT2022/map_pre_sub_o1.pdf", width =6, height = 6)
# plot(map, xlim = lims[[1]], ylim = lims[[2]])
# dev.off()
# 
# pdf("figures/06OCT2022/map_pre_sub_o56.pdf", width =6, height = 6)
# plot(map, optimization_number = 56, xlim = lims[[1]], ylim = lims[[2]])
# dev.off()



map <- keepSingleOptimization(map, optimization_number = opti_nr)
lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
# adjust limits for pre resubmission map
lims$xlim[2] <- lims$xlim[2] + 1
lims$ylim[1] <- lims$ylim[1] + 1

if(opti_nr != 1 & length(map$optimizations) > 1){
  map <- keepSingleOptimization(map, optimization_number = opti_nr)
}

#agNames(map)[agNames(map) == "BA.4/BA.5"] <- "BA.4/5"
#titers_adjusted <- readr::read_csv("./data/titer_data/monogram_titrations_w_montefiori_ba45_ba2121_adjusted.csv")


# angle for html page
angle <- list(
  rotation = c(-1.4647, -0.0093, -0.171),
  translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.45
)

# straight orientation
angle <- list(
  rotation = c(-1.5032, 0.0083, -0.0068),
  translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.45
)

# old D29
angle <- list(
  rotation = c(-1.3971,0.0037,-0.0190),# rotation = c(-1.4370, 0.0062, -0.5350),
  translation = c(0, 0.05,0.1), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.45
)


# old D29
angle <- list(
  rotation = c(-1.4681,0.004,-0.0162),# rotation = c(-1.4370, 0.0062, -0.5350),
  translation = c(0, 0.05,0.1), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.45
)
lndscp_list <- list()
#----------------------- Titer data

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


ag_order <- seq(1, 2*length(plot_antigens), 2)
names(ag_order) <- plot_antigens

for(d_adj in BO_D_ajudstment){
  
  # now drop B+O, B+O arm day 29 make copy of B+O arm day 
  data <- BO_format_visit_code(data_b, d_adj) %>%
    filter(visit_code %in% day_visno)
  
  
  # filter out non responders
  if(lab != "Montefiori"){
    non_responders <- data %>% filter(D614G == 20) %>% filter(B.1.617.2 == 20) %>%
      filter(BA.1 == 20) %>% filter(B.1.351 == 20) %>% filter(`BA.4/5` == 20) %>% pull(subjid)
    
    data <- data %>% filter(!(subjid %in% non_responders))
    warning(paste0(length(unique(non_responders)), " subjects were filtered out because they were non responders"))
    
  }
  
  if(remove_all_breakthrough_sub){
    data <- subset_to_visit_subjects(data, day_visno)
  }
  
  # now combine v1, v3 but group by arm and infection
  sr_group_data <- group_data_sr_group(data, by_arm = T, by_visit = T, by_infection = T, by_age = F) %>%
    arm_code_wo_vacc()
  
  if(only_gmt_landscapes & d_adj %in% c("D57D85")){
    sr_group_data <- sr_group_data %>%
      filter(!(arm_code == "B+O, B+O" & visit_code == "D1"))
  }
  
  # Highlight certain antigens
  highlighted_ags <- c(plot_antigens, "BA.4/BA.5(2)")
  highlighted_ags <- gsub("BA.4/5", "BA.4/BA.5", highlighted_ags)
  arm_cols <- sr_group_colors[c(unique(sr_group_data$arm_code), "M", "Pf", "S"), "Color"]
  names(arm_cols) <- c(unique(sr_group_data$arm_code), "M", "Pf", "S")
  
  sr_group_data_b <- sr_group_data
  
  for(inf_stat in inf_names){
    
    figure_dir <- file.path("figures",date, paste0("stage", stage), inf_stat, "landscapes_imprinting", paste0(d_adj,"_b+o_adj"), ifelse(remove_all_breakthrough_sub, "wo_breakthrough", "w_breakthrough"), paste0("optimization_", opti_nr),
                            ifelse(pre_opacity, "full_pre", "."))
    suppressWarnings(dir.create(figure_dir, recursive = T))
    
    
    
    sr_group_data <- sr_group_data_b %>% 
      filter(inf_code == inf_stat)
    
    titers_adjusted <- sr_group_data
    colnames(titers_adjusted) <- gsub("BA.4/5", "BA.4/BA.5", colnames(titers_adjusted))
    titers_adjusted["BA.4/BA.5(2)"] <- titers_adjusted[,"BA.4/BA.5"]
    
    titerdata <- titers_adjusted
    titerdata %>%
      pivot_longer(
        cols = highlighted_ags,
        names_to = "variant",
        values_to = "titer"
      ) -> titerdata
    
    titerdata %>%
      group_by(
        arm_code,
        visit_code,
        v_manuf_code
      ) -> titerdata
    
    #titerdata <- titerdata %>% filter(!(variant %in% c("BA.4/BA.5", "BA.4/BA.5(2)")))
    titerdata %>%
      group_map(
        get_titertable
      ) -> titertables
    
    lndscp_fits <- lapply(
      titertables,
      function(titertable) {
        
        ablandscape.fit(
          titers = titertable[,c("BA.4/BA.5", "BA.1", "B.1.617.2", "B.1.351", "D614G"), drop = FALSE],
          # titers = titertable[,c("BA.1", "D614G"), drop = FALSE],
          bandwidth = 1,
          degree = 1,
          method = "cone",
          error.sd = 1,
          acmap = map,
          control = list(
            optimise.cone.slope = TRUE
          )
        )
        
      }
    )
    print("after fits")
    
    titertables_groups <- group_data(titerdata)
    
    # Add impulses
    
    if(file.exists(file.path(figure_dir, "titertools_gmt_data.csv")) & !recalculate_gmt){
      gmt_data <- read.csv(file.path(figure_dir, "titertools_gmt_data.csv"))
    } else {
      
      titerdata %>%
        group_by(
          visit_code,
          arm_code,
          v_manuf_code,
          variant
        ) %>%
        summarize(gmt = titertools::gmt(titer, dilution_stepsize = 0)["mean", "estimate"])-> gmt_data
      
      write.csv(x = gmt_data, file = file.path(figure_dir, "titertools_gmt_data.csv"), row.names=FALSE)
    }
    
    print("after gmts")
    # Add landscapes
    data3js_base <- base_plot_data3js(map, lndscp_fits, highlighted_ags, lims, ag_plot_names)
  #  for(v_manuf in unique(titertables_groups$v_manuf_code)){
    
    
    for(v_manuf in "M"){
      
      for(v1 in c("D1")){
        
        
        
        visit_rows <- c(grep("TRUE", v1 == titertables_groups$visit_code))
        manuf_rows <- grep(v_manuf, titertables_groups$v_manuf_code)
        
        #remove B+O, B+O D1 landscape
        
        lndscp_fits_t <- lndscp_fits[intersect(manuf_rows, visit_rows)]
        titertables_groups_t <- titertables_groups[intersect(manuf_rows, visit_rows),]
        
       
        if(v1 != last(day_visno)){
          
          for(v2 in day_visno[c((grep("TRUE", v1 == day_visno)+1):length(day_visno))]) {
            
            # create different data3js bases with gmts and then do plot for each data 3js base  
            data3js_list <- list()
            for(gmt_name in names(gmts)){
              
              gmt <- gmts[[gmt_name]]
             
              data3js_list[[paste0(gmt_name, "_2x")]] <- add_plane_to_data3js(data3js_base, gmt %>% 
                                                   filter(arm_code == "2x Vax") %>%
                                                   filter(ag_name == "D614G") %>%
                                                   filter(visit_code == v2) %>%
                                                    pull(gmt),
                                                 plane_color = sr_group_colors["P", "Color"])
              
              
              data3js_list[[paste0(gmt_name, "_3x")]] <- add_plane_to_data3js(data3js_base, gmt %>% 
                                                                                filter(arm_code == "3x Vax") %>%
                                                                                filter(ag_name == "D614G") %>%
                                                                                filter(visit_code == v2) %>%
                                                                                pull(gmt),
                                                                              plane_color = "grey20")
              
              data3js_list[[paste0(gmt_name, "_BA1")]] <- add_plane_to_data3js(data3js_base, gmt %>% 
                                                                                filter(arm_code == "BA.1 conv") %>%
                                                                                filter(ag_name == "BA.1") %>%
                                                                                filter(visit_code == v2) %>%
                                                                                pull(gmt),
                                                                              plane_color = sr_group_colors["O", "Color"])
              
              
              data3js_list[[paste0(gmt_name, "_2x3x")]] <- add_plane_to_data3js(data3js_list[[paste0(gmt_name, "_2x")]], gmt %>% 
                                                                                  filter(arm_code == "3x Vax") %>%
                                                                                  filter(ag_name == "D614G") %>%
                                                                                  filter(visit_code == v2) %>%
                                                                                  pull(gmt),
                                                                                plane_color = "grey20")
              
              
              data3js_list[[paste0(gmt_name, "_2x3xBA1")]] <- add_plane_to_data3js(data3js_list[[paste0(gmt_name, "_2x3x")]], gmt %>% 
                                                                                     filter(arm_code == "BA.1 conv") %>%
                                                                                     filter(ag_name == "BA.1") %>%
                                                                                     filter(visit_code == v2) %>%
                                                                                     pull(gmt),
                                                                                   plane_color = "darkred")
              
              
              
            }
            
            
            # get matching lndscp_fits
            visit_rows <- c(grep("TRUE", v1 == titertables_groups$visit_code), grep("TRUE", v2 == titertables_groups$visit_code))
            manuf_rows <- grep(v_manuf, titertables_groups$v_manuf_code)
            
            #remove B+O, B+O D1 landscape
            
            lndscp_fits_t <- lndscp_fits[intersect(manuf_rows, visit_rows)]
            titertables_groups_t <- titertables_groups[intersect(manuf_rows, visit_rows),]
            
            for(data_name in names(data3js_list)){
              
              data3js <- data3js_list[[data_name]]
              lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups_t, lndscp_fits_t, map, gmt_data, highlighted_ags, ag_plot_names)
              
              
              lndscp <-r3js(
                lndscp_3js,
                rotation = angle$rotation,
                zoom = angle$zoom
              )
              
              #lndscp_list[[paste0(inf_stat, "_",v1, "_",v2, "_", v_manuf)]] <- lndscp
              save_name <- file.path(figure_dir,paste0(data_name, "_", v1, "_", v2, "_", v_manuf,"_derek_gmt_landscapes_pre_styling"))
              screenshot_html_landscape_to_png(lndscp, save_name)
              
              
              
            }
            
           
          
          }
          
        }
        
      }
      
    }
    
   
    
    
    
  }
}

