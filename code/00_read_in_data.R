# all specifications to read in the data. Are used in several files

correct_magnitude <- F
remove_all_breakthrough_sub <- FALSE
do_age_split_plots <- FALSE
only_p_arm_visit_comparison_plots <- TRUE
two_visit_comparison_plots <- FALSE

# combine inf stat for stage 4
stage <- 1
combine_inf_stat <- ifelse(stage == 1, FALSE, TRUE)
date <- "26SEP2024"
day <- "366"

inf_names <- c(non_inf = "non_inf", inf = "inf")


BO_D_ajudstment <- c("wo_2dose") #c("wo_2dose", "D85D147", "D85D147D237")
day_visno <- c("D1","D15", "D29", "D91")

lab <- "Monogram"
titer_threshold <- 40

# # # uncomment below for stage 4
# stage <- 4
# BO_D_ajudstment <-  c("wo") #c("wo_2dose", "D85D147")
# day_visno <- c("D1","D29", "D91")
# combine_inf_stat <- TRUE
# day <- "91"
# lab <- "Montefiori"


# read the data

# This is date with all time points, should be the most complete data set
if(date == "26SEP2024"){
  
  undetected_bt <- rbind(read_csv(file.path("data", "titer_data", "18APR2023", "non_inf_undetected_BT.csv")),
                         read_csv(file.path("data", "titer_data", "18APR2023", "inf_undetected_BT.csv"))) %>%
    filter(subjid != "COV.00545") # filter out delta breakthrough sample, could be an assay issue at 6m
  
  plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5") #BA.2.12.1 only for D1, D15 and Moderna arm
  if(stage == 1){
    
    # Monogram read in
    if(lab == "Monogram"){
      
     
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Moderna.csv")), sep= ",") %>%
        plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Pfizer.csv")), sep= ",")) %>%
        plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Sanofi.csv")), sep= ",")) 
      
      
    } else {
      # Duke lab here
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Duke_Moderna.csv")), sep= ",") %>%
        plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Duke_Sanofi.csv")), sep= ",")) 
      
      titer_threshold <- 10
      
      plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.2.12.1", "BA.2.75", "BA.2.75.2", "BA.4/5", "BA.4.6", "BA.4+R346T", "BQ.1.1", "XBB.1", "XBB.1.5")
    }
  
  } else {
    # Stage 4 here
    
    if(lab == "Monogram"){
      
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Stage4.csv")), sep= ",")
      
    } else {
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Duke_Stage4.csv")), sep= ",") 
      
      titer_threshold <- 10
      
      plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.2.12.1", "BA.2.75", "BA.2.75.2", "BA.4/5", "BA.4.6", "BA.4+R346T", "BQ.1.1", "XBB.1", "XBB.1.5")
      
    }
    
  }
 
  
} else if(date == "29JUN2022"){
  data <- read.csv(file.path("data","titer_data", date, "COVAIL_Landscape_data_Moderna_D15.csv"), sep= ",")
  plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1")
  
} else if(date %in% c("19SEP2022", "14DEC2022", "22MAR2023", "18APR2023")) {
  
  if(stage == 1){
    if(lab == "Montefiori"){
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Duke_Moderna_D", day,".csv")), sep= ",")
      
      if(date == "18APR2023"){
        data <- data %>%
          plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Duke_Sanofi_D15_D", day,".csv")), sep= ",")) 
      }
      
      titer_threshold <- 10
      
      plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.2.12.1", "BA.2.75", "BA.2.75.2", "BA.4/5", "BA.4.6", "BA.4+R346T", "BQ.1.1", "XBB.1", "XBB.1.5")
      
    } else {
      
      if(day == "271"){
        data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Moderna_D", day,".csv")), sep= ",") 
        
        undetected_bt <- rbind(read_csv(file.path("data", "titer_data", date, "non_inf_undetected_BT.csv")),
                               read_csv(file.path("data", "titer_data", date, "inf_undetected_BT.csv"))) %>%
          filter(subjid != "COV.00545") # filter out delta breakthrough sample, could be an assay issue at 6m
                               
        
      } else {
        data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Moderna_D", ifelse(date == "18APR2023", "271", day),".csv")), sep= ",") %>%
          plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Pfizer_D", day,".csv")), sep= ",")) %>%
          plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Sanofi_D", day,".csv")), sep= ",")) 
        
        data %>%
          filter(VISNO != 7) -> data
        
        undetected_bt <- rbind(read_csv(file.path("data", "titer_data", date, "non_inf_undetected_BT.csv")),
                               read_csv(file.path("data", "titer_data", date, "inf_undetected_BT.csv"))) %>%
          filter(subjid != "COV.00545") # filter out delta breakthrough sample, could be an assay issue at 6m
        
      }
      
      
      plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5", "BA.2.12.1")
     # plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5")
      
    }
  } else {
    
    if(lab == "Montefiori"){
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Duke_Stage4_D", day,".csv")), sep= ",")
      titer_threshold <- 10
      plot_antigens <- c("D614G", "BA.1", "BA.4/5", "BQ.1.1", "XBB.1")
      
    } else {
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Stage4_D", day,".csv")), sep= ",")
      plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5")
    }
    
  }
  
  
} else {
  if(stage == 1){
    data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Moderna_D", day,".csv")), sep= ",")
    
    plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5")
  } else {
    if(lab == "Montefiori"){
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Duke_Stage4_D15.csv")), sep= ",") 
      day_visno <- c("D1", "D15")
      plot_antigens <- c("D614G", "BA.1", "BA.4/5", "BQ.1.1", "XBB.1")
      titer_threshold <- 10
    } else {
      data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Stage4_D29.csv")), sep= ",")
      plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5")
    }
    
  }
  
}
