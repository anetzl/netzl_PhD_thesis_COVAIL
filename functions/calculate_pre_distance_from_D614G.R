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

day <- "91"
data <- read.csv(file.path("data","titer_data", "06OCT2022", paste0("COVAIL_Landscape_data_Monogram_Moderna_D", day,".csv")), sep= ",")
plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5")

data <- format_data(data, sr_group_code, combine_b_o_arm = FALSE)

data_b <- adjust_titers_to_mean_pre(data, plot_antigens)

data_b %>%
  group_data_sr_group(by_arm = T, by_visit = T, by_infection = T, by_age = F) %>%
  filter(inf_code == "non_inf") %>%
  filter(arm_code == "P:M") %>%
  filter(visit_code == "D1") %>%
  data_to_long(., 40, plot_antigens)-> p_data

p_data %>% sr_group_gmt_calc(., 40) %>%
  calculate_fold_change_from_ag() -> p_gmt

scaled_d <- log2(abs(as.numeric(p_gmt$fold_change)))
scaled_d <- scaled_d/max(scaled_d)
names(scaled_d) <- p_gmt$ag_name

write.csv(scaled_d, file.path("data", "metadata", "scaled_distance_from_D614G_pre.csv"))

# for all antigens except BA.4/5
# p_gmt <- p_gmt %>% filter(ag_name != "BA.4/5")

# now calculate matrix of 
fold_changes <- abs(outer(p_gmt$logtiter, p_gmt$logtiter, `-`))
rownames(fold_changes) <- p_gmt$ag_name
colnames(fold_changes) <- p_gmt$ag_name

scaled_fc <- fold_changes/max(fold_changes)

write.csv(scaled_fc, file.path("data", "metadata", "scaled_fc_from_D614G_pre_all_ags.csv"))
# n


# now calculate matrix of 
fold_changes <- outer(p_gmt$logtiter, p_gmt$logtiter, `-`)
rownames(fold_changes) <- p_gmt$ag_name
colnames(fold_changes) <- p_gmt$ag_name

scaled_fc <- fold_changes/max(fold_changes)

write.csv(scaled_fc, file.path("data", "metadata", "scaled_fc_from_D614G_pre_all_ags_lower_higher.csv"))
# n


