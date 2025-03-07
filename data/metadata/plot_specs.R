inf_code_shapes <- c("inf" = 24, "non_inf" = 21, "breakthrough" = 21,
                     "inf_breakthrough" = 24, "non_inf_breakthrough" = 21,
                     "inf_oosboost" = 24, "non_inf_oosboost" = 21)

ag_order <- c("D614G" = 1, "B.1.617.2" = 3, "B.1.351" = 5, "BA.1" = 7, "BA.4/5" = 9, "BA.2.12.1" = 11, "BQ.1.1" = 11, "XBB.1" = 13, "XBB.1.5" = 15)

ag_plot_names <- c("D614G" = "D614G", "B.1.617.2" = "Delta", "B.1.351" = "Beta", "BA.1" = "BA.1", "BA.2.12.1" = "BA.2.12.1", 
    "BA.2.75" = "BA.2.75",
    "BA.2.75.2" = "BA.2.75.2", "BA.4/5" = "BA.4/BA.5",
    "BA.4.6" = "BA.4.6", "BA.4+R346T" = "BA.4+R346T", "BQ.1.1" = "BQ.1.1", "XBB.1" = "XBB.1", "XBB.1.5" = "XBB.1.5")

# below for their arm order
# arm_code_order_no_vacc <- c("P", "P+O","O", "D+O", "B+O", "B+O, B+O",
#                             "P+B", "B")
# arm_code_order <- c("P:M","P+O:M", "O:M", "D+O:M", "B+O:M", "B+O, B+O:M",
#                     "P:Pf","P+O:Pf", "O:Pf", "P+B:Pf", "B+O:Pf", "B:Pf",
#                     "P:S", "P+B:S", "B:S",
#                     arm_code_order_no_vacc)

# arm_plot_names <- c("Prototype", "Prototype + Omicron BA.1","Omicron BA.1", "Delta + Omicron BA.1", "Beta + Omicron BA.1", "2x(Beta + Omicron BA.1)", #Moderna
#                     "Prototype", "Prototype + Omicron BA.1", "Omicron BA.1", "Prototype + Beta", "Beta + Omicron BA.1", "Beta", #Pfizer
#                     "Prototype", "Prototype + Beta", "Beta", #Sanofi
#                     "Prototype", "Prototype + Omicron BA.1","Omicron BA.1", "Delta + Omicron BA.1", "Beta + Omicron BA.1", "2x(Beta + Omicron BA.1)", "Prototype + Beta", "Beta",
#                     "D29-D1", "D91-D1", "D91-D29", "D1", "D29", "D91", "D15", "D57", "Moderna", "Pfizer", "Sanofi") #no vacc manuf
# names(arm_plot_names) <- c(arm_code_order, "D29  D1", "D91  D1", "D91  D29", "D1", "D29", "D91", "D15", "D57", "M", "Pf", "S")

# below for our arm order
arm_code_order_no_vacc <- c("P","O","P+O", "D+O", "B+O", "B+O, B+O",
                            "P+B", "B","P+BA.1", "P+BA.4")
arm_code_order <- c("P:M", "O:M","P+O:M", "D+O:M", "B+O:M", "B+O, B+O:M",
                    "P:Pf", "O:Pf","P+O:Pf", "P+B:Pf", "B+O:Pf", "B:Pf",
                    "P:S", "P+B:S", "B:S","P+BA.1:Pf", "P+BA.4:Pf",
                    arm_code_order_no_vacc)

arm_plot_names <- c("Prototype", "Omicron BA.1","Prototype + Omicron BA.1", "Delta + Omicron BA.1", "Beta + Omicron BA.1", "2x(Beta + Omicron BA.1)", #Moderna
                    "Prototype", "Omicron BA.1","Prototype + Omicron BA.1", "Prototype + Beta", "Beta + Omicron BA.1", "Beta", #Pfizer
                    "Prototype", "Prototype + Beta", "Beta", #Sanofi
                    "Wildtype + Omicron BA.1", "Wildtype + Omicron BA.4/5", # Stage 4
                    "Prototype","Omicron BA.1","Prototype + Omicron BA.1", "Delta + Omicron BA.1", "Beta + Omicron BA.1", "2x(Beta + Omicron BA.1)", "Prototype + Beta", "Beta",
                    "Wildtype + Omicron BA.1", "Wildtype + Omicron BA.4/5", # Stage 4
                    "Pre to 1m", "Pre to 3m", "1m to 3m","Pre to 6m", "1m to 6m", "3m to 6m", "2w to 1m", "6m to 9m", "Pre to 9m", 
                    "Pre to 12m", "6m to 12m", "9m to 12m",
                    "Pre", "1m", "3m", "2w", "D57", "D85", "D147", "6m", "9m", "12m",
                    "Moderna", "Pfizer", "Sanofi",
                    "D614G", "BA.1") #no vacc manuf
names(arm_plot_names) <- c(arm_code_order, "D29  D1", "D91  D1", "D91  D29", "D181  D1", "D181  D29", "D181  D91", "D29  D15", "D271  D181","D271  D1", 
                           "D366  D1", "D366  D181", "D366  D271",
                            "D1", "D29", "D91", "D15", "D57","D85", "D147", "D181", "D271", "D366",
                            "M", "Pf", "S",
                           "D614G", "BA.1")


inf_code_order <- rev(c("inf", "non_inf", "inf_breakthrough", "non_inf_breakthrough", "breakthrough", "inf_oosboost", "non_inf_oosboost", "oosboost"))
age_code_order <- c("<65", ">65", "combined")
visit_code_order <- c("D1","D15", "D29", "D57","D85","D91","D147", "D181", "D271", "D366")

visit_diff_order <- c("D15 - D1", "D29 - D1", "D57 - D1", "D91 - D1", "D181 - D1", "D271 - D1", "D91 - D29", "D181 - D29", "D181 - D91", "D29 - D15",
                      "D271 - D181")

line_vals <- c("<65" = "solid", ">65" = "solid", "combined" = "solid", "NA" = "solid",
"inf" = "solid", "non_inf" = "solid", "breakthrough" = "solid", "non_inf_breakthrough" = "solid", "inf_breakthrough" = "solid",
"non_inf_oosboost" = "solid", "inf_oosboost" = "solid",
"V3" = "solid", "V4" = "solid", "V1" = "solid",
"M" = "solid", "Pf" = "dashed", "S" = "dotted",
"D1" = "solid", "D29"= "solid", "D15"= "solid", "D91"= "solid", "D85"= "solid", "D57"= "solid", "D147" = "solid", 
"D181" = "solid", "D271" = "solid", "D366" = "solid")

point_spacing <- 0.2

x_axis_title <- c("map_distance" = "Map distance from D614G by vaccine composition", 
"pre_distance" = "Pre-titer fold drop from D614G by vaccine composition",
"antigen" = "Pre-titer fold change from variant by vaccine composition",
"map_antigen" = "Scaled map distance from variant by vaccine composition",
"titer_fc" = "Pre-titer 2-fold change from variant by vaccine composition",
"min_map_antigen" = "Scaled map distance from variant by closest vaccine composition strain")

x_axis_tick_size <- 9

visit_plot_order <- seq(1, 2*length(visit_code_order), 2)
names(visit_plot_order) <- visit_code_order

