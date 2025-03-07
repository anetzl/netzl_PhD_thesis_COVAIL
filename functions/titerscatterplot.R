source(file.path("data", "metadata", "plot_specs.R"))

`+.uneval` <- function(a,b) {
  `class<-`(modifyList(a,b), "uneval")
}

map_distance_spacing <- function(data){
  
  map <- read.acmap(file.path("data", "maps", "map_ndsubset_no_outliers_slope_adjusted.ace"))
  
  map_antigens <- agNames(map)[agNames(map) %in% unique(c(data$ag_name,"BA.4/BA.5"))]
 
  map <- subsetMap(map, antigens = map_antigens)
  
  # get Euclidean distance from D614G
  coords <- agCoords(map)
  euc_dist <- dist(coords)
  
  dist_to_D614G <- c(0,dist(coords)[1:4])
  names(dist_to_D614G) <- map_antigens
  scaled_d <- dist_to_D614G/max(dist_to_D614G)
  
  names(scaled_d) <- ag_plot_names[names(scaled_d)]
  names(scaled_d) <- gsub("D614G", "Prototype", names(scaled_d))
  
  data$arm_composition <- arm_plot_names[data$arm_code]
  
  arm_composition <- lapply(data$arm_composition, function(x){
    if(grepl("2x", x)){
      return(c("Beta", "BA.1"))
    } else {
      strsplit(x,"\\+")[[1]]
    }
  })
  
  data$x <- unlist(lapply(arm_composition, function(x){
    mean(scaled_d[x])
  }))
  
  return(data)
  
}


map_distance_spacing_by_ag <- function(data){
  map <- read.acmap(file.path("data", "maps", "map_ndsubset_no_outliers_slope_adjusted.ace"))
  
  map_antigens <- agNames(map)[agNames(map) %in% unique(c(data$ag_name,"BA.4/BA.5"))]
 
  map <- subsetMap(map, antigens = map_antigens)
  
  # get Euclidean distance from D614G
  coords <- agCoords(map)

  euc_dist <- dist(coords)
  euc_dist <- as.matrix(euc_dist)
  scaled_d <- euc_dist/(max(euc_dist))
  
  colnames(scaled_d) <- c(ag_plot_names[c(colnames(scaled_d)[1:4])], "BA.4/5")
  colnames(scaled_d) <- gsub("D614G", "Prototype", colnames(scaled_d))
  
  rownames(scaled_d) <- c(rownames(scaled_d)[1:4], "BA.4/5")

  data$arm_composition <- arm_plot_names[data$arm_code]
  
  arm_composition <- lapply(data$arm_composition, function(x){
    if(grepl("2x", x)){
      return(c("Beta", "BA.1"))
    } else {
      strsplit(x,"\\+")[[1]]
    }
  })
  
  names(arm_composition) <- 1:nrow(data)

  data$x <- unlist(lapply(1:nrow(data), function(x){
    mean(as.numeric(scaled_d[data[[x,"ag_name"]], arm_composition[[x]]]))
  }))
  
  data$x <- data$x/max(data$x)
  
  return(data)
  
}


map_distance_spacing_by_min_ag <- function(data){
  map <- read.acmap(file.path("data", "maps", "map_ndsubset_no_outliers_slope_adjusted.ace"))
  
  map_antigens <- agNames(map)[agNames(map) %in% unique(c(data$ag_name,"BA.4/BA.5"))]
 
  map <- subsetMap(map, antigens = map_antigens)
  
  # get Euclidean distance from D614G
  coords <- agCoords(map)

  euc_dist <- dist(coords)
  euc_dist <- as.matrix(euc_dist)
  scaled_d <- euc_dist/(max(euc_dist))
  
  colnames(scaled_d) <- c(ag_plot_names[c(colnames(scaled_d)[1:4])], "BA.4/5")
  colnames(scaled_d) <- gsub("D614G", "Prototype", colnames(scaled_d))
  
  rownames(scaled_d) <- c(rownames(scaled_d)[1:4], "BA.4/5")

  data$arm_composition <- arm_plot_names[data$arm_code]
  
  arm_composition <- lapply(data$arm_composition, function(x){
    if(grepl("2x", x)){
      return(c("Beta", "BA.1"))
    } else {
      strsplit(x,"\\+")[[1]]
    }
  })
  
  names(arm_composition) <- 1:nrow(data)

  data$x <- unlist(lapply(1:nrow(data), function(x){
    min(as.numeric(scaled_d[data[[x,"ag_name"]], arm_composition[[x]]]))
  }))
  
  data$x <- data$x/max(data$x)
  
  return(data)
  
}


pre_distance_spacing <- function(data){
  

  scaled_D <- read.csv(file.path("data", "metadata", "scaled_distance_from_D614G_pre.csv"))
  scaled_d <- scaled_D[,2]
  names(scaled_d) <- ag_plot_names[scaled_D[,1]]
  names(scaled_d) <- gsub("D614G", "Prototype", names(scaled_d))
  
  data$arm_composition <- arm_plot_names[data$arm_code]
  
  arm_composition <- lapply(data$arm_composition, function(x){
    if(grepl("2x", x)){
      return(c("Beta", "BA.1"))
    } else {
      strsplit(x,"\\+")[[1]]
    }
  })
  
  data$x <- unlist(lapply(arm_composition, function(x){
    mean(scaled_d[x])
  }))
  
  data$x <- data$x/max(data$x)
  
  return(data)
  
}


pre_distance_spacing_by_ag <- function(data){
  

  scaled_D <- read.csv(file.path("data", "metadata", "scaled_fc_from_D614G_pre_all_ags.csv"))
  rownames(scaled_D) <- scaled_D$X
  colnames(scaled_D) <- c("X", rownames(scaled_D))
  scaled_d <- scaled_D[,rownames(scaled_D)]
  
  colnames(scaled_d) <- ag_plot_names[colnames(scaled_d)]
  colnames(scaled_d) <- gsub("D614G", "Prototype", colnames(scaled_d))
  
 # rownames(scaled_d) <- ag_plot_names[rownames(scaled_d)]
 # rownames(scaled_d) <- gsub("D614G", "Prototype", rownames(scaled_d))
  
  data$arm_composition <- arm_plot_names[data$arm_code]
  
  arm_composition <- lapply(data$arm_composition, function(x){
    if(grepl("2x", x)){
      return(c("Beta", "BA.1"))
    } else {
      strsplit(x,"\\+")[[1]]
    }
  })
 
 names(arm_composition) <- 1:nrow(data)

  data$x <- unlist(lapply(1:nrow(data), function(x){
   # print(scaled_d[[x["ag_name"]])
    mean(as.numeric(scaled_d[data[[x,"ag_name"]], arm_composition[[x]]]))
  }))
  
  data$x <- data$x/max(data$x)
  
  return(data)
  
}

pre_distance_spacing_by_ag_lower_higher <- function(data){
  

  scaled_D <- read.csv(file.path("data", "metadata", "scaled_fc_from_D614G_pre_all_ags_lower_higher.csv"))
  rownames(scaled_D) <- scaled_D$X
  colnames(scaled_D) <- c("X", rownames(scaled_D))
  scaled_d <- scaled_D[,rownames(scaled_D)]
  
  colnames(scaled_d) <- ag_plot_names[colnames(scaled_d)]
  colnames(scaled_d) <- gsub("D614G", "Prototype", colnames(scaled_d))
  
 # rownames(scaled_d) <- ag_plot_names[rownames(scaled_d)]
 # rownames(scaled_d) <- gsub("D614G", "Prototype", rownames(scaled_d))

  data$arm_composition <- arm_plot_names[data$arm_code]
  
  arm_composition <- lapply(data$arm_composition, function(x){
    if(grepl("2x", x)){
      return(c("Beta", "BA.1"))
    } else {
      strsplit(x,"\\+")[[1]]
    }
  })
 
 names(arm_composition) <- 1:nrow(data)

  data$x <- unlist(lapply(1:nrow(data), function(x){
   # print(scaled_d[[x["ag_name"]])
    mean(as.numeric(scaled_d[data[[x,"ag_name"]], arm_composition[[x]]]))
  }))
  
   data$x <- data$x*(-1)
  
  return(data)
  
}

x_position_spacing <- function(data, position_by = "arm_code"){
  
  if(position_by == "map_distance"){
    return(map_distance_spacing(data))
  }
  
  if(position_by == "pre_distance"){
    return(pre_distance_spacing(data))
  }

  if(position_by == "antigen"){
    return(pre_distance_spacing_by_ag(data))
  }

  if(position_by == "titer_fc"){
    return(pre_distance_spacing_by_ag_lower_higher(data))
  }

  if(position_by == "map_antigen"){
    return(map_distance_spacing_by_ag(data))
  }

  if(position_by == "min_map_antigen"){
    return(map_distance_spacing_by_min_ag(data))
  }
  
  unique_labels <- unique(data %>% pull(position_by))
  
  unique_pos <- length(unique_labels)
  total_len <- point_spacing*(unique_pos-1)
  
  seq_spaces <- seq(from = -total_len/2, to = total_len/2, by = point_spacing)
  names(seq_spaces) <- unique_labels
  
  data$x <- unlist(lapply(1:nrow(data), function(x){
    target <- as.character(data[[x, position_by]])
    ag_order[data$ag_name[x]] + seq_spaces[target]
  }))
  
  
  return(data)
  
}


titer_differences_scatterplot_dodge <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                                       facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, color_by = "sr_group",
                                       difference_by = "visit_code", x_position_by = "visit_code",  shape_by = "inf_code", line_by = "age_code",
                                       y_label = "Log2 titer difference", cols_to_keep = c(), show_gmt_label = F, show_group_count = F,
                                       show_mean_line = F, mean_line_color = "black",
                                       gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code",
                                       p_norm = FALSE, dodge_group = "arm_code", stat_lm = FALSE) {
  
  
  ag_order <- ag_order[antigens]
  if(TRUE %in% grepl(antigens[1], colnames(data))) {
    # format titer to longer
    data %>%
      pivot_longer(cols = antigens, names_to = "ag_name", values_to = "titer") %>%
      mutate(logtiter = log2(titer/10),
             logtiter = ifelse(titer < titer_thresh, log2((titer_thresh/20)), logtiter)) -> data_long
    
  } else {
    data_long <- data 
  }
  
  # create plot colors
  sr_groups <- unique(as.character(data_long$sr_group))
  
  plot_colors <- unlist(lapply(sr_groups, function(x) {
    blend_sr_group_colors(x, sr_group_colors, split_by = color_by)
  }))
  
  sr_group_colors[sr_groups, "Color"] <- plot_colors
  
  plot_colors <- sr_group_colors$Color
  
  names(plot_colors)<- rownames(sr_group_colors)
  
  # set sr_order
  if(!is.null(sr_group_order)) {
    data_long$sr_group <- factor(as.character(data_long$sr_group), levels =sr_group_order)
  }
  
  
  if(y_label == "Titer"){
    if(p_norm){
      max_y <- 2
    } else {
      max_y <- 12 #12
    }
  } else {
    if(p_norm){
      max_y <- 2
    } else {
      max_y <- 5 #max(ceiling(data_long$logtiter), na.rm = T) +0.5
    }
    
  }
  
  min_y_long <- ifelse(min(data_long$logtiter) < 0, floor(min(data_long$logtiter)), 0) 
  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep) %>%
    calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(titer, 0), "\n", fold_change))
  
  min_y <- ifelse(min(sr_group_gmt_plotdata$logtiter) < 0, ifelse(y_label == "Log2 titer difference", -1.5, -12), 0) #was -10 
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  # add x position
  data_long <- x_position_spacing(data_long, position_by = x_position_by)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by)
  
  sr_group_gmt_plotdata$ag_name <- factor(ag_plot_names[sr_group_gmt_plotdata$ag_name], levels = unique(ag_plot_names))
 
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", "", x)
    x <- gsub("B\\+O, B\\+O", "2x(B+O)", x)
    x
  }
  
  if(length(color_by) > 1){
    color_by <- "sr_group"
  }
  
  # plot gmts of different serum groups together
  
  # max_y <- max(ceiling(sr_group_gmt_plotdata$logtiter), na.rm = T) + 1
  
  spacing <- 0.4
  if(y_label != "Log2 titer difference") {
    spacing <- 0.6
  }
  
  if(length(unique(sr_group_gmt_plotdata$visit_diff)) >1){
    sr_group_gmt_plotdata$visit_diff <- factor(sr_group_gmt_plotdata$visit_diff, levels =  c("D15 - D1", "D29 - D1", "D91 - D1","D29 - D15", "D91 - D15", "D91 - D29"))
  }
  
  sr_group_gmt_plotdata$y_label_pos <- unlist(lapply(1:nrow(sr_group_gmt_plotdata), function(x){
    if(sr_group_gmt_plotdata$inf_code[x] == "non_inf"){
      if(sr_group_gmt_plotdata$age_code[x] == "<65") {
        max_y + spacing/2
      } else {
        max_y - spacing/2
      }
    } else {
      if(sr_group_gmt_plotdata$age_code[x] == "<65") {
        max_y - spacing + spacing/2
      } else {
        max_y - spacing + spacing/2
      }
    }
  }))
  
  if(min(sr_group_gmt_plotdata$x) < 0){
    symm_val <- max(abs(min(sr_group_gmt_plotdata$x)), max(sr_group_gmt_plotdata$x))
    lims <- c(-symm_val - 0.1, symm_val + 0.1)
  }  else {
    lims <- c(-0.1, 1.1)
  }
  
  sr_group_gmt_plotdata %>%
    ggplot(
      aes_string(
        x = "x",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) + 
    geom_hline(yintercept = 0, color = "black") -> gp_gmt

    if(x_position_by == "titer_fc"){
      gp_gmt <- gp_gmt + 
        geom_vline(xintercept = 0, color = "grey50", alpha = 0.3)
    }

    # pointrange if I want them dodged, x position by. If I just want it dodged then I can add position = position_dodge(width = 1)
    # geom_pointrange(
    #   aes(ymin = lower,
    #       ymax = upper)
    # ) +
    # errorbar if I want shifted bars
    gp_gmt + geom_errorbar(
      aes_string(ymin = "lower",
          ymax = "upper",
          group = dodge_group),
      width = 0,
      size = 0.7,
      position = position_dodge(width = 0.1)
    ) +
    geom_point(
      size = 3
    ) +
    geom_point(
      size = 2,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_x_continuous(
      limits = lims,
      name = x_axis_title[x_position_by]
    ) + 
    labs(
      y = "Boost magnitude (fold change)"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(min_y-0.1, max_y+0.1), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      plot.title = element_text(size = 6),
      strip.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) -> gp_gmt
  
  # do the facets
  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
  } else {
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
  }
  
   if(stat_lm){

    gp_gmt <- gp_gmt  +
      stat_smooth(method = "lm") +
      stat_regline_equation(label.x.npc = "left")

      return(gp_gmt)
  }

  if(show_gmt_label){
    gp_gmt <- gp_gmt + 
      geom_text(aes_string(color = color_by) + aes(x = ifelse(age_code == "<65", x + point_spacing/2, x-point_spacing/2), y = y_label_pos,
                                                   label = fold_change),
                color = "black",
                size = 3) 
  }
  
  if(y_label != "Log2 titer difference") {
    gp_gmt <- gp_gmt + scale_y_continuous(
      breaks = seq(floor(min_y), max_y, by = 2),
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x))*10), paste0(round(2^x)*10)),
      name = "Boost magnitude (linear)",
      limits = c(floor(min_y), max_y)
    ) 
  } else {
    
    gp_gmt <- gp_gmt + scale_y_continuous(
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x)), "x"), paste0(round(2^x), "x")),
      breaks= c(round(min_y):round(max_y)),
      limits = c(round(min_y), max_y)
    )  + 
      ylab("Boost magnitude (fold change)")
    
  }
  
  if(shape_by == "inf_code"){
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
          g <- ggplotGrob(gp_gmt)
g$respect <- TRUE
require(grid)
grid.draw(g)
  return(g)
  
}

plot_to_square <- function(gp_gmt){
  g = ggplotGrob(gp_gmt)
g$respect = TRUE
require(grid)
grid.draw(g)

return(g)
}
