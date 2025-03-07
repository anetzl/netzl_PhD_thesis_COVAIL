source(file.path("data", "metadata", "plot_specs.R"))

PLOTLY <- FALSE

`+.uneval` <- function(a,b) {
  `class<-`(modifyList(a,b), "uneval")
}

facet_labeller <- function(x){
  x <- gsub("NA-", "", x)
  x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
  x <- gsub("-", " ", x)
  x <- arm_plot_names[x]
  x
}

x_by_visit_diff <- function(data){
  
  visits <- unlist(lapply(data$visit_diff, function(x){
    d <- strsplit(x, " - ")[[1]][1]
    d <- gsub("D", "")
    as.numeric(d)
  }))
  
  data$x <- visits

  return(data)
}

x_position_spacing <- function(data, position_by = "arm_code", x_position_var = "ag_name"){
  
  if(position_by == "visit_diff"){
    return(x_by_visit_diff(data))
  }
  unique_labels <- unique(data %>% pull(position_by))
  
  unique_pos <- length(unique_labels)
  total_len <- point_spacing*(unique_pos-1)
  
  seq_spaces <- seq(from = -total_len/2, to = total_len/2, by = point_spacing)
  names(seq_spaces) <- unique_labels
#  print(visit_plot_order)
  if(x_position_var == "visit_code"){
    # data$x <- unlist(lapply(1:nrow(data), function(x){
    #   target <- as.character(data[[x, position_by]])
    #   
    #   #print(visit_plot_order[as.character(data$visit_code[x])])
    #   visit_plot_order[as.character(data$visit_code[x])] + seq_spaces[target]
    # }))
    
    data$x <- as.numeric(gsub("D", "", as.character(data$visit_code)))
  } else {
    data$x <- unlist(lapply(1:nrow(data), function(x){
      target <- as.character(data[[x, position_by]])
      ag_order[data$ag_name[x]] + seq_spaces[target]
    }))
  }
  
  
  
  return(data)
  
}

titerlineplot <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                          facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, 
                          color_by = "sr_group", shape_by = "inf_code", line_by = "age_code", x_position_by = "age_code",
                          cols_to_keep = c("arm_code", "age_code", "inf_code"), to_long = T, 
                          show_gmt_label = F, show_group_count = F, show_mean_line = F, mean_line_color = "black",
                          gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code",
                          pre_titer_adj = FALSE, single_visit = NULL, highlight_samples = c()) {
  
  ag_order <- ag_order[antigens]
  # format titer to longer
  if(to_long){
    data %>%
      pivot_longer(cols = antigens, names_to = "ag_name", values_to = "titer") %>%
      mutate(logtiter = log2(titer/10),
             logtiter = ifelse(titer < titer_thresh, log2((titer_thresh/20)), logtiter)) -> data_long
  } else{
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

  
 # plot_colors <-  plot_colors[as.character(unique(data_long[, color_by]))]
 
  # set sr_order
  if(!is.null(sr_group_order)) {
    data_long$sr_group <- factor(as.character(data_long$sr_group), levels =sr_group_order)
  }
  

#  max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
  max_y <- 15

  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep, pre_adj = pre_titer_adj) %>%
    filter(!all_below_thresh) -> sr_group_gmt_plotdata 

  if(!is.null(single_visit)){
    data_long <- data_long %>% filter(visit_code %in% single_visit)
    sr_group_gmt_plotdata %>% filter(visit_code %in% single_visit)
  }
   
   sr_group_gmt_plotdata <- sr_group_gmt_plotdata %>% calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(as.numeric(titer), 0), "\n", fold_change))
  
  
 # sr_group_gmt_plotdata <- sr_group_to_column(sr_group_gmt_plotdata)
  
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  data_long <- x_position_spacing(data_long, position_by = x_position_by)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by)
  
  # fill  NA values of lower and upper CI with mean value
  sr_group_gmt_plotdata$lower[is.na(sr_group_gmt_plotdata$lower)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$lower)]
  sr_group_gmt_plotdata$upper[is.na(sr_group_gmt_plotdata$upper)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$upper)]
  
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", " ", x)
    x <- arm_plot_names[x]
    x
  }

  # facet_labeller <- function(x){
  #   strsplit(x, "-")[[1]][2]
  # }

  sub_visit <- FALSE
  if(length(color_by) > 1){
    if("visit_code" %in% color_by | length(unique(data$visit_code)) > 1){
      sub_visit <- TRUE
    }
    color_by <- "sr_group"
  }
  
  data_long$visit_code <- factor(data_long$visit_code, levels = visit_code_order)
  sr_group_gmt_plotdata$visit_code <- factor(sr_group_gmt_plotdata$visit_code, levels = visit_code_order)

  data_long %>%
    ggplot(
      aes_string(
        x = "ag_name",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) + 
    geom_line(
      aes(group = sr_name),
      linetype = "solid",
      alpha = 0.4
    ) + 
    geom_point(
      aes(group = sr_name),
      alpha = 0.4
    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
      aes_string(group = "sr_group",
          linetype = line_by),
      size = 1.3
    ) +
    geom_pointrange(
      data = sr_group_gmt_plotdata,
      aes(ymin = lower,
            ymax = upper)
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = ag_name, y = y, label = label),
              color = "black",
              size = 2.5
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
    scale_x_discrete(
      limits = antigens,
      labels = ag_plot_names[antigens],
      name = "Antigen variant"
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+1), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 10)
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    ) -> gp
  
  if(length(highlight_samples) > 0){
    
    data_long %>%
      filter(subjid %in% highlight_samples) -> data_highlight
    
    gp <- gp + 
      geom_line(data = data_highlight, aes(group = sr_name), color = "red", alpha = 1)
    
  }
  
  visits <- unique(sr_group_gmt_plotdata$visit_code)
  
  if(line_by != "v_manuf_code" & "V3" %in% visits & "V4" %in% visits){
    line_by <- "visit_code"
  }
  # plot gmts of different serum groups together
  max_y <- 13
  sr_group_gmt_plotdata$y <- max_y #- 0.5
  
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
    geom_line(
      aes_string(group = "sr_group",
          linetype = line_by),
      size = 1
    ) +
    geom_pointrange(
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
    scale_x_continuous(
      breaks = ag_order,
      labels = ag_plot_names[names(ag_order)],
      limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(ag_order, na.rm = T))+0.5),
      name = "Antigen variant"
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+0.5), 0.05)
    ) +
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
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    ) -> gp_gmt
  

  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
    gp <- gp + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"sr_group", "arm_code", "v_manuf_code", "count", "y", "age_code", "visit_code", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    if(sub_visit){
      visits <- unique(sr_group_gmt_plotdata$visit_code)
      
      target_visit <- visits[visits != "D15"]
     
     #   count_label <- count_label[grepl(target_visit[1], count_label$sr_group),]
      count_label <- count_label[target_visit[1] == count_label$visit_code,]
    }
    
    count_label <- count_label %>%
      group_by(sr_group,arm_code, v_manuf_code, age_code) %>%
      summarize(total_count = sum(count),
                label = paste0("n (", age_code, ") = ", total_count),
                x = x,
                y = ifelse(age_code == "combined", y, ifelse(age_code == "<65", y - 1, y-2)),
                inf_code = inf_code)
    
    count_label <- count_label %>%
      mutate(y = ifelse(inf_code == "inf", y-0.5, y-1))
    
    count_label$label <- gsub("\\(combined\\)", "", count_label$label)
      
    
  #  count_label <- count_label[!duplicated(count_label$arm_code),]
    
    
  } else {
    gp <- gp + 
      facet_wrap(
        as.formula(paste("~", gmt_facetter)),
        labeller = as_labeller(facet_labeller),
        nrow = facet_n_row
      ) 
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"sr_group", "arm_code", "count", "y", "age_code", "visit_code", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    if(sub_visit){
      visits <- unique(sr_group_gmt_plotdata$visit_code)
    
      target_visit <- visits[!("D15" == visits)]
   
      #   count_label <- count_label[grepl(target_visit[1], count_label$sr_group),]
      count_label <- count_label[target_visit[1] == count_label$visit_code,]
    }
    
    count_label <- count_label %>%
      group_by(sr_group, arm_code, age_code) %>%
      summarize(total_count = sum(count),
                label = paste0("n (", age_code, ") = ", total_count),
                x = x,
                y = ifelse(age_code == "combined", y, ifelse(age_code == "<65", y - 1, y-2)),
                inf_code = inf_code)
    
    count_label <- count_label %>%
      mutate(y = ifelse(inf_code == "inf", y, y-1))
    
    count_label$label <- gsub("\\(combined\\)", "", count_label$label)
    
   # count_label <- count_label[!duplicated(count_label$arm_code),]
  }
  
  if(shape_by == "inf_code"){
    gp <- gp + 
      scale_shape_manual(values = inf_code_shapes)
    
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt +
      geom_text(data = sr_group_gmt_plotdata,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = ag_order[ag_name], y = y, label = label),
                color = "black",
                size = 2.5
      )
  }
  
  if(show_group_count){
   
    gp_gmt <- gp_gmt +
      geom_text(data = count_label,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = x, y = y+0.5, label = label, hjust = 1),
                color = "black",
                size = 4
      )
  }
  
  if(show_mean_line){
    gp <- gp +
      geom_line(aes_string(y = "gmt_arm", group = "sr_group", linetype = line_by), color = mean_line_color, 
                size = 1.2, alpha = 0.3)
    
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.2, alpha = 0.3)
    
  }
  
  return(list("all" = gp, "gmt" = gp_gmt))
  
  
}


titer_differences_lineplot <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                          facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, color_by = "sr_group",
                          difference_by = "visit_code", x_position_by = "visit_code",  shape_by = "inf_code", line_by = "age_code",
                          y_label = "Log2 titer difference", cols_to_keep = c(), show_gmt_label = F, show_group_count = F,
                          show_mean_line = F, mean_line_color = "black",
                          gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code",
                          p_norm = FALSE) {
  
 
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
    max_y <- 13
  } else {
    if(p_norm){
      max_y <- 3.5
    } else {
      max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
    }
    
  }

  min_y_long <- ifelse(min(data_long$logtiter) < 0, floor(min(data_long$logtiter)), 0) 
  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep) %>%
    calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(titer, 0), "\n", fold_change))
 
  min_y <- ifelse(min(sr_group_gmt_plotdata$logtiter) < 0, floor(min(sr_group_gmt_plotdata$lower)), 0) 
  
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  # add x position
  data_long <- x_position_spacing(data_long, position_by = x_position_by)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by)
  
  # fill  NA values of lower and upper CI with mean value
  sr_group_gmt_plotdata$lower[is.na(sr_group_gmt_plotdata$lower)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$lower)]
  sr_group_gmt_plotdata$upper[is.na(sr_group_gmt_plotdata$upper)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$upper)]
  
  
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", "", x)
    x <- arm_plot_names[x]
    x
  }
  
  if(length(color_by) > 1){
    color_by <- "sr_group"
  }
  
  data_long %>%
    ggplot(
      aes_string(
        x = "ag_name",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = "inf_code"
      )
    ) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_line(
      aes(group = sr_name),
      alpha = 0.4
    ) + 
    geom_point(
      alpha = 0.4
    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
      aes(group = sr_group), 
    #      linetype = line_by),
      size = 1.3
    ) +
    geom_pointrange(
      data = sr_group_gmt_plotdata,
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = ag_name, y = y, label = label),
              color = "black",
              size = 3
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_x_discrete(
      limits = antigens,
      labels = ag_plot_names[antigens],
      name = "Antigen variant"
    ) +
    labs(
      y = "Boost magnitude (fold_change)"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(min_y_long-0.5, max_y+1), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    ) -> gp
  
  
  
  # plot gmts of different serum groups together
 
 # max_y <- max(ceiling(sr_group_gmt_plotdata$logtiter), na.rm = T) + 1
  
  if(y_label == "Titer"){
    max_y <- 13
  } else {
    if(p_norm){
      max_y <- 3.5
    } else {
      max_y <- 4 +0.5
      min_y <- -3
    }
  }
  
  spacing <- 0.4
  if(y_label != "Log2 titer difference") {
    spacing <- 0.6
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
 
  sr_group_gmt_plotdata %>%
    ggplot(
      aes_string(
        x = "x",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = "inf_code"
      )
    ) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_line(
      aes_string(group = "sr_group", 
          linetype = line_by),
      size = 1
    ) +
    # pointrange if I want them dodged, x position by. If I just want it dodged then I can add position = position_dodge(width = 1)
    geom_pointrange(
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype= "none") + 
    scale_x_continuous(
      breaks = ag_order,
      labels = ag_plot_names[names(ag_order)],
      limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(ag_order, na.rm = T))+0.5),
      name = "Antigen variant"
    ) + 
    labs(
      y = "Boost magnitude (fold change)"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(min_y-0.5, max_y+0.5), 0.05)
    ) +
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
    )  -> gp_gmt
  
  # do the facets
  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
    gp <- gp + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"v_manuf_code", "count", "y_label_pos", color_by, shape_by, gmt_grid_col, gmt_grid_row))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    dots = sapply(c(gmt_grid_row, gmt_grid_col), . %>% {as.formula(paste0('~', .))})
    
    count_label <- count_label %>%
      group_by(.dots = dots)%>%
      summarize(total_count = mean(count),
                label = paste0("n = ", total_count),
                x = x,
                y = mean(y_label_pos),# -0.5,
                inf_code = inf_code)

    #  count_label <- count_label[!duplicated(count_label$label),]
    
    
  } else {
    gp <- gp + 
      facet_wrap(
        vars(sr_group),
        labeller = as_labeller(facet_labeller),
        nrow = facet_n_row
      ) 
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter, "count", "y_label_pos", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    count_label <- count_label %>%
      group_by(arm_code) %>%
      summarize(total_count = mean(count),
                label = paste0("n = ", total_count),
                x = x,
                y = mean(y_label_pos)-0.5,
                inf_code = inf_code)
    
    
  }
  
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt + 
      geom_text(aes_string(color = color_by) + aes(x = ifelse(age_code == "<65", x + point_spacing/2, x-point_spacing/2), y = y_label_pos,
                                                   label = fold_change),
                color = "black",
                size = 3) 
  }
  
  if(show_group_count){
   
    gp_gmt <- gp_gmt +
      geom_text(data = count_label,
                mapping = aes(x = x, y = y, label = label, hjust = 1),
                color = "black",
                size = 4
      )
  }
  
  if(y_label != "Log2 titer difference") {

    gp <- gp + scale_y_continuous(
      breaks = seq(round(min_y), max_y, by = 2),
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x))*10), paste0(round(2^x)*10)),
      name = "Boost magnitude (linear)"
    ) 
    
    gp_gmt <- gp_gmt + scale_y_continuous(
      breaks = seq(round(min_y), max_y, by = 2),
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x))*10), paste0(round(2^x)*10)),
      name = "Boost magnitude (linear)"
    ) 
  } else {
    
    gp <- gp + scale_y_continuous(
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x)), "x"), paste0(round(2^x), "x")),
      breaks= c(round(min_y_long):5)
    )  + 
      ylab("Boost magnitude (fold change)")
    
    gp_gmt <- gp_gmt + scale_y_continuous(
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x)), "x"), paste0(round(2^x), "x")),
      breaks= c(round(min_y):round(max_y)),
      limits = c(round(min_y), max_y)
    )  + 
      ylab("Boost magnitude (fold change)")
    
  }
  
  if(shape_by == "inf_code"){
    gp <- gp + 
      scale_shape_manual(values = inf_code_shapes)
    
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
  if(show_mean_line){
    gp <- gp +
      geom_line(aes_string(y = "gmt_arm", group = "arm_code", linetype = line_by), color = mean_line_color, size = 1.2,
                alpha = 0.3)
    
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.2, alpha = 0.3)
    
  }
  
  return(list("all" = gp, "gmt" = gp_gmt))
  
  
}

titerlineplot_over_time <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                          facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, 
                          color_by = "sr_group", shape_by = "inf_code", line_by = "age_code", x_position_by = "age_code",
                          cols_to_keep = c("arm_code", "age_code", "inf_code"), to_long = T, 
                          show_gmt_label = F, show_group_count = F, show_mean_line = F, mean_line_color = "black",
                          gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code") {
  
  ag_order <- ag_order[antigens]
  # format titer to longer
  if(to_long){
    data %>%
      pivot_longer(cols = antigens, names_to = "ag_name", values_to = "titer") %>%
      mutate(logtiter = log2(titer/10),
             logtiter = ifelse(titer < titer_thresh, log2((titer_thresh/20)), logtiter)) -> data_long
  } else{
    data_long <- data 
  }
  
  # set sr_order
  if(!is.null(sr_group_order)) {
    data_long$sr_group <- factor(as.character(data_long$sr_group), levels =sr_group_order)
  }
  
  
  #  max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
  max_y <- 15
  
  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep) %>%
    filter(!all_below_thresh) %>%
    calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(as.numeric(titer), 0), "\n", fold_change))
  
  # sr_group_gmt_plotdata <- sr_group_to_column(sr_group_gmt_plotdata)
  
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  data_long <- x_position_spacing(data_long, position_by = x_position_by)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by)
  # 
  
#  visit_order <- c("D1", "D15", "D29", "D57", "D91")
  sr_group_gmt_plotdata$visit_code <- as.numeric(gsub("D", "", sr_group_gmt_plotdata$visit_code))
  
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", " ", x)
    x <- arm_plot_names[x]
    x
  }
  
  # facet_labeller <- function(x){
  #   strsplit(x, "-")[[1]][2]
  # }
  sub_visit <- FALSE
  if(length(color_by) > 1){
    if("visit_code" %in% color_by | length(unique(data$visit_code)) > 1){
      sub_visit <- TRUE
    }
    color_by <- "sr_group"
  }
  
  
  # plot gmts of different serum groups together
  max_y <- 13
  sr_group_gmt_plotdata$y <- max_y #- 0.5
  
  sr_group_gmt_plotdata %>%
    ggplot(
      aes_string(
        x = "visit_code",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) + 
    geom_line(
      aes_string(group = "sr_group",
                 linetype = line_by),
      size = 1
    ) +
    geom_pointrange(
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+0.5), 0.05)
    ) +
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
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    ) -> gp_gmt
  
  
  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
   
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col))
      ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"sr_group", "arm_code", "v_manuf_code", "count", "y", "age_code", "visit_code", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    if(sub_visit){
      visits <- unique(sr_group_gmt_plotdata$visit_code)
      
      target_visit <- visits[visits != "D15"]
      
      #   count_label <- count_label[grepl(target_visit[1], count_label$sr_group),]
      count_label <- count_label[target_visit[1] == count_label$visit_code,]
    }
    
    count_label <- count_label %>%
      group_by(sr_group,arm_code, v_manuf_code, age_code) %>%
      summarize(total_count = sum(count),
                label = paste0("n (", age_code, ") = ", total_count),
                x = x,
                y = ifelse(age_code == "combined", y, ifelse(age_code == "<65", y - 1, y-2)),
                inf_code = inf_code)
    
    count_label <- count_label %>%
      mutate(y = ifelse(inf_code == "inf", y, y-1))
    
    count_label$label <- gsub("\\(combined\\)", "", count_label$label)
    
    
    #  count_label <- count_label[!duplicated(count_label$arm_code),]
    
    
  } else {
    
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"sr_group", "arm_code", "count", "y", "age_code", "visit_code", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    if(sub_visit){
      visits <- unique(sr_group_gmt_plotdata$visit_code)
      
      target_visit <- visits[!("D15" == visits)]
      
      #   count_label <- count_label[grepl(target_visit[1], count_label$sr_group),]
      count_label <- count_label[target_visit[1] == count_label$visit_code,]
    }
    
    count_label <- count_label %>%
      group_by(sr_group, arm_code, age_code) %>%
      summarize(total_count = sum(count),
                label = paste0("n (", age_code, ") = ", total_count),
                x = x,
                y = ifelse(age_code == "combined", y, ifelse(age_code == "<65", y - 1, y-2)),
                inf_code = inf_code)
    
    
    
    count_label$label <- gsub("\\(combined\\)", "", count_label$label)
    
    # count_label <- count_label[!duplicated(count_label$arm_code),]
  }
  
  if(shape_by == "inf_code"){
    
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt +
      geom_text(data = sr_group_gmt_plotdata,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = ag_order[ag_name], y = y, label = label),
                color = "black",
                size = 2.5
      )
  }
  
  if(show_group_count){
    
    gp_gmt <- gp_gmt +
      geom_text(data = count_label,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = x, y = y+0.5, label = label, hjust = 1),
                color = "black",
                size = 4
      )
  }
  
  if(show_mean_line){
    
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.2)
    
  }
  
  return(gp_gmt)
  
  
}



titerlineplot_dodge <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                          facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, 
                          color_by = "sr_group", shape_by = "inf_code", line_by = "age_code", x_position_by = "age_code",
                          cols_to_keep = c("arm_code", "age_code", "inf_code"), to_long = T, 
                          show_gmt_label = F, show_group_count = F, show_mean_line = F, mean_line_color = "black",
                          gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code", dodge_group = "visit_code", pre_titer_adj = FALSE,
                          flu_plot = FALSE, thinner_pre = FALSE,
                          dodge_width = 0.4,
                          add_group_count_to_legend = FALSE,
                          x_axis_var = "ag_name",
                          plot_grouping_var = "sr_group") {
  
 # ag_order <- ag_order[antigens]
  # format titer to longer
  if(to_long){
    data %>%
      pivot_longer(cols = antigens, names_to = "ag_name", values_to = "titer") %>%
      mutate(logtiter = log2(titer/10),
             logtiter = ifelse(titer < titer_thresh, log2((titer_thresh/20)), logtiter)) -> data_long
  } else{
    data_long <- data 
  }
 # data_long <- data_long %>%
#    filter(!is.na(titer))
  # create plot colors
  
  sr_groups <- unique(as.character(data_long$sr_group))
  
  if(length(color_by) == 1){
    if(color_by == "lab_code"){
    plot_colors <- c("Monogram" = "#ed94a7", "Montefiori" = "#6C992E", "Adjusted" = "#bb94f1",
      "Monogram-PNEUT50" = "#ed94a7", "Monogram-PNEUT80" = "#ed94a7", "Montefiori-PNEUT50" = "#6C992E", "Montefiori-PNEUT80" = "#9bd748")
    } else {

      plot_colors <- unlist(lapply(sr_groups, function(x) {
        blend_sr_group_colors(x, sr_group_colors, split_by = color_by)
      }))
  
      sr_group_colors[sr_groups, "Color"] <- plot_colors
  
      plot_colors <- sr_group_colors$Color
  
      names(plot_colors)<- rownames(sr_group_colors)

    }
  } else {
      plot_colors <- unlist(lapply(sr_groups, function(x) {
        blend_sr_group_colors(x, sr_group_colors, split_by = color_by)
      }))
  
      sr_group_colors[sr_groups, "Color"] <- plot_colors
  
      plot_colors <- sr_group_colors$Color
  
      names(plot_colors)<- rownames(sr_group_colors)

  }

  
  
  # plot_colors <-  plot_colors[as.character(unique(data_long[, color_by]))]
  
  # set sr_order
  if(!is.null(sr_group_order)) {
    data_long$sr_group <- factor(as.character(data_long$sr_group), levels =sr_group_order)
  }
  
  
  #  max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
  max_y <- 15
  
  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep, pre_adj = pre_titer_adj) %>%
    filter(!all_below_thresh) %>%
    calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(as.numeric(titer), 0), "\n", fold_change))

  # sr_group_gmt_plotdata <- sr_group_to_column(sr_group_gmt_plotdata)
  
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  data_long <- x_position_spacing(data_long, position_by = x_position_by, x_position_var = x_axis_var)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by, x_position_var = x_axis_var)
  # 
  # fill  NA values of lower and upper CI with mean value
  sr_group_gmt_plotdata$lower[is.na(sr_group_gmt_plotdata$lower)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$lower)]
  sr_group_gmt_plotdata$upper[is.na(sr_group_gmt_plotdata$upper)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$upper)]
  
  facet_labeller <- function(x){
  
    if(length(x) == 2 & "inf" %in% x){
      x <- c("inf" = "Prior infection", "non_inf" = "Naive")[x]
      return(x)
    }
    
    if(length(x) == 4 & "inf" %in% x){
      x <- c("inf" = "Prior infection", "non_inf" = "Naive", "inf_breakthrough" = "Prior infection + BT", "non_inf_breakthrough" = "Naive + BT")[x]
      return(x)
    }
    
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", " ", x)
    
   # if(x %in% names(arm_plot_names)){
      x <- arm_plot_names[x]
  #  }
    x
  }
  
  # facet_labeller <- function(x){
  #   strsplit(x, "-")[[1]][2]
  # }
  sub_visit <- FALSE
  if(length(color_by) > 1){
    if("visit_code" %in% color_by | length(unique(data$visit_code)) > 1){
      sub_visit <- TRUE
    }
    color_by <- "sr_group"
  }
  
  data_long$visit_code <- factor(data_long$visit_code, levels = visit_code_order)
  sr_group_gmt_plotdata$visit_code <- factor(sr_group_gmt_plotdata$visit_code, levels = visit_code_order)

  
  plot_color_labels <- gsub("NA|non_inf|inf|-|P|O|\\+|\\:|M", "", names(plot_colors))
  plot_color_labels <- arm_plot_names[plot_color_labels]
  
  names(plot_color_labels) <- names(plot_colors)
  
  data_long %>%
    ggplot(
      aes_string(
        x = "ag_name",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) + 
    geom_line(
      aes_string(group = "sr_name", 
          color = ifelse(PLOTLY, "sr_name", color_by)),
      linetype = "solid",
      alpha = 0.4
    ) + 
 #   geom_point(
#      alpha = 0.4
#    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
      aes_string(group = "sr_group",
                 linetype = line_by),
      size = 1.3
    ) +
    geom_pointrange(
      data = sr_group_gmt_plotdata,
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors,
                       name = "",
                       labels = plot_color_labels) +
    scale_fill_manual(values = plot_colors,
                      labels = plot_color_labels,
                      name = "") +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      threshold = paste0("<", titer_thresh),
      ymin = log2(titer_thresh/10)
    ) + 
    scale_x_discrete(
      limits = antigens,
      labels = ag_plot_names[antigens],
      name = "Antigen variant"
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+1), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 10)
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.2
    ) -> gp
  
  
  visits <- unique(sr_group_gmt_plotdata$visit_code)
  
  if(line_by != "v_manuf_code" & "V3" %in% visits & "V4" %in% visits){
    line_by <- "visit_code"
  }
  # plot gmts of different serum groups together
  max_y <- 13
  sr_group_gmt_plotdata$y <- max_y #- 0.5
  
  if(length(unique(sr_group_gmt_plotdata$visit_code)) <=4 & dodge_group == "visit_code"){
    sr_group_gmt_plotdata$visit_code <- factor(sr_group_gmt_plotdata$visit_code, levels = c("D91","D15", "D1", "D29","D181","D271", "D57", "D85", "D147"))
  }

  if(dodge_group == "arm_code"){
    sr_group_gmt_plotdata$arm_code <- factor(as.character(sr_group_gmt_plotdata$arm_code), levels = arm_code_order)
  }

  p_arm <- sr_group_gmt_plotdata %>% filter(arm_code == "P") %>% filter(visit_code == "D1")
  
  sr_group_gmt_plotdata %>%
    ggplot(
      aes_string(
        x = "x",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) -> gp_gmt
  
  if(thinner_pre){
    gp_gmt + geom_line(data = 
                         sr_group_gmt_plotdata %>%
                         filter(visit_code %in% c("D1", "pD1")),
      aes_string(group = plot_grouping_var,
                 linetype = line_by),
      size = 0.5
    ) + geom_line(data = 
                    sr_group_gmt_plotdata %>%
                    filter(!(visit_code %in% c("D1", "pD1"))),
                  aes_string(group = plot_grouping_var,
                             linetype = line_by),
                  size = 1
    ) -> gp_gmt
    } else {
    gp_gmt + geom_line(
      aes_string(group =  plot_grouping_var,
                 linetype = line_by),
      size = 1
    ) -> gp_gmt
    
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
          group = dodge_group), # usually visit code
      width = 0,
      position = position_dodge(width = dodge_width)
    ) +
    geom_point(
      size = 2
    ) +
    geom_point(
      size = 1.5,
      fill = "white",
      shape = 21
    ) -> gp_gmt
    
    if(show_mean_line){
      gp_gmt + geom_line(
        data = p_arm,
        aes_string(group = "sr_group",
                   linetype = line_by),
        size = 1
      ) +
        geom_point(
          data = p_arm,
          size = 2
        ) +
        geom_point(
          data = p_arm,
          size = 1.5,
          fill = "white",
          shape = 21
        ) -> gp_gmt
    }
    
    gp_gmt + scale_color_manual(values = plot_colors,
                                name = "",
                                labels = plot_color_labels) +
      scale_fill_manual(values = plot_colors,
                        labels = plot_color_labels,
                        name = "") +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
      coord_cartesian(
        ylim = shrinkrange(c(-0.5, max_y+0.5), 0.05)
      ) +
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
      ) -> gp_gmt
    
    if(x_axis_var == "visit_code"){
      
      gp_gmt <- gp_gmt + 
        scale_x_continuous(
          breaks = as.numeric(gsub("D", "", as.character(unique(sr_group_gmt_plotdata$visit_code)))), #visit_plot_order,
          # labels = ag_plot_names[names(ag_order)],
        #  limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(visit_plot_order, na.rm = T))+0.5),
          name = "Time since exposure (days)"
        )
      
    } else {
      
      gp_gmt <- gp_gmt + 
        scale_x_continuous(
          breaks = ag_order,
          labels = ag_plot_names[names(ag_order)],
          limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(ag_order, na.rm = T))+0.5),
          name = "Antigen variant"
        )
    }
     
    if(!flu_plot){
      gp_gmt <- gp_gmt +
        annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.2
    ) 
    }
 
  
  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
    gp <- gp + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
  
    
  } else {
    gp <- gp + 
      facet_wrap(
        vars(sr_group),
        labeller = as_labeller(facet_labeller),
        nrow = facet_n_row
      ) 
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
  
  }
  
  if(shape_by == "inf_code"){
    gp <- gp + 
      scale_shape_manual(values = inf_code_shapes)
    
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt +
      geom_text(data = sr_group_gmt_plotdata,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = ag_order[ag_name], y = y, label = label),
                color = "black",
                size = 2.5
      )

    gp <- gp +
      geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = ag_name, y = y+2, label = label),
              color = "black",
              size = 2.5
    )
  }
  
  if(show_group_count){
    
    count_label <- count_label_function(sr_group_gmt_plotdata, max_y+1, gmt_facetter, color_by, shape_by, gmt_grid_row,
                                        cols_to_keep = cols_to_keep)
  
    
    if(add_group_count_to_legend){
      
      count_label %>%
        mutate(label_split = str_split(label, " = ", 2),
               label_count = sapply(label_split, function(x){
                 x[2]
               }),
               new_label = paste0(" (n = ", label_count, ")")) -> count_label
      
      plot_color_labels <- paste0(plot_color_labels, count_label$new_label[match(names(plot_color_labels), count_label$sr_group)])
  
      names(plot_color_labels) <- names(plot_colors)
      
      gp_gmt <- gp_gmt + scale_color_manual(values = plot_colors,
                         name = "",
                         labels = plot_color_labels) +
        scale_fill_manual(values = plot_colors,
                           name = "",
                           labels = plot_color_labels)
      
    } else {
      
      gp_gmt <- gp_gmt +
        geom_text(data = count_label,
                  #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                  mapping = aes_string(x = "x", y = "y + 1", label = "label", hjust = "1", color = color_by),
                 # color = "black",
                  size = 2.5
        )
      
    }
    
  }
  
  if(show_mean_line){
   
    gp <- gp +
      geom_line(aes_string(y = "gmt_arm", group = "sr_group", linetype = line_by), color = mean_line_color, 
                size = 1.2, alpha = 0.3)
    
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.2, alpha = 0.3)
    
  }
  
  return(list("all" = gp, "gmt" = gp_gmt))
  
  
}



titerlineplot_dodge_alpha_pre <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                                facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, 
                                color_by = "sr_group", shape_by = "inf_code", line_by = "age_code", x_position_by = "age_code",
                                cols_to_keep = c("arm_code", "age_code", "inf_code"), to_long = T, 
                                show_gmt_label = F, show_group_count = F, show_mean_line = F, mean_line_color = "black",
                                gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code", dodge_group = "visit_code", pre_titer_adj = FALSE,
                                flu_plot = FALSE, thinner_pre = FALSE,
                                dodge_width = 0.4,
                                add_group_count_to_legend = FALSE,
                                x_axis_var = "ag_name",
                                plot_grouping_var = "sr_group") {
  
  # ag_order <- ag_order[antigens]
  # format titer to longer
  if(to_long){
    data %>%
      pivot_longer(cols = antigens, names_to = "ag_name", values_to = "titer") %>%
      mutate(logtiter = log2(titer/10),
             logtiter = ifelse(titer < titer_thresh, log2((titer_thresh/20)), logtiter)) -> data_long
  } else{
    data_long <- data 
  }
  # data_long <- data_long %>%
  #    filter(!is.na(titer))
  # create plot colors
  
  sr_groups <- unique(as.character(data_long$sr_group))
  
  if(length(color_by) == 1){
    if(color_by == "lab_code"){
      plot_colors <- c("Monogram" = "#ed94a7", "Montefiori" = "#6C992E", "Adjusted" = "#bb94f1",
                       "Monogram-PNEUT50" = "#ed94a7", "Monogram-PNEUT80" = "#ed94a7", "Montefiori-PNEUT50" = "#6C992E", "Montefiori-PNEUT80" = "#9bd748")
    } else {
      
      plot_colors <- unlist(lapply(sr_groups, function(x) {
        blend_sr_group_colors(x, sr_group_colors, split_by = color_by)
      }))
      
      sr_group_colors[sr_groups, "Color"] <- plot_colors
      
      plot_colors <- sr_group_colors$Color
      
      names(plot_colors)<- rownames(sr_group_colors)
      
    }
  } else {
    plot_colors <- unlist(lapply(sr_groups, function(x) {
      blend_sr_group_colors(x, sr_group_colors, split_by = color_by)
    }))
    
    sr_group_colors[sr_groups, "Color"] <- plot_colors
    
    plot_colors <- sr_group_colors$Color
    
    names(plot_colors)<- rownames(sr_group_colors)
    
  }
  
  
  
  # plot_colors <-  plot_colors[as.character(unique(data_long[, color_by]))]
  
  # set sr_order
  if(!is.null(sr_group_order)) {
    data_long$sr_group <- factor(as.character(data_long$sr_group), levels =sr_group_order)
  }
  
  
  #  max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
  max_y <- 15
  
  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep, pre_adj = pre_titer_adj) %>%
    filter(!all_below_thresh) %>%
    calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(as.numeric(titer), 0), "\n", fold_change))
  
  # sr_group_gmt_plotdata <- sr_group_to_column(sr_group_gmt_plotdata)
  
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  data_long <- x_position_spacing(data_long, position_by = x_position_by, x_position_var = x_axis_var)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by, x_position_var = x_axis_var)
  # 
  # fill  NA values of lower and upper CI with mean value
  sr_group_gmt_plotdata$lower[is.na(sr_group_gmt_plotdata$lower)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$lower)]
  sr_group_gmt_plotdata$upper[is.na(sr_group_gmt_plotdata$upper)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$upper)]
  
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", " ", x)
    # if(x %in% names(arm_plot_names)){
    x <- arm_plot_names[x]
    #  }
    x
  }
  
  # facet_labeller <- function(x){
  #   strsplit(x, "-")[[1]][2]
  # }
  sub_visit <- FALSE
  if(length(color_by) > 1){
    if("visit_code" %in% color_by | length(unique(data$visit_code)) > 1){
      sub_visit <- TRUE
    }
    color_by <- "sr_group"
  }
  
  data_long$visit_code <- factor(data_long$visit_code, levels = visit_code_order)
  sr_group_gmt_plotdata$visit_code <- factor(sr_group_gmt_plotdata$visit_code, levels = visit_code_order)
  
  
  plot_color_labels <- gsub("NA|non_inf|inf|-|P|O|\\+|\\:|M", "", names(plot_colors))
  plot_color_labels <- arm_plot_names[plot_color_labels]
  
  names(plot_color_labels) <- names(plot_colors)
  
  data_long %>%
    ggplot(
      aes_string(
        x = "ag_name",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) + 
    geom_line(
      aes_string(group = "sr_name", 
                 color = ifelse(PLOTLY, "sr_name", color_by)),
      linetype = "solid",
      alpha = 0.4
    ) + 
    #   geom_point(
    #      alpha = 0.4
    #    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
      aes_string(group = "sr_group",
                 linetype = line_by),
      size = 1.3
    ) +
    geom_pointrange(
      data = sr_group_gmt_plotdata,
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors,
                       name = "",
                       labels = plot_color_labels) +
    scale_fill_manual(values = plot_colors,
                      labels = plot_color_labels,
                      name = "") +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      threshold = paste0("<", titer_thresh),
      ymin = log2(titer_thresh/10)
    ) + 
    scale_x_discrete(
      limits = antigens,
      labels = ag_plot_names[antigens],
      name = "Antigen variant"
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+1), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 10)
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    ) -> gp
  
  
  visits <- unique(sr_group_gmt_plotdata$visit_code)
  
  if(line_by != "v_manuf_code" & "V3" %in% visits & "V4" %in% visits){
    line_by <- "visit_code"
  }
  # plot gmts of different serum groups together
  max_y <- 13
  sr_group_gmt_plotdata$y <- max_y #- 0.5
  
  if(length(unique(sr_group_gmt_plotdata$visit_code)) <=4 & dodge_group == "visit_code"){
    sr_group_gmt_plotdata$visit_code <- factor(sr_group_gmt_plotdata$visit_code, levels = c("D91","D15", "D1", "D29","D181","D271", "D57", "D85", "D147"))
  }
  
  if(dodge_group == "arm_code"){
    sr_group_gmt_plotdata$arm_code <- factor(as.character(sr_group_gmt_plotdata$arm_code), levels = arm_code_order)
  }
  
  p_arm <- sr_group_gmt_plotdata %>% filter(arm_code == "P") %>% filter(visit_code == "D1")
  
  arm_codes <- unique(as.character(sr_group_gmt_plotdata$arm_code))
  alpha_vals <- rep(1, length(arm_codes))
  names(alpha_vals) <- arm_codes
  alpha_vals["pre"] <- 0.4
  
  sr_group_gmt_plotdata$alpha_val <- sr_group_gmt_plotdata$arm_code

  sr_group_gmt_plotdata %>%
    arrange(arm_code) %>%
    ggplot(
      aes_string(
        x = "x",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) -> gp_gmt
  
  if(thinner_pre){
    gp_gmt + geom_line(data = 
                         sr_group_gmt_plotdata %>%
                         filter(visit_code %in% c("D1", "pD1")),
                       aes_string(group = plot_grouping_var,
                                  linetype = line_by),
                       size = 0.5
    ) + geom_line(data = 
                    sr_group_gmt_plotdata %>%
                    filter(!(visit_code %in% c("D1", "pD1"))),
                  aes_string(group = plot_grouping_var,
                             linetype = line_by),
                  size = 1
    ) -> gp_gmt
  } else {
    gp_gmt + geom_line(
      aes_string(group =  plot_grouping_var,
                 linetype = line_by,
                 alpha = "alpha_val"),
      size = 1
    ) -> gp_gmt
    
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
               group = dodge_group,
               alpha = "alpha_val"), # usually visit code
    width = 0,
    position = position_dodge(width = dodge_width)
  ) +
    geom_point(
      aes_string(alpha = "alpha_val"),
      size = 2
    ) +
    geom_point(
      aes_string(alpha = "alpha_val"),
      size = 1.5,
      fill = "white",
      shape = 21
    ) -> gp_gmt
  
  if(show_mean_line){
    gp_gmt + geom_line(
      data = p_arm,
      aes_string(group = "sr_group",
                 linetype = line_by,
                 alpha = "alpha_val"),
      size = 1
    ) +
      geom_point(
        data = p_arm,
        size = 2
      ) +
      geom_point(
        data = p_arm,
        size = 1.5,
        fill = "white",
        shape = 21
      ) -> gp_gmt
  }
  
  gp_gmt + scale_color_manual(values = plot_colors,
                              name = "",
                              labels = plot_color_labels) +
    scale_fill_manual(values = plot_colors,
                      labels = plot_color_labels,
                      name = "") +
    scale_linetype_manual(values= line_vals) +
    scale_alpha_manual(values = alpha_vals, guide = "none") + 
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+0.5), 0.05)
    ) +
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
    ) -> gp_gmt
  
  if(x_axis_var == "visit_code"){
    
    gp_gmt <- gp_gmt + 
      scale_x_continuous(
        breaks = as.numeric(gsub("D", "", as.character(unique(sr_group_gmt_plotdata$visit_code)))), #visit_plot_order,
        # labels = ag_plot_names[names(ag_order)],
        #  limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(visit_plot_order, na.rm = T))+0.5),
        name = "Time since exposure (days)"
      )
    
  } else {
    
    gp_gmt <- gp_gmt + 
      scale_x_continuous(
        breaks = ag_order,
        labels = ag_plot_names[names(ag_order)],
        limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(ag_order, na.rm = T))+0.5),
        name = "Antigen variant"
      )
  }
  
  if(!flu_plot){
    gp_gmt <- gp_gmt +
      annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = -1,
        ymax = log2(titer_thresh/10),
        fill = "grey50",
        color = NA,
        alpha = 0.3
      ) 
  }
  
  
  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
    gp <- gp + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    
    
  } else {
    gp <- gp + 
      facet_wrap(
        vars(sr_group),
        labeller = as_labeller(facet_labeller),
        nrow = facet_n_row
      ) 
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
  }
  
  if(shape_by == "inf_code"){
    gp <- gp + 
      scale_shape_manual(values = inf_code_shapes)
    
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt +
      geom_text(data = sr_group_gmt_plotdata,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = ag_order[ag_name], y = y, label = label),
                color = "black",
                size = 2.5
      )
    
    gp <- gp +
      geom_text(data = sr_group_gmt_plotdata,
                mapping = aes(x = ag_name, y = y+2, label = label),
                color = "black",
                size = 2.5
      )
  }
  
  if(show_group_count){
    
    count_label <- count_label_function(sr_group_gmt_plotdata %>%
                                          filter(arm_code != "pre"), max_y, gmt_facetter, color_by, shape_by, gmt_grid_row,
                                        cols_to_keep = cols_to_keep)
    
    if(add_group_count_to_legend){
      
      count_label %>%
        mutate(label_split = str_split(label, " = ", 2),
               label_count = sapply(label_split, function(x){
                 x[2]
               }),
               new_label = paste0(" (n = ", label_count, ")")) -> count_label
      
      plot_color_labels <- paste0(plot_color_labels, count_label$new_label[match(names(plot_color_labels), count_label$sr_group)])
      
      names(plot_color_labels) <- names(plot_colors)
      
      gp_gmt <- gp_gmt + scale_color_manual(values = plot_colors,
                                            name = "",
                                            labels = plot_color_labels) +
        scale_fill_manual(values = plot_colors,
                          name = "",
                          labels = plot_color_labels)
      
    } else {
      
      gp_gmt <- gp_gmt +
        geom_text(data = count_label,
                  #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                  mapping = aes(x = x, y = y+1, label = label, hjust = 1, color = arm_code),
                  # color = "black",
                  size = 2.5
        )
      
    }
    
  }
  
  if(show_mean_line){
    
    gp <- gp +
      geom_line(aes_string(y = "gmt_arm", group = "sr_group", linetype = line_by), color = mean_line_color, 
                size = 1.2, alpha = 0.3)
    
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.2, alpha = 0.3)
    
  }
  
  return(list("all" = gp, "gmt" = gp_gmt))
  
  
}


titer_differences_lineplot_dodge <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                                       facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, color_by = "sr_group",
                                       difference_by = "visit_code", x_position_by = "visit_code",  shape_by = "inf_code", line_by = "age_code",
                                       y_label = "Log2 titer difference", cols_to_keep = c(), show_gmt_label = F, show_group_count = F,
                                       show_mean_line = F, mean_line_color = "black",
                                       gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code",
                                       p_norm = FALSE, dodge_group = "visit_code", highlight_samples = c()) {
  
  
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
    max_y <- 13
  } else {
    if(p_norm){
      max_y <- 3.5
    } else {
      max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
    }
    
  }
  
  min_y_long <- ifelse(min(data_long$logtiter) < 0, floor(min(data_long$logtiter)), 0) 
  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep) %>%
    calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(titer, 0), "\n", fold_change))
  
  min_y <- ifelse(min(sr_group_gmt_plotdata$logtiter) < 0, floor(min(sr_group_gmt_plotdata$lower)), 0) 
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  # add x position
  data_long <- x_position_spacing(data_long, position_by = x_position_by)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by)
  # fill  NA values of lower and upper CI with mean value
  sr_group_gmt_plotdata$lower[is.na(sr_group_gmt_plotdata$lower)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$lower)]
  sr_group_gmt_plotdata$upper[is.na(sr_group_gmt_plotdata$upper)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$upper)]
  
  
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", "", x)
    x <- arm_plot_names[x]
    x[is.na(x)] <- ""
    x
  }
  
  if(length(color_by) > 1){
    color_by <- "sr_group"
  }
  
  data_long %>%
    ggplot(
      aes_string(
        x = "ag_name",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = "inf_code"
      )
    ) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_line(
      aes(group = sr_name),
      alpha = 0.4
    ) + 
    geom_point(
      alpha = 0.4
    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
      aes(group = sr_group), 
      #      linetype = line_by),
      size = 1.3
    ) +
    geom_pointrange(
      data = sr_group_gmt_plotdata,
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) -> gp
  if(show_gmt_label){
    gp + geom_text(data = sr_group_gmt_plotdata,
                   mapping = aes(x = ag_name, y = y, label = label),
                   color = "black",
                   size = 3
    ) -> gp
  }
    gp + scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_x_discrete(
      limits = antigens,
      labels = ag_plot_names[antigens],
      name = "Antigen variant"
    ) +
    labs(
      y = "Boost magnitude (fold_change)"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(min_y_long-0.5, max_y+0.5), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    ) -> gp
  
  
    if(length(highlight_samples) > 0){
      
      data_long %>%
        filter(subjid %in% highlight_samples) -> data_highlight
      
      gp <- gp + 
        geom_line(data = data_highlight, aes(group = sr_name), color = ifelse("non_inf" %in% unique(data_long$inf_code), "red", "blue"))
      
    }
  
  # plot gmts of different serum groups together
  
  # max_y <- max(ceiling(sr_group_gmt_plotdata$logtiter), na.rm = T) + 1
  
  if(y_label == "Titer"){
    max_y <- 13
  } else {
    if(p_norm){
      max_y <- 3.5
    } else {
      max_y <- ifelse("visit_diff" %in% cols_to_keep, 8+0.5, 4 +0.5)
      min_y <- -3
    }
  }
  
  spacing <- 0.4
  if(y_label != "Log2 titer difference") {
    spacing <- 0.6
  }
  
  if(length(unique(sr_group_gmt_plotdata$visit_diff)) >1){
    sr_group_gmt_plotdata$visit_diff <- factor(sr_group_gmt_plotdata$visit_diff, levels = visit_diff_order)
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
  

  if(dodge_group == "arm_code"){
    sr_group_gmt_plotdata$arm_code <- factor(as.character(sr_group_gmt_plotdata$arm_code), levels = arm_code_order)
  }


  sr_group_gmt_plotdata %>%
    ggplot(
      aes_string(
        x = "x",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = "inf_code"
      )
    ) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_line(
      aes_string(group = "sr_group", 
                 linetype = line_by),
      size = 1
    ) +
    # pointrange if I want them dodged, x position by. If I just want it dodged then I can add position = position_dodge(width = 1)
    # geom_pointrange(
    #   aes(ymin = lower,
    #       ymax = upper)
    # ) +
    # errorbar if I want shifted bars
    geom_errorbar(
      aes_string(ymin = "lower",
          ymax ="upper",
          group = dodge_group),
      width = 0,
      position = position_dodge(width = 0.4)
    ) +
    geom_point(
      size = 2
    ) +
    geom_point(
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype= "none") + 
    scale_x_continuous(
      breaks = ag_order,
      labels = ag_plot_names[names(ag_order)],
      limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(ag_order, na.rm = T))+0.5),
      name = "Antigen variant"
    ) + 
    labs(
      y = "Boost magnitude (fold change)"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(min_y-0.5, max_y+0.5), 0.05)
    ) +
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
    )  -> gp_gmt
  
  # do the facets
  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
    gp <- gp + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"v_manuf_code", "count", "y_label_pos", color_by, shape_by, gmt_grid_col, gmt_grid_row))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    dots = sapply(c(gmt_grid_row, gmt_grid_col), . %>% {as.formula(paste0('~', .))})
    
    count_label <- count_label %>%
      group_by(.dots = dots)%>%
      summarize(total_count = mean(count),
                label = paste0("n = ", total_count),
                x = x,
                y = mean(y_label_pos), #-0.5,
                inf_code = inf_code)
    
    #  count_label <- count_label[!duplicated(count_label$label),]
    
    
  } else {
    gp <- gp + 
      facet_wrap(
        vars(sr_group),
        labeller = as_labeller(facet_labeller),
        nrow = facet_n_row
      ) 
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter, "count", "y_label_pos", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    count_label <- count_label %>%
      group_by(arm_code) %>%
      summarize(total_count = mean(count),
                label = paste0("n = ", total_count),
                x = x,
                y = mean(y_label_pos)-0.5,
                inf_code = inf_code)
    
    
  }
  
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt + 
      geom_text(aes_string(color = color_by) + aes(x = ifelse(age_code == "<65", x + point_spacing/2, x-point_spacing/2), y = y_label_pos,
                                                   label = fold_change),
                color = "black",
                size = 3) 
  }
  
  if(show_group_count){
    
    count_label <- count_label_function(sr_group_gmt_plotdata, max_y+2, gmt_facetter, color_by, shape_by, gmt_grid_row,
                                        cols_to_keep = cols_to_keep)
  
    gp_gmt <- gp_gmt +
      geom_text(data = count_label,
                mapping = aes(x = x, y = y, label = label, hjust = 1),
                color = "black",
                size = 3
      )
    count_label %>%
      mutate(ag_name = last(plot_antigens)) -> count_label
    
    gp <- gp +
      geom_text(data = count_label,
                mapping = aes(x = ag_name, y = y, label = label, hjust = 1),
                color = "black",
                size = 3
      )
  }
  
  if(y_label != "Log2 titer difference") {
    
    gp <- gp + scale_y_continuous(
      breaks = seq(round(min_y), max_y, by = 2),
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x))*10), paste0(round(2^x)*10)),
      name = "Boost magnitude (linear)"
    ) 
    
    gp_gmt <- gp_gmt + scale_y_continuous(
      breaks = seq(round(min_y), max_y, by = 2),
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x))*10), paste0(round(2^x)*10)),
      name = "Boost magnitude (linear)"
    ) 
  } else {
    
    gp <- gp + scale_y_continuous(
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x)), "x"), paste0(round(2^x), "x")),
      breaks= c(round(min_y_long):round(max_y)),
      limits = c(round(min_y), max_y),
      oob= scales::oob_squish
    )  + 
      coord_cartesian(ylim = c(round(min_y), max_y)) + 
      ylab("Boost magnitude (fold change)")
    
    gp_gmt <- gp_gmt + scale_y_continuous(
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x)), "x"), paste0(round(2^x), "x")),
      breaks= c(round(min_y):round(max_y)),
      limits = c(round(min_y), max_y),
      oob= scales::oob_squish
    )  + 
      ylab("Boost magnitude (fold change)")
    
  }
  
  if(shape_by == "inf_code"){
    gp <- gp + 
      scale_shape_manual(values = inf_code_shapes)
    
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
  if(show_mean_line){
    gp <- gp +
      geom_line(aes_string(y = "gmt_arm", group = "arm_code", linetype = line_by), color = mean_line_color, size = 1.2,
                alpha = 0.3)
    
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.2, alpha = 0.3)
    
  }
  
  return(list("all" = gp, "gmt" = gp_gmt))
  
  
}



titer_differences_time_lineplot_dodge <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                                             facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, color_by = "sr_group",
                                             difference_by = "visit_code", x_position_by = "visit_diff",  shape_by = "inf_code", line_by = "age_code",
                                             y_label = "Log2 titer difference", cols_to_keep = c(), show_gmt_label = F, show_group_count = F,
                                             show_mean_line = F, mean_line_color = "black",
                                             gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code",
                                             p_norm = FALSE) {
  
  
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
    max_y <- 13
  } else {
    if(p_norm){
      max_y <- 3.5
    } else {
      max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
    }
    
  }
  
  min_y_long <- ifelse(min(data_long$logtiter) < 0, floor(min(data_long$logtiter)), 0) 
  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep) %>%
    calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(titer, 0), "\n", fold_change))
  
  min_y <- ifelse(min(sr_group_gmt_plotdata$logtiter) < 0, floor(min(sr_group_gmt_plotdata$lower)), 0) 
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  # add x position
  data_long <- x_position_spacing(data_long, position_by = x_position_by)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by)
  # fill  NA values of lower and upper CI with mean value
  sr_group_gmt_plotdata$lower[is.na(sr_group_gmt_plotdata$lower)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$lower)]
  sr_group_gmt_plotdata$upper[is.na(sr_group_gmt_plotdata$upper)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$upper)]
  
  
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", "", x)
    x <- gsub("B\\+O, B\\+O", "2x(B+O)", x)
    x <- arm_plot_names[x]
    x
  }
  
  if(length(color_by) > 1){
    color_by <- "sr_group"
  }
  
 
  
  # plot gmts of different serum groups together
  
  # max_y <- max(ceiling(sr_group_gmt_plotdata$logtiter), na.rm = T) + 1
  
  if(y_label == "Titer"){
    max_y <- 13
  } else {
    if(p_norm){
      max_y <- 3.5
    } else {
      max_y <- 4 +0.5
    }
  }
  
  spacing <- 0.4
  if(y_label != "Log2 titer difference") {
    spacing <- 0.6
  }
  
  if(length(unique(sr_group_gmt_plotdata$visit_diff)) >1){
    sr_group_gmt_plotdata$visit_diff <- factor(sr_group_gmt_plotdata$visit_diff, levels = visit_diff_order)
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
  
  sr_group_gmt_plotdata %>%
    ggplot(
      aes_string(
        x = "x",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = "inf_code"
      )
    ) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_line(
      aes_string(group = "sr_group", 
                 linetype = line_by),
      size = 1
    ) +
    # pointrange if I want them dodged, x position by. If I just want it dodged then I can add position = position_dodge(width = 1)
    # geom_pointrange(
    #   aes(ymin = lower,
    #       ymax = upper)
    # ) +
    # errorbar if I want shifted bars
    geom_errorbar(
      aes(ymin = lower,
          ymax = upper,
          group = visit_diff),
      width = 0,
      position = position_dodge(width = 0.4)
    ) +
    geom_point(
      size = 2
    ) +
    geom_point(
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype= "none") + 
    scale_x_continuous(
      limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(ag_order, na.rm = T))+0.5),
      name = "Days since vaccine"
    ) + 
    labs(
      y = "Boost magnitude (fold change)"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(min_y-0.5, max_y+0.5), 0.05)
    ) +
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
    )  -> gp_gmt
  
  # do the facets
  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
   
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"v_manuf_code", "count", "y_label_pos", color_by, shape_by, gmt_grid_col, gmt_grid_row))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    dots = sapply(c(gmt_grid_row, gmt_grid_col), . %>% {as.formula(paste0('~', .))})
    
    count_label <- count_label %>%
      group_by(.dots = dots)%>%
      summarize(total_count = mean(count),
                label = paste0("n = ", total_count),
                x = x,
                y = mean(y_label_pos)-0.5,
                inf_code = inf_code)
    
    #  count_label <- count_label[!duplicated(count_label$label),]
    
    
  } else {
   
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter, "count", "y_label_pos", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    count_label <- count_label %>%
      group_by(arm_code) %>%
      summarize(total_count = mean(count),
                label = paste0("n = ", total_count),
                x = x,
                y = mean(y_label_pos)-0.5,
                inf_code = inf_code)
    
    
  }
  
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt + 
      geom_text(aes_string(color = color_by) + aes(x = ifelse(age_code == "<65", x + point_spacing/2, x-point_spacing/2), y = y_label_pos,
                                                   label = fold_change),
                color = "black",
                size = 3) 
  }
  
  if(show_group_count){
    
    gp_gmt <- gp_gmt +
      geom_text(data = count_label,
                mapping = aes(x = x, y = y, label = label, hjust = 1),
                color = "black",
                size = 4
      )
  }
  
  if(y_label != "Log2 titer difference") {
  
    gp_gmt <- gp_gmt + scale_y_continuous(
      breaks = seq(round(min_y), max_y, by = 2),
      labels = function(x) ifelse(x <0, paste0("-", round(2^abs(x))*10), paste0(round(2^x)*10)),
      name = "Boost magnitude (linear)"
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
  
  if(show_mean_line){
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.2, alpha = 0.3)
    
  }
  
  return(list("gmt" = gp_gmt))
  
  
}


plot_for_plotly <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                          facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, 
                          color_by = "sr_group", shape_by = "inf_code", line_by = "age_code", x_position_by = "age_code",
                          cols_to_keep = c("arm_code", "age_code", "inf_code"), to_long = T, 
                          show_gmt_label = F, show_group_count = F, show_mean_line = F, mean_line_color = "black",
                          gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code",
                          pre_titer_adj = FALSE, single_visit = NULL) {
  
  ag_order <- ag_order[antigens]
  # format titer to longer
  if(to_long){
    data %>%
      pivot_longer(cols = antigens, names_to = "ag_name", values_to = "titer") %>%
      mutate(logtiter = log2(titer/10),
             logtiter = ifelse(titer < titer_thresh, log2((titer_thresh/20)), logtiter)) -> data_long
  } else{
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

  
 # plot_colors <-  plot_colors[as.character(unique(data_long[, color_by]))]
 
  # set sr_order
  if(!is.null(sr_group_order)) {
    data_long$sr_group <- factor(as.character(data_long$sr_group), levels =sr_group_order)
  }
  

#  max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
  max_y <- 15

  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep, pre_adj = pre_titer_adj) %>%
    filter(!all_below_thresh) -> sr_group_gmt_plotdata 

  if(!is.null(single_visit)){
    data_long <- data_long %>% filter(visit_code %in% single_visit)
    sr_group_gmt_plotdata %>% filter(visit_code %in% single_visit)
  }
   
   sr_group_gmt_plotdata <- sr_group_gmt_plotdata %>% calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(as.numeric(titer), 0), "\n", fold_change))
  
  
 # sr_group_gmt_plotdata <- sr_group_to_column(sr_group_gmt_plotdata)
  
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  data_long <- x_position_spacing(data_long, position_by = x_position_by)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by)
  
  # fill  NA values of lower and upper CI with mean value
  sr_group_gmt_plotdata$lower[is.na(sr_group_gmt_plotdata$lower)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$lower)]
  sr_group_gmt_plotdata$upper[is.na(sr_group_gmt_plotdata$upper)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$upper)]
  
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", " ", x)
    x <- arm_plot_names[x]
    x
  }

  # facet_labeller <- function(x){
  #   strsplit(x, "-")[[1]][2]
  # }
  sub_visit <- FALSE
  if(length(color_by) > 1){
    if("visit_code" %in% color_by | length(unique(data$visit_code)) > 1){
      sub_visit <- TRUE
    }
    color_by <- "sr_group"
  }

  data_long$color <- plot_col
  
  data_long %>%
    ggplot(
      aes_string(
        x = "ag_name",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) + 
    geom_line(
      aes(group = sr_name),
      linetype = "solid",
      alpha = 0.4
    ) + 
    geom_point(
      aes(group = sr_name),
      alpha = 0.4
    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
      aes_string(group = "sr_group",
          linetype = line_by),
      size = 1.3
    ) +
    geom_pointrange(
      data = sr_group_gmt_plotdata,
      aes(ymin = lower,
            ymax = upper)
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = ag_name, y = y, label = label),
              color = "black",
              size = 2.5
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
    scale_x_discrete(
      limits = antigens,
      labels = ag_plot_names[antigens],
      name = "Antigen variant"
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+1), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 10)
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    ) -> gp
  
  return(gp)
  visits <- unique(sr_group_gmt_plotdata$visit_code)
  
  if(line_by != "v_manuf_code" & "V3" %in% visits & "V4" %in% visits){
    line_by <- "visit_code"
  }
  # plot gmts of different serum groups together
  max_y <- 13
  sr_group_gmt_plotdata$y <- max_y #- 0.5
  
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
    geom_line(
      aes_string(group = "sr_group",
          linetype = line_by),
      size = 1
    ) +
    geom_pointrange(
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
    scale_x_continuous(
      breaks = ag_order,
      labels = ag_plot_names[names(ag_order)],
      limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(ag_order, na.rm = T))+0.5),
      name = "Antigen variant"
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+0.5), 0.05)
    ) +
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
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    ) -> gp_gmt
  

  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
    gp <- gp + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"sr_group", "arm_code", "v_manuf_code", "count", "y", "age_code","vistic_code", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    if(sub_visit){
      visits <- unique(sr_group_gmt_plotdata$visit_code)
      
      target_visit <- visits[visits != "D15"]
     
      #   count_label <- count_label[grepl(target_visit[1], count_label$sr_group),]
      count_label <- count_label[target_visit[1] == count_label$visit_code,]
    }
    
    count_label <- count_label %>%
      group_by(sr_group,arm_code, v_manuf_code, age_code) %>%
      summarize(total_count = sum(count),
                label = paste0("n (", age_code, ") = ", total_count),
                x = x,
                y = ifelse(age_code == "combined", y, ifelse(age_code == "<65", y - 1, y-2)),
                inf_code = inf_code)
    
    count_label <- count_label %>%
      mutate(y = ifelse(inf_code == "inf", y, y-1))
    
    count_label$label <- gsub("\\(combined\\)", "", count_label$label)
      
    
  #  count_label <- count_label[!duplicated(count_label$arm_code),]
    
    
  } else {
    gp <- gp + 
      facet_wrap(
        vars(sr_group),
        labeller = as_labeller(facet_labeller),
        nrow = facet_n_row
      ) 
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"sr_group", "arm_code", "count", "y", "age_code", "visit_code", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    if(sub_visit){
      visits <- unique(sr_group_gmt_plotdata$visit_code)
    
      target_visit <- visits[!("D15" == visits)]
   #   count_label <- count_label[grepl(target_visit[1], count_label$sr_group),]
      count_label <- count_label[target_visit[1] == count_label$visit_code,]
    }
    
    count_label <- count_label %>%
      group_by(sr_group, arm_code, age_code) %>%
      summarize(total_count = sum(count),
                label = paste0("n (", age_code, ") = ", total_count),
                x = x,
                y = ifelse(age_code == "combined", y, ifelse(age_code == "<65", y - 1, y-2)),
                inf_code = inf_code)
    
    count_label <- count_label %>%
      mutate(y = ifelse(inf_code == "inf", y, y-1))
    
    count_label$label <- gsub("\\(combined\\)", "", count_label$label)
    
   # count_label <- count_label[!duplicated(count_label$arm_code),]
  }
  
  if(shape_by == "inf_code"){
    gp <- gp + 
      scale_shape_manual(values = inf_code_shapes)
    
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt +
      geom_text(data = sr_group_gmt_plotdata,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = ag_order[ag_name], y = y, label = label),
                color = "black",
                size = 2.5
      )
  }
  
  if(show_group_count){
   
    gp_gmt <- gp_gmt +
      geom_text(data = count_label,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = x, y = y+0.5, label = label, hjust = 1),
                color = "black",
                size = 4
      )
  }
  
  if(show_mean_line){
    gp <- gp +
      geom_line(aes_string(y = "gmt_arm", group = "sr_group", linetype = line_by), color = mean_line_color, 
                size = 1.2, alpha = 0.3)
    
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.2, alpha = 0.3)
    
  }
  
  return(list("all" = gp, "gmt" = gp_gmt))
  
  
}



count_label_function <- function(sr_group_gmt_plotdata, ymax, gmt_facetter, color_by, shape_by, gmt_grid_row, cols_to_keep = c()){
  
  
  visits_in_data <- visit_code_order[visit_code_order %in% sr_group_gmt_plotdata$visit_code]
  inf_code_in_data <- inf_code_order[inf_code_order %in% sr_group_gmt_plotdata$inf_code]
  arms_in_data <- arm_code_order[arm_code_order %in% sr_group_gmt_plotdata$arm_code]
  
  count_label <- sr_group_gmt_plotdata %>% 
    ungroup() %>%
    select(all_of(unique(c(gmt_facetter,"sr_group", "arm_code", "v_manuf_code", "count", "y", "age_code", "visit_code", color_by, shape_by, "inf_code"), cols_to_keep))) -> sub
  
  if("visit_diff" %in% cols_to_keep){
    
    sub %>%
      group_by(sr_group,arm_code, v_manuf_code, age_code, visit_code, visit_diff) %>%
      summarize(x = max(sr_group_gmt_plotdata$x),
                pretty_visit = arm_plot_names[as.character(visit_code)],
                label = paste0("n(", pretty_visit, ifelse(grepl("breakthrough", as.character(inf_code)), " BT", ""), ") = ", count),
                y = ymax - (match(visit_code, visits_in_data)+1)/2,
                y = y - (match(inf_code, inf_code_in_data)+1)/2,
                inf_code = inf_code) %>%
      unique() -> count_label
    
  } else if(color_by == "arm_code"){
   
    sub %>%
      group_by(sr_group,arm_code, v_manuf_code, age_code, visit_code) %>%
      summarize(x = max(sr_group_gmt_plotdata$x),
                pretty_visit = arm_plot_names[as.character(visit_code)],
                label = paste0("n(", arm_code, ifelse(grepl("breakthrough", as.character(inf_code)), " BT", ""), ") = ", count),
                y = ymax - (match(arm_code, arms_in_data)+1)/2,
                y = y - (match(inf_code, inf_code_in_data)+1)/2,
                inf_code = inf_code) %>%
      unique() -> count_label
    
    
  } else {
    
    sub %>%
    group_by(sr_group,arm_code, v_manuf_code, age_code, visit_code) %>%
      summarize(x = max(sr_group_gmt_plotdata$x),
                pretty_visit = arm_plot_names[as.character(visit_code)],
                label = paste0("n(", pretty_visit, ifelse(grepl("breakthrough", as.character(inf_code)), " BT", ""), ") = ", count),
                y = ymax - (match(visit_code, visits_in_data)+1)/2,
                y = y - (match(inf_code, inf_code_in_data)+1)/2,
                inf_code = inf_code) %>%
      unique() -> count_label
  }
    
  
  if(gmt_grid_row == "visit_code"){
    count_label$y <- ymax
  }
  return(count_label)
}
