padding <- 0.25
landscape_opacity <- 1
pre_opacity <- 0.5
# Functions to remove buttons
addObject3js <- function(
    data3js,
    object,
    number_of_ids = 1
){
  
  # Generate an object ID
  if(is.null(data3js$lastID)){ data3js$lastID <- 0 }
  object$ID <- max(data3js$lastID) + seq_len(number_of_ids)
  
  # If object is interactive and highlighted add a reference to itself to
  # it's highlight group by default
  if(!is.null(object$properties$interactive)){
    object$group <- object$ID
  }
  
  # Add the object to the plot data
  data3js$plot[[length(data3js$plot)+1]] <- object
  
  # Update the ID of the last object added
  data3js$lastID <- object$ID
  
  # Return the new data
  data3js
  
}

remove_buttons <- function(data3js){
  
  new_data3js = data3js
  
  new_data3js = data3js
  
  new_data3js[['lastID']] = 0
  new_data3js[['plot']] = list()
  
  N = data3js[['lastID']] 
  
  
  
  
  for (i in 1:N)
  {
    obj = data3js[['plot']][[i]]
    
    
    
    if ('toggle' %in% names(obj[['properties']])){
      obj[['properties']][['toggle']] <- NULL
    }
    
    new_data3js = addObject3js(new_data3js,obj)
    
    
  }
  
  
  
  return (new_data3js)
  
}





# sams landscape functions to add landscape from lndscp fits list
get_titertable <- function(data, group) {
  
  data %>% 
    select(
      variant,
      subjid,
      titer
    ) %>%
    mutate(
      titer = replace(titer, is.na(titer), "*")
    ) %>%
    pivot_wider(
      id_cols = subjid,
      names_from = variant,
      values_from = titer
    ) %>% 
    as.matrix() -> titermatrix
  
  attr(titermatrix, "arm") <- group$arm_code
  attr(titermatrix, "visit") <- group$visit_code
  rownames(titermatrix) <- titermatrix[,"subjid"]
  titermatrix <- titermatrix[,-1]
  
#  print(titermatrix)
#  print(titermatrix[titermatrix[,"BA.4/BA.5"] != "*",,drop=F])
#  titermatrix[titermatrix[,"BA.4/BA.5"] != "*",,drop=F]
#  titermatrix[titermatrix[,"BA.4/BA.5(2)"] != "*",,drop=F]
  
  return(titermatrix)
  
}

base_plot_data3js <- function(map, lndscp_fits, highlighted_ags, lims, ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1,
                              add_border = TRUE, add_axis = TRUE, max_z = 12){
  
  if(alternative_ba5){
    # set points and coordinates of highlighted ags
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags)), 12) / 5
    ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags], agOutline(map)[agNames(map) =="BA.4/BA.5"])
    ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags], agFill(map)[agNames(map) =="BA.4/BA.5"])
    labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags], "BA.4/BA.5(2)")
    
  } else if("BA.4/BA.5" %in% agNames(map) ) {
    
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    x_coords["BA.4/BA.5"]<- agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1]
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    y_coords["BA.4/BA.5"]<- agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2]
    z_coords <- rep(0.02, length(highlighted_ags)-1)
    ag_point_size <- c(rep(14, length(highlighted_ags)-1)) / 5
    ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
    ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
    labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
    
  } else {
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
    ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
    ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
    labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
  }
  
  border_col <- "grey50"
  
  z_lims <- c(0,max_z)
  axis_at <- seq(z_lims[1], z_lims[2],ifelse(max_z >9, 2, 1))
  # Setup plot
  data3js <- ablandscapes:::lndscp3d_setup(
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    aspect.z = 0.5,
    options = list(
      lwd.grid =  0.05,
      sidegrid.lwd = 1,
      sidegrid.col = border_col,
      sidegrid.at = list("z" = axis_at),
      zaxt = "log"
    ),
    show.axis = FALSE
  )
  
  if(add_axis){

    axis_labels <- 2^axis_at*10
    
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
    # labeloffset = 0.11,
      cornerside = "f",
      size = 20,
      alignment = "right"
    )
  }

  # Add basemap
  data3js <- lndscp3d_map(
    data3js = data3js,
    fit = lndscp_fits[[1]],
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = c(0, 10),
    show.map.sera = FALSE,
    options = list(
      opacity.basemap = 0.3
    )
  )
  
  data3js <- r3js::points3js(
    data3js,
    x          = x_coords,
    y          = y_coords,
    z          = z_coords,
    size       = ag_point_size,
    col        = ag_col,
    fill       = ag_fill,
    lwd        = 0.5,
    opacity    = 1,
    highlight  = list(col = "red"),
    label      = labels,
    toggle     = "Basepoints",
    depthWrite = FALSE,
    shape      = "circle filled"
  )
  
  if(add_border){
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[1]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[2],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    # y border
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[1]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[2], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)

    data3js <- r3js::box3js(
      data3js,
      col   = border_col
    )
    
  }
  
  return(data3js)
}

plot_idvl_landscapes_from_list <- function(data3js, idvl_landscapes, sr_colors){
 
  for (x in seq_along(idvl_landscapes)) {
    
    surface_options <- list()
    surface_options$col.surface = sr_colors[x]
    surface_options$col.surface.grid = adjustcolor(
      "grey80",
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    )
    surface_options$opacity.surface = 0.2
    
    data3js <- lndscp3d_surface(
      data3js = data3js,
      object = idvl_landscapes[[x]],
      toggle = x,
      options = surface_options,
      crop2chull = FALSE,
      grid_spacing = 0.5,
      padding = padding
    )
    
  }
  
  return(data3js)
}


plot_landscapes_from_list <- function(data3js, titertables_groups, lndscp_fits,map, gmt_data, highlighted_ags,
                                      ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1, hide_buttons = TRUE,
                                      color_by = "arm_code"){
  
  if(alternative_ba5){
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags)), 12) / 5
   # text_x <- c(c(x_coords[1:4] + ag_point_size[1:4]*0.15),c(x_coords[5:6] - ag_point_size[5:6]*0.25))
  #  text_y <- c(y_coords[1:4], c(y_coords[5:6] - ag_point_size[5:6]*0.2))
    text_x <- c(x_coords[1:6] - ag_point_size[1:6]*0.2)
    text_y <- c(y_coords[1:6] - ag_point_size[1:6]*0.2)
    text_plot <- c(ag_plot_names[agNames(map)[agNames(map) %in% highlighted_ags]], "BA.4/BA.5(2)")
    
  } else if("BA.4/BA.5" %in% agNames(map) ) {
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    x_coords["BA.4/BA.5"]<- agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1]
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    y_coords["BA.4/BA.5"]<- agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2]
    z_coords <- rep(0.02, length(highlighted_ags)-1)
    ag_point_size <- c(rep(14, length(highlighted_ags)-1)) / 5
  #  text_x <- c(c(x_coords[1:3] + ag_point_size[1:3]*0.15),c(x_coords[4:5] - ag_point_size[4:5]*0.2))
  #  text_y <- c(y_coords[1:3], c(y_coords[4:5] - ag_point_size[4:5]*0.2))
    text_x <- c(x_coords[1:5] - ag_point_size[1:5]*0.2)
    text_y <- c(y_coords[1:5] - ag_point_size[1:5]*0.2)
    text_plot <- c(ag_plot_names[agNames(map)[agNames(map) %in% highlighted_ags]])
    
  } else {
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
    ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
    ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
    labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
  }
  
  if(length(titertables_groups) >1){
     min_offset <- -0.1
    max_offset <- 0.1
    offset <- seq(from = min_offset, to = max_offset, by = (max_offset - min_offset)/length(lndscp_fits))
    
  } else {
    offset <- 0
  }
   
  for (i in seq_along(lndscp_fits)) {
    
   # message(i)
    arm <- titertables_groups$arm_code[i]
    visit <- titertables_groups$visit_code[i]
    v_manuf <- titertables_groups$v_manuf_code[i]
    lndscp_fit <- lndscp_fits[[i]]
    bt <- FALSE
  
    coords <- cbind(x_coords, y_coords)
    
    # Add titers
    
    gmts <- filter(gmt_data, visit_code == visit, arm_code == arm, v_manuf_code == v_manuf)

    # if (is.na(filter(gmts, variant == "BA.4/BA.5")$gmt)){
    #   gmts <- gmts[!is.na(gmts$gmt),]
    #   coords <- coords[rownames(coords) != "BA.4/BA.5",]
    # } # next
    
    bt_in_gmt <- TRUE %in% grepl("breakthrough", colnames(gmt_data))
    if(bt_in_gmt){
      bt_gmts <- gmts %>% 
        filter(breakthrough)
      bt_gmts <- bt_gmts[match(rownames(coords), bt_gmts$variant),]
      bt <- titertables_groups$breakthrough[i]
    }
    
    gmts <- gmts[match(rownames(coords), gmts$variant),]
    
    
    for (j in seq_len(nrow(coords))) {
      
      data3js <- r3js::lines3js(
        data3js,
        x = rep(coords[j, 1], 2),
        y = rep(coords[j, 2], 2),
        z = c(0, ifelse(bt_in_gmt, max(bt_gmts$gmt[j], gmts$gmt[j]), gmts$gmt[j])),
        col = "grey50",
        highlight = list(col = "red"),
        interactive = FALSE,
     #   toggle = sprintf("Titers, %s, %s", arm, visit),
        geometry = TRUE,
        lwd = 0.4 #was 1
      )

     # print(coords[j,])
     # print(offset[j])
     # print(gmts$gmt[j])
      
      data3js <- r3js::points3js(
        data3js,
        x         = coords[j, 1] + offset[i],
        y         = coords[j, 2],
        z         = gmts$gmt[j],
      #  size      = 0.7, #was 2
        size      = 0.9, #was 2
      #  col       = "grey50",
      col  = ifelse(color_by == "arm_code", ifelse(visit == "D1", "grey80", arm_cols[arm]), arm_cols[v_manuf]), 
      opacity   = ifelse(color_by == "arm_code", ifelse(visit == "D1", pre_opacity, landscape_opacity), landscape_opacity) # was 1
      
      )
      
      if(bt_in_gmt){
        data3js <- r3js::points3js(
          data3js,
          x         = coords[j, 1] + offset[i],
          y         = coords[j, 2],
          z         = bt_gmts$gmt[j],
          #  size      = 0.7, #was 2
          size      = 0.9, #was 2
          #  col       = "grey50",
          col  = ifelse(color_by == "arm_code", ifelse(visit == "D1", "grey80", arm_cols[arm]), arm_cols[v_manuf]), 
          opacity   = ifelse(color_by == "arm_code", ifelse(visit == "D1", pre_opacity, landscape_opacity), landscape_opacity) # was 1
          
        )
      }
      
     # text_x <- c(c(agCoords(map)[agNames(map) %in% highlighted_ags[1:4], 1] + agSize(map)[agNames(map) %in% highlighted_ags[1:4]]*0.025), c(agCoords(map)[agNames(map) == "BA.4/BA.5", 1]- agSize(map)[agNames(map) == "BA.4/BA.5"]*0.07))
    #  text_y <- c(c(agCoords(map)[agNames(map) %in% highlighted_ags[1:4], 2]), c(agCoords(map)[agNames(map) == "BA.4/BA.5", 2] - agSize(map)[agNames(map) == "BA.4/BA.5"]*0.07))
      
      # set points and coordinates of highlighted ags

      
      data3js <- r3js::text3js(
        data3js,
        x          = text_x,
        y          = text_y,
        z          = z_coords,
        text       = text_plot,
        toggle     = "Labels",
        size       = c(rep(16*0.02, length(text_x))), #agSize(map)[agNames(map) %in% highlighted_ags]*0.02,
        alignment  = "right"
      )
      
    }
    
    # Add landscapes
    data3js <- lndscp3d_surface(
      data3js = data3js,
      object = lndscp_fit,
      # zlim = c(0, 10),
      crop2chull = FALSE,
      # crop2base = TRUE,
      toggle = sprintf("Landscape, %s, %s", arm, visit),
      grid_spacing = 0.5,
      padding = padding,
      options = list(
        opacity.surface = ifelse(color_by != "arm_code", 1, ifelse(visit == "D1", pre_opacity, landscape_opacity)), # ifelse(visit == "D1", 0.3, 1),
        col.surface = ifelse(color_by == "arm_code", ifelse(visit == "D1", "grey", arm_cols[arm]), arm_cols[v_manuf]) #ifelse(visit == "D1", "grey80", arm_cols[arm])
       # opacity.surface = 0.5
      )
    )
    
  }
  
  if(hide_buttons){
    data3js <- remove_buttons(data3js)
  }
 
  
  
  return(data3js)
}



plot_flu_landscapes_from_list <- function(data3js, titertables_groups, lndscp_fits,map, gmt_data, highlighted_ags,
                                      ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1, hide_buttons = TRUE, sr_group_cols){
  
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
    ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
    ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
    labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
  
  if(length(lndscp_fits) >1){
     min_offset <- -0.1
    max_offset <- 0.1
    offset <- seq(from = min_offset, to = max_offset, by = (max_offset - min_offset)/length(lndscp_fits))
    
  } else {
    offset <- 0
  }
   
  for (i in seq_along(lndscp_fits)) {
    
   # message(i)
    sr_group_temp <- titertables_groups$sr_group[i]
    
    lndscp_fit <- lndscp_fits[[i]]
    
    
    
    coords <- cbind(x_coords, y_coords)
    
    # Add titers
    gmts <- filter(gmt_data, sr_group == sr_group_temp)
    
    gmts <- gmts[match(rownames(coords), gmts$variant),]
    
 
    for (j in seq_len(nrow(coords))) {
      
      data3js <- r3js::lines3js(
        data3js,
        x = rep(coords[j, 1], 2),
        y = rep(coords[j, 2], 2),
        z = c(0, gmts$gmt[j]),
        col = "grey50",
        highlight = list(col = "red"),
        interactive = FALSE,
       # toggle = sprintf("Titers, %s", sr_group_temp),
        geometry = TRUE,
        lwd = 0.4 #was 1
      )
      
      data3js <- r3js::points3js(
        data3js,
        x         = coords[j, 1] + offset[i],
        y         = coords[j, 2],
        z         = gmts$gmt[j],
      #  size      = 0.7, #was 2
        size      = 0.9, #was 2
      #  col       = "grey50",
        toggle = sprintf("Titers, %s", sr_group_temp),
        col  = sr_group_cols[sr_group_temp],
        opacity   = 1 # was 1
     
      )
      
    }
   
    # Add landscapes
    data3js <- lndscp3d_surface(
      data3js = data3js,
      object = lndscp_fit,
      # zlim = c(0, 10),
      crop2chull = FALSE,
      # crop2base = TRUE,
      grid_spacing = 0.5,
      padding = padding,
      toggle = sprintf("Landscape, %s", sr_group_temp),
      options = list(
        opacity.surface = 0.5,
        col.surface = sr_group_cols[sr_group_temp]
       # opacity.surface = 1
      )
    )
   
  }
  
  if(hide_buttons){
    data3js <- remove_buttons(data3js)
  }
 
  
  
  return(data3js)
}


# GMT landscape plotting
plot_D1_v2_gmt_landcapes <- function(lndscp_fits, titertables_groups, day_visno, map,
                                     highlighted_ags, lims, ag_plot_names, gmt_data, angle, figure_dir,
                                     save_html_landscape = TRUE, inf_stat, lndscp_name = "_gmt_landscapes"){
  
  lndscp_list <- list()
  # Add landscapes
  for(v_manuf in unique(titertables_groups$v_manuf_code)){
    
    for(v1 in c("D1")){
      
      visit_rows <- c(grep("TRUE", v1 == titertables_groups$visit_code))
      manuf_rows <- grep(v_manuf, titertables_groups$v_manuf_code)
      
      lndscp_fits_t <- lndscp_fits[intersect(manuf_rows, visit_rows)]
      titertables_groups_t <- titertables_groups[intersect(manuf_rows, visit_rows),]
      
      
      if(v1 != last(day_visno)){
        
        for(v2 in day_visno[c((grep("TRUE", v1 == day_visno)+1):length(day_visno))]) {
          # get matching lndscp_fits
          visit_rows <- c(grep("TRUE", v1 == titertables_groups$visit_code), grep("TRUE", v2 == titertables_groups$visit_code))
          manuf_rows <- grep(v_manuf, titertables_groups$v_manuf_code)
          
          #remove B+O, B+O D1 landscape
          
          lndscp_fits_t <- lndscp_fits[intersect(manuf_rows, visit_rows)]
          titertables_groups_t <- titertables_groups[intersect(manuf_rows, visit_rows),]
          
          data3js <- base_plot_data3js(map, lndscp_fits, highlighted_ags, lims, ag_plot_names)
          lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups_t, lndscp_fits_t, map, gmt_data, highlighted_ags, ag_plot_names)
          
          
          lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups_t, lndscp_fits_t, map, gmt_data, highlighted_ags, ag_plot_names)
          
          
          lndscp <-r3js(
            lndscp_3js,
            rotation = angle$rotation,
            zoom = angle$zoom
          )
          
          lndscp_list[[paste0(inf_stat, "_",v1, "_",v2, "_", v_manuf)]] <- lndscp
          if(save_html_landscape){
            save_name <- file.path(figure_dir,paste0(v1, "_", v2, "_", v_manuf,lndscp_name))
            screenshot_html_landscape_to_png(lndscp, save_name)
          }
         
    
         
        }
        
      }
      
    }
    
  }
  
  
  return(lndscp_list)
}


# plot individual landscapes
plot_d1_v2_idvl_arm_landscapes <- function(lndscp_fits, titertables_groups, day_visno, map,
                                       highlighted_ags, lims, ag_plot_names, gmt_data, angle, figure_dir,
                                       individual_lndscp_fits, individual_sr_colors,
                                       arm_code_order_no_vacc, save_html_landscape = TRUE, inf_stat,
                                       lndscp_name = "_idvl_landscapes"){
  
  lndscp_list <- list()
  # now do landscapes of arms, both visits
  for(v_manuf in unique(titertables_groups$v_manuf_code)){
    
    for(v1 in c("D1")){
      
      if(v1 != last(day_visno)){
        
        for(v2 in day_visno[c((grep("TRUE", v1 == day_visno)+1):length(day_visno))]) {
          # get matching lndscp_fits
          visit_rows <- c(grep("TRUE", v1 == titertables_groups$visit_code), grep("TRUE", v2 == titertables_groups$visit_code))
          manuf_rows <- grep(v_manuf, titertables_groups$v_manuf_code)
          
          titertables_groups_temp <- titertables_groups[intersect(manuf_rows, visit_rows),]
          lndscp_fits_temp <- lndscp_fits[intersect(manuf_rows, visit_rows)]
          
          lndscp_fits_idvl_t <- individual_lndscp_fits[intersect(manuf_rows, visit_rows)]
          colors_idvl_t <- individual_sr_colors[intersect(manuf_rows, visit_rows)]
          
          individual_landscape_plots <- list()
          
          for(arm in arm_code_order_no_vacc){
            
            arm_rows <- arm==titertables_groups_temp$arm_code
            
            if(!(TRUE %in% arm_rows)){
              next()
            } else {
              
              lndscp_fits_t <- lndscp_fits_temp[arm_rows]
              titertables_groups_t <- titertables_groups_temp[arm_rows,]
              sr_colors <- colors_idvl_t[arm_rows][[1]]
              lndscp_fits_idvl <- lndscp_fits_idvl_t[arm_rows] #[[1]]
              
              null_lndscps <- lapply(lndscp_fits_idvl, function(x){
                is.null(x)
              })
              
              lndscp_fits_idvl <- lndscp_fits_idvl[!unlist(null_lndscps)][[1]]
             
              data3js <- base_plot_data3js(map, lndscp_fits, highlighted_ags, lims, ag_plot_names)
              data3js <- plot_idvl_landscapes_from_list(data3js, lndscp_fits_idvl, sr_colors)
              lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups_t, lndscp_fits_t, map, gmt_data, highlighted_ags, ag_plot_names)
              
              lndscp <-r3js(
                lndscp_3js,
                rotation = angle$rotation,
                zoom = angle$zoom
              )
              
              lndscp_list[[paste0(inf_stat, "_",v1, "_",v2, "_", v_manuf, "_", arm)]] <- lndscp
              
              if(save_html_landscape){
                save_name <- file.path(figure_dir,paste0(v1, "_", v2, "_", v_manuf,"_",arm, lndscp_name))
                screenshot_html_landscape_to_png(lndscp, save_name)
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  return(lndscp_list)
}

fit_landscapes <- function(sr_group_data, highlighted_ags, map, group_vars = c("visit_code","arm_code", "v_manuf_code")){
  
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
    group_by_at(vars(all_of(group_vars))) -> titerdata
  
  #titerdata <- titerdata %>% filter(!(variant %in% c("BA.4/BA.5", "BA.4/BA.5(2)")))
  titerdata %>%
    group_map(
      get_titertable
    ) -> titertables
  
  lndscp_fits <- lapply(
    titertables,
    function(titertable) {
     
      if(is.null(dim(titertable))){
        titert <- matrix(titertable, nrow = 1, byrow = TRUE)
        colnames(titert) <- names(titertable)
        titertable <- titert
      }
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
  
  titertables_groups <- group_data(titerdata)
  
  # Add impulses
  titerdata %>%
    group_by_at(vars(all_of(c(group_vars, "variant")))) %>%
    summarize(gmt = titertools::gmt(titer, dilution_stepsize = 0)["mean", "estimate"])-> gmt_data
  
  
  return(list("titerdata" = titerdata, "titertables_groups" = titertables_groups, 
              "titertables" = titertables, "gmt_data" = gmt_data,
              "lndscp_fits" = lndscp_fits))
}


fit_idvl_landscapes <- function(titertables, titertable_groups, map, only_breakthroughs = TRUE){
  individual_lndscp_fits <- list()
  
  for(table_nr in 1:length(titertables)){
    
    if(only_breakthroughs){
      
      if(titertable_groups$breakthrough[table_nr]){
        
        if(is.null(dim(titertables[[table_nr]]))){
          titert <- matrix(titertables[[table_nr]], nrow = 1, byrow = TRUE)
          colnames(titert) <- names(titertables[[table_nr]])
          titertables[[table_nr]] <- titert
        }
        
        individual_lndscp_fits[[table_nr]] <- lapply(1:nrow(titertables[[table_nr]]), function(sr){
          
          ablandscape.fit(
            titers    = titertables[[table_nr]][sr, c("BA.1", "B.1.617.2", "B.1.351", "D614G", "BA.4/BA.5"), drop = F],
            acmap     = map,
            bandwidth = 1,
            degree    = 1,
            method = "cone",
            error.sd = 1,
            control = list(
              optimise.cone.slope = TRUE
            )
          )
        })
        
      }
      
    } else {
      individual_lndscp_fits[[table_nr]] <- lapply(1:nrow(titertables[[table_nr]]), function(sr){
        
        ablandscape.fit(
          titers    = titertables[[table_nr]][sr, c("BA.1", "B.1.617.2", "B.1.351", "D614G", "BA.4/BA.5"), drop = F],
          acmap     = map,
          bandwidth = 1,
          degree    = 1,
          method = "cone",
          error.sd = 1,
          control = list(
            optimise.cone.slope = TRUE
          )
        )
      })
    }
    
    
  }
  
  return(individual_lndscp_fits)
}


lndscp_combo_fit <- function(sr_group_data, lndscp_combo, highlighted_ags, map){
  
  if(lndscp_combo == "_all_pre_post_BT"){
    sr_group_data <- sr_group_data %>%
      mutate(breakthrough = VISNO >= btfl_visno)
    
  } else if(lndscp_combo == "_all_pre_pre_BT"){
    
    visno_order <- sort(unique(sr_group_data$VISNO))
    
    sr_group_data <- sr_group_data %>%
      mutate(breakthrough = VISNO == btfl_visno,
            btfl_visno = ifelse(breakthrough, visno_order[match(btfl_visno, visno_order)-1], btfl_visno),
             breakthrough = VISNO == btfl_visno) %>%
      filter(VISNO <= btfl_visno)
  }
  
  fit_result <- fit_landscapes(sr_group_data, highlighted_ags, map, group_vars = c("visit_code", "arm_code", "v_manuf_code", "breakthrough"))
  
  return(fit_result)
}
