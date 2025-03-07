library(r3js)

source("./functions/remove_reactivity_bias.R")

# set plotting parameters
plot_zlim <- c(-1, 17)

plot_xlim <-  read.csv("./data/metadata/xlim_zoom.csv")$x
plot_xlim[2] <- plot_xlim[2] + 1
plot_ylim <-  read.csv("./data/metadata/ylim_zoom.csv")$x

plot_xlim <- c(-3.435233, 7.564767)
plot_ylim <- c(-3.844888, 3.155112)


# for rotation towards left rotation = c(-1.537, 0, -0.8350)
# Set viewing angles
angle <- list(
  rotation = c(-1.437, 0, -0.5350),# rotation = c(-1.4370, 0.0062, -0.5350),
  translation = c(0, 0.05,0.1), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.3342
 # zoom = 1.1646 # higher is more zoomed out
)

# function to fit slope and serum coordinates for multi exposure landscapes
fit_cone_all_sera <- function(
    pars,
    ag_coords,
    log_titers,
    colbase
) {
  
  if(length(colbase)>1) {
    pred_titers_all <- unlist(lapply(1:length(colbase), function(ser) {
      coords <- c(pars[paste0("x",ser)], pars[paste0("y",ser)])
      ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
      predicted_logtiters <- colbase[ser] - ag_distances*pars["slope"]
      return(predicted_logtiters)
    }))
  } else {
    coords <- c(pars["x"], pars["y"])
    ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
    pred_titers_all <- colbase - ag_distances*pars["slope"]
  }
  sum((log_titers - pred_titers_all)^2, na.rm = T)
  
}

# fit for each serum group level and remove reactivity bias
ablandscape_coordinates_slope_fit <- function(ag_coords, log_titers, remove_bias = TRUE){
  
  log <-as.matrix(log_titers)
  
  unadj_log <- log
  unadj_colb <- apply(log, 2, function(x) max(x, na.rm = T))
  
  # remove reactivity bias for each serum
  if(remove_bias){
    log <- remove_reactivity_bias_logtiter(log)
    colb <- apply(log, 2, function(x) max(x, na.rm = T))
  }
  
  
  fit <-  optim(
    par = c(
      x = unname(ag_coords[rownames(log)[apply(log, 2, which.max)],1]),
      y = unname(ag_coords[rownames(log)[apply(log, 2, which.max)],2]),
      slope = 1
    ),
    fn = fit_cone_all_sera,
    method = "L-BFGS-B",
    ag_coords = ag_coords,
    log_titers = log,
    colbase = colb,
    control = list(maxit = 500)
  )
  
  pars <- fit[["par"]]
  fit$slope <- pars["slope"]
  
  fit$sr_cone_coords <- as.matrix(cbind(pars[grepl("x", names(pars))], pars[grepl("y", names(pars))]))
  
  fit$colbases_adj <- colb
  
  fit$log_titers_adj <- log
  
  fit$colbases <- unadj_colb
  
  fit$log_titers <- unadj_log
  
  fit$bias_removed <- remove_bias
  
  return(fit)
}


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


# fit the z value per landscape
fit_lndscp_val <- function(x, y, sr_cone_coords, sr_colbases, sr_cone_slopes) {
  
  sr_distances <- as.matrix(dist(rbind(c(x, y), sr_cone_coords)))[1, -1]
  mean(sr_colbases - sr_distances*sr_cone_slopes)
  
}

obtain_xyz_matrix <- function(lndscp_xlim, lndscp_ylim, sr_group_cone_coords, sr_group_colbases, sr_group_cone_slopes){
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  
  grid_z_matrix[] <- vapply(
    seq_len(length(grid_z_matrix)),
    \(n) {
      fit_lndscp_val(
        grid_x_matrix[n], grid_y_matrix[n],
        sr_group_cone_coords,
        sr_group_colbases,
        sr_group_cone_slopes
      )
    }, numeric(1)
  )
  
  return(list("grid_x_matrix" = grid_x_matrix, "grid_y_matrix" = grid_y_matrix, "grid_z_matrix" = grid_z_matrix))
}

plot_base_map <- function(map, sr_group){
  
  # # Work out which sera are out of bounds
  sr_x_coords <- srCoords(map)[,1]
  sr_y_coords <- srCoords(map)[,2]
  margin <- 0.2
  sr_out_of_bound <- sr_x_coords < plot_xlim[1] + margin |
    sr_x_coords > plot_xlim[2] - margin |
    sr_y_coords < plot_ylim[1] + margin |
    sr_y_coords > plot_ylim[2] - margin
  
  # Set map subset
  map_subset <- map
  srShown(map_subset)[sr_out_of_bound] <- FALSE
  srShown(map_subset)[srGroups(map_subset) != sr_group] <- FALSE
  
  srShown(map_subset) <- FALSE
  
  # Plot the base plot
  data3js <- lndscp_r3jsplot(
    fit = list(acmap = map_subset),
    aspect.z = 0.5,
    show.surface = FALSE,
    show.titers = FALSE,
    output = "data3js",
    xlim = plot_xlim,
    ylim = plot_ylim,
    zlim = plot_zlim,
    show.sidegrid = TRUE,
    show.axis = FALSE,
    options = list(
      opacity.basemap.ags = 1,
      cex.basemap.ags = 3,
      cex.basemap.sr = 1.5,
      # lwd.grid             = 1.5,
      opacity.basemap.sr = 1
    )
  )
  
  return(data3js)
  
}

add_landscape_plot <- function(data3js, grid_x_matrix, grid_y_matrix, grid_z_matrix, color, landscape_name = "Mean landscape", opacity_val= 0.8){
  
  # Add the surface
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = color,
    opacity = opacity_val,
    toggle = landscape_name,
    wireframe = FALSE,
    doubleSide = TRUE
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = adjustcolor(
      color,
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    ),
    opacity = opacity_val,
    toggle = landscape_name,
    wireframe = TRUE,
    doubleSide = TRUE
  )
  
  return(data3js)
  
}

add_titer_50_plane <- function(data3js, z_height = 5, color = "grey80", opacity_val = 0.2, name = "Titer 50"){
  # Add the titer 50 plane
  x_grid <- seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.5)
  y_grid <- seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.5)
  z_grid <- matrix(log2(z_height), length(x_grid), length(y_grid))
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = color,
    opacity = opacity_val,
    toggle = name
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = color,
    opacity = opacity_val*2,
    toggle = name,
    wireframe = TRUE
  )
  
  # Draw border
  data3js <- r3js::lines3js(
    data3js,
    x = c(plot_xlim[1], plot_xlim[1], plot_xlim[2], plot_xlim[2], plot_xlim[1]),
    y = c(plot_ylim[1], plot_ylim[2], plot_ylim[2], plot_ylim[1], plot_ylim[1]),
    z = rep(plot_zlim[1], 5),
    lwd = 2,
    col = "grey70"
  )
  
  return(data3js)
}

render_widget <- function(data3js, remove_buttons = TRUE){
  if(remove_buttons) {
    data3js <- remove_buttons(data3js)
    
  }
  
  # Create html widget
  widget <- r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom
  )
  
  htmlwidgets::onRender(
    widget,
    jsCode = paste0("function(el, x, data){
    el.style.outline = 'solid 2px #eeeeee';
    }")
  )
}

# # Function to plot gmt and individual landscapes for specific serum group for fitted landscape
plot_sr_group_lndscp <- function(map, sr_group_logtiters, sr_group_cone_coords, sr_group_colbases, sr_group_cone_slopes, sr_group, 
                                 remove_buttons = TRUE, adjust_reactivity_bias = TRUE) {
  
  
  lndscp_xlim <- range(agCoords(map)[,1], na.rm = T)
  lndscp_ylim <- range(agCoords(map)[,2], na.rm = T)
  
  lndscp_ylim[1] <- lndscp_ylim[1] -2.5
  
  sr_group_logtiters_adj <- remove_reactivity_bias_logtiter(sr_group_logtiters)
  sr_group_colb_adj <- apply(sr_group_logtiters_adj, 2, function(x) max(x, na.rm = T))
  
  if(adjust_reactivity_bias) {
    
    sr_group_mean_logtiters <- rowMeans(sr_group_logtiters_adj, na.rm = T)
    sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
    
    
  } else {
    sr_group_mean_logtiters <- rowMeans(sr_group_logtiters, na.rm = T)
    sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
    
  }
  
  
  # Get fitted surface
  #xyz_matrix_gmt <- obtain_xyz_matrix(lndscp_xlim, lndscp_ylim, sr_group_cone_coords, sr_group_colbases, sr_group_cone_slopes)
  xyz_matrix_gmt <- obtain_xyz_matrix(lndscp_xlim, lndscp_ylim, sr_group_cone_coords, sr_group_colb_adj, sr_group_cone_slopes)
  grid_x_matrix <- xyz_matrix_gmt$grid_x_matrix
  grid_y_matrix <- xyz_matrix_gmt$grid_y_matrix
  grid_z_matrix <- xyz_matrix_gmt$grid_z_matrix
  
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  
 
  # Plot the base plot
  data3js <- plot_base_map(map, sr_group)
  
  # add gmt surface
  data3js <- add_landscape_plot(data3js, grid_x_matrix =  grid_x_matrix, grid_y_matrix = grid_y_matrix,
                                grid_z_matrix = grid_z_matrix, color = blend_sr_group_colors(sr_group, sr_group_colors), 
                                landscape_name = "GMT landscape")
  
  
  # Add individual landscapes
  for (i in seq_along(sr_group_colbases)) {
    
    ## Calculate the individual surface
    sr_cone_coord <- sr_group_cone_coords[i,]
    sr_colbase <- sr_group_colbases[i]
    sr_cone_slope <- sr_group_cone_slopes#[i]
    
    grid_dists <- as.matrix(dist(rbind(
      sr_cone_coord,
      cbind(
        as.vector(grid_x_matrix),
        as.vector(grid_y_matrix)
      )
    )))[1, -1]
    
    grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
    grid_z_matrix[] <- sr_colbase - grid_dists*sr_cone_slope
    
    # Add the individual surface
    data3js <- add_landscape_plot(data3js, grid_x_matrix =  grid_x_matrix, grid_y_matrix = grid_y_matrix,
                                  grid_z_matrix = grid_z_matrix, color = "grey70", 
                                  landscape_name = "GMT landscape", opacity_val = 0.2)
    
  }
  #
  
  ag_coords_titer <- agCoords(map)[names(sr_group_mean_logtiters),]
  # # Add the titers
  data3js <- lndscp3d_titers(
    data3js = data3js,
    object = list(
      coords = ag_coords_titer[!is.na(sr_group_mean_logtiters),],
      logtiters = sr_group_mean_logtiters[!is.na(sr_group_mean_logtiters)],
      indices = which(!is.na(sr_group_mean_logtiters)),
      acmap = map
    ),
    zlim = plot_zlim,
    options = list(
      cex.titer = 1,
      col.impulse = "grey60"
    )
  )
  
  # Add the titer 50 plane
  data3js <- add_titer_50_plane(data3js)
 
  #render widget
  render_widget(data3js, remove_buttons)
}



# Function to plot serum group gmt landscapes for fitted parameters
plot_sr_group_lndscp_gmt <- function(map, landscape_pars, sr_groups, remove_buttons =TRUE, adjust_reactivity_bias = TRUE) {
  
  lndscp_xlim <- range(agCoords(map)[,1], na.rm = T)
  lndscp_ylim <- range(agCoords(map)[,2], na.rm = T)
  
  sr_group <- sr_groups[1]
  
  if(adjust_reactivity_bias) {
    
  
    landscape_pars[[sr_group]]$log_titers <- remove_reactivity_bias_logtiter(landscape_pars[[sr_group]]$log_titers)
    landscape_pars[[sr_group]]$sr_colbases <- apply(landscape_pars[[sr_group]]$log_titers, 2, function(x) max(x, na.rm = T))
    
  }

  
  # Get fitted surface
  xyz_matrix_gmt <- obtain_xyz_matrix(lndscp_xlim, lndscp_ylim, landscape_pars[[sr_group]]$sr_cone_coords, landscape_pars[[sr_group]]$colbases, landscape_pars[[sr_group]]$slope)
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- xyz_matrix_gmt$grid_x_matrix
  grid_y_matrix <- xyz_matrix_gmt$grid_y_matrix
  grid_z_matrix <- xyz_matrix_gmt$grid_z_matrix
  
  # do base plot
  data3js <- plot_base_map(map, sr_group)
  
  # add first gmt surface
  data3js <- add_landscape_plot(data3js, grid_x_matrix =  grid_x_matrix, grid_y_matrix = grid_y_matrix,
                     grid_z_matrix = grid_z_matrix, color = blend_sr_group_colors(sr_group, sr_group_colors), 
                     landscape_name = sr_group)

  # Add other gmt_surfaces
  for (i in 2:length(sr_groups)) {
    
    sr_group <- sr_groups[[i]]
    
    if(adjust_reactivity_bias) {
      
      landscape_pars[[sr_group]]$log_titers <- remove_reactivity_bias_logtiter(landscape_pars[[sr_group]]$log_titers)
      landscape_pars[[sr_group]]$sr_colbases <- apply(landscape_pars[[sr_group]]$log_titers, 2, function(x) max(x, na.rm = T))
      
    }
    
    grid_z_matrix[] <- vapply(
      seq_len(length(grid_z_matrix)),
      \(n) {
        fit_lndscp_val(
          grid_x_matrix[n], grid_y_matrix[n],
          landscape_pars[[sr_group]]$sr_cone_coords,
          landscape_pars[[sr_group]]$colbases,
          landscape_pars[[sr_group]]$slope
        )
      }, numeric(1)
    )
    
    # Add the surface
    data3js <- add_landscape_plot(data3js, grid_x_matrix =  grid_x_matrix, grid_y_matrix = grid_y_matrix,
                                  grid_z_matrix = grid_z_matrix, color = blend_sr_group_colors(sr_group, sr_group_colors), 
                                  landscape_name = sr_group)
    
  }
  
  
  # Add the titer 50 plane
  # Add the titer 50 plane
  data3js <- add_titer_50_plane(data3js)
  
  #render widget
  render_widget(data3js, remove_buttons)
  
}


plot_subtract_landscapes <- function(map, sr_groups, landscape_pars, adjust_reactivity_bias = TRUE, remove_buttons = TRUE, parallel_lift_val = 0){
  
  lndscp_xlim <- range(agCoords(map)[,1], na.rm = T)
  lndscp_ylim <- range(agCoords(map)[,2], na.rm = T)
  
  sr_group <- sr_groups[1]
  
  if(adjust_reactivity_bias) {
    
    for(sg in sr_groups){
      landscape_pars[[sg]]$log_titers <- remove_reactivity_bias_logtiter(landscape_pars[[sg]]$log_titers)
      landscape_pars[[sg]]$sr_colbases <- apply(landscape_pars[[sg]]$log_titers, 2, function(x) max(x, na.rm = T))
    }
  }
  
  # do base plot
  data3js <- plot_base_map(map, sr_group)
  
  # Get fitted surface for first sg group
  xyz_matrix_gmt <- obtain_xyz_matrix(lndscp_xlim, lndscp_ylim, landscape_pars[[sr_groups[1]]]$sr_cone_coords, landscape_pars[[sr_groups[1]]]$colbases, landscape_pars[[sr_groups[1]]]$slope)
  grid_x_matrix <- xyz_matrix_gmt$grid_x_matrix
  grid_y_matrix <- xyz_matrix_gmt$grid_y_matrix
  grid_z_matrix_1 <- parallel_lift_val + xyz_matrix_gmt$grid_z_matrix
  
  grid_z_matrix_2 <- parallel_lift_val + obtain_xyz_matrix(lndscp_xlim, lndscp_ylim, landscape_pars[[sr_groups[2]]]$sr_cone_coords, landscape_pars[[sr_groups[2]]]$colbases, landscape_pars[[sr_groups[2]]]$slope)$grid_z_matrix
  
  # subtract the two 
  grid_z_matrix <- parallel_lift_val + grid_z_matrix_1 - grid_z_matrix_2
  
  
  
  data3js <- add_landscape_plot(data3js, grid_x_matrix, grid_y_matrix, grid_z_matrix_1, blend_sr_group_colors(sr_groups[1], sr_group_colors), sr_groups[1], opacity_val = 0.3)
  data3js <- add_landscape_plot(data3js, grid_x_matrix, grid_y_matrix, grid_z_matrix_2, color = blend_sr_group_colors(sr_groups[2], sr_group_colors), sr_groups[2], opacity_val = 0.3)
  data3js <- add_landscape_plot(data3js, grid_x_matrix, grid_y_matrix, grid_z_matrix, "#90dcff", "Difference", opacity_val = 0.6)
  
  sub_matrix_neg <- function(matrix, ref_mat, val = 0) {
    sub_mat <- matrix
    sub_mat[ref_mat >= val] <- NA
    
    return(sub_mat)
  }
  
  # try it with subsetting landscape to larger 0, smaller and color "#d15247"
  neg_z <- sub_matrix_neg(grid_z_matrix, grid_z_matrix)

  data3js <- add_landscape_plot(data3js, grid_x_matrix , grid_y_matrix, neg_z, "#d15247", "Difference", opacity_val = 0.6)
  
  
  # Add the titer 50 plane
 # data3js <- add_titer_50_plane(data3js)
  
 # Add plane at fold change of 1 
  data3js <- add_titer_50_plane(data3js, z_height = 1)
  
  #render widget
  render_widget(data3js, remove_buttons)
}


plot_single_landscape_panel <- function(landscape, label, label_size = 10, label_x_pos = 2, label_y_pos = 9,
                                        sr_group_label = "", sr_group_y_pos = 0, sr_group_size = 3, show_border = FALSE){
  
  
  to_save <- file.path("temp.html")
  png_save <- gsub(".html", ".png", to_save)
  saveWidget(landscape, to_save, selfcontained = FALSE)
  webshot(url=to_save,file = png_save)
  temp_plot <- readPNG(png_save)
  
  qplot(c(1:10),c(1:10), geom="blank") +
    annotation_custom(rasterGrob(temp_plot, height = unit(0.7, "npc")), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    annotate(geom="text", x=label_x_pos, y=label_y_pos, label=label,size= label_size, hjust = 0) + 
    annotate(geom="text", x=label_x_pos, y=sr_group_y_pos, label=sr_group_label,size= sr_group_size, hjust = 0) +
    theme_void() -> p

if(show_border) {
    p + theme(panel.border = element_rect(color = "grey50",
                                      fill = NA,
                                      size = 0.5))-> p
}
  
  if (file.exists(to_save)) {
    #Delete file if it exists
    file.remove(to_save)
  }
  if (file.exists(png_save)) {
    #Delete file if it exists
    file.remove(png_save)
  }
  
 return(p) 
}




screenshot_html_landscape_to_png <- function(landscape, save_name ){

  to_save <- file.path(paste0(save_name, ".html"))
  saveWidget(landscape, to_save, selfcontained = FALSE)

}
