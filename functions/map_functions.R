makeMap <- function(table, baseMap = NULL, dimensions = 2,
                    nOptimisations = 500, mcb = "none", dilution_stepsize = 1, options = list()) {
  # Take a titer table, make a map. Optionally apply plotspec and re-align
  # map to an already existing map.
  m <- acmap(titer_table = table)
  
  dilutionStepsize(m) <- dilution_stepsize
  
  m <- optimizeMap(
    map                     = m,
    number_of_dimensions    = dimensions,
    number_of_optimizations = nOptimisations,
    minimum_column_basis    = mcb,
    options = options
  )
  
  agNames(m) <- rownames(table)
  srNames(m) <- colnames(table)
  
  if(!(is.null(baseMap))) {
    m <- applyPlotspec(m, baseMap)
    m <- realignMap(m, baseMap)
  }
  
  ptDrawingOrder(m) <- rev(seq_len(numPoints(m)))
  
  return(m)
}



apply_style <- function(map){
  
  nag <- numAntigens(map)
  nsr <- numSera(map)
  N <- nsr + nag
  ptDrawingOrder(map) <- c(N:(nag+1), 1:nag)
  
  srOutlineWidth(map) <- 1
  srSize(map) <- 9
  agSize(map) <- 18
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.7)
  
  agSize(map)[agNames(map) == "B.1.1.7+E484K"] <- 16
  
  return (map)
  
}