# Split titers into pre and post titers
longTiterTableToPrePost <- function(
  longtitertable,
  vaccine_trial_year = NULL
){
  
  longtitertable %>%
    mutate(
      sr_parent   = collate(ag.applyFunction(sr_records, function(sr){ sr$parent_id })),
      sr_parent_records = ag.applyFunction(sr_records, function(sr){ parent.env(sr) }),
      sample_type = collate(ag.applyFunction(sr_records, function(sr){ sr$meta$sample_type })),
      logtiter    = acutils:::as.logtiter(titer, dilution_stepsize = 1),
      sample_type = recode(
        sample_type,
        "pre-vaccination"  = "pre",
        "post-vaccination" = "post"
      )
    ) %>%
    select(
      -sr,
      -sr_records,
      -srag_records
    ) %>%
    pivot_wider(
      names_from  = sample_type,
      values_from = c(titer, logtiter)
    )
  
}


titerlndscp_longTiterTablePrePost <- function(
  plotdata,
  title = waiver(),
  subtitle = waiver(),
  ylim = c(-1,9),
  boost.col = "red"
){
  
  plotdata %>%
    arrange(ag_plotid) %>%
    ggplot() +
    # Create colored rectangles to indicate clusters
    geom_tile(
      aes(
        x = ag_plotid,
        y = -1.5+(ylim[2]+1.5)/2,
        fill = ag_cluster
      ),
      alpha = 0.05,
      width  = 1,
      height = ylim[2]+1.5,
      show.legend = FALSE
    ) +
    geom_tile(
      aes(
        x = ag_plotid,
        y = -1.25,
        fill = ag_cluster
      ),
      width  = 1,
      height = 0.5,
      show.legend = FALSE
    ) +
    # Add line for pre titers
    geom_line(
      aes(
        x = ag_plotid,
        y = logtiter_pre,
        group = sr_parent
      ),
      linetype = "dashed"
    ) +
    # Add body polygon
    geom_ribbon(
      aes(
        x     = ag_plotid,
        ymax  = apply(cbind(logtiter_pre, logtiter_post), 1, min),
        group = sr_parent
      ),
      ymin = ylim[1],
      fill = "#333333"
    ) +
    # Add boost polygon
    geom_ribbon(
      aes(
        x = ag_plotid,
        ymin = apply(cbind(logtiter_pre, logtiter_post), 1, min),
        ymax = apply(cbind(logtiter_pre, logtiter_post), 1, max),
        group = sr_parent
      ),
      fill = boost.col
    ) +
    # Add decay polygon
    geom_ribbon(
      aes(
        x = ag_plotid,
        ymin = apply(cbind(logtiter_pre, logtiter_post), 1, max),
        ymax = logtiter_post,
        group = sr_parent
      ),
      fill = "#ffa64d"
    ) +
    # Add line for post titers
    geom_line(
      aes(
        x = ag_plotid,
        y = logtiter_post,
        group = sr_parent
      )
    ) +
    # Add points for pre titers
    geom_point(
      aes(
        x = ag_plotid,
        y = logtiter_pre
      ),
      size = 2,
      shape = "circle filled",
      fill = "#ffffff"
    ) +
    # Add points for post titers
    geom_point(
      aes(
        x = ag_plotid,
        y = logtiter_post
      ),
      size = 1
    ) +
    coord_cartesian(
      expand = FALSE,
      ylim = c(ylim[1]-0.5, ylim[2])
    ) +
    scale_fill_manual(
      values = h3ClusterColors
    ) +
    scale_x_discrete(
      labels = function(x){
        acdb.agShort(agdb.getIDs(x, plotdata$ag_records))
      }
    ) +
    scale_y_continuous(
      breaks = function(lims){ ceiling(lims[1]):floor(lims[2]) },
      labels = function(x){
        labs <- 2^x*10
        labs[x == -1] <- "<10"
        labs
      }
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "",
      y = ""
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size  = 10
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank()
    ) -> gp
  
  # Mark range of potential exposure
  if ("sr_parent_records" %in% colnames(plotdata)) {
    age      <- as.numeric(plotdata$sr_parent_records[[1]]$meta$age)
    vac_year <- as.numeric(plotdata$sr_parent_records[[1]]$meta$vaccine_trial_year)
    ag_years <- acdb.agYear(plotdata$ag_records)
    ag_lims  <- range(as.numeric(plotdata$ag_plotid[
      ag_years >= vac_year - age & ag_years <= vac_year
    ]))
    gp <- gp + annotate(
      "rect",
      xmin = ag_lims[1] - 0.5,
      xmax = ag_lims[2] + 0.5,
      ymin = ylim[1],
      ymax = ylim[2],
      fill = "blue",
      alpha = 0.1,
      color = NA
    )
  }
  
  # Return the plot
  gp
  
}
