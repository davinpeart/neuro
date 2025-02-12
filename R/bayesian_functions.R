# graphical predictive check for bayesian regression models
pp_stat_dens <- function(yrep, y, ynu1 = NULL, ynu2 = NULL, nbin = 4, stat1 = mean, stat2 = sd,
                         fill_limits = c("#FEF8F9", "#FA819E"), y_col = "teal") {
  # check fill_cols is length-2 vector of colours
  is.colour <- function(x) {
    !("try-error" %in% class(try(col2rgb(x), silent = TRUE)))
  }
  if(!is.colour(fill_limits) | length(fill_limits) != 2)
    stop("fill_cols must be a length-2 vector of colours")
  
  # check y_col is a single colour in nature palette
  ops <- c("background", "grey", "olive", "green", "teal", 
           "blue", "purple", "red", "orange", "yellow")
  if(!(y_col %in% ops) | length(y_col) != 1)
    stop(paste("y_col must be one of", 
               paste0(ops, collapse = ", ")
    ))
  
  # set palette
  ramp_palette <- grDevices::colorRampPalette(fill_limits)
  
  # calculate statistics
  yrep_enf <- enframe_descriptives(yrep, stat1 = stat1, stat2 = stat2)
  y_enf <- enframe_descriptives(y, stat1 = stat1, stat2 = stat2)
  
  # plot values
  graph <- 
    ggplot2::ggplot(yrep_enf, ggplot2::aes(x = stat1, y = stat2)) +
    ggplot2::geom_density2d_filled(contour_var = "ndensity", bins = nbin,
                                   linewidth = 0, alpha = 1, adjust = 1) +
    ggplot2::scale_fill_manual(values = c(ramp_palette(nbin)), 
                               name = "Prediction") +
    ggplot2::geom_segment(
      data = y_enf, ggplot2::aes(x = stat1, xend = stat1, y = -Inf,  yend = stat2),
      colour = nature_palette(y_col, 4), linewidth = .7, linetype = "dotted") +
    ggplot2::geom_segment(
      data = y_enf, ggplot2::aes(x = -Inf, xend = stat1, y = stat2,  yend = stat2),
      colour = nature_palette(y_col, 4), linewidth = .7, linetype = "dotted") +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(data = y_enf, aes(fill = ""), colour = nature_palette(y_col, 5),
                        size = 4, alpha = 1, shape = 21) +
    ggplot2::scale_fill_manual(values = nature_palette(y_col, 3), name = "Sample") +
    ggplot2::labs(x = as.character(substitute(stat1)), y = as.character(substitute(stat2))) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_line(colour = "grey50", linewidth = .7),
      axis.ticks.y = ggplot2::element_line(colour = "grey50", linewidth = .7),
      axis.ticks.x = ggplot2::element_line(colour = "grey50", linewidth = .5),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.ticks.length.x = ggplot2::unit(-.1, "cm"),
      text = ggplot2::element_text(family = "Arial", colour = "grey20"),
      axis.title = ggplot2::element_text(face = "bold", size = rel(.9)),
      legend.title = ggplot2::element_text(face = "bold", size = rel(.9))) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  return(graph)
}