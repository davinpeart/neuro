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

# extract decriptives from vector or matrix of observations
enframe_descriptives <- function(y, stat1 = mean, stat2 = var) {
  if(is.vector(y)) {
    opv <- data.frame(s1 = stat1(y), s2 = stat2(y))
    colnames(opv) <- c(as.character(substitute(stat1)), as.character(substitute(stat2)))
    return(opv)
  }
  if(is.matrix(y)) {
    I <- dim(y)[1]
    s1 <- vector("numeric", I)
    s2 <- vector("numeric", I)
    for(i in 1:I) {
      s1[i] <- stat1(y[i, ])
      s2[i] <- stat2(y[i, ])
    }
    op <- data.frame(s1, s2)
    colnames(op) <- c(as.character(substitute(stat1)), as.character(substitute(stat2)))
    return(op)
  }
}

# enframe proportion of integers
enframe_prop_integer <- function(y) {
  if(is.vector(y)) {
    x <- table(y)
    z <- as.vector(unname(x))
    return(data.frame(integer = as.integer(names(x)), freq = z, prop = z/length(y)))
  }
  if(is.array(y)) {
    I <- dim(y)[1]
    mx <- vector("numeric", I)
    mn <- vector("numeric", I)
    x <- vector("list", I)
    nm <- vector("list", I)
    for(i in 1:I) {
      x[[i]] <- table(y[i, ])/dim(y)[2]
      nm[[i]] <- as.integer(names(x[[i]]))
      mx[i] <- max(nm[[i]])
      mn[i] <- min(nm[[i]])
    }
    z <- seq.int(min(mn), max(mx))
    J <- length(z)
    w <- vector("list", J)
    names(w) <- as.character(z)
    for(j in 1:J) {
      for(i in 1:I) {
        w[[j]][i] <- x[[i]][names(w)[j]]
      }
    }
    op <- data.frame(integer = z)
    op$mean <- 0
    op$lower <- 0
    op$upper <- 0
    for(k in 1:J) {
      op[k, "mean"] <- mean(w[[as.character(op[k, "integer"])]], na.rm = T)
      op[k, "lower"] <- unname(quantile(w[[as.character(op[k, "integer"])]], .05, na.rm = T))
      op[k, "upper"] <- unname(quantile(w[[as.character(op[k, "integer"])]], .95, na.rm = T))
    }
    return(op)
  }
}

# calculate bayes factor using savage-dickey method
savage.dickey.bf <- function(posterior.samples, point.null = 0, plot = T,
                             prior.dist = "normal", prior.location = 0, prior.scale = 2) {
  if(prior.dist == "normal") {
    prior.fn <- function(q, mu = prior.location, sigma = prior.scale) {
      dnorm(x = q, mean = mu, sd = sigma)
    }
  }
  if(prior.dist == "cauchy") {
    prior.fn <- function(q, mu = prior.loaction, gamma = prior.scale) {
      dcauchy(x = q, location = mu, scale = gamma)
    }
  }
  prior.density <- prior.fn(q = point.null)
  posterior.spline <- polspline::logspline(posterior.samples)
  posterior.density <- polspline::dlogspline(point.null, posterior.spline)
  if(plot) {
    if(point.null < posterior.spline$range[1] & point.null < posterior.spline$range[2]) {
      interval <- seq(from = point.null, to = posterior.spline$range[2],
                      by = (posterior.spline$range[2] - point.null)/1000)
    }
    if(point.null > posterior.spline$range[1] & point.null < posterior.spline$range[2]) {
      interval <- seq(from = posterior.spline$range[1], to = posterior.spline$range[2],
                      by = (posterior.spline$range[2] - posterior.spline$range[1])/1000)
    }
    if(point.null > posterior.spline$range[1] & point.null > posterior.spline$range[2]) {
      interval <- seq(from = posterior.spline$range[1], to = point.null,
                      by = (point.null - posterior.spline$range[1])/1000)
    }
    graph <-
      ggplot2::ggplot(data = rbind(data.frame(Distribution = rep("Prior", times = length(interval)),
                                              Theta = interval,
                                              Density = prior.fn(q = interval)),
                                   data.frame(Distribution = rep("Posterior", times = length(interval)),
                                              Theta = interval,
                                              Density = polspline::dlogspline(interval, posterior.spline))),
                      mapping = ggplot2::aes(x = Theta, y = Density, ymax = Density,
                                             group = Distribution, fill = Distribution, colour = Distribution)) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = point.null), colour = "lightgrey") +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = 0), linewidth = 0, alpha = .5) +
      ggplot2::geom_point(data = rbind(data.frame(Distribution = "Prior", Theta = 0, Density = prior.density),
                                       data.frame(Distribution = "Posterior", Theta = 0,
                                                  Density = posterior.density)),
                          size = 4, alpha = .5) +
      ggplot2::scale_colour_manual(values = c(neuro::nature_palette("blue")[3], neuro::nature_palette("blue")[2])) +
      ggplot2::scale_fill_manual(values = c(neuro::nature_palette("blue")[3], neuro::nature_palette("blue")[1])) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     text = ggplot2::element_text(family = "Helvetica Neue"),
                     axis.line = ggplot2::element_line(colour = "black", linewidth = .2),
                     axis.ticks = ggplot2::element_line(colour = "black", linewidth = .2),
                     axis.title = ggplot2::element_text(size = 11),
                     axis.text = ggplot2::element_text(size = 10),
                     legend.title = ggplot2::element_text(size = 11),
                     legend.text = ggplot2::element_text(size = 10)) +
      ggplot2::geom_text(data = data.frame(Theta = quantile(interval, .8),
                                           Density = quantile(polspline::dlogspline(interval, posterior.spline), .85)),
                         inherit.aes = FALSE, mapping = ggplot2::aes(x = Theta, y = Density), size = 3.5,
                         label = paste("BF[10] : ", round(prior.density / posterior.density, digits = 2), sep = ""),
                         parse = TRUE, colour = "black") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
  }
  cat(paste0("Approximate BF (Savage-Dickey) in favour of null hypothesis Theta = ",
             point.null, " : ", round(posterior.density/prior.density, digits = 10), "\n"))
  if(plot) {
    return(graph)
  }
}
