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
savage.dickey.bf <- function(brmsfit, ..., point.null = 0, plot = T, gamma = 1000,
                             colour_scheme = "blue") {

  # function to parse terms from ...
  extract_expression_terms <- function(string) {
    terms <- strsplit(string, "[+-/*()]")[[1]]
    for(i in 1:length(terms)) {
      terms[i] <- paste0(str_extract_all(terms[i], "[^ ]")[[1]], collapse = "")
    }
    return(terms[terms != "" & !grepl(pattern = "exp|log", x = terms)])
  }

  # posterior
  code = deparse(substitute(...))
  terms <- unique(extract_expression_terms(string = code))
  if(length(terms) == 1) {
    posterior.code <- gsub(pattern = terms, x = code,
                           replacement = paste0("as.data.frame(brmsfit$fit)$", terms))
  } else {
    posterior.code <- gsub(pattern = terms[1], x = code,
                           replacement = paste0("as.data.frame(brmsfit$fit)$", terms[1]))
    for(i in 2:length(terms)) {
      posterior.code <- gsub(pattern = terms[i], x = posterior.code,
                             replacement = paste0("as.data.frame(brmsfit$fit)$", terms[i]))
    }
  }
  posterior.samples <- eval(parse(text = posterior.code))
  posterior.spline <- polspline::logspline(posterior.samples)
  posterior.density <- polspline::dlogspline(point.null, posterior.spline)

  # prior
  prior.terms <- terms
  for(i in 1:length(prior.terms)) {
    if(paste0("prior_", prior.terms[i]) %in% names(brmsfit$fit)) {
      prior.terms[i] <- paste0("prior_", prior.terms[i])
    } else {
      prior.terms[i] <- paste0("prior_", strsplit(prior.terms[i], "_")[[1]][1])
    }
  }
  if(length(prior.terms) == 1) {
    prior.code <- paste0("as.data.frame(brmsfit$fit)$", prior.terms)
  } else {
    prior.code <- gsub(pattern = terms[1], x = code,
                       replacement = paste0("as.data.frame(brmsfit$fit)$", prior.terms[1]))
    for(i in 2:length(prior.terms)) {
      prior.code <- gsub(pattern = terms[i], x = prior.code,
                         replacement = paste0("as.data.frame(brmsfit$fit)$", prior.terms[i]))
    }
  }
  prior.samples <- eval(parse(text = prior.code))
  prior.spline <- polspline::logspline(prior.samples)
  prior.fn <- function(q, spline = prior.spline) {
    polspline::dlogspline(q = q, fit = spline)
  }
  prior.density <- prior.fn(q = point.null)

  # plot densities
  if(plot) {
    colour_scheme <- strsplit(x = colour_scheme, split = "-")[[1]]
    if(length(colour_scheme) == 1) {
      colour_scheme[2] <- colour_scheme[1]
    }
    interval <- seq(from = min(point.null, posterior.spline$range[1], prior.spline$range[1]),
                    to = max(point.null, posterior.spline$range[2], prior.spline$range[2]),
                    by = (max(point.null, posterior.spline$range[2], prior.spline$range[2]) -
                            min(point.null, posterior.spline$range[1], prior.spline$range[1])) / gamma)
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
      ggplot2::geom_point(data = rbind(data.frame(Distribution = "Prior", Theta = point.null, Density = prior.density),
                                       data.frame(Distribution = "Posterior", Theta = point.null,
                                                  Density = posterior.density)),
                          size = 4, alpha = .5) +
      ggplot2::scale_colour_manual(values = c(neuro::nature_palette(colour_scheme[1])[3],
                                              neuro::nature_palette(colour_scheme[2])[2])) +
      ggplot2::scale_fill_manual(values = c(neuro::nature_palette(colour_scheme[1])[3],
                                            neuro::nature_palette(colour_scheme[2])[1])) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     text = ggplot2::element_text(family = "Helvetica Neue"),
                     axis.line = ggplot2::element_line(colour = "black", linewidth = .2),
                     axis.ticks = ggplot2::element_line(colour = "black", linewidth = .2),
                     axis.ticks.length = ggplot2::unit(1.75, "mm"),
                     axis.title = ggplot2::element_text(size = 11),
                     axis.text = ggplot2::element_text(size = 10),
                     legend.title = ggplot2::element_text(size = 11),
                     legend.text = ggplot2::element_text(size = 10)) +
      ggplot2::geom_text(data = data.frame(Theta = quantile(interval, .8),
                                           Density = quantile(polspline::dlogspline(interval, posterior.spline), .95)),
                         inherit.aes = FALSE, mapping = ggplot2::aes(x = Theta, y = Density), size = 3.5,
                         label = paste("BF[10] : ", signif(prior.density / posterior.density, digits = 3), sep = ""),
                         parse = TRUE, colour = "black") +
      ggplot2::scale_x_continuous(expand = c(.04, .04)) +
      ggplot2::scale_y_continuous(expand = c(.02, .02))
  }

  # return BF01
  cat(paste0("Approximate BF (Savage-Dickey) in favour of null hypothesis Theta = ",
             point.null, " : ", signif(posterior.density/prior.density, digits = 3), "\n"))

  # return plot
  if(plot) {
    return(graph)
  } else {
    return(posterior.density/prior.density)
  }
}

# visualize posterior of the expected value
ribbon_epreds <- function(brmsfit, x, y, method = mean, probs = c(.8, .99), plot = "preds",
                          wrap, points = F, hline = F, marginalize = F) {

  # obtain posterior of the expected value
  if(plot == "preds") {
    epreds <- as.data.frame(t(brms::posterior_epred(brmsfit)))
    colnames(epreds) <- gsub("V", "", colnames(epreds))
  }

  if(plot == "reps") {
    epreds <- as.data.frame(t(brms::posterior_predict(brmsfit)))
    colnames(epreds) <- gsub("V", "", colnames(epreds))
  }

  # integrate with information used from the observed data
  epreds <- cbind(brmsfit$data, epreds)

  # make tidy for summary and plotting
  epreds <- tidyr::pivot_longer(
    data = epreds,
    cols = colnames(epreds)[!(colnames(epreds) %in% colnames(brmsfit$data))],
    names_to = "Draw",
    values_to = "epred")

  # summarise plotting statistics
  if(!missing(wrap)) {
    epred_sum <- epreds %>%
      dplyr::group_by({{x}}, {{wrap}})
  }

  if(missing(wrap) & marginalize) {
    xname <- deparse(substitute(x))

    epred_sum <- epreds %>%
      dplyr::group_by({{x}}, Draw) %>%
      dplyr::reframe(xname = unique({{x}}),
                     epred = mean(epred)) %>%
      dplyr::ungroup() %>%
      group_by({{x}})
  }

  if(missing(wrap) & !marginalize) {
    epred_sum <- epreds %>%
      dplyr::group_by({{x}})
  }

  epred_sum <- epred_sum %>%
    dplyr::summarize(pred = method(epred),
                     back_lower = bayestestR::hdi(epred, ci = max(probs))[[2]],
                     back_upper = bayestestR::hdi(epred, ci = max(probs))[[3]],
                     front_lower = bayestestR::hdi(epred, ci = min(probs))[[2]],
                     front_upper = bayestestR::hdi(epred, ci = min(probs))[[3]])

  # plot
  plot <-
    ggplot2::ggplot(data = epred_sum, ggplot2::aes(x = {{x}})) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = back_lower, ymax = back_upper), fill = "#EEF5F9", colour = "#B9CCE0", linewidth = .4) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = front_lower, ymax = front_upper), fill = "#CCE1F2", colour = "#B8D0E8", linewidth = .4) +
    ggplot2::geom_line(ggplot2::aes(y = pred), colour = "#618BC3", linewidth = .6) +
    ggplot2::theme(axis.line = ggplot2::element_line(),
                   axis.title = ggplot2::element_text(face = "bold"),
                   axis.ticks.length = ggplot2::unit(1, "mm"),
                   panel.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(),
                   axis.text.y = ggplot2::element_text(),
                   strip.text = ggplot2::element_blank()) +
    ggplot2::labs(y = "Posterior") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = .5),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 8, face = "bold"),
                   text = ggplot2::element_text(colour = "black", family = "Helvetica"),
                   plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = .5),
                   axis.ticks = ggplot2::element_line(colour = "black", linewidth = .25),
                   axis.ticks.length = ggplot2::unit(.75, "mm"),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 8),
                   axis.title = ggplot2::element_text(size = 9, face = "bold"))

  # add wrap
  if(!missing(wrap)) {
    plot <- plot +
      ggplot2::facet_wrap(vars({{wrap}}), scales = "free_y", nrow = 2, axes = "all_y")
  }

  # add hline at 0
  if(hline) {
    plot <-
      plot +
      ggplot2::geom_hline(yintercept = 0, colour = "grey30", linetype = "dotted")
  }

  # add mean line
  if(points & marginalize) {
    if(!missing(wrap)) {
      plot <-
        plot +
        ggplot2::geom_point(data = brmsfit$data %>%
                              dplyr::group_by({{x}}, {{wrap}}) %>%
                              dplyr::summarize(y = mean({{y}})), aes(y = y),
                            colour = neuro::nature_palette("blue")[6], linewidth = .6)
    }
    if(missing(wrap)) {
      plot <-
        plot +
        ggplot2::geom_point(data = brmsfit$data %>%
                              dplyr::group_by({{x}}, {{wrap}}) %>%
                              dplyr::summarize(y = mean({{y}})), aes(y = y),
                            colour = neuro::nature_palette("blue")[6], stroke = .6)
    }
  }

  # add points
  if(points & !marginalize) {
    plot <- plot +
      ggplot2::geom_point(data = brmsfit$data, mapping = ggplot2::aes(x = {{x}}, y = {{y}}),
                          shape = 16, colour = neuro::nature_palette("blue")[6], fill = NA, alpha = 1, size = 1.3)
  }

  return(list(plot, epreds, epred_sum))
}

# plot posterior predictive distribution for negative binomial variables
plot_posterior_pred <- function(data, col, fac, draws,
                                pt.size = 3, pt.stroke = .2, pt.shape = 21,
                                lwd.bar = .2, lwd.den = .8, adj = 1,
                                xlab = NULL, xlim = NULL) {

  # extract distributional parameters
  mu1 <- draws[draws["Group"] == "E", ][["mu"]]
  mu2 <- draws[draws["Group"] == "V", ][["mu"]]
  pi1 <- draws[draws["Group"] == "E", ][["shape"]]
  pi2 <- draws[draws["Group"] == "V", ][["shape"]]

  # draw from posterior predictive distribution
  preds <- data.frame(x = c(rnbinom(10000, mu = exp(mu1), size = exp(pi1)),
                            rnbinom(10000, mu = exp(mu2), size = exp(pi2))))

  # graph posterior predictions
  p <-
    ggplot(
      enframe_prop_integer(preds$x),
      aes(x = integer, y = prop)) +
    geom_bar(stat = "identity", fill = "#98A9C0", colour = "#55587D", linewidth = lwd.bar) +
    geom_density(inherit.aes = F, data = preds, aes(x = x), adjust = adj,
                 colour = "#55587D", linewidth = lwd.den, bounds = c(0, Inf)) +
    theme(panel.background = element_rect(fill = "#EBEEF2"),
          axis.line = element_line(colour = "grey40"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = xlim) +
    guides(x = "prism_offset", y = "prism_offset")

  # plot distribution of sample data
  s <-
    ggplot(
      data, aes(x = {{col}}, y = {{fac}})) +
    geom_point(shape = pt.shape, size = pt.size, fill = "#98A9C0", colour = "#55587D",
               stroke = pt.stroke, position = position_jitter(height= 0.01)) +
    theme(text = element_text(colour = "grey20"),
          panel.background = element_rect(fill = "#EBEEF2"),
          axis.ticks = element_line(linewidth = .3),
          axis.line = element_line(colour = "grey40", linewidth = .3),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face = "bold")) +
    labs(x = xlab) +
    scale_y_discrete() +
    scale_x_continuous(limits = xlim) +
    guides(x = "prism_offset", y = "prism_offset")

  p / s + plot_layout(axes = "collect", heights = c(1, .2))
}
