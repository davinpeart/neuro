# graphical predictive check for bayesian regression models
#' Plot density of two descriptive statistics of draws against observed values
#'
#' @param yrep Matrix of posterior predictive distributions with the same number of columns as rows in y.
#' @param y Vector of observed values.
#' @param ynu1 (Optional) Vector of previously observed values.
#' @param ynu2 (Optional) Vector of previously observed values.
#' @param nbin Number of bins of density cutoffs.
#' @param stat1 Function to be performed on y and each draw of yrep.
#' @param stat2 Function to be performed on y and each draw of yrep.
#' @param fill Colour of the gradient mapped to the density of descriptives. See ?nature_palette or ?springer.
#' @param y_col Colour of plotted values of y. See ?nature_palette or ?springer.
#' @param adjust Bandwidth of the density boundaries passed to geom_density2d_filled.
#' @param linetype Linetype of the plot of y passed to geom_segment.
#' @param linewidth Linewidth of the plot of y passed to geom_segment.
#' @param stroke Linewidth of the shape of y passed to geom_point.
#'
#' @returns A ggplot.
#' @export
#'
#' @examples
#' data <- data.frame(y = rnorm(n = 100))
#' draws <- matrix(data = rnorm(n = 100 * 2000), nrow = 100)
#' pp_stat_dens(yrep = draws, y = data$y)
pp_stat_dens <- function(yrep, y, nbin = 4, stat1 = mean, stat2 = sd,
                         fill = "yellow", y_col = "red",
                         linetype = "dotted", adjust = 1, linewidth = .5, stroke = .5,
                         ynu1 = NULL, ynu2 = NULL) {

  # check fill_cols is length-2 vector of colours
  is.colour <- function(x) {
    !("try-error" %in% class(try(col2rgb(x), silent = TRUE)))
  }

  # check y_col is a single colour in nature palette
  ops <- c("background", "grey", "olive", "green", "teal",
           "blue", "purple", "red", "orange", "yellow", "peach", "pink",
           "darkpurple", "navy")
  if(!(y_col %in% ops) | !(fill %in% ops) | length(y_col) != 1 | length(fill) != 1)
    stop(paste("y_col and fill must each be one of:",
               paste0(ops, collapse = ", ")
    ))

  # set palette
  ramp_palette <- grDevices::colorRampPalette(c("white", nature_palette(fill, 3)))

  # calculate statistics
  yrep_enf <- enframe_descriptives(yrep, stat1 = stat1, stat2 = stat2)
  y_enf <- enframe_descriptives(y, stat1 = stat1, stat2 = stat2)

  # plot values
  graph <-
    ggplot2::ggplot(yrep_enf, ggplot2::aes(x = stat1, y = stat2)) +
    ggplot2::geom_density2d_filled(contour_var = "ndensity", bins = nbin,
                                   linewidth = 0, alpha = 1, adjust = adjust) +
    ggplot2::scale_fill_manual(values = c(ramp_palette(nbin)),
                               name = "Replication") +
    ggplot2::geom_segment(
      data = y_enf, ggplot2::aes(x = stat1, xend = stat1, y = -Inf,  yend = stat2),
      colour = nature_palette(y_col, 4), linewidth = linewidth, linetype = linetype) +
    ggplot2::geom_segment(
      data = y_enf, ggplot2::aes(x = -Inf, xend = stat1, y = stat2,  yend = stat2),
      colour = nature_palette(y_col, 4), linewidth = linewidth, linetype = linetype) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(data = y_enf, aes(fill = ""), colour = nature_palette(y_col, 5),
                        size = 4, alpha = 1, shape = 21, stroke = stroke) +
    ggplot2::scale_fill_manual(values = nature_palette(y_col, 3), name = "Sample") +
    ggplot2::labs(x = as.character(substitute(stat1)), y = as.character(substitute(stat2))) +
    theme_nature()
  return(graph)
}

# extract decriptives from vector or matrix of observations
#' Calculate desriptives of a vector of observations or matrix of posterior predicive distributions
#'
#' @param y Vector or matrix of values to be described.
#' @param stat1 Function to calculate descriptive 1.
#' @param stat2 Function to calculate descriptive 2.
#'
#' @returns Dataframe containing stat1 and stat2 of y.
#' @export
#'
#' @examples
#' data <- data.frame(y = rnorm(n = 100))
#' enframe_descriptives(y = data$y)
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
#' Obtain the proportion of each integer in a set of values
#'
#' @param y Vector of observed values or array of posterior predictive distributions.
#'
#' @returns Dataframe containing proportion of each integer contained in y.
#' @export
#'
#' @examples
#' data <- data.frame(y = rpois(n = 1000, lambda = 3))
#' enframe_prop_integer(y = data$y)
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

# function to check if strings are function names
exists_vec <- function(x) {
  v <- c()
  for(i in 1:length(x)) {
    v[i] <- exists(x[i], envir = .GlobalEnv)
  }
  v
}

# function to parse terms from ...
extract_expression_terms <- function(string) {
  terms <- strsplit(string, "[+-/*()]")[[1]]
  for(i in 1:length(terms)) {
    terms[i] <- paste0(stringr::str_extract_all(terms[i], "[^ ]")[[1]], collapse = "")
  }
  terms <- terms[terms != "" & !grepl(pattern = "^", x = terms, fixed = TRUE)]
  terms <- terms[!exists_vec(terms)]
  terms <- terms[is.na(suppressWarnings(as.double(terms)))]
  return(terms)
}

# calculate bayes factor using savage-dickey method
#' Bayesian hypothesis testing using Savage-Dickey density ratio as Bayes Factor.
#' Prior samples and saving of parameters when fitting required. See example code
#'
#' @param brmsfit An object of class brmsfit.
#' @param ... An expression indicating the posterior to be tested against the null.
#' @param point.null The value of the point null hypothesis.
#' @param plot Logical indicating whether to return a plot of the hypothesis. If F, returns a Bayes Factor in favour of the null.
#' @param gamma Resolution of sampling from the prior and posterior distributions for plotting.
#' @param colour_scheme String indicating a single value to be supplied to nature_pallete(). Optionally
#' a string containing two values separated by "-". See ?nature_pallete for acceptable inputs.
#'
#' @returns See argument plot.
#' @export
#'
#' @examples
#' data <- data.frame(a = rep(c("a", "b"), each = 30), b = c(rnorm(n = 30), rnorm(n = 30, mean = 3)))
#' fit <- brms::brm(b ~ 1 + a, data = data, sample_prior = "yes", save_pars = save_pars(all = TRUE))
#' savage.dickey(brmsfit = fit, b_a)
savage.dickey <- function(brmsfit, ..., point.null = 0, plot = T, gamma = 1000,
                          colour_scheme = "blue", lab) {

  # posterior
  code = paste0(deparse(substitute(...)), collapse = "")
  terms <- unique(extract_expression_terms(string = code))
  if(length(terms) == 1) {
    posterior.code <- gsub(pattern = terms, x = code,
                           replacement = paste0("as.data.frame(brmsfit$fit)$", paste0("`", terms, "`")))
  } else {
    posterior.code <- gsub(pattern = terms[1], x = code,
                           replacement = paste0("as.data.frame(brmsfit$fit)$", paste0("`", terms[1], "`")))
    for(i in 2:length(terms)) {
      posterior.code <- gsub(pattern = terms[i], x = posterior.code,
                             replacement = paste0("as.data.frame(brmsfit$fit)$", paste0("`", terms[i], "`")))
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
    prior.code <- gsub(pattern = terms[1], x = code,
                       replacement = paste0("as.data.frame(brmsfit$fit)$", paste0("`", prior.terms[1], "`")))
  } else {
    prior.code <- gsub(pattern = terms[1], x = code,
                       replacement = paste0("as.data.frame(brmsfit$fit)$", paste0("`", prior.terms[1], "`")))
    for(i in 2:length(prior.terms)) {
      prior.code <- gsub(pattern = terms[i], x = prior.code,
                         replacement = paste0("as.data.frame(brmsfit$fit)$", paste0("`", prior.terms[i], "`")))
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

  if(!missing(lab)) {
    graph <-
      graph + ggplot2::labs(x = lab)
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
#' Visualization of the posterior of the expected value by group.
#'
#' @param brmsfit Fit of class brmsfit.
#' @param x Expression of variable to be mapped to the x axis.
#' @param y Expression of variable to be mapped to the y axis (the dependent variable).
#' @param wrap Expression of variable to be passed to facet_wrap.
#' @param grid Two-sided formula to be passed to facet_grid. The LHS should contain the variable to be
#' displayed in the rows of the grid and the RHS should contain the variable to be mapped to the columns.
#' @param method A function indicating the method of calculating the point estimate of the posterior.
#' @param probs A numeric vector of length 2 indicating the length of the CIs.
#' @param plot A string indicating whether to plot "preds" or "reps".
#' @param points Logical indicating whether to plot points.
#' @param hline Logical indicating whether to draw a line at x = 0.
#' @param marginalize Logical indicating whether samples should be marginalized within draws. If there are multiple
#' draws within each component of the specified design (e.g., repeated measures within subjects), then marginalize should be TRUE.
#' @param nrow Number of rows passed to facet_wrap.
#'
#' @returns A list of length 3 containing the specified plot, a dataframe containing samples bound to the original observations,
#' and a summarized dataframe containing point estimates and CIs for the specified design.
#' @export
#'
#' @examples
#' data("iris")
#' sepal <- brms::brm(Sepal.Length ~ 1 + Sepal.Width + (1|Species), data = iris,
#' save_pars = brms::save_pars(all = TRUE), sample_prior = TRUE)
#' ribbon_epreds(brmsfit = sepal, x = Sepal.Width, y = Sepal.Length, wrap = Species, nrow = 1)
#'
#'
ribbon_epreds <- function(brmsfit, x, y, wrap, grid, method = median, probs = c(.5, .89), plot = "preds",
                          points = T, hline = F, marginalize = F, nrow = NULL) {

  # obtain variable names as strings
    yname <- deparse(substitute(y))
    xname <- deparse(substitute(x))

    if(!missing(wrap)) {
      wrapvar <- deparse(substitute(wrap))
    }

    if(!missing(grid)) {
      gridvars <- strsplit(extract_expression_terms(deparse(substitute(grid))), split = "~")[[1]]
    }

  # obtain posterior of the expected value
  if(plot == "preds") {
    epreds <- as.data.frame(t(brms::posterior_epred(brmsfit, resp = yname)))
    colnames(epreds) <- gsub("V", "", colnames(epreds))
  }

  if(plot == "reps") {
    epreds <- as.data.frame(t(brms::posterior_predict(brmsfit, resp = yname)))
    colnames(epreds) <- gsub("V", "", colnames(epreds))
  }

  # bind to information used from the observed data
  epreds <- cbind(brmsfit$data, epreds)

  # make tidy for summary and plotting
  epreds <- tidyr::pivot_longer(
    data = epreds,
    cols = colnames(epreds)[!(colnames(epreds) %in% colnames(brmsfit$data))],
    names_to = "Draw",
    values_to = "epred")

  # summarize plotting statistics
  if(!missing(wrap) & !marginalize) {
    epred_sum <- epreds %>%
      dplyr::group_by({{x}}, {{wrap}})
  }

  if(!missing(wrap) & marginalize) {
    epred_sum <- epreds %>%
      dplyr::group_by({{x}}, {{wrap}}, Draw) %>%
      dplyr::reframe(xname = unique({{x}}),
                     wrapnamename = unique({{wrap}}),
                     epred = mean(epred)) %>%
      dplyr::group_by({{x}}, {{wrap}})
  }

  if(!missing(grid) & marginalize) {
    gridvar1 <- gridvars[1]
    gridvar2 <- gridvars[2]
    epred_sum <- epreds %>%
      dplyr::group_by({{x}}, !!sym(gridvars[1]), !!sym(gridvars[2]), Draw) %>%
      dplyr::reframe(xname = unique({{x}}),
                     gridvar1 = unique(!!sym(gridvars[1])),
                     gridvar2 = unique(!!sym(gridvars[2])),
                     epred = mean(epred)) %>%
      dplyr::group_by({{x}}, !!sym(gridvars[1]), !!sym(gridvars[2]))
  }

  if(!missing(grid) & !marginalize) {
    epred_sum <- epreds %>%
      dplyr::group_by({{x}}, !!sym(gridvars[1]), !!sym(gridvars[2]))
  }

  if(missing(wrap) & missing(grid) & marginalize) {
    xname <- deparse(substitute(x))

    epred_sum <- epreds %>%
      dplyr::group_by({{x}}, Draw) %>%
      dplyr::reframe(xname = unique({{x}}),
                     epred = mean(epred)) %>%
      dplyr::ungroup() %>%
      group_by({{x}})
  }

  if(missing(wrap) & missing(grid) & !marginalize) {
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
      ggplot2::facet_wrap(ggplot2::vars({{wrap}}), scales = "free_y", nrow = nrow, axes = "all_y")
  }

  # add grid
  if(!missing(grid)) {
    plot <- plot +
      ggplot2::facet_grid(grid)
  }

  # add hline at 0
  if(hline) {
    plot <-
      plot +
      ggplot2::geom_hline(yintercept = 0, colour = "grey30", linetype = "dotted")
  }

  # add mean points
  if(points & marginalize) {
    if(!missing(wrap)) {
      plot <-
        plot +
        ggplot2::geom_point(data = brmsfit$data %>%
                              dplyr::group_by({{x}}, {{wrap}}) %>%
                              dplyr::summarize(y = mean({{y}})), aes(y = y),
                            colour = neuro::nature_palette("blue")[6], stroke = .6)
    }
    if(!missing(grid)) {
      plot <-
        plot +
        ggplot2::geom_point(data = brmsfit$data %>%
                              dplyr::group_by({{x}}, !!sym(gridvars[1]), !!sym(gridvars[2])) %>%
                              dplyr::summarize(y = mean({{y}})), aes(y = y),
                            colour = neuro::nature_palette("blue")[6], stroke = .6)
    }
    if(missing(wrap) & missing(grid)) {
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

# function for vectorized brms priors
priors <- function(..., data, prior, exception = "Intercept", class = "b", dpar = "", nlpar = "",
                   lb = NA, ub = NA, resp = "") {
  pars <- grep(pattern = exception, invert = T, value = T, colnames(model.matrix(..., data)))
  priors <- brms::set_prior(prior = prior, class = class, coef = pars[1], nlpar = nlpar, dpar = dpar,
                            lb = lb, ub = ub, resp = resp)
  for(i in 2:length(pars)) {
    priors <- c(priors,
                brms::set_prior(prior = prior, class = class, coef = pars[i], nlpar = nlpar, dpar = dpar,
                                lb = lb, ub = ub, resp = resp))
  }
  return(priors)
}

# logit function
#' Logit
#'
#' @param p A probability between 0 and 1 to be tranformed to unbounded parameter space.
#'
#' @returns A numeric value.
#' @export
#'
#' @examples
#' logit(p = .75)
logit <- function(p) {
  log(p / (1 - p))
}

# double exponential curve
#' Exponential change
#'
#' @param x Numeric value(s) of time.
#' @param intercept Numeric intercept at t = 0.
#' @param asymptote Numeric level of asymptote of change.
#' @param alpha Numeric rate of change.
#' @param beta Numeric rate of change.
#'
#' @returns Numeric value of change at t = x.
#' @export
#'
#' @examples
#' double exponential(x = 3, intercept = 1, asymptote = 5, alpha = 5, beta = 7)
double_exponential <- function(x, intercept, asymptote, alpha, beta) {
  intercept - asymptote*(exp(-x / alpha) - exp(-x / beta))
}

# exponential curve
#' Exponential change from intercept toward asymptote.
#'
#' @param x Numeric value of time.
#' @param intercept Numeric value of intercept at t = 0.
#' @param asymptote Numeric asymptote of change.
#' @param beta Numeric rate of change.
#'
#' @returns Numeric value fo change at t = x.
#' @export
#'
#' @examples
#' asymptotic_exponential(x = .3, intercept = 2, asymptote = 1, beta = 1)
asymptotic_exponential <- function(x, intercept, asymptote, beta) {
  asymptote - (asymptote - intercept)*exp(-x / beta)
}

# generalized logistic
#' Logistic change
#'
#' @param A The left asymptote.
#' @param K The right asymptote.
#' @param B Maximum rate of change.
#' @param Q Translation along the x axis (useful in specification of y at t = 0).
#' @param M Alternative method of translation along the x axis; starting time.
#' @param t Time.
#' @param v Asymptotic skew of change.
#'
#' @returns Numeric value of y at t = t.
#' @export
#'
#' @examples
#' generalized_logistic(A = 1, K = 10, B = 1, Q = 10, M = 0, t = seq(0, 10, .1), v = .5)
generalized_logistic <- function(A, K, B, Q, M, t, v) {
  A + ((K - A) / ((1 + Q*exp(-1*B*(t-M)))^(1/v)))
}

#' Formulae for generalized non-linear regression
#'
#' @title
#' Generate a formula for use in a negative binomial, poisson, or gamma regression function
#'
#' @description
#' These formulae are parameterized such that the mean is constrained to non-negative reals.
#' Each non-linear parameter will be modeled in logarithmic space.
#'
#' @param dv Expression containing the name of the variable on the LHS of the formula.
#' @param t Expression containing the name of the variable on the RHS of the formula.
#' @param curve String indicating the type of curve to be generated. One of "exp",
#' "double exp", or "logistic".
#' @param t0 (Optional) The starting time.
#'
#' @returns An object of class formula. See functions asymptotic_exponential(),
#' double_exponential(), and generalized_logistic for more details on the non-linear parameters.
#'
#' When curve = "exp", returns the following formula:
#' y ~ exp(L) - (exp(L) - exp(M)) * exp(-(t)/exp(o))
#' L = intercept
#' M = asymptote
#' o = rate of change
#'
#' When curve = "double exp", returns the following formula:
#' y ~ exp(C) - exp(D) * (exp(-(t)/exp(r)) - exp(-(t)/(exp(r) + exp(s))))
#' C = intercept
#' D = asymptote
#' r = rate of change 1
#' s = rate of change 2 (above r, to constrain to growth or shrinkage)
#'
#' When curve = "logistic", returns the following formula:
#' y ~exp(A) + (exp(K)/((1 + exp(Q) * exp(-exp(B) * (t)))^(1/exp(v))))
#' A = left asymptote
#' K = right asymptote (above r, to constrain to growth)
#' Q = horizontal shift
#' B = maximum rate of change
#' v = asymptotic asymmetry
#' This curve is symmetrical about the origin when A is negative, K = -2*A, and all other parameters are set to 1.
#'
#' @export
#'
#' @examples
#' x <- runif(n = 200, min = 0, max = 5)
#' data <- data.frame(x = x, y = rpois(n = 200, lambda = 6 - (6 - 1) * exp(-x / 1) ))
#' fit <-
#' brms::brm(brms::bf(
#'   nlnb_formula(dv = y, t = x, curve = "exp", t0 = 0),
#'   L + M + o ~ 1,
#'   nl = TRUE),
#'   data = data,
#'   family = poisson(link = "identity"),
#'   control = list(adapt_delta = .9, max_treedepth = 12))
nlnb_formula <- function(dv, t, curve, t0) {

  dv_name <- deparse(substitute(dv))
  t_name <- deparse(substitute(t))
  if(!missing(t0)) {
    t0_text <- paste0(" - ", t0)
  } else {
    t0_text <- ""
  }

  if(curve == "exp") {
    x <-
      paste0(dv_name, " ~ exp(L) - (exp(L) - exp(M))*exp(-(", t_name,
             t0_text, ") / exp(o))")
  }

  if(curve == "double exp") {
    x <-
      paste0(dv_name, " ~ exp(C) - exp(D)*(exp(-(", t_name, t0_text,
             ") / exp(r)) - exp(-(", t_name, t0_text, ") / (exp(r) + exp(s))))")
  }

  if(curve == "logistic") {
    x <-
      paste0(dv_name, " ~ exp(A) + (exp(K) / ((1 + exp(Q)*exp(-exp(B)*(",
             t_name, t0_text, ")))^(1/exp(v))))")
  }
  output <- as.formula(x)
  return(output)
}

