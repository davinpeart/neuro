# matrix of p values for correlation coeficients
#' Matrix of p-values for hypothesis testing of correlation coefficients
#'
#' @param x A wide dataframe with columns containing numeric variables to be tested
#' @param method Method of calculating correlation coefficient supplied to cor.test.
#' One of: "pearson", "kendall", "spearman".
#'
#' @returns A symmetric numeric matrix containing p-values for the correlation coefficients calculated for
#' each pair of variables in x.
#' @export
#'
#' @examples
#' x <- data.frame(a = rnorm(n = 10), b = rnorm(n = 10), c = rnorm(n = 10))
#' cor_test_mat(x = x, method = "spearman")
cor_test_mat <- function(x, method = "pearson") {
  N <- ncol(x)
  vars <- colnames(x)
  y <- matrix(nrow = N, ncol = N)
  colnames(y) <- vars
  rownames(y) <- vars
  for(n in 1:N) {
    for(i in 1:N) {
      y[vars[n], vars[i]] <-
        cor.test(x = x[[vars[n]]], y = x[[vars[i]]], conf.level = 0.95,
                 alternative = "two.sided", method = method,
                 exact = FALSE)$p.value
    }
  }
  return(y)
}

# plot correlation matrix given r and p matrices
#' Plot of correlation matrix showing coefficients and significance
#'
#' @param rmat Numeric matrix of correlation coefficients.
#' @param pmat Numeric matrix of identical structure to rmat containing correlation coefficient p-values.
#' @param sepvars Character vector of dimnames of rmat and pmat to be moved to the edge of the plot.
#' @param show_legend Logical indicating whether to plot a legend for colour map of correlation coefficient.
#' @param lwd Linewidth of squares containing correlation coefficients.
#' @param mid_fill String indicating colour for a correlation coefficient of 0.
#' @param axis_rel Relative size of axis text.
#' @param legend_rel Relative size of legend text.
#' @param guide_rel Relative size of legend title.
#' @param fill_cols Character vector of length 2 containing colours mapped to correlation coefficients of -1 and 1.
#' @param line_col Colour of the lines separating the squares containing each correlation coefficient.
#' @param rsize Size of the text indicating the values of the correlation coefficients.
#' @param psize Size of the text indicating significance.
#'
#' @returns A ggplot depicting the correlation matrix
#' @export
#'
#' @examples
#' x <- data.frame(a = rnorm(n = 10), b = rnorm(n = 10), c = rnorm(n = 10))
#' rmat <- cor(x = x)
#' pmat <- cor_test_mat(x = x, method = "spearman")
#' plot_cor_pmat(rmat = rmat, pmat = pmat)
plot_cor_mat <- function(rmat, pmat, sepvars = NULL, show_legend = T, lwd = 0, mid_fill = "white",
                         axis_rel = .7, legend_rel = .7, guide_rel = .9,
                         fill_cols = c("#657EA1", "#FA819E"), line_col = "white", rsize = 3, psize = 2) {
  # check if matrices were constructed using the same variables
  if(!identical(colnames(rmat), colnames(pmat)) |
     !identical(rownames(rmat), rownames(pmat)))
    stop("rmat and pmat must be correlation matrices with identical dimensions")

  # check show_legend is logical
  if(!is.logical(show_legend))
    stop("show_legend must be true or false")

  # check fill_cols is length-2 vector of colours
  is.colour <- function(x) {
    !("try-error" %in% class(try(col2rgb(x), silent = TRUE)))
  }
  if(!is.colour(fill_cols) | length(fill_cols) != 2)
    stop("fill_cols must be a length-2 vector of colours")

  # check line_col is a single colour
  if(!is.colour(line_col) | length(line_col) != 1)
    stop("line_col must be a single colour")

  # enframe matrices in tidy dataframe
  d <- data.frame(var1 = rep(rownames(rmat), times = ncol(rmat)),
                  var2 = rep(colnames(rmat), each = nrow(rmat)),
                  r = NA, p = NA, sig = NA, sep = "no")
  for(n in 1:nrow(d)) {
    d$r[n] <- rmat[d$var1[n], d$var2[n]]
    d$p[n] <- pmat[d$var1[n], d$var2[n]]
    if(d$p[n] < .05)
      d$sig[n] <- "*"
    else
      d$sig[n] <- ""
  }

  # relevel factors containing variables by matching
  if(is.character(sepvars)) {
    d$var1 <- factor(d$var1)
    d$var2 <- factor(d$var2)
    fac_relev_byvar <- function(f, v) {
      c(levels(f)[!(levels(f) %in% v)], levels(f)[levels(f) %in% v])
    }
    limits <- fac_relev_byvar(d$var1, sepvars)
  } else
    limits <- NULL

  # make plot
  ggplot2::ggplot(d, ggplot2::aes(x = var1, y = var2, fill = r)) +
    ggplot2::geom_tile(colour = line_col, linewidth = lwd) +
    ggplot2::scale_fill_gradient2(
      name = "Pearson's r",
      low = fill_cols[1],
      mid = mid_fill,
      high = fill_cols[2],
      midpoint = 0,
      breaks = seq(-1, 1, 0.5),
      limits = c(-1, 1),
      guide = "colourbar",
      aesthetics = "fill"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::geom_text(ggplot2::aes(label = round(r, digits = 2)), size = psize,
                       colour = "grey20") +
    ggplot2::geom_text(ggplot2::aes(label = sig), colour = "grey20",
                       nudge_y = 0.25, nudge_x = .35, size = rsize) +
    ggplot2::scale_x_discrete(expand = c(0, 0), limits = limits) +
    ggplot2::scale_y_discrete(expand = c(0, 0), limits = limits) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(0, "cm"),
      panel.background = ggplot2::element_blank(),
      text = element_text(family = "Arial", colour = "grey20"),
      axis.text = ggplot2::element_text(size = rel(axis_rel)),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1),
      legend.title = ggplot2::element_text(face = "bold", size = rel(guide_rel)),
      legend.text = ggplot2::element_text(size = rel(legend_rel)),
      legend.position = if(!show_legend) "none" else "right"
    )
}

# function for significance testing of correlation coefficients between groups
#' Significance testing of the difference of correlation coefficients to a null hypothesis of 0
#'
#' @param r1 Numeric value of a correlation coefficient.
#' @param r2 Numeric value of another correlation coefficient.
#' @param n1 The number of samples used to calculate r1.
#' @param n2 The number of samples used to calculate r2.
#' @param return A string indicating the desired output. "p" returns a numeric p-value.
#' "k" returns a length-2 numeric vector containing the z-scores of each correlation coefficient.
#' "z" returns the test statistic used in the difference between the correlation coefficients.
#'
#' @returns See return parameter.
#' @export
#'
#' @examples
#' r1 <- runif(n = 1, min = -1, max = 1)
#' r2 <- runif(n = 1, min = -1, max = 1)
#' fisher_r_to_z(r1 = r1, r2 = r2, n1 = 10, n2 = 10, return = "p")
fisher_r_to_z <- function(r1, r2, n1, n2, return = "p") {
  if(r1 == r2 & r1 == 1) {
    return(if(return == "p") 1 else if(return == "k") NA else if(return == "z") 0)
  }
  r <- c(r1, r2)
  k <- vector("numeric", 2)
  for(i in 1:2) {
    k[i] <- 0.5 * (log(1 + r[i]) - log(1 - r[i]))
  }
  z <- (k[1] - k[2]) / sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))
  p <- pnorm(q = z, mean = 0, sd = 1, lower.tail = if(z > 0) FALSE)*2
  return(if(return == "p") p else if(return == "k") k else if(return == "z") z)
}

# volcano plot of comparisons of correlation coefficients between groups
#' Volcano plot of significance testing of differences between two matrices of correlation coefficients
#'
#' @param cor1 Numeric matrix of correlation coefficients.
#' @param cor2 Numieric matrix of correlation coefficients with identical dimensions to cor1.
#' @param plab Numeric value of -log10(p) above which points will be labeled.
#' @param return A string indicating whether to return a "graph" or "data".
#' @param colours A character vector of length 3 containing colours to fill negative, small, and positive fold-change values.
#'
#' @returns See parameter return.
#' @export
#'
#' @examples
#' cor1 <- cor(x = data.frame(a = rep(c("a", "b"), each = 10), b = rnorm(n = 20), c = rnorm(n = 20), d = rnorm(n = 20))[1:10, 2:4])
#' cor2 <- cor(x = data.frame(a = rep(c("a", "b"), each = 10), b = rnorm(n = 20), c = rnorm(n = 20), d = rnorm(n = 20))[11:20, 2:4])
#' plot_volcano(cor1 = cor1, cor2 = cor2)
plot_volcano <- function(cor1, cor2, plab = .1, return = "graph",
                         colours = c("#FA819E", "#5f83ea", "grey70")) {
  # check if matrices have the same variables in the same order
  if(!identical(colnames(cor1), colnames(cor2)) |
     !identical(rownames(cor1), rownames(cor2)))
    stop("rmat and pmat must be correlation matrices with identical dimensions")

  # matrix for p-values of comparison of correlation coefficients
  p_fc <- matrix(nrow = nrow(cor1), ncol = ncol(cor1), dimnames = dimnames(cor1))
  for (i in 1:nrow(p_fc)) {
    for(k in 1:ncol(p_fc)) {
      p_fc[i,k] <- fisher_r_to_z(r1 = cor1[i,k], r2 = cor2[i,k], n = 7, return = "p")
    }
  }
  p_fc[lower.tri(p_fc, diag = T)] <- NA

  # matrix for fold change of correlation coefficients
  r_fc <- cor1 / cor2
  r_fc[lower.tri(r_fc, diag = T)] <- NA

  # store matrices in tidy data frame
  mats_to_tidy <- function(fc, p) {
    d <- data.frame(var1 = rep(rownames(fc), times = ncol(fc)),
                    var2 = rep(colnames(fc), each = nrow(fc)),
                    r = NA, p = NA)
    for(n in 1:nrow(d)) {
      d$r[n] <- fc[d$var1[n], d$var2[n]]
      d$p[n] <- p[d$var1[n], d$var2[n]]
    }
    return(d)
  }
  vol_df <- na.omit(mats_to_tidy(r_fc, p_fc))

  # store pair names, log fc, -log p, change type, and label names
  vol_df$neg_lp <- -log(vol_df$p, base = 10)
  vol_df$lfc <- log(abs(vol_df$r), base = 2)
  vol_df$pair <- NA
  vol_df$label <- NA
  vol_df$class <- NA
  for(n in 1:nrow(vol_df)) {
    vol_df$pair[n] <- paste(vol_df$var1[n], vol_df$var2[n], sep = "-")
    if(vol_df$p[n] < plab) {
      vol_df$label[n] <- vol_df$pair[n]
    } else {
      vol_df$label[n] <- NA
    }
    if(vol_df$p[n] < plab & vol_df$lfc[n] > 1) {
      vol_df$class[n] <- "Increase"
    } else {
      if(vol_df$p[n] < plab & vol_df$lfc[n] < -1) {
        vol_df$class[n] <- "Decrease"
      } else {
        vol_df$class[n] <- "Small"
      }
    }
  }

  if(return == "data")
    return(vol_df)

  # calculate axis limits
  ylim <- c(0, max(ceiling(vol_df$neg_lp / .5) * .5))
  xmin <- -ceiling(abs(min(vol_df$lfc)))
  xlim <- c(xmin, ceiling(max(vol_df$lfc)))
  if(min(vol_df$lfc) > 0)
    xlim[1] <- -xlim[1]

  # plot values
  graph <- ggplot2::ggplot(
    data = vol_df,
    mapping = ggplot2::aes(
      x = lfc, y = neg_lp, colour = class, fill = class, label = label)) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1), colour = "gray70", linetype = "dotted", linewidth = .4) +
    ggplot2::geom_hline(
      yintercept = -log(c(0.10, 0.05), base = 10),
      colour = "grey70", linetype = "dotted", linewidth = .4) +
    ggplot2::geom_point(shape = 16, alpha = 1, size = 2) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0), limits = xlim, breaks = (xlim[1] + 1):(xlim[2] - 1)) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0), limits = ylim, breaks = seq(0, ylim[2], .5)) +
    ggplot2::labs(fill = "Change Type", colour = "Change Type",
                  x = expression("Log"[2]*"FC"), y = expression("-Log"[10]*"p-value")) +
    ggrepel::geom_text_repel(max.overlaps = Inf, size = 2.5) +
    theme_cell(background = F, xticks = "down")
}

# function for location of k-means elbow
#' Scree plot of SS of k-means up to a specified value of k
#'
#' @param matrix Numeric matrix of correlation coefficients.
#' @param max.k Maximum value of k up to which SS will be displayed.
#'
#' @returns A ggplot scree plot of SS or k-means outputs up to max.k.
#' @export
#'
#' @examples
#' matrix <- cor(x = data.frame(a = rnorm(n = 20), b = rnorm(n = 20), c = rnorm(n = 20), d = rnorm(n = 20)))
#' k_means_scree_plot(matrix = matrix, max.k = 7)
k_means_scree_plot <- function(matrix, max.k) {
  list.k <- vector("list", max.k)
  within.ss <- vector("numeric", max.k)
  for(i in 1:max.k) {
    list.k[[i]] <- kmeans(matrix, i)
    within.ss[i] <- sum(list.k[[i]]$withinss)
  }
  plot <- ggplot2::ggplot(data.frame(x = 1:max.k, y = within.ss), aes(x = x, y = y)) +
    ggplot2::geom_line(colour = "black") +
    ggplot2::geom_point(shape = 21, fill = nature_palette("blue", 3), colour = nature_palette("blue", 5), size = 3) +
    ggplot2::scale_x_continuous(breaks = 1:max.k) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    ggplot2::labs(y = "Within-Cluster SS", x = "Number of Custers", title = "K-Means Elbow") +
    neuro::theme_nature()
  return(list(plot = plot, clusters = list.k))
}


















