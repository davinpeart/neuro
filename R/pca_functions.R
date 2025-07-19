# scree plot of variance explained
#' Scree plot of explained variance from prcomp object
#'
#' @param pc PCA computed using prcomp.
#' @param ncomp Number of components up to which to display.
#' @param nret Number of components to retain.
#' @param return A string, either "data" or "plot".
#'
#' @returns See return.
#' @export
#'
#' @examples
#' data <- data.frame(a = rnorm(n = 20), b = rnorm(n = 20), c = rnorm(n = 20), d = rnorm(n = 20))
#' pc <- prcomp(x = data)
plot_importance <- function(pc, ncomp = 9, nret = 3, return = NULL) {
  require(ggplot2)
  imp <- as.data.frame(summary(pc)$importance)[c("Proportion of Variance",
                                                "Cumulative Proportion"),]
  imp$var <- rownames(imp)
  df <- data.frame(dim = rep(colnames(imp)[1:ncol(imp)-1], each = 2))
  df$var = rep(imp[["var"]], times = nrow(df)/2)
  df$val <- vector("numeric", nrow(df))
  for(n in 1:nrow(df)) {
    df[["val"]][n] <- imp[[df[["dim"]][n]]][grep(df[["var"]][n], imp[["var"]])]
  }
  df$cpr <- vector("numeric", nrow(df))
  dims <- unique(df$dim)
  for(k in 1:length(dims)) {
    df[df["dim"] == dims[k], ]$cpr <-
      df[df["dim"] == dims[k] & df["var"] == "Cumulative Proportion", ]$val
  }
  df <- df[order(df$var, decreasing = F), ]
  df$dim <- fct_reorder(as.factor(df$dim), df$cpr)
  df$ret <- vector("character", nrow(df))
  df$ret[df["var"] == "Proportion of Variance"] <- "Individual"
  df$ret[df["var"] == "Cumulative Proportion" &
           as.numeric(str_extract(df$dim, "[:digit:][:digit:]?")) <= nret] <- "Cumulative"
  df$ret[df["var"] == "Cumulative Proportion" &
           as.numeric(str_extract(df$dim, "[:digit:][:digit:]?")) > nret] <- "Not Retained"
  df$ret <- as.factor(df$ret)

  if(return == "data") return(df)

  ggplot2::ggplot(data = df,
                  ggplot2::aes( x = dim,  y = val, group = ret, colour = ret,
                                fill = ret, alpha = ret)) +
    ggplot2::geom_bar(stat = "identity", position = "identity") +
    ggplot2::geom_point(data = df[df["var"] == "Proportion of Variance", ],
                        shape = 16, size = 3, colour = "#495576",
                        alpha = .6) +
    ggplot2::geom_line(data = df[df["var"] == "Proportion of Variance", ],
                       linewidth = 1, colour = "#495576",
                       alpha = .6) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = c(0, 0)) +
    ggplot2::scale_fill_discrete(labels = c("Cumulative", "Individual", "Not Retained"),
                                 type = c("#5A6D9F", "#7E92C7",
                                          "grey75"), name = "Type") +
    ggplot2::scale_colour_discrete(labels = c("Cumulative", "Individual", "Not Retained"),
                                   type = c(nature_palette("background", 4),
                                            "grey40",
                                            nature_palette("background", 4)), name = "Type") +
    scale_alpha_manual(values = c("Cumulative" = .6,
                                  "Individual" = 1,
                                  "Not Retained" = .6), name = "Type") +
    ggplot2::labs(x = "Principal Component", y = "Proportion of Variance") +
    ggplot2::scale_x_discrete(labels = str_extract(levels(df$dim), "[:digit:][:digit:]?"),
                              limits = levels(df$dim)[1:ncomp]) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_line(colour = "grey50", linewidth = .7),
      axis.ticks.y = ggplot2::element_line(colour = "grey50", linewidth = .7),
      axis.ticks.x = ggplot2::element_line(colour = "grey50", linewidth = .5),
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid = ggplot2::element_line(colour = NA),
      panel.grid.major.y = ggplot2::element_line(colour = "white"),
      axis.ticks.length.x = ggplot2::unit(.05, "cm"),
      text = ggplot2::element_text(family = "Arial", colour = "grey20"),
      axis.title = ggplot2::element_text(face = "bold", size = rel(.9)),
      legend.title = ggplot2::element_text(face = "bold", size = rel(.9)))
}

# plot r-squared between original variables and principal components
plot_rsquared <- function(pc, data, nvar = 12, ncomp = 3, return = NULL,
                          pcvars = TRUE, pcobs = TRUE, bars = "#7E92C7") {
  if(!is.null(nvar)) {
    if(nvar %% 1 != 0 | nvar < 1 | nvar > ncol(data)) {
      stop("nvar must be a positive integer less than the number of columns in data")
    }
  }
  if(!is.null(ncomp)) {
    if(ncomp %% 1 != 0 | ncomp < 1 | ncomp > ncol(pc$x)) {
      stop("nvar must be a positive integer")
    }
  }
  require(ggplot2)
  require(patchwork)
  temp <- cor(cbind(pc$x[pcobs, 1:ncomp], data[pcobs, pcvars]), method = "pearson")^2
  temp[lower.tri(temp, diag = T)] <- NA
  temp <- as.data.frame(na.omit(t(temp[,-c(1:ncomp)])[, 1:ncomp]))
  if(ncomp == 1)
    names(temp) <- "PC1"
  ulim <- max(ceiling(temp/.1)*.1)
  df <- vector("list", ncomp)
  for(i in 1:ncomp) {
    df[[i]] <- temp[i]
    df[[i]]$dim <- colnames(df[[i]])
    colnames(df[[i]])[1] <- "prop"
    df[[i]]$var <- as.factor(rownames(df[[i]]))
    df[[i]]$var <- fct_reorder(df[[i]]$var, df[[i]]$prop, .desc = T)
  }
  if(return == "data") {
    for(i in 1:length(df)) {
      names(df)[i] <- unique(df[[i]]$dim)
    }
    return(df)
  }
  plots <- vector("list", ncomp)
  for(n in 1:ncomp) {
    plots[[n]] <-
      ggplot2::ggplot(df[[n]], ggplot2::aes(y = var, x = prop)) +
      ggplot2::geom_bar(stat = "identity", fill = bars, alpha = .7) +
      ggplot2::geom_text(ggplot2::aes(label = var, x = 0), hjust = 0, vjust = .5,
                         size = 2.5, colour = "grey20", nudge_x = 0.005) +
      ggplot2::labs(x = "Proportion of Variance", y = "Variable", title = unique(df[[n]]$dim)) +
      ggplot2::theme(
        axis.line.x = ggplot2::element_line(colour = "grey50", linewidth = .5),
        axis.ticks.x = ggplot2::element_line(colour = "grey50", linewidth = .5),
        axis.ticks.y = ggplot2::element_line(colour = "grey50", linewidth = .5),
        panel.background = ggplot2::element_rect(fill = "grey95"),
        panel.grid = ggplot2::element_line(colour = NA),
        panel.grid.major.x = ggplot2::element_line(colour = "white"),
        axis.ticks.length.y = ggplot2::unit(0, "cm"),
        axis.text.y = ggplot2::element_blank(),
        text = ggplot2::element_text(family = "Arial", colour = "grey20"),
        axis.title = ggplot2::element_text(face = "bold", size = rel(.9)),
        plot.title = ggplot2::element_text(face = "bold", size = rel(.9),
                                           hjust = .5),
        plot.margin = ggplot2::unit(c(.3,.5,.3,.3), "cm")) +
      ggplot2::scale_y_discrete(limits = levels(df[[n]]$var)[1:nvar],
                                expand = c(0, .6)) +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .25),
                                  expand = c(0, 0))
  }
  if(ncomp == 1)
    return(plots[[1]]) else
  cplot <- plots[[1]] + plots[[2]]
  if(ncomp == 2)
    return(cplot + patchwork::plot_layout(axis_titles = "collect")) else
  for(k in 3:ncomp) {
    cplot <- cplot + plots[[k]]
  }
  return(cplot + patchwork::plot_layout(axis_titles = "collect"))
}

