# cell-inspired theme for ggplot2
#' Theme for ggplot2 inspired by journals of the cell press
#'
#' @param background a logical object of length 1 indicating where to show or remove the grey background
#' @param xticks a string of length 1 (either "up" or "down") indicating the direction of the ticks on the x axis
#'
#' @returns a ggplot2 object
#' @export
#'
#' @examples
#' data <- data.frame(a = LETTERS[1:3], b = 1:3)
#' plot <- ggplot(data = data, aes(x = a, y = b)) +
#' geom_point()
#' plot
#' plot + theme_cell()
theme_cell <- function(background = T, xticks = "up") {
  if(xticks != "up" & xticks != "down") {
    stop("xticks must be set to either up or down")
  }
  require(ggplot2)
  ggplot2::theme(
    axis.line.y = ggplot2::element_line(colour = "grey50", linewidth = .7),
    axis.ticks.y = ggplot2::element_line(colour = "grey50", linewidth = .7),
    axis.ticks.x = ggplot2::element_line(colour = "grey50", linewidth = .5),
    panel.background = ggplot2::element_rect(
      fill = if(background) "grey95" else "white"),
    panel.grid = ggplot2::element_line(colour = NA),
    panel.grid.major.y = ggplot2::element_line(colour = "white"),
    axis.ticks.length.x = ggplot2::unit(
      if(xticks == "up") -.1 else .1, "cm"),
    text = ggplot2::element_text(family = "Arial", colour = "grey20"),
    axis.title = ggplot2::element_text(face = "bold", size = rel(.9)),
    legend.title = ggplot2::element_text(face = "bold", size = rel(.9)))
}



# nature-inspired theme for ggplot2
#' Theme for ggplot2 inspired by springer nature journals
#'
#' @returns a ggplot2 object
#' @export
#'
#' @examples
#' data <- data.frame(a = LETTERS[1:3], b = 1:3)
#' plot <- ggplot(data = data, aes(x = a, y = b)) +
#' geom_point()
#' plot
#' plot + theme_nature()
theme_nature <- function() {
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 panel.grid = ggplot2::element_blank(),
                 strip.background = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(size = 10, face = "bold"),
                 text = ggplot2::element_text(family = "Helvetica Neue"),
                 axis.line = ggplot2::element_line(colour = "black", linewidth = .2),
                 axis.ticks = ggplot2::element_line(colour = "black", linewidth = .2),
                 axis.ticks.length = ggplot2::unit(1.75, "mm"),
                 axis.title = ggplot2::element_text(size = 11),
                 axis.text = ggplot2::element_text(size = 10),
                 plot.title = ggplot2::element_text(size = 11, hjust = .5),
                 legend.title = ggplot2::element_text(size = 11),
                 legend.text = ggplot2::element_text(size = 10))
}














