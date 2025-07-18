# gradients of 6
#' Colour palettes used in nature journals
#'
#' @param colour A string indicating the hue of the desired pallette. Must be one of:
#' "background", "grey", "olive", "green", "teal", "blue", "purple", "red", "orange", "yellow".
#' If NULL (the default), returns all hues.
#' @param n an integer between 1 and 6 inclusive. Indicates the index of the palette desired.
#'
#' @returns If neither colour nor n supplied, a list of legnth 10 with each element containing
#' a character vector of hex codes of length 6. If colour but not n supplied, a character vector of
#' hex codes of length 6.
#' @export
#'
#' @examples
#' data <- data.frame(a = LETTERS[1:3], b = 1:3)
#' plot <- ggplot(data = data, aes(x = a, fill = a, y = b)) +
#' geom_point()
#' plot
#' plot + scale_fill_manual(values = c(nature_palette(colour = "blue", n = 3),
#' nature_palette(colour = "purple", n = 3),
#' nature_palette(colour = "teal", n = 3)))
nature_palette <- function(colour = NULL, n = NULL) {
  if(!is.null(n)) {
    if(n %% 1 != 0 | n < 1 | 6 < n) {
      stop("n must be an integer in the range of 1 to 6")
    }
  }
  colours <- list(
    background = c("#f6f2ec", "#e3ddca", "#c6bfa3",
                   "#a7a083", "#868062", "#605c47"),
    grey = c("#e5e3ea", "#c7ccd8", "#97a1b0",
             "#6e788c", "#46536a", "#202e43"),
    olive = c("#efedb0", "#d8d465", "#c0bf31",
              "#959f38", "#637331", "#344319"),
    green = c("#d7e3c3", "#9dc375", "#62aa44",
              "#468e41", "#256d37", "#18381b"),
    teal = c("#cbe2e9", "#95cbcd", "#50b3b3",
             "#14939b", "#006475", "#01354a"),
    blue = c("#c7e2f4", "#9cc4e2", "#5790c4",
             "#086ca8", "#0b4b84", "#192b57"),
    purple = c("#e7d3e6", "#cea5c8", "#b376aa",
               "#a04a8c", "#75256d", "#46174d"),
    red = c("#f6ccc9", "#e6a0a2", "#d46261",
            "#bd3b3b", "#952422", "#701311"),
    orange = c("#f8dcbb", "#f5bc7d", "#ed9645",
               "#e06c24", "#af4d25", "#7f3218"),
    yellow = c("#fcedc0", "#f0d685", "#e3c24e",
               "#c59c2b", "#99752a", "#675221")
  )
  if(!is.null(colour)) {
    if(!grepl(paste0(names(colours), collapse = "|"), colour)) {
      stop(paste("colour must be one of:", paste0(names(colours), collapse = ", ")))
    }
  }
  if(is.null(colour) & is.null(n)) {
    return(colours)
  } else {
    if(!is.null(colours) & is.null(n)) {
      return(colours[[colour]])
    } else {
      if(is.null(colour) & !is.null(n)) {
        n.colours <- vector("list", length(colours))
        names(n.colours) <- names(colours)
        for(k in 1:length(colours)) {
          n.colours[[k]] <- colours[[k]][n]
        }
        return(n.colours)
      } else {
        return(colours[[colour]][n])
      }
    }
  }
}

# five saturated colours
#' Saturated colours
#'
#' @param colour a string containing the desired colour. One of: "purple", "green", "blue", "yellow", or "red".
#' NULL returns all colours.
#'
#' @returns A character vector of hex codes of length 1 if colour supplied, or length 5 if not.
#' @export
#'
#' @examples
#' data <- data.frame(a = LETTERS[1:3], b = 1:3)
#' plot <- ggplot(data = data, aes(x = a, fill = a, y = b)) +
#' geom_point()
#' plot
#' plot + scale_fill_manual(values = c(colours5(purple = "blue",
#' colours5("yellow"),
#' colours5("red")))
colours5 <- function(colour = NULL) {
  colours <-
    list(purple = "#5f83ea",
         green = "#68c7a3",
         blue = "#5bc3e8",
         yellow = "#f7ce67",
         red = "#FA819E")
  if(!is.null(colour)) {
    if(!grepl(paste0(names(colours), collapse = "|"), colour)) {
      stop(paste("colour must be one of:", paste0(names(colours), collapse = ", ")))
    } else {
      return(colours[[colour]])
    }
  } else {
    return(colours)
  }
}
