#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
filterWidget <- function(colName = "GridCell", dataset = data.frame(GridCell = c("A", "B", "C", "D"),
                                           Total = c(5,5,4,1)), width = NULL, height = NULL) {

  # forward options using x
  x = list(
    dataset = dataset,
    colName = colName
  )

  # create widget
  htmlwidgets::createWidget(
    name = 'filterWidget',
    x,
    width = width,
    height = height,
    package = 'filterWidget'
  )
}

#' Shiny bindings for filterWidget
#'
#' Output and render functions for using filterWidget within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a filterWidget
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name filterWidget-shiny
#'
#' @export
filterWidgetOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'filterWidget', width, height, package = 'filterWidget')
}

#' @rdname filterWidget-shiny
#' @export
renderfilterWidget <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, filterWidgetOutput, env, quoted = TRUE)
}
