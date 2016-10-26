#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
barwidget <- function(dataset = data.frame(GridCell = c("A", "B", "C", "D"),
                                           Total = c(5,5,4,1)), width = NULL, height = NULL) {

  # forward options using x
  x = list(
    dataset = dataset
  )

  # create widget
  htmlwidgets::createWidget(
    name = 'barwidget',
    x,
    width = width,
    height = height,
    package = 'barwidget'
  )
}

#' Shiny bindings for barwidget
#'
#' Output and render functions for using barwidget within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a barwidget
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name barwidget-shiny
#'
#' @export
barwidgetOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'barwidget', width, height, package = 'barwidget')
}

#' @rdname barwidget-shiny
#' @export
renderBarwidget <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, barwidgetOutput, env, quoted = TRUE)
}
