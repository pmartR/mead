#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
filterBarWidget <-function(colName = "GridCell", dataset = data.frame(GridCell = c("A", "B", "C", "D"),
                                                                      Total = c(5,5,4,1)), width = NULL, height = NULL) {
  
  # forward options using x
  x = list(
    dataset = dataset,
    colName = colName
  )
  
  # create widget
  htmlwidgets::createWidget(
    name = 'filterBarWidget',
    x,
    width = width,
    height = height,
    package = 'filterBarWidget'
  )
}

#' Shiny bindings for filterBarWidget
#'
#' Output and render functions for using filterBarWidget within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a filterBarWidget
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name filterBarWidget-shiny
#'
#' @export
filterBarWidgetOutput <- function(outputId, width = '300px', height = '150px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'filterBarWidget', width, height, package = 'filterBarWidget')
}

#' @rdname filterBarWidget-shiny
#' @export
renderfilterBarWidget <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, filterBarWidgetOutput, env, quoted = TRUE)
}
