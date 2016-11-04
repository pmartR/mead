sums_hist = function(thesums=NULL, xlab="", ylab=""){
  if(is.null(thesums)){
    p = qplot(0)
  } else {
    p = ggplot(data.frame(sums=thesums), aes(x=sums))
    p = p + geom_histogram()
    p = p + xlab(xlab) + ylab(ylab) 
    p = p + scale_x_log10(labels = scales::comma)
  }
  return(p)
}

output_phyloseq_print_html = function(physeq){
  HTML(
    paste(
      '<p class="phyloseq-print">',
      paste0(capture.output(print(physeq)), collapse=" <br/> "),
      "</p>"
    )
  )
}

numericInputRow <- function(inputId, label, value, min = NA, max = NA, step = NA, class="form-control", ...){
  inputTag <- tags$input(id = inputId, type = "number", value = value, class=class, ...)
  if (!is.na(min)) 
    inputTag$attribs$min = min
  if (!is.na(max)) 
    inputTag$attribs$max = max
  if (!is.na(step)) 
    inputTag$attribs$step = step
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      inputTag)
}

parse_date_columns <- function(..., date.formats = c("%m/%d/%Y", "%Y/%m/%d")) {
  dat <- read.table(...)
  for (col.idx in seq_len(ncol(dat))) {
    x <- dat[, col.idx]
    if (!is.character(x) | is.factor(x)) next
    if (all(is.na(x))) next
    for (format in date.formats) {
      complete.x <- !(is.na(x))
      d <- as.Date(parse_date_time(as.character(x), format, quiet = TRUE))
      d.na <- d[complete.x]
      if (any(is.na(d.na))) next
      dat[, col.idx] <- d         
    }
  }
  dat
  
}