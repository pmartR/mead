require(shiny)
get_plot_output_list <- function(meta_data) {
  # try to coerce columns to Date
  date_formats <- c("%d/%m/%Y","%m/%d/%y")
  date_column <- which(sapply(meta_data, function(x) !all(is.na(as.Date(as.character(x),format = date_formats)))))
  if (length(date_column == 1)) {
    #date_obj <- meta_data[,date_column]
    meta_data <- meta_data[,-date_column] 
    # TODO: Render Date as a bar chart
  }
  # get the class of each remaining column in the data frame
  column_class <- sapply(meta_data, class)
  numerical <- meta_data[, which(column_class %in% c("numeric", "integer"))]
  input_n_numeric <- ncol(numerical)
  
  # numerical plots
  numerical_list <- lapply(1:input_n_numeric, function(i) {
    plotname <- paste("plot", i, sep = "")
    plot_output_object <- filterWidget::filterWidgetOutput(plotname, height = 10, width = 20)
    plot_output_object <- renderfilterWidget({
      filterWidget(as.character(names(numerical)[i]), meta_data)
    })
  })
  # TODO: combine bar charts, dates, and histograms into one
  # plot_output_list <- c(numerical_list, categorical_list)
  plot_output_list <- numerical_list
  do.call(tagList, plot_output_list) # needed to display properly.
  return(plot_output_list)
}

get_checkbox_output_list <- function(meta_data) {
  column_class <- sapply(meta_data, class)
  categorical <- meta_data[, which(column_class %in% c("character", "factor", "logical"))]
  # If categorical check for non-uniqueness
  non_unique_columns <- which(unlist(lapply(categorical, function(x) length(unique(x)) != nrow(categorical) & length(unique(x)) != 1)))
  categorical <- categorical[, non_unique_columns]
  input_n_categorical <- ncol(categorical)
  # categorical boxes
  categorical_list <- lapply(1:input_n_categorical, function(i) {
    #think of a better way to capture the columnname, but for now map it back to raw
    column_num <- which(names(meta_data) == (names(categorical)[i]))
    #paste the column num 
    boxname = names(categorical)[i]
    labelname <- names(categorical)[i]
    box_output_object <- checkboxGroupInput(inputId = boxname, label = labelname, choices = unique(categorical[,i]), selected = unique(categorical[,i]))
  })
  do.call(tagList, categorical_list)
  return(categorical_list)
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