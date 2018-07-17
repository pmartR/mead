require(shiny)
get_plot_output_list <- function(meta_data) {
  # get the class of each remaining column in the data frame
  column_class <- sapply(meta_data, class)
  numerical <- meta_data[, which(column_class %in% c("numeric"))]
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

# thank you stack overflow! http://grokbase.com/t/r/r-help/127vrkcj7e/r-add-leading-zeros
add.dt0 <- function(x, width=8 ){
  sprintf(paste('%0', width, switch(is.numeric(x)+1, 's', 'i'),
                sep=''), x)
}

get_checkbox_output_list <- function(meta_data) {
  # check to see if any dates parse
  #browser()
  date.formats = c("%m/%d/%Y", "%Y/%m/%d")
  for (col.idx in seq_len(ncol(meta_data))) {
    x <- meta_data[, col.idx]
    if (!(is.character(x) | is.integer(x)) | is.factor(x)) next #epoch time could be numeric
    if (all(is.na(x))) next
    for (format in date.formats) {
      complete.x <- !(is.na(x))
      if (is.integer(x)) {
        if (nchar(x[1]) < 8) {
          d <- as.Date(lubridate::parse_date_time(add.dt0(as.character(x)), format, quiet = TRUE))
        }else{
          d <- as.Date(lubridate::as_datetime(x))
        }
      } else{
        d <- as.Date(lubridate::parse_date_time(as.character(x), format, quiet = TRUE))
      }
      d.na <- d[complete.x]
      if (any(is.na(d.na))) next
      meta_data[, col.idx] <- d 
    }
  }
  # find factor columns
  column_class <- sapply(meta_data, class)
  categorical <- meta_data[, which(column_class %in% c("character", "factor", "logical", "Date", "integer"))]
  # If categorical check for non-uniqueness
  non_unique_columns <- which(unlist(lapply(categorical, function(x) length(unique(x)) != nrow(categorical) & !all(is.na(x)))))
  categorical <- categorical[, non_unique_columns]
  input_n_categorical <- ncol(categorical)
  # categorical boxes
  categorical_list <- lapply(1:input_n_categorical, function(i) {
    #think of a better way to capture the columnname, but for now map it back to raw
    column_num <- which(names(meta_data) == (names(categorical)[i]))
    #paste the column num 
    boxname = names(categorical)[i]
    labelname <- names(categorical)[i]
    box_output_object <- checkboxGroupInput(inputId = boxname, label = labelname, choices = unique(categorical[,i]), selected = unique(categorical[,i]), inline = TRUE)
  })
  #browser()
  do.call(tagList, categorical_list)
  return(categorical_list)
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
