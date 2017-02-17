#---------------------------------------- Error Handling ----------------------------------------# 
#validate tests a condition and
#returns a validation error if the test fails. Validation errors are designed
#to interact with the Shiny framework in a pleasing way

biom_not_matching_metadata <- function(biom_input, metadata_input) {
  if (any(!(biom_input %in% metadata_input))) {
    "Mismatched sample names between biom and QIIME"
  } else if (any(sample_names(biom_input) %in% sample_names(metadata_input))) {
    NULL
  } else {
    FALSE
  }
}
