#---------------------------------------- Error Handling ----------------------------------------# 
#validate tests a condition and
#returns a validation error if the test fails. Validation errors are designed
#to interact with the Shiny framework in a pleasing way

biom_not_matching_metadata <- function(biom_input, metadata_input) {
  if (all(!(biom_input %in% metadata_input))) {
    "Samples in biom do not match samples in QIIME"
  } else if (all(sample_names(biom_input) %in% sample_names(metadata_input))) {
    NULL
  } else {
    FALSE
  }
}

biom_mismatching_metadata <- function(biom_input, metadata_input) {
  if (any(!(biom_input %in% metadata_input))) {
    TRUE
  } else if (any(biom_input %in% metadata_input)) {
    FALSE
  } else {
    FALSE
  }
}

metadata_not_matching_biom <- function(biom_input, metadata_input) {
  if (all(!(metadata_input %in% biom_input))) {
    "Samples in QIIME do not match samples in biom"
  } else if (any(biom_input %in% metadata_input)) {
    NULL
  } else {
    FALSE
  }
}

metadata_mismatching_biom <- function(biom_input, metadata_input) {
  if (any(!(metadata_input %in% biom_input))) {
    TRUE
  } else if (any(biom_input %in% metadata_input)) {
    FALSE
  } else {
    FALSE
  }
}