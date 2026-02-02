#' Custom Name Repair
#'
#' Used in the survFMM function when only two subgroups are specified
#'
#' @param names
#'
#' @returns Vector of repaired names
#'
#' @examples
custom_name_repair <- function(names) {
  if (length(names) == 1 && names %in% c("value")) {
    names <- "2"
  } else {
    names
  }
  return(names)
}
