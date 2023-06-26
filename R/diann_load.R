#' diann_load
#'
#' load data from DIAnn
#'
#' @param file Path for data file
#'
#' @return A dataframe
#'
#' @export
diann_load <- function (file){
  as.data.frame(data.table::fread(file, stringsAsFactors = FALSE))
}
