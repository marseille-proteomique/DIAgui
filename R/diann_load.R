#' diann_load
#'
#' load data from DIAnn
#'
#' @param file Path for data file
#'
#' @return A dataframe
#'
#' @export
diann_load <- function(file){
  if(grepl("\\.parquet$", file)){
    df <- as.data.frame(arrow::read_parquet(file))
  }
  else{
    df <- as.data.frame(data.table::fread(file, stringsAsFactors = FALSE))
  }

  return(df)
}
