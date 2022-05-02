#' Run the dia app
#'
#' \code{runDIAgui} works exactly the same as runExample from \code{\link{shiny}} package.
#'
#' @export
runDIAgui <- function(){
  appDir <- system.file("shiny-examples", "myapp", package = "DIAgui")
  if(appDir == ""){
    stop("Couldn't find example directory. Try re-installing 'DIAgui'.", call = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
