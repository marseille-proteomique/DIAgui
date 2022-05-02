#check if limma is installed
.onAttach <- function(libname, pkgname){
  if(!("limma" %in% rownames(installed.packages()))){
    packageStartupMessage(
      paste0(
        "Please install 'limma' package by runing",
        " 'BiocManager::install('limma')'",
        " Otherwise, you won't be able to see MDS plot."
      )
    )
  }
}

.onLoad <- function(...){
  WD <<- getwd()
  message("A variable named 'WD' has been created in order to facilitate file exploration when you are using the app.
           You can change it but it has to be a file path.")
  message("")
  message("Welcome to DIAgui package ! To run the app, run runDIAgui() function.")
}
