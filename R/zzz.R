#check if limma is installed
.onAttach <- function(libname, pkgname){
  if (!requireNamespace("limma", quietly = TRUE)) {
    packageStartupMessage(
      "Please install 'limma' package by runing",
      " 'BiocManager::install('limma')'",
      " Otherwise, you won't be able to see MDS plot."
    )
  }
  packageStartupMessage(
    "\n",
    "Welcome to DIAgui package! To launch the app, run runDIAgui() function.\n",
    "To access the documentation, run browseVignettes(package = 'DIAgui').")

  WD <<- getwd()
  packageStartupMessage(
    "A variable named 'WD' has been created to make it easier to browse files\n",
    "when using the app. You can change its value but it must be a valid file path.")
}
