.onAttach <- function(libname, pkgname){
  packageStartupMessage(
    "\n",
    "Welcome to DIAgui package! To launch the app, run runDIAgui() function.\n",
    "To access the documentation, run browseVignettes(package = 'DIAgui').")

  WD <<- getwd()
  packageStartupMessage(
    "A variable named 'WD' has been created to make it easier to browse files\n",
    "when using the app. You can change its value but it must be a valid file path.")
}
