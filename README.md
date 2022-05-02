# DIAgui
DIAgui is a r package that contains a user-friendly shiny app to process output from DIA-nn proteomics software. Since shiny can be quite slow to process big data file,
a function that does the complete workflow and saves all results in a directory is also present. 
The process is as follow : filtering data according q.values, applying MaxLFQ algorithm for quantification, getting some other informations, plot some graphs to check data quality.

This package is based on [diann-rpackage](https://github.com/vdemichev/diann-rpackage) from Vadim Demichev. The c++ source code is exactly the same as the one from V.Demichev,
it is used for MaxLFQ algorithm. Since [iq](https://cran.r-project.org/web/packages/iq/index.html) r package is way faster, this code will be removed.

# How to install and use DIAgui ?
First, go to Rstudio. Before installing DIAgui, you will need to install [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
package from bioconductor in order to use all functionnalities from the app.
Run this commands :

* if(!requireNamespace("BiocManager", quietly = TRUE)){

* install.packages("BiocManager")  

* }

* BiocManager::install("limma")  

You can now install DIAgui from github : 

* if(!requireNamespace("devtools", quietly = TRUE)){

* install.packages("devtools") 

* }

* devtools::install_github("mgerault/DIAgui")

You can now load it and run the app with this commands : 

* library(DIAgui)

* runDIAgui()
