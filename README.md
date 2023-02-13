# DIAgui
DIAgui is an R package that contains a user-friendly shiny app to process output from DIA-nn proteomics software. Since shiny can be quite slow to process big data file,
a function that does the complete workflow and saves all results in a directory is also present. 
The process is as follow: filtering data according q.values, applying MaxLFQ algorithm for quantification, getting some other informations, plot some graphs to check data quality.

This package is based on [diann-rpackage](https://github.com/vdemichev/diann-rpackage) from Vadim Demichev. The c++ source code is exactly the same as the one from V.Demichev,
it is used for MaxLFQ algorithm. You can also use [iq](https://cran.r-project.org/web/packages/iq/index.html) R package which is way faster to run MaxLFQ algorithm.

In addition to filter data and calculate MaxLFQ intensities, you can also get the Top3 and iBAQ quantification within the app.

# How to install and use DIAgui ?
First, go to Rstudio. Install DIAgui from github: 

```c
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools") 
}
devtools::install_github("mgerault/DIAgui")
```

You can now load it and run the app with this commands: 

```c
library(DIAgui)
runDIAgui()
```

If you want to learn more about the app, check the [tutorial](https://www.youtube.com/watch?v=vfvh15Q93eU) !
