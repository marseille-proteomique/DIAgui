# DIAgui
DIAgui is an R package that contains a user-friendly shiny app to process output from DIA-nn proteomics software. Since shiny can be quite slow to process big data file,
a function that does the complete workflow and saves all results in a directory is also present. 
The process is as follow: filtering data according q.values, applying MaxLFQ algorithm for quantification, getting some other informations, plot some graphs to check data quality.

This package is based on [diann-rpackage](https://github.com/vdemichev/diann-rpackage) from Vadim Demichev. The c++ source code is exactly the same as the one from V.Demichev,
it is used for MaxLFQ algorithm. You can also use [iq](https://cran.r-project.org/web/packages/iq/index.html) R package which is way faster to run MaxLFQ algorithm.

In addition to filter data and calculate MaxLFQ intensities, you can also get the Top3 and iBAQ quantification within the app.

### Publication

DOI: https://doi.org/10.1093/bioadv/vbae001 

# How to install and use DIAgui ?

### Requirement

* R version >= 4.0
* Rstudio 

### Installation
Go to Rstudio. Install DIAgui from github: 

```c
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools") 
}
devtools::install_github("marseille-proteomique/DIAgui")
```

You can now load it and run the app with this commands: 

```c
library(DIAgui)
runDIAgui()
```

If you want to learn more about the app, check the [tutorial](https://www.youtube.com/watch?v=vfvh15Q93eU) !

### Toy dataset

If you want to test the application, you can download a toy dataset [here](https://drive.google.com/file/d/1BVAGqKIkdqIqhebunM7K_FdSrOhXudA9/view?usp=sharing) which contain a small report from saccharomyces cerevisiae with its corresponding FASTA file. For information about the data, you access access its documentation in R via ``` ?small_report ```
