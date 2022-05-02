#load small example data, this data were obtained thanks to diann-rpackage from Vadim Demichev
small_report <- as.data.frame(data.table::fread("data-raw/diann_report.tsv", stringsAsFactors = FALSE))

usethis::use_data(small_report, overwrite = TRUE)
