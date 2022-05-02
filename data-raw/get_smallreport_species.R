#load small example data, this data were obtained thanks to diann-rpackage from Vadim Demichev
small_report <- as.data.frame(data.table::fread("data-raw/diann_report.tsv", stringsAsFactors = FALSE))


### See all species available in seqinr swissprot bank
library(seqinr)
seqinr::choosebank("swissprot")
seqsp <- seqinr::query("seqsp","ST=REVIEWED", virtual=T) #get all sequences reviewed
sp <- seqinr::query("sp","PS seqsp")  #get species from those sequences
all_species <- seqinr::getName(sp)    #get name of this species
seqinr::closebank()

usethis::use_data(all_species, overwrite = TRUE)
usethis::use_data(small_report, overwrite = TRUE)
