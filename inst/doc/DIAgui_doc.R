## ---- eval=FALSE--------------------------------------------------------------
#  if(!requireNamespace("devtools", quietly = TRUE)){ #check if you already have the devtools package
#   install.packages("devtools")  #if not, install it
#  }
#  devtools::install_github("marseille-proteomique/DIAgui")

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  library("DIAgui")

## ---- eval=FALSE--------------------------------------------------------------
#  runDIAgui()      #this function will directly start the app

## ---- eval=FALSE--------------------------------------------------------------
#  report <- small_report
#  View(report) #take a look at the data

## ---- eval=FALSE--------------------------------------------------------------
#  readr::write_tsv(report, "small_report.tsv") # save the data as a tsv file

## ---- eval=FALSE--------------------------------------------------------------
#  report <- diann_load("small_report.tsv")

## ---- eval=FALSE--------------------------------------------------------------
#  precursor <- diann_matrix(report, proteotypic.only = TRUE, method = "max")
#  
#  head(precursor)  # check the results

## ---- eval=FALSE--------------------------------------------------------------
#  peptide <- diann_matrix(report, id.header = "Modified.Sequence",
#                          proteotypic.only = TRUE, method = "max")
#  
#  head(peptide)  # check the results

## ---- eval=FALSE--------------------------------------------------------------
#  peptide_maxlfq <- report %>% dplyr::filter(Q.Value <= 0.01 & PG.Q.Value <= 0.01 & Protein.Q.Value <= 1 & GG.Q.Value <= 1)
#  peptide_maxlfq <- iq::preprocess(peptide_maxlfq,
#                              intensity_col = "Precursor.Normalised",
#                              primary_id = "Modified.Sequence.",
#                              sample_id  = "File.Name",
#                              secondary_id = "Precursor.Id",
#                              median_normalization = FALSE,
#                              pdf_out = NULL)
#  peptide_maxlfq <- iq::fast_MaxLFQ(peptide_maxlfq)
#  peptide_maxlfq <- peptide_maxlfq$estimate
#  peptide_maxlfq <- as.data.frame(peptide_maxlfq)
#  
#  head(peptide_maxlfq)  # check the results

## ---- eval=FALSE--------------------------------------------------------------
#  peptide_maxlfq <- report %>% dplyr::filter(Q.Value <= 0.01 & PG.Q.Value <= 0.01 & Protein.Q.Value <= 1 & GG.Q.Value <= 1)
#  peptide_maxlfq <- diann_maxlfq(peptide_maxlfq,
#                            group.header = "Modified.Sequence",
#                            id.header = "Precursor.Id",
#                            quantity.header = "Precursor.Normalised",
#                            count_pep = FALSE
#          )
#  
#  head(peptide_maxlfq)  # check the results

## ---- eval=FALSE--------------------------------------------------------------
#  protein <- re %>% report %>% dplyr::filter(Q.Value <= 0.01 & PG.Q.Value <= 0.01 & Protein.Q.Value <= 1 & GG.Q.Value <= 1)
#  n_cond <- length(unique(df$File.Name))
#  
#  protein_maxlfq <- diann_maxlfq(protein,
#                    group.header="Protein.Group",
#                    id.header = "Precursor.Id",
#                    quantity.header = "Precursor.Normalised",
#                    only_countsall = FALSE,
#                    Top3 = TRUE
#                    )
#  
#  # extract useful information from report
#  nc <- ncol(protein_maxlfq)
#  protein_maxlfq$Protein.Group <- rownames(protein_maxlfq)
#  rownames(protein_maxlfq) <- 1:nrow(protein_maxlfq)
#  protein <- protein[(protein$Protein.Group %in% protein_maxlfq$Protein.Group),]
#  protein <- protein[order(protein$Protein.Group),]
#  protein_maxlfq$Protein.Names <- unique(protein[,c("Protein.Group", "Protein.Names")])$Protein.Names
#  protein_maxlfq$First.Protein.Description <- unique(protein[,c("Protein.Group", "First.Protein.Description")])$First.Protein.Description
#  protein_maxlfq$Genes <- unique(protein[,c("Protein.Group", "Genes")])$Genes
#  protein_maxlfq <- protein_maxlfq[,c((nc+1):ncol(protein_maxlfq), 1:nc)]
#  
#  # get iBAQ quantification
#  protein_seq <- getallseq(pr_id = protein_maxlfq$Protein.Group,
#                           spec = "SACCHAROMYCES CEREVISIAE")
#  # here, the function makes a query to swissprot to get the amino-acid sequence from each protein, so it can be long
#  # you can also put one or several FASTA files (go check getallseq documentation)
#  
#  # to compute iBAQ quantification you'll the raw intensities
#  raw <- diann_matrix(report, id.header = "Protein.Group",
#                               quantity.header = "Precursor.Quantity",
#                               method = "sum")
#  raw$Genes <- NULL
#  raw$Protein.Names <- NULL
#  # compute iBAQ quantification
#  raw <- get_iBAQ(raw, proteinDB = protein_seq,
#                  id_name = "Protein.Group",
#                  ecol = 2:(n_cond+1),
#                  peptideLength = c(5,36),
#                  proteaseRegExp = DIAgui:::getProtease("trypsin"),
#                  log2_transformed = FALSE)
#  
#  raw <- raw[,-c(2:(n_cond+1))]
#  
#  protein_maxlfq <- dplyr::left_join(protein_maxlfq, raw, by = "Protein.Group")
#  
#  
#  head(protein_maxlfq)  # check the results

## ---- eval=FALSE--------------------------------------------------------------
#  genes <- diann_matrix(report,
#                     id.header="Genes",
#                     quantity.header="Genes.MaxLFQ.Unique",
#                     proteotypic.only = TRUE,
#                     get_pep = TRUE, only_pepall = TRUE,
#                     Top3 = TRUE,
#                     method = "max")
#  
#  head(genes)  # check the results

## ---- eval=FALSE--------------------------------------------------------------
#  report_process("small_report.tsv", # needs to be path to a report file
#                 get_iBAQ = TRUE, get_Top3 = FALSE,
#                 species = "SACCHAROMYCES CEREVISIAE")

## ---- eval=FALSE--------------------------------------------------------------
#  names(genes)[1] <- "id"

## ---- eval=FALSE--------------------------------------------------------------
#  # density plot
#  densityDIA(genes, transformation =  "log2", area = TRUE, data_type = "intensity", tit = "My plot")
#  
#  # MDS plot
#  MDS_DIA(genes, transformation =  "log2", data_type = "intensity", tit = "My plot")
#  
#  # interactive heatmap
#  heatmapDIA(genes, transformation = "log2", print_val = FALSE, data_type = "intensity")

## ---- eval=FALSE--------------------------------------------------------------
#  validDIA(genes, "log2", data_type = "intensity", prop_cut = 0.75) # we keep approximately 90% of the genes
#  
#  validDIA(genes, "log2", data_type = "intensity", prop_cut = 0.3) # we keep all genes

