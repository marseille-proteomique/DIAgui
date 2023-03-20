#' report_process
#'
#' Process data from DIAnn in a single workflow, save it and also save some quality check plots.
#'
#' @param data Path to data (report from DIAnn)
#' @param header.id Id column name (protein, gene, ...)
#' @param sample.id Id column that contains the fractions names
#' @param quantity.id Quantity column name (Precursor.Normalised)
#' @param secondary.id Id column that contains precursor (usefull when get_pep is TRUE)
#' @param id_to_add Some secondary id you want to keep in your final data
#' @param qv Precursor q-value threshold
#' @param pg.qv Protein group q-value threshold
#' @param p.qv Uniquely identified protein q-value threshold
#' @param gg.qv Gene group q-value threshold
#' @param only_proteotypic If TRUE, will only keep proteotypic
#' @param get_pep logical; get peptide count ?
#' @param only_pepall logical; should only keep peptide counts all or also peptide counts for each fractions ?
#' @param get_Top3 If TRUE, will calculate Top3 quantification
#'                 (only use when \code{header.id} is set to 'Protein.group' or 'Genes)
#' @param get_iBAQ If TRUE, will calculate iBAQ quantification
#'                 (only use when \code{header.id} is set to 'Protein.group' )
#' @param fasta When iBAQ is set to TRUE, you can add path to fasta files.
#'              If NULL, will get theoretical peptides from swissprot bank (quite long).
#' @param species When fasta is NULL, you need to specify the specie from your data
#' @param peptide_length The minimum and maximum peptide length
#' @param format Format from your saved file (only xlsx, csv and txt are supported)
#'
#' @importFrom utils write.csv write.table
#'
#' @export

report_process <- function(data, header.id = "Protein.Group", sample.id = "File.Name",
                           quantity.id = "Precursor.Normalised", secondary.id = "Precursor.Id",
                           id_to_add = c("Protein.Names", "First.Protein.Description", "Genes"),
                           qv = 0.01, pg.qv = 0.01, p.qv = 1, gg.qv = 1, only_proteotypic = TRUE,
                           get_pep = TRUE, only_pepall = FALSE, get_Top3 = FALSE, get_iBAQ = FALSE,
                           fasta = NULL, species = NULL, peptide_length = c(5,36),
                           format = c("xlsx", "csv", "txt")){

  if(inherits(data, "character")){
    if(file.exists(data)){
      extension <- stringr::str_split(data, "\\.")[[1]]
      if(length(extension) == 1){
        stop("Don't forget to type the format of your file !")
      }

      dname <- extension[length(extension)-1]
      dname <- stringr::str_split(dname, "/")[[1]]
      dname <- dname[length(dname)]
      dname <- paste(dname, header.id, sep="_")

      extension <- extension[length(extension)]
    }
    else{
      stop("This file doesn't exists ! Check your file name")
    }
  }
  else {
    stop("Please, provide a file name")
  }


  if(!file.exists(dname)){
    dir.create(dname)
  }

  report <- diann_load(data)  #load your report file
  brut <- report

  if(header.id == "Protein.Group" | header.id == "Genes"){
    if(only_proteotypic){
      report <- report[which(report[["Proteotypic"]] != 0), ]
    }
    ### Filtering on q values
    report <- report %>% dplyr::filter(Q.Value <= qv & PG.Q.Value <= pg.qv &
                                         Protein.Q.Value <= p.qv & GG.Q.Value <= gg.qv)

    ### iq process

    #preprocess the data in order to use MaxLFQ (see documentation from iq : https://cran.r-project.org/web/packages/iq/index.html  --> see the vignettes)
    iq_report <- iq::preprocess(report,intensity_col = quantity.id,
                                primary_id = header.id,
                                sample_id  = sample.id,
                                secondary_id = secondary.id,
                                median_normalization = FALSE,
                                pdf_out = file.path(dname, "qc-plots.pdf")
                                )
  }
  else{
    iq_report <- diann_matrix(report, id.header = header.id,
                              proteotypic.only = only_proteotypic,
                              q = qv,
                              protein.q = p.qv,
                              pg.q = pg.qv,
                              gg.q = gg.qv,
                              method = "max")
  }


  if(get_Top3 & (header.id == "Protein.Group" | header.id == "Genes")){
    Top3 <- iq::preprocess(report,
                           intensity_col = "Precursor.Quantity",
                           primary_id = header.id,
                           sample_id  = sample.id,
                           secondary_id = secondary.id,
                           median_normalization = FALSE,
                           pdf_out = NULL)
    Top3 <- Top3 %>% dplyr::group_by(protein_list) %>%
      tidyr::spread(sample_list, quant)
    Top3$id <- NULL
    top3_f <- function(x){
      if(sum(!is.na(x)) < 3){
        x <- NA
      }
      else{
        x <- x[order(x, decreasing = TRUE)][1:3]
        x <- mean(x)
      }
    }
    Top3 <- Top3 %>% dplyr::group_by(protein_list) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), top3_f))

    Top3 <- as.data.frame(Top3)
    rownames(Top3) <- Top3$protein_list
    Top3$protein_list <- NULL
    colnames(Top3) <- paste0("Top3_", colnames(Top3))
    Top3 <- Top3[order(rownames(Top3)),]
  }

  if(get_pep & (header.id == "Protein.Group" | header.id == "Genes")){
    ### Get peptides counts from preprocess
    pc <- iq_report %>% dplyr::group_by(protein_list, sample_list) %>%
      dplyr::mutate("countpep" = length(unique(id)))
    pc <- unique(pc[,c("protein_list", "sample_list", "countpep")])
    pc <- tidyr::spread(pc, sample_list, countpep)
    pc[is.na(pc)] <- 0
    pc <- as.data.frame(pc)
    rownames(pc) <- pc$protein_list
    pc$protein_list <- NULL
    pc <- pc[order(rownames(pc)),]
    colnames(pc) <- paste0("pep_count_", colnames(pc))
    pc$peptides_counts_all <- unname(apply(pc, 1, max))
    pc <- pc[,c(ncol(pc), 1:(ncol(pc)-1))]
  }

  if(header.id == "Protein.Group" | header.id == "Genes"){
    ### Get the protein group
    ## Perform fast MaxLFQ function from iq
    iq_report <- iq::fast_MaxLFQ(iq_report) #return a list, check it, call element with '$'

    iq_report <- iq_report$estimate #estimate is the dataset
    #other element is annotation, depends on your data


    # quick handle of data  --> add column protein group and put it at the beginning, remove rownames
    iq_report <- as.data.frame(iq_report)
    iq_report <- iq_report[order(rownames(iq_report)),]

    if(get_Top3){
      iq_report <- cbind(iq_report, Top3)
    }
    ## Add peptide counts
    if(get_pep & only_pepall){ #only keep peptide counts all
      iq_report$peptides_counts_all <- pc$peptides_counts_all
    }
    else if(get_pep){
      iq_report <- cbind(iq_report, pc)
    }
  }

  ### Add gene names and other informations; reshape data frame
  nc <- ncol(iq_report)
  iq_report[[header.id]] <- rownames(iq_report)
  rownames(iq_report) <- 1:nrow(iq_report)

  report <- report[(report[[header.id]] %in% iq_report[[header.id]]),]
  report <- report[order(report[[header.id]]),]
  col_n <- colnames(report)
  for(i in id_to_add){
    if(i %in% col_n){
      iq_report[[i]] <- unique(report[,c(header.id, i)])[[i]]
    }
    else{
      message(paste(i, "is not in your colnames data. Check the id you want to add."))
    }
  }
  iq_report <- iq_report[,c((nc+1):ncol(iq_report), 1:nc)]  # reorder columns

  if(get_iBAQ & header.id == "Protein.Group"){
    if(!is.null(fasta)){
      d_seq <- getallseq(pr_id = iq_report$Protein.Group,
                         fasta_file = TRUE,
                         bank_name = fasta)
    }
    else{
      d_seq <- getallseq(pr_id = iq_report$Protein.Group,
                         spec = species)
    }
    n_cond <- length(unique(report$File.Name))
    n_info <- length((nc+1):ncol(iq_report))
    brut <- diann_matrix(brut, id.header = "Protein.Group",
                         quantity.header = "Precursor.Quantity",
                         proteotypic.only = only_proteotypic,
                         q = qv, protein.q = p.qv,
                         pg.q = pg.qv, gg.q = gg.qv,
                         method = "sum")
    brut$Genes <- NULL
    brut$Protein.Names <- NULL
    brut <- get_iBAQ(brut, proteinDB = d_seq,
                     id_name = "Protein.Group",
                     ecol = n_info:(n_cond+1),
                     peptideLength = peptide_length,
                     proteaseRegExp = getProtease("trypsin"),
                     log2_transformed = TRUE)
    brut <- brut[,-c(n_info:(n_cond+1))]
    iq_report <- dplyr::left_join(iq_report, brut, by = "Protein.Group")
  }

  ### SAVE YOUR FILE
  filename <- paste(header.id, format(Sys.time(), "%H%M_%e-%m-%y"), sep = "_")
  filename <- file.path(dname, filename)

  format <- match.arg(format)
  if(format == "xlsx"){
    openxlsx::write.xlsx(iq_report, file = paste0(filename, ".xlsx"))
  }
  else if(format == "csv"){
    write.csv(iq_report, file = paste0(filename, ".csv"), row.names = FALSE)
  }
  else if(format == "txt"){
    write.table(iq_report, file = paste0(filename, ".txt"),  row.names = FALSE)
  }

  ### SAVE SOME PLOTS
  names(iq_report)[1] <- "id"

  ggplot2::ggsave(file.path(dname, "density.png"), DIAgui::densityDIA(iq_report, "none", TRUE, header.id),
                 height = 8, width = 16)

  ggplot2::ggsave(file.path(dname, "MDS.png"), DIAgui::MDS_DIA(iq_report, "none",  header.id))

  g <- ggplot2::ggplot(report, ggplot2::aes(iRT, RT, color = PG.Q.Value)) +
    ggplot2::geom_point() + ggplot2::facet_wrap(vars({{sample.id}})) +
    ggplot2::labs(title = "Report data filtered") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(file.path(dname, "RT.png"), g, height = 8, width = 16)

  message(paste("All your results have been saved ! Go check the directory", dname))
}
