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
#' @param get_pep logical; get peptide count ?
#' @param only_pepall logical; should only keep peptide counts all or also peptide counts for each fractions ?
#' @param format Format from your saved file (only xlsx, csv and txt are supported)
#'
#'
#' @export

report_process <- function(data, header.id = "Protein.Group", sample.id = "File.Name",
                           quantity.id = "Precursor.Normalised", secondary.id = "Precursor.Id",
                           id_to_add = c("Protein.Names", "First.Protein.Description", "Genes"),
                           qv = 0.01, pg.qv = 0.01, p.qv = 1, gg.qv = 1,
                           get_pep = TRUE, only_pepall = FALSE, format = c("xlsx", "csv", "txt")){

  if(length(stringr::str_subset(class(data), "^character$")) >= 1){
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



  if(get_pep){
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

  ### Get the protein group
  ## Perform fast MaxLFQ function from iq
  iq_report <- iq::fast_MaxLFQ(iq_report) #return a list, check it, call element with '$'

  iq_report <- iq_report$estimate #estimate is the dataset
  #other element is annotation, depends on your data


  # quick handle of data  --> add column protein group and put it at the beginning, remove rownames
  iq_report <- as.data.frame(iq_report)
  iq_report <- iq_report[order(rownames(iq_report)),]

  ## Add peptide counts
  if(get_pep & only_pepall){ #only keep peptide counts all
    iq_report$peptides_counts_all <- pc$peptides_counts_all
  }
  else if(get_pep){
    iq_report <- cbind(iq_report, pc)
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


  ### SAVE YOUR FILE
  filename <- paste(header.id, format(Sys.time(), "%H%M_%e-%m-%y"), sep = "_")
  filename <- file.path(dname, filename)

  if(format == "xlsx"){
    openxlsx::write.xlsx(iq_report, file = paste0(filename, ".xlsx"))
  }
  else if(format == "csv"){
    write.csv(iq_report, file = paste0(filename, ".csv"), row.names = FALSE)
  }
  else if(format == "txt"){
    write.table(iq_report, file = paste0(filename, ".txt"),  row.names = FALSE)
  }

  ### SAVE SOME PPLOTS
  names(iq_report)[1] <- "id"

  ggplot2::ggsave(file.path(dname, "density.png"), DIAgui::densityDIA(iq_report, "none", TRUE, header.id),
                 height = 8, width = 16)

  ggplot2::ggsave(file.path(dname, "MDS.png"), DIAgui::MDS_DIA(iq_report, "none",  header.id))

  g <- ggplot2::ggplot(report, ggplot2::aes(iRT, RT, color = PG.Q.Value)) +
    ggplot2::geom_point() + ggplot2::facet_wrap(vars({{sample.id}})) +
    ggplot2::labs(title = "Report data filtered") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(file.path(dname, "RT.png"), g, height = 8, width = 16)

  g <- report[, c(sample.id, "Proteotypic")]
  g$Proteotypic <- as.character(g$Proteotypic)
  g <- ggplot2::ggplot(g, ggplot2::aes(Proteotypic, fill = Proteotypic)) +
    ggplot2::geom_bar() +
    ggplot2::facet_wrap(vars({{sample.id}})) +
    ggplot2::labs(title = "Report data filtered") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(file.path(dname, "Proteotypic.png"), g,  height = 8, width = 10)

  message(paste("All your results have been saved ! Go check the directory", dname))
}


