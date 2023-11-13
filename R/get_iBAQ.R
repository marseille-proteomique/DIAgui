#' get_iBAQ
#'
#' Function to calculate iBAQ quantities from a dataset.
#'
#' @details
#' This function use three function from \pkg{SafeQuant} package from github eahrne/SafeQuant. Those three function
#' are simple and small, so to facilitate dependencies and lower the size of this package, I chose to copy-paste those three.
#'
#' Erik Ahrne (2021). SafeQuant: A Toolbox for the Analysis of Proteomics Data. R package version 2.3.4.
#'
#'
#' @param dataset A data frame that you want to get the iBAQ quantities from.
#' @param proteinDB A list that contains for each protein the theoretical peptide sequence.
#'                  It doesn't have to be the same length as the number of proteins from the dataset;
#'                  if it's the case, value from proteins that are not in this list will be set to NA.
#' @param id_name The name of the column that contains the protein IDs.
#' @param ecol The indexes of the intensity columns.
#' @param peptideLength A vector that contains the minimum and maximum peptide length.
#' @param nbMiscleavages The number of miss cleavages.
#' @param proteaseRegExp The enzyme you used to digest proteins.
#' @param log2_transformed Logical to tell if you want to log2 transformed the data.
#' @param keep_original Logical to tell if you want to keep the original intensity values.
#'
#' @return The dataset with the iBAQ quantities and number of theoretical peptides.
#'
#' @examples
#' library(DIAgui)
#' data <- small_report %>% dplyr::filter(Q.Value <= 0.01 & PG.Q.Value <= 0.01)
#' data <- diann_maxlfq(data,
#'                   group.header="Protein.Group",
#'                   id.header = "Precursor.Id",
#'                   quantity.header = "Precursor.Normalised"
#'                   )
#' nc <- ncol(data)
#' data$Protein.Group <- rownames(data)
#' rownames(data) <- 1:nrow(data)
#' data <- data[order(data$Protein.Group),]
#' data <- data[,c((nc+1):ncol(data), 1:nc)]
#' sequence_list <- getallseq(pr_id = data$Protein.Group)
#' iBAQ <- get_iBAQ(data, sequence_list, ecol = 2:4)
#'
#' @export

get_iBAQ <- function (dataset, proteinDB = list(),
                      id_name = "Protein.Group",
                      ecol = 5:7,
                      peptideLength = c(5, 36), nbMiscleavages = 0,
                      proteaseRegExp = getProtease("trypsin"),
                      log2_transformed = TRUE,
                      keep_original = TRUE){
  check_numeric <- names(dataset[,ecol])[!unlist(lapply(lapply(dataset[,ecol], class), function(x) x == "numeric"))]
  if(length(check_numeric)){
    stop(paste("The column(s)", paste(check_numeric, collapse = ", "), "is/are not of class 'numeric'. Please check your values in ecol."))
  }
  nbPeptides <- vector(length = nrow(dataset))
  ID <- dataset[[id_name]]
  if(is.null(ID)){
    stop(paste("The column name", id_name, "didn't return any value. Please check it."))
  }
  idx <- 0
  for (i in ID) {
    idx <- idx + 1
    ac <- as.character(i)
    nbPep <- NA
    if (length(proteinDB[[ac]])) {
      if(!is.na(proteinDB[[ac]])){
        nbPep <- getNbDetectablePeptides(getPeptides(proteinDB[[ac]],
                                                     proteaseRegExp = proteaseRegExp,
                                                     nbMiscleavages = nbMiscleavages),
                                         peptideLength = peptideLength)
      }
      else {
        warning("WARN: ", ac, " NOT FOUND IN PROTEIN DATABASE")
      }
    }
    else {
      warning("WARN: ", ac, " NOT FOUND IN PROTEIN DATABASE")
    }
    nbPeptides[idx] <- nbPep
  }
  dataset_IBAQ <- dataset[,ecol]
  dataset_IBAQ <- dataset_IBAQ/nbPeptides
  if(log2_transformed){
    dataset_IBAQ <- log2(dataset_IBAQ)
    dataset[,ecol] <- log2(dataset[,ecol])
  }
  names(dataset_IBAQ) <- paste0("iBAQ_", names(dataset_IBAQ))
  dataset_IBAQ$nbTrypticPeptides <- nbPeptides

  last_ecol = ecol[length(ecol)]
  dataset_IBAQ <- as.data.frame(cbind(dataset[,1:last_ecol], dataset_IBAQ))
  if(last_ecol < ncol(dataset)){
    dataset_IBAQ <- as.data.frame(cbind(dataset_IBAQ, dataset[,(last_ecol + 1):ncol(dataset)]))
  }
  if(!keep_original){
    dataset_IBAQ <- dataset_IBAQ[,-ecol]
  }
  return(dataset_IBAQ)
}


#' getProtease
#'
#' Returns a regular corresponding to the given enzyme in order to digest a
#' proteic sequence.
#'
#' @param protease name of the enzyme; currently "trypsin" or "lys-c"
#'
#' @return a regular expression
#' @export
#'
#' @note function from SafeQuant R package by Erik Ahrne (2021). SafeQuant: A
#'   Toolbox for the Analysis of Proteomics Data. R package version 2.3.4.
#'
#' @examples
#' getProtease()

getProtease <- function (protease = "trypsin") {
  if (protease == "trypsin") {
    return("[KR](?!P)")
  }
  else if (protease == "lys-c") {
    return("K(?!P)")
  }
  else {
    stop("Unknown Protease:", protease)
  }
}


getNbDetectablePeptides <- function (peptides, peptideLength = c(5, 36)) {
  okLength <- (nchar(peptides) >= peptideLength[1]) & (nchar(peptides) <=
                                                         peptideLength[2])
  return(sum(okLength))
}


getPeptides <- function (proteinSeq, proteaseRegExp = getProtease("trypsin"),
                         nbMiscleavages = 0, minLength = 0, maxLength = Inf) {
  allAA <- as.vector(unlist(strsplit(proteinSeq, "")))
  fcPeptides <- strsplit(proteinSeq, proteaseRegExp, perl = TRUE)[[1]]
  matchPos <- gregexpr(proteaseRegExp, proteinSeq, perl = TRUE)[[1]]
  separator <- allAA[matchPos]
  if (length(separator) < length(fcPeptides))
    separator <- c(separator, "")
  fcPeptides <- paste(fcPeptides, separator, sep = "")
  fcPeptides <- fcPeptides[nchar(fcPeptides) > 0]
  if (nbMiscleavages == 0) {
    allPeptides = fcPeptides
  }
  else {
    allPeptides = vector()
    k = 1
    for (i in 1:length(fcPeptides)) {
      for (j in i:(i + nbMiscleavages)) {
        if (j <= length(fcPeptides)) {
          pept <- paste(fcPeptides[i:j], collapse = "")
          kNew = k + length(pept) - 1
          allPeptides[k:kNew] = pept
          k = kNew + 1
        }
      }
    }
  }
  nbPeptides = nchar(allPeptides)
  return(allPeptides[nbPeptides >= minLength & nbPeptides <= maxLength] %>% unique)
}
