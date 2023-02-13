#' diann_matrix
#'
#' Generate matrix with quantities (zero will be replaced by NA values) and
#' eventually peptides count.
#'
#' @param x data, output from diann_load
#' @param id.header Id column name (protein, gene, ...)
#' @param quantity.header Quantity column name
#' @param proteotypic.only logical; Only proteotypic peptides and the
#'   respective proteins should be considered?
#' @param q Precursor q-value threshold
#' @param protein.q Uniquely identified protein q-value threshold
#' @param pg.q Protein group q-value threshold
#' @param gg.q Gene group q-value threshold
#' @param get_pep logical; get peptide count?
#' @param only_pepall logical; should only keep peptide counts all or also
#'   peptide counts for each fractions?
#' @param margin quantities below exp(margin) might be treated as NA
#' @param Top3 logical; get Top3 absolute quantification
#' @param method When for one identifier there are several values, take either
#'   maximum of these or sum them all. When using 'sum', any additional
#'   information than the id.header should be ignored. 'max' reports for only
#'   one id, 'sum' get the sum of the intensities of all ids.
#'
#' @return A dataframe containing the quantities from the id you selected
#'
#' @export

diann_matrix <- function (x, id.header = "Precursor.Id", quantity.header = "Precursor.Normalised",
                          proteotypic.only = FALSE, q = 0.01, protein.q = 1, pg.q = 0.01,
                          gg.q = 1, get_pep = FALSE, only_pepall = FALSE, margin = -10, Top3 = FALSE,
                          method = c("max", "sum")){
  df <- data.table::as.data.table(x)
  if(proteotypic.only){
    df <- df[which(df[["Proteotypic"]] != 0), ]
  }
  dft <- df[which(df[[id.header]] != "" & df[["Q.Value"]] <= q & df[["Protein.Q.Value"]] <=
                    protein.q & df[["PG.Q.Value"]] <= pg.q & df[["GG.Q.Value"]] <=
                    gg.q),]

  info <- c()
  if(id.header == "Precursor.Id" | id.header == "Stripped.Sequence" | id.header == "Modified.Sequence"){
    info <- c("Stripped.Sequence", "Modified.Sequence")
  }
  info <- c(info, "Protein.Group", "Protein.Names", "Genes")
  info <- info[!(info %in% id.header)]
  dft$add_info <- apply(as.data.frame(dft)[,info], 1, function(x) paste(x, collapse = " "))

  df <- unique(dft[, c("File.Name", id.header, quantity.header, "add_info"), with = FALSE]) # should remove unique since it will be handled
  is_duplicated = any(duplicated(paste0(df[["File.Name"]],
                                        ":", df[[id.header]])))
  if (is_duplicated) {
    warning("Multiple quantities per id: the maximum of these will be calculated")
    out <- pivot_aggregate(df, "File.Name", id.header, quantity.header, method = method)
    out <- tidyr::separate(out, add_info, into = info, sep = " ")
  }
  else {
    out <- pivot(df, "File.Name", id.header, quantity.header)
    out <- tidyr::separate(out, add_info, into = info, sep = " ")
  }
  dft$add_info <- NULL
  if(Top3 & id.header == "Genes"){
    x <- unique(dft[which(dft[["Genes"]] != ""), c("File.Name",
                                                   "Genes", "Precursor.Id", "Precursor.Quantity"), with = FALSE])
    x[["File.Name"]] <- as.character(x[["File.Name"]])
    x[["Genes"]] <- as.character(x[["Genes"]])
    x[["Precursor.Id"]] <- as.character(x[["Precursor.Id"]])
    x[["Precursor.Quantity"]] <- as.numeric(x[["Precursor.Quantity"]])
    if (any(x[["Precursor.Quantity"]] < 0, na.rm = T))
      stop("Only non-negative quantities accepted")
    is_duplicated = any(duplicated(paste0(x[["File.Name"]],
                                          ":", x[["Genes"]], ":", x[["Precursor.Id"]])))
    if (is_duplicated)
      warning("Multiple quantities per id: the maximum of these will be calculated")
    x[["Precursor.Quantity"]][which(x[["Precursor.Quantity"]] == 0)] <- NA
    x[["Precursor.Quantity"]] <- log(x[["Precursor.Quantity"]])
    x[["Precursor.Quantity"]][which(x[["Precursor.Quantity"]] <= margin)] <- NA
    x <- x[!is.na(x[["Precursor.Quantity"]]), ]
    genes <- unique(x[["Genes"]])
    samples <- unique(x[["File.Name"]])
    top3_res <- list()
    for(i in genes){
      top3 <- x[which(x[["Genes"]] == i),]
      top3[["Genes"]] <- NULL
      top3 <- tidyr::spread(top3, File.Name, Precursor.Quantity)
      top3[["Precursor.Id"]] <- NULL
      top3_res[[i]] <- top3
    }
    top3_res <- lapply(top3_res, function(x){
      x <- apply(x, 2, function(y){
        if(sum(!is.na(y)) < 3){
          y <- NA
        }
        else{
          y <- y[order(y, decreasing = TRUE)][1:3]
          y <- mean(y)
        };
        y
      })
      x <- as.data.frame(t(x))
      n <- samples[!(samples %in% colnames(x))]
      if(!purrr::is_empty(n)){
        for(i in n){
          x[[i]] <- NA
        }
      };
      x
    })
    p = names(top3_res)
    top3_res <- Reduce(rbind, top3_res)
    colnames(top3_res) <- paste0("Top3_", colnames(top3_res))
    top3_res <- exp(top3_res)
    top3_res[[id.header]] <- p
    top3_res <- top3_res[order(top3_res[[id.header]]),]

    out <- merge(out, top3_res, by=id.header)
  }
  if(get_pep){
    x <- dft[,c("File.Name", id.header, "Genes.MaxLFQ.Unique", "Precursor.Id"), with = FALSE]
    pep <- as.data.frame(matrix(0, nrow = nrow(out), ncol = 2 + length(unique(x$File.Name))))
    colnames(pep) <- c(id.header, "peptides_counts_all", paste0("pep_count_", unique(x$File.Name)))
    for (i in unique(x$File.Name)){
      frac <- which(x$File.Name == i)
      a <- x[frac,]
      r = 1
      for (k in unique(df[[id.header]])){
        gene <- which(a[[id.header]] == k)
        b <- a[gene,]
        np <- length(unique(b[["Precursor.Id"]]))
        pep[r, c(id.header, paste0("pep_count_", i))] <- c(k, np)
        r = r + 1
      }
    }
    g <- pep[[id.header]]
    pep[[id.header]] <- NULL
    pep <- as.data.frame(apply(pep, 2, as.numeric))
    pep$peptides_counts_all <- apply(pep, 1, max)
    pep[[id.header]] <- g

    if(only_pepall){
      pep <- pep[,c(1, ncol(pep))]
    }

    out <- merge(out, pep, by=id.header)
  }

  return(out)
}


### interns functions from diann-rpackage from vdemichev

pivot_aggregate <- function (df, sample.header, id.header, quantity.header, method = c("max", "sum")){
  x <- data.table::melt.data.table(df, id.vars = c(sample.header, id.header, "add_info"),
                                   measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  if(method == "max"){
    piv <- x %>% dplyr::group_by(!!dplyr::sym(sample.header), !!dplyr::sym(id.header)) %>%
      dplyr::summarise("results" = max(value, na.rm = TRUE),
                       "add_info" = add_info[which(value == max(value, na.rm = TRUE))])
  }
  else if(method == "sum"){
    piv <- x %>% dplyr::group_by(!!dplyr::sym(sample.header), !!dplyr::sym(id.header)) %>%
      dplyr::summarise("results" = sum(value, na.rm = TRUE),
                       "add_info" = add_info[which(value == max(value, na.rm = TRUE))])
  }
  else{
    stop("method can only be max or sum")
  }
  piv <- tidyr::spread(piv, !!dplyr::sym(sample.header), results)
  piv <- piv[order(piv[[id.header]]),]

  return(piv)
}

pivot <- function (df, sample.header, id.header, quantity.header){
  x <- data.table::melt.data.table(df, id.vars = c(sample.header, id.header, "add_info"),
                                   measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  x$variable <- NULL
  piv <- tidyr::spread(x, !!dplyr::sym(sample.header), value)
  piv <- piv[order(piv[[id.header]]),]

  return(piv)
}



