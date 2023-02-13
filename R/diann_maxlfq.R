#' diann_maxlfq
#'
#' Quantify using the MaxLFQ algorithm https://doi.org/10.1074/mcp.M113.031591.
#' This function can be used to calculate protein quantities from peptide/precursor/fragment quantities or, for example,
#' to calculate peptide quantities from precursor/fragment quantities.
#' Totally based on function from diann-rpackage from Vadim Demichev
#'
#' @param x data, output from diann_load
#' @param sample.header Sample id column name
#' @param group.header Colulmn name corresponding to the group Id, e.g. protein Id
#' @param id.header Id column name (protein, gene, ...)
#' @param quantity.header Quantity column name
#' @param margin quantities below exp(margin) might be treated as NA
#' @param count_pep logical; get peptide count?
#' @param only_countsall logical; should only keep peptide counts all or also peptide counts for each fractions?
#' @param Top3 logical; get Top3 absolute quantification
#'
#' @return A dataframe containing the quantities from the id you selected
#'
#' @export

diann_maxlfq <- function (x, sample.header = "File.Name", group.header = "Protein.Names",
                          id.header = "Precursor.Id", quantity.header = "Precursor.Normalised",
                          margin = -10, count_pep = TRUE, only_countsall = FALSE, Top3 = FALSE){
  df <- data.table::as.data.table(x)
  df <- unique(df[which(df[[group.header]] != ""), c(sample.header,
                                                     group.header, id.header, quantity.header), with = FALSE])
  if(Top3){
    tp3 <- x[, c(id.header, "Precursor.Quantity")]
    df <- dplyr::left_join(df, tp3, by = id.header)
    df[["Precursor.Quantity"]] <- as.numeric(df[["Precursor.Quantity"]])
  }
  df[[sample.header]] <- as.character(df[[sample.header]])
  df[[group.header]] <- as.character(df[[group.header]])
  df[[id.header]] <- as.character(df[[id.header]])
  df[[quantity.header]] <- as.numeric(df[[quantity.header]])
  if (any(df[[quantity.header]] < 0, na.rm = T))
    stop("Only non-negative quantities accepted")
  is_duplicated = any(duplicated(paste0(df[[sample.header]],
                                        ":", df[[group.header]], ":", df[[id.header]])))
  if (is_duplicated)
    warning("Multiple quantities per id: the maximum of these will be calculated")
  if (margin > -2) {
    margin <- -2
    warning("margin reset to -2.0")
  }
  if (margin < -20) {
    margin <- -20
    warning("margin reset to -20.0")
  }
  df[[quantity.header]][which(df[[quantity.header]] == 0)] <- NA
  df[[quantity.header]] <- log(df[[quantity.header]])
  df[[quantity.header]][which(df[[quantity.header]] <= margin)] <- NA
  df <- df[!is.na(df[[quantity.header]]), ]
  if(Top3){
    df[["Precursor.Quantity"]] <- log(df[["Precursor.Quantity"]])
  }
  proteins <- unique(df[[group.header]])
  m <- length(proteins)
  samples <- unique(df[[sample.header]])
  n <- length(samples)
  result <- matrix(NA, nrow = m, ncol = n)
  rownames(result) = proteins
  colnames(result) = samples
  nbPep <- list()
  all_iden <- list()
  all_piv <- list()
  if(Top3){
    all_iden_brut <- list()
    all_piv_brut <- list()
  }
  for (i in 1:length(proteins)) {
    if (is_duplicated) {
      piv <- cast_aggregate(df[which(df[[group.header]] ==
                                       proteins[i]), ], sample.header, id.header, quantity.header)
      if(Top3){
        piv_brut <- cast_aggregate(df[which(df[[group.header]] ==
                                              proteins[i]), ], sample.header, id.header, "Precursor.Quantity")
      }
    }
    else {
      piv <- cast(df[which(df[[group.header]] == proteins[i]),
      ], sample.header, id.header, quantity.header)
      if(Top3){
        piv_brut <- cast(df[which(df[[group.header]] == proteins[i]),
        ], sample.header, id.header, "Precursor.Quantity")
      }
    }
    if (nrow(piv) == 1 | ncol(piv) == 1) {
      res <- col_max(as.vector(piv), nrow(piv), ncol(piv))
      all_piv[[proteins[i]]] <- piv
      nbPep[[proteins[i]]] <- nrow(piv)
      if(Top3){
        all_piv_brut[[proteins[i]]] <- piv_brut
      }
    }
    else {
      piv[is.na(piv)] <- -1e+06
      ref = col_max(as.vector(piv), nrow(piv), ncol(piv))
      columns <- which(ref > margin)
      identified <- piv[, columns]
      nbPep[[proteins[i]]] <- nrow(identified)
      all_iden[[proteins[i]]] <- identified
      if(Top3){
        piv_brut[is.na(piv_brut)] <- -1e+06
        ref_brut = col_max(as.vector(piv_brut), nrow(piv_brut), ncol(piv_brut))
        columns_brut <- which(ref_brut > margin)
        identified_brut <- piv_brut[, columns_brut]
        all_iden_brut[[proteins[i]]] <- identified_brut
      }
      if (ncol(identified) >= 2) {
        res <- maxlfq_solve(as.vector(identified), nrow(identified),
                            ncol(identified), margin * 1.001)
        ref[columns] <- res
      }
      else res <- ref
      res[which(res <= margin)] <- NA
    }
    result[i, match(colnames(piv), samples)] <- res
  }
  result <- exp(result)
  result <- as.data.frame(result)
  result <- result[order(rownames(result)),]

  if(Top3){
    all_iden_brut <- lapply(all_iden_brut, function(x){
      x[which(x <= margin)] <- NA;
      x
    })

    top3_iden <- lapply(all_iden_brut, function(x){
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
    p = names(top3_iden)
    top3_iden <- Reduce(rbind, top3_iden)
    rownames(top3_iden) <- p

    top3_piv <- lapply(all_piv_brut, function(x){
      x <- apply(x, 2, function(y){
        y <- NA;
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
    p = names(top3_piv)
    top3_piv <- Reduce(rbind, top3_piv)
    rownames(top3_piv) <- p

    sumup <- rbind(top3_iden, top3_piv)
    sumup <- sumup[order(rownames(sumup)),]
    colnames(sumup) <- paste0("Top3_", colnames(sumup))
    sumup <- exp(sumup)

    result <- cbind(result, sumup)
  }

  if(count_pep){
    nbPep <- nbPep[order(names(nbPep))]
    nbPep <- unname(unlist(nbPep))
    result$peptide_counts_all <- nbPep

    if(!only_countsall){
      all_iden <- lapply(all_iden, function(x){
          x[which(x <= margin)] <- NA;
          x
        })
      all_iden <- lapply(all_iden, function(x){
        apply(x, 2, function(k){
          sum(!is.na(k))
        })
      })
      all_iden <- lapply(all_iden, function(x){x <- as.data.frame(t(x));x})
      all_iden <- lapply(all_iden, function(x){
        n <- samples[!(samples %in% colnames(x))]
        if(!purrr::is_empty(n)){
          for(i in n){
            x[[i]] <- 0
          }
        };
        x
      })
      p = names(all_iden)
      all_iden <- Reduce(rbind, all_iden)
      rownames(all_iden) <- p

      all_piv <- lapply(all_piv, function(x){
        apply(x, 2, function(k){
          sum(!is.na(k))
        })
      })
      all_piv <- lapply(all_piv, function(x){x <- as.data.frame(t(x));x})
      all_piv <- lapply(all_piv, function(x){
        n <- samples[!(samples %in% colnames(x))]
        if(!purrr::is_empty(n)){
          for(i in n){
            x[[i]] <- 0
          }
        };
        x
      })
      p = names(all_piv)
      all_piv <- Reduce(rbind, all_piv)
      rownames(all_piv) <- p

      sumup <- rbind(all_iden, all_piv)
      sumup <- sumup[order(rownames(sumup)),]
      colnames(sumup) <- paste0("pep_count_", colnames(sumup))

      result <- cbind(result, sumup)
    }
  }

  return(result)
}


### interns functions from diann r package from V. Demichev

cast_aggregate <- function (df, sample.header, id.header, quantity.header){
  x <- data.table::melt.data.table(df, id.vars = c(sample.header, id.header),
                                   measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- data.table::dcast.data.table(x, as.formula(paste0(id.header,
                                                           "~", sample.header)),
                                      value.var = "value", fun.aggregate = function(x) max(x,
                                                                                           na.rm = TRUE))
  piv[[1]] <- NULL
  piv = as.matrix(piv)
  piv[is.infinite(piv)] <- NA
  piv
}
cast <- function (df, sample.header, id.header, quantity.header){
  x <- data.table::melt.data.table(df, id.vars = c(sample.header, id.header),
                                   measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- data.table::dcast.data.table(x, as.formula(paste0(id.header,
                                                           "~", sample.header)),
                                      value.var = "value")
  piv[[1]] <- NULL
  as.matrix(piv)
}

