#' getallseq
#'
#' Get theoretical peptide sequence from proteins (from data bank or fasta files), thanks to \pkg{seqinr} package.
#'
#' @details
#' Charif D, Lobry J (2007). “SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to
#' biological sequences retrieval and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.),
#' Structural approaches to sequence evolution: Molecules, networks, populations, series Biological and Medical Physics,
#' Biomedical Engineering, 207-232. Springer Verlag, New York. ISBN : 978-3-540-35305-8.
#'
#' If the protein is grouped, it will only take the first protein.
#'
#' @param spec The species name corresponding to the proteins from \code{pr_id}.
#' @param pr_id The proteins that you want to get the theoretical peptide sequence from.
#' @param bank_name The name of the data bank you want to use to get theoretical peptide.
#'                  If you want to use you're own fasta files, you can put the path (or several in a double)
#'                  to your fasta files and set fasta_file to TRUE.
#' @param fasta_file Logical to tell if you use fasta file in \code{bank_name}.
#' @param verbose Logical to tell if you want to print advancement.
#'
#' @return A list containing the peptide sequence from each protein.
#'
#' @export

getallseq <- function(spec = "Saccharomyces cerevisiae",
                      pr_id = "P09938",
                      bank_name = "swissprot",
                      fasta_file = FALSE,
                      verbose = TRUE){
  seq_res <- list()
  n = length(pr_id)
  ni = 1

  if(fasta_file){
    n_bk <- length(bank_name)
    if(n_bk == 1){
      fasta_base <- seqinr::read.fasta(file = bank_name, as.string = TRUE, forceDNAtolower = FALSE)
    }
    else if(n_bk > 1){
      fasta_base <- list()
      for (i in 1:n_bk){
        fasta_base[[i]] <- seqinr::read.fasta(file = bank_name[i], as.string = TRUE, forceDNAtolower = FALSE)
      }
      fasta_base <- do.call(c, fasta_base)
    }
    names(fasta_base) <- unlist(lapply(stringr::str_split(names(fasta_base), "\\|"), function(x) x[2]))

    for(i in pr_id){
      pr <- i
      if(stringr::str_detect(pr, "\\;")){
        pr <- stringr::str_split(pr, "\\;")[[1]][1]  #take only the first protein if grouped
        seq_res[[i]] <- fasta_base[[pr]][1]
      }
      else{
        seq_res[[i]] <- fasta_base[[pr]][1]
      }

      if(verbose){
        message(paste0(pr, "; ", ni, "/", n))
      }
      ni = ni + 1
    }
  }
  else{
    mybank <- seqinr::choosebank(bank_name)

    for(i in pr_id){
      pr <- i
      if(stringr::str_detect(pr, "\\;")){
        pr <- stringr::str_split(pr, "\\;")[[1]] #take only the first protein if grouped
        get_quer <- TRUE
        k <- 1
        while(get_quer){
          quer <- paste0("SP=", spec, " AND AC=", pr[k])
          search <- try(seqinr::query("search", quer), silent = TRUE)
          get_quer <- inherits(search, "try-error") & k < length(pr)
          k <- k+1
        }
        if(!inherits(search, "try-error")){
          info <- stringr::str_split(quer, "=")[[1]]
          info <- info[length(info)]
          message(paste("Only", info, "is taken from the group", i))
        }
      }
      else{
        quer <- paste0("SP=", spec, " AND AC=", pr[1])
        search <- try(seqinr::query("search", quer), silent = TRUE)
      }

      if(!inherits(search, "try-error")){
        if(length(search$req)){  # query didn't throw error but returned nothing; for example if protein is not from same species
          if(verbose){
            message(paste0(quer, "; ", ni, "/", n))
          }
          ni = ni + 1
          sequ <- seqinr::getSequence(search$req[[1]])
          sequ <- paste(sequ, collapse = "")
          seq_res[[i]] <- sequ
          if (sequ == "NA"){
            mybank <- seqinr::choosebank(bank_name)
          }
          if(stringr::str_detect(seqinr::getName(search$req[[1]]), "^\\d{1}")){
            message("Reopening bank. Last protein name starts with a digit.")
            mybank <- seqinr::choosebank(bank_name)
          }
        }
        else{
          message(paste("Couldn't find sequence in bank", bank_name, "for protein", i))
          seq_res[[i]] <- NA
        }
      }
      else{
        message(paste("Couldn't find sequence in bank", bank_name, "for protein", i))
        seq_res[[i]] <- NA
      }
    }
    seqinr::closebank()
  }
  return(seq_res)
}
