#' imputationDIA
#'
#' Perform data imputation for your DIA data
#'
#' @param data Processed data from DIAnn (from iq processing or diann_matrix for example)
#' @param transformation Which transformation do you want to apply (log2 or none)
#' @param iteration Number of iterations for the imputation. See \code{\link[mice]{mice}} and \code{\link[imputeLCMD]{imputeLCMD}} for more details.
#' @param method Method used for the imputation. See \code{\link[mice]{mice}} and \code{\link[imputeLCMD]{imputeLCMD}} for more details.
#'
#' @details If your data contains iBAQ and or Top3 absolute quantification, then their column
#'   names should start with 'iBAQ_' and 'Top3_' respectively.
#'   Other numeric columns like the number of peptides used for quantification or the
#'   number of theoritical peptides will not be imputed, please keep the default name obtained
#'   with the app or the functions \code{diann_matrix}  or \code{diann_maxlfq}.
#'
#' @return The imputed data
#'
#' @export

imputationDIA <- function(data, transformation = c("none", "log2"),
                          iteration = 3, method = c("norm.predict", "zeros", "pmm", "midastouch",
                                                    "sample", "cart", "rf", "mean", "norm", "norm.nob",
                                                    "norm.boot", "lasso.norm", "lasso.select.norm",
                                                    "knn", "MinDet", "MinProb", "QRILC"
                                                    )
                          ){
  method <- match.arg(method)
  LCMD <- method %in% c("knn", "MinDet", "MinProb", "QRILC")
  if(LCMD){
    if(!("impute" %in% installed.packages())){
      message("Installing impute package from bioconductor")
      if(!("BiocManager" %in% installed.packages())){
        install.packages("BiocManager")
      }
      if(!("pcaMethods" %in% installed.packages())){
        BiocManager::install("pcaMethods")
      }
      BiocManager::install("impute")
    }
  }

  check_num <- mapply(class, data) == "numeric"
  if(all(!check_num)){
    stop("Your data doesn't contain any numeric value !")
  }

  transformation <- match.arg(transformation)
  tit2 <- ""
  if(transformation == "log2"){
    data[,which(check_num)] <- log2(data[,which(check_num)])
  }

  to_impute <- names(data)[which(check_num)]
  no_imp <- grep("peptides_counts_all|nbTrypticPeptides|^pep_count_", to_impute)
  if(length(no_imp)){
    to_impute <- to_impute[-no_imp]
  }

  if(method == "zeros"){
    data[,to_impute] <- apply(data[,to_impute], 2,
                              function(y){
                                if(any(is.na(y))){
                                  y[which(is.na(y))] <- 0
                                };
                                y
                              })
  }
  else{
    if(any(grepl("^Top3_", to_impute))){
      imp <- data[,to_impute[grep("^Top3_", to_impute)]]
      is_imp <- any(apply(imp, 2, function(x) any(is.na(x))))
      if(is_imp){
        org_names <- colnames(imp)
        imp <- data.frame(imp, check.names = TRUE)

        if(LCMD){
          if(method == "knn"){
            imp <- impute::impute.knn(as.matrix(imp))$data
          }
          else if(method == "QRILC"){
            imp <- imputeLCMD::impute.QRILC(as.matrix(imp))[[1]]
          }
          else if(method == "MinDet"){
            imp <- imputeLCMD::impute.MinDet(as.matrix(imp))
          }
          else if(method == "MinProb"){
            imp <- imputeLCMD::impute.MinProb(as.matrix(imp))
          }
        }
        else{
          imp <- mice::mice(as.matrix(imp), m = iteration, maxit = 5,
                            method = method, printFlag = FALSE, seed = 123)
          imp <- mice::complete(imp)
        }

        colnames(imp) <- org_names
        data[,to_impute[grep("^Top3_", to_impute)]] <- imp
      }

      to_impute <- to_impute[-grep("^Top3_", to_impute)]
    }
    if(any(grepl("^iBAQ_", to_impute))){
      imp <- data[,to_impute[grep("^iBAQ_", to_impute)]]
      is_imp <- any(apply(imp, 2, function(x) any(is.na(x))))
      if(is_imp){
        org_names <- colnames(imp)
        imp <- data.frame(imp, check.names = TRUE)

        if(LCMD){
          if(method == "knn"){
            imp <- impute::impute.knn(as.matrix(imp))$data
          }
          else if(method == "QRILC"){
            imp <- imputeLCMD::impute.QRILC(as.matrix(imp))[[1]]
          }
          else if(method == "MinDet"){
            imp <- imputeLCMD::impute.MinDet(as.matrix(imp))
          }
          else if(method == "MinProb"){
            imp <- imputeLCMD::impute.MinProb(as.matrix(imp))
          }
        }
        else{
          imp <- mice::mice(as.matrix(imp), m = iteration, maxit = 5,
                            method = method, printFlag = FALSE, seed = 123)
          imp <- mice::complete(imp)
        }

        colnames(imp) <- org_names
        data[,to_impute[grep("^iBAQ_", to_impute)]] <- imp
      }

      to_impute <- to_impute[-grep("^iBAQ_", to_impute)]
    }

    if(length(to_impute)){
      imp <- data[,to_impute]
      is_imp <- any(apply(imp, 2, function(x) any(is.na(x))))
      if(is_imp){
        org_names <- colnames(imp)
        imp <- data.frame(imp, check.names = TRUE)

        if(LCMD){
          if(method == "knn"){
            imp <- impute::impute.knn(as.matrix(imp))$data
          }
          else if(method == "QRILC"){
            imp <- imputeLCMD::impute.QRILC(as.matrix(imp))[[1]]
          }
          else if(method == "MinDet"){
            imp <- imputeLCMD::impute.MinDet(as.matrix(imp))
          }
          else if(method == "MinProb"){
            imp <- imputeLCMD::impute.MinProb(as.matrix(imp))
          }
        }
        else{
          imp <- mice::mice(as.matrix(imp), m = iteration, maxit = 5,
                            method = method, printFlag = FALSE, seed = 123)
          imp <- mice::complete(imp)
        }

        colnames(imp) <- org_names
        data[,to_impute] <- imp
      }
    }
  }

  data <- as.data.frame(data)
  return(data)
}

