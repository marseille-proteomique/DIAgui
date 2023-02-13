#' densityDIA
#'
#' Print density from each fraction from processed data from DIAnn
#'
#' @param data Processed data from DIAnn (from iq processing or diann_matrix for example)
#' @param transformation Which transformation do you want to apply (log2 or none)
#' @param area logical; print the area under curve
#' @param tit Title of your plot
#' @param data_type The type of data you want to visualize; either 'intensity', 'Top3', 'iBAQ' or 'all'.
#'
#' @return ggplot2 density graph
#'
#' @export

densityDIA <- function(data, transformation = c("log2", "none"), area = FALSE,
                       tit = "", data_type = c("intensity", "Top3", "iBAQ", "all")){
  d <- as.data.frame(data)
  cl <- lapply(d, class)
  cl <- cl == "numeric"
  if(sum(cl) == 0){
    message("No numeric data !")
    return(NULL)
  }
  else if(sum(cl) == 1){
    frac <- names(d)[cl]
    rw <- rownames(d)
    d <- as.data.frame(d[,cl])
    names(d) <- frac
    rownames(d) <- rw
  }
  else{
    d <- d[,cl]
    to_rm <- stringr::str_which(colnames(d), "nbTrypticPeptides|peptides_counts_all|^pep_count_")
    if(length(to_rm) > 0){
      if(length(to_rm) == ncol(d)){
        message("No numeric data !")
        return(NULL)
      }
      else{
        d <- d[,-to_rm]
      }
    }
    data_type <- match.arg(data_type)
    if(data_type == "Top3"){
      d <- d[,stringr::str_which(colnames(d), "^Top3_")]
    }
    else if(data_type == "iBAQ"){
      d <- d[,stringr::str_which(colnames(d), "^iBAQ_")]
    }
    else if(data_type == "intensity"){
      idx <- stringr::str_which(colnames(d), "^iBAQ_|^Top3_")
      if(length(idx) > 0){
        d <- d[,-idx]
      }
    }
    else if(data_type == "all"){
      message("You chose to keep all numeric data, they may differ completly.")
    }
    else{
      stop("data_type can only be 'intensity', 'Top3', 'iBAQ' or 'all' .")
    }
  }
  tit2 <- ""
  transformation <- match.arg(transformation)
  if(transformation == "log2"){
    if(stringr::str_length(tit) == 0){
      tit2 <- ", Log2 transformed"
    }
    d <- log2(d)
  }
  else if(transformation != "none"){
    transformation <- "none"
  }
  d <- tidyr::gather(d, fraction, intensity)
  if(stringr::str_length(tit) == 0){
    tit <- paste0("Density plot", tit2)
  }

  if(area){
    g <- ggplot2::ggplot(d, ggplot2::aes(intensity, color = fraction, fill = fraction)) + ggplot2::geom_density(alpha = 0.1) +
      ggplot2::labs(title = tit,
           y = "Density",
           x = "Intensity",
           color = "Fraction",
           fill = "Fraction") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  else{
    g <- ggplot2::ggplot(d, ggplot2::aes(intensity, color = fraction)) + ggplot2::geom_density() +
      ggplot2::labs(title = paste0("Density plot", tit),
           y = "Density",
           x = "Intensity",
           color = "Fraction") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  return(g)
}
