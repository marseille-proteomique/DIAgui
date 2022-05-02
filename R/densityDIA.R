#' densityDIA
#'
#' Print density from each fraction from processed data from DIAnn
#'
#' @param data Processed data from DIAnn (from iq processing or diann_matrix for example)
#' @param transformation Which transformation do you want to apply (log2 or none)
#' @param area logical; print the area under curve
#' @param tit Title of your plot
#'
#' @return ggplot2 density graph
#'
#' @export

densityDIA <- function(data, transformation = c("log2", "none"), area = FALSE,
                       tit = ""){
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
  }
  tit2 <- ""
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
