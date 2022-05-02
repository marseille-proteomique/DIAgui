#' MDS_DIA
#'
#' Print MDS of processed data from DIAnn
#'
#' @param data Processed data from DIAnn (from iq processing or diann_matrix for example)
#' @param transformation Which transformation do you want to apply (log2 or none)
#' @param tit Title of your plot
#'
#' @return ggplot2 MDS graph
#'
#' @export

MDS_DIA <- function(data, transformation = c("log2", "none"), tit = ""){
  m <- as.data.frame(data)
  cl <- lapply(m, class)
  cl <- cl == "numeric"
  m <- m[,cl]
  tit2 <- ""
  if(transformation == "log2"){
    if(stringr::str_length(tit) == 0){
      tit2 <- ", Log2 transformed"
    }
    m <- log2(m)
  }
  else if(transformation != "none"){
    transformation <- "none"
  }
  if(stringr::str_length(tit) == 0){
    tit <- paste0(deparse(substitute(data)), tit2)
  }

  bad <- rowSums(is.finite(as.matrix(m))) < ncol(m)
  if(nrow(as.matrix(m)[!bad, , drop = FALSE]) == 0){
    p <- ggplot2::ggplot(data.frame(x = c(0,1), y = c(0,1)), ggplot2::aes(x,y, label = "s")) +
      ggplot2::geom_text(x=0.5, y=0.5, label = "No rows with only finite value", size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank())
  }
  else{
    m <- limma::plotMDS(m, plot = FALSE)
    m <- data.frame(x = m$x, y = m$y, row.names = rownames(m$distance.matrix.squared))
    m$name <- rownames(m)


    p <- ggplot2::ggplot(m, ggplot2::aes(x, y, color = name)) +
      ggplot2::geom_point(size = 1.8) +
      ggrepel::geom_label_repel(ggplot2::aes(label = name), max.overlaps = Inf) +
      ggplot2::labs(title = tit,
           subtitle = "MDS plot",
           x = "Dimension 1", y = "Dimension 2") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "none") +
      ggplot2::xlim(-max(abs(m$x)), max(abs(m$x))) +
      ggplot2::ylim(-max(abs(m$y)), max(abs(m$y)))
  }

  return(p)
}
