#' corrDIA
#'
#' Plot correlation plot from each fraction from processed data from DIAnn
#'
#' @param data Processed data from DIAnn (from iq processing or diann_matrix for example)
#' @param transformation Which transformation do you want to apply (log2 or none)
#' @param tit Title of your plot
#' @param data_type The type of data you want to visualize; either 'intensity', 'Top3', 'iBAQ' or 'all'.
#' @param gradient_color Three colors for the printed gradient color on the correlation plot
#' @param plot_pairs Plot pairwise correlation plot. If FALSE, will plot correlation matrix plot.
#'
#' @return ggplot2 correlation plot
#'
#' @export

corrDIA <- function(data, transformation = c("none", "log2"),
                    tit = "", data_type = c("intensity", "Top3", "iBAQ", "all"),
                    gradient_color = c("#09009D", "#ffffff", "#BE0010"),
                    plot_pairs = FALSE){
  if(length(gradient_color) != 3){
    stop("Please provide exactly three colors for the gradient color.")
  }
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

  if(stringr::str_length(tit) == 0){
    tit <- paste0("Correlation plot", tit2)
  }

  if(plot_pairs){
    gn <- unname(unlist(mapply(function(x, y){
      n <- paste(x, "vs", names(d)[y:ncol(d)]);
      n
    },
    names(d), 1:ncol(d), SIMPLIFY = FALSE)))

    g <- lapply(gn,
                function(x){
                  x <- strsplit(x, " vs ")[[1]]
                  if(x[1] == x[2]){
                    g <- ggplot2::ggplot(d) +
                      ggplot2::geom_density(ggplot2::aes(.data[[x[1]]])) +
                      ggplot2::theme_bw()

                    if(x[1] != names(d)[ncol(d)]){
                      g <- g +
                        ggplot2::theme(axis.title.x = ggplot2::element_blank())
                    }
                    if(x[1] != names(d)[1]){
                      g <- g +
                        ggplot2::theme(axis.title.y = ggplot2::element_blank())
                    }
                    else{
                      g <- g + ggplot2::labs(y = names(d)[1])
                    }
                  }
                  else{
                    cf <- coef(lm(d[[x[2]]] ~ d[[x[1]]]))
                    g <- ggplot2::ggplot(d) +
                      ggplot2::geom_point(ggplot2::aes(.data[[x[1]]], .data[[x[2]]])) +
                      ggplot2::geom_abline(intercept = cf[[1]],
                                           slope = cf[[2]],
                                           color = "red") +
                      ggplot2::theme_bw()

                    if(x[2] != names(d)[ncol(d)]){
                      g <- g +
                        ggplot2::theme(axis.title.x = ggplot2::element_blank())
                    }
                    if(x[1] != names(d)[1]){
                      g <- g +
                        ggplot2::theme(axis.title.y = ggplot2::element_blank())
                    }
                  };
                  g
                })
    names(g) <- gn

    gn <- expand.grid(names(d), names(d))
    gn <- paste(gn[[1]], "vs", gn[[2]])
    g <- lapply(gn, function(x) g[[x]])

    g <- cowplot::plot_grid(plotlist = g, nrow = ncol(d), ncol = ncol(d),
                            labels = tit, hjust = -2)
  }
  else{
    corr <- d %>%
      cor(use = "pairwise.complete.obs") %>%
      as.data.frame() %>%
      tibble::rownames_to_column("with") %>%
      tidyr::gather("fraction", "corr", -with)

    g <- ggplot2::ggplot(corr, ggplot2::aes(with, fraction, fill = corr)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = round(corr, 3)),
                         color = "black", size = 6) +
      ggplot2::scale_fill_gradientn(colors = gradient_color,
                                    values = c(0,0.5,1),
                                    breaks = c(-1,0,1),
                                    limits = c(-1,1)) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::labs(title = tit) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     axis.text = ggplot2::element_text(size = 12),
                     axis.title = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.title = ggplot2::element_blank())
  }

  return(g)
}


