#' validDIA
#'
#' Print plot showing proportion of valid values for each protein (or else) among the experiment
#'
#' @details If you have replicates or several conditions and you want to see the proportion of valid values
#'          for each replicate for example. Then, the id column in your data must be the first one;
#'          as in the output of the shiny app of DIAgui.
#'
#' @param data Processed data from DIAnn (from iq processing or diann_matrix for example)
#' @param transformation Which transformation do you want to apply (log2 or none)
#' @param tit Title of your plot
#' @param data_type The type of data you want to visualize; either 'intensity', 'Top3', 'iBAQ' or 'all'.
#' @param design The design of your experiment that match your columns names.
#'               For example, if you have columns named '37_B1_Treatment', put \code{c("temperature", "replicate", "condition")}.
#'               If NULL, will take all values.
#' @param to_check The grouping variable (same as in your design; like 'replicate' for example)
#' @param prop_cut The minimum proportion to show on the graph
#'
#' @return ggplot2 graph showing proportion of valid values
#'
#' @export

validDIA <- function(data, transformation = c("log2", "none"),
                     tit = "", data_type = c("intensity", "Top3", "iBAQ", "all"),
                     design = NULL, to_check = NULL, prop_cut = 0.75){
  if(!inherits(prop_cut, "numeric")){
    stop("'prop_cut' needs to be a numeric between 0 and 1")
  }
  else if(prop_cut > 1){
    prop_cut <- 1
  }
  else if(prop_cut < 0){
    prop_cut <- 0
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
  if(transformation == "log2"){
    if(stringr::str_length(tit) == 0){
      tit2 <- "\nLog2 transformed"
    }
    d <- log2(d)
    d <- as.data.frame(apply(d, 2, function(x) {x[is.infinite(x)] <- NA; x}))
  }
  else if(transformation != "none"){
    transformation <- "none"
  }
  if(stringr::str_length(tit) == 0){
    tit <- paste0("Valid value proportion", tit2)
  }

  if(is.null(design)){
    d <- apply(d, 1, function(x) sum(!is.na(x))/length(x))
    d <- d[order(d, decreasing = TRUE)]
    d <- as.data.frame(d)
    colnames(d) <- "prop_valid"

    g <- ggplot2::ggplot(d, ggplot2::aes(100*(1:nrow(d))/nrow(d), prop_valid)) +
      ggplot2::geom_point() +
      ggplot2::ylim(c(prop_cut,1)) +
      ggplot2::labs(x = "Percentage of IDs",
                    y = "Proportion of valid values",
                    title = tit) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  else if(inherits(design, "character")){
    if(length(to_check) == 1 & to_check %in% design){
      d$id <- data[[1]]
      d <- d %>% tidyr::gather("cond", "value", -id) %>%
        tidyr::separate("cond", into = design) %>%
        group_by_at(c("id", to_check)) %>%
        mutate(prop_valid =  sum(!is.na(value))/length(value))

      d <- unique(d[,c("id", to_check, "prop_valid")])
      d <- d %>% group_by_at(to_check) %>% mutate(ord = 100*(order(order(prop_valid, decreasing = TRUE))/nrow(data)))

      g <- ggplot2::ggplot(d, ggplot2::aes_string("ord", "prop_valid", color = to_check)) +
        ggplot2::geom_point() +
        ggplot2::ylim(c(prop_cut,1)) +
        ggplot2::facet_wrap(~d[[to_check]]) +
        ggplot2::labs(x = "Percentage of IDs",
                      y = "Proportion of valid values",
                      title = tit) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    else{
      stop("'to_check' needs to be one word with the exact same spelling as you put it in 'design'.")
    }
  }
  else{
    stop("The 'design' argument needs to be a character vector.")
  }

  return(g)
}

