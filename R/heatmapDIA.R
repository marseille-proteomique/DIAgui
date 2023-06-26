#' heatmapDIA
#'
#' Print interactive heatmap of processed data from DIAnn
#'
#' @param data Processed data from DIAnn (from iq processing or diann_matrix for example)
#' @param transformation Which transformation do you want to apply (log2, z.score on proteins or on fractions, none)
#' @param maxna How many NAs do you authorize per proteins
#' @param print_val logical; do you want to print values on heatmpa ?
#' @param nm_id The name of column that contains the IDs (if NULL, will take it automatically)
#' @param data_type The type of data you want to visualize; either 'intensity', 'Top3', 'iBAQ' or 'all'.
#' @param gradient_color Three colors for the printed gradient color on the heatmap
#' @param static Logical to return a static plot (ggplot object) or an interactive one (plotly object)
#'
#' @return Plotly heatmap
#'
#' @export

heatmapDIA <- function(data, transformation = c("none", "log2", "z.score_proteins", "z.score_fraction"),
                       maxna = 0, print_val = TRUE, nm_id = NULL,
                       data_type = c("intensity", "Top3", "iBAQ", "all"),
                       gradient_color = c("#09009D", "#ffffff", "#BE0010"),
                       static = FALSE){
  transformation <- match.arg(transformation)
  data_type <- match.arg(data_type)
  if(length(gradient_color) != 3){
    stop("Please provide exactly three colors for the gradient color.")
  }

  H <- as.data.frame(data)
  if(length(names(H)) == 1){
    frac <- names(H)
  }
  if(is.null(nm_id)){
    if(!("id" %in% names(H))){
      if(sum(is.na(as.numeric(rownames(H)))) == 0){
        stop("Your data doesn't have any rownames (only numeric). Please provide the name of your column id via 'nm_id'.")
      }
      got_id <- apply(H, 2, function(x) sum(x == rownames(H), na.rm = TRUE) == length(x))
      if(sum(got_id) == 1){
        names(H)[got_id] <- "id"
      }
      else{
        H$id <- rownames(H)
      }
    }
  }
  else{
    idx <- stringr::str_which(names(H), paste0("^", nm_id, "$"))
    if(length(idx)){
      names(H)[idx] <- "id"
    }
    else{
      stop(paste("Please provide a valid column name. You enter :", nm_id,
                 "and the column names are :", paste(names(H), collapse = ", "), "."))
    }
  }
  cl <- lapply(H, class)
  cl <- cl == "numeric"
  if(sum(cl) == 0){
    message("No numeric data !")
    return(NULL)
  }
  else if(sum(cl) == 1){
    frac <- names(H)[cl]
    rw <- rownames(H)
    H <- as.data.frame(cbind(H[,"id"], H[,cl]))
    H[,2] <- as.numeric(H[,2])
    names(H) <- c("id", frac)
    rownames(H) <- rw
  }
  else{
    H <- as.data.frame(cbind(H[,"id"], H[,cl]))
    names(H)[1] <- "id"
    to_rm <- stringr::str_which(colnames(H), "nbTrypticPeptides|peptides_counts_all|^pep_count_")
    if(length(to_rm) > 0){
      if(length(to_rm) == ncol(H)){
        message("No numeric data !")
        return(NULL)
      }
      else{
        H <- H[,-to_rm]
      }
    }

    if(data_type == "Top3"){
      H <- H[,c(1, stringr::str_which(colnames(H), "^Top3_"))]
    }
    else if(data_type == "iBAQ"){
      H <- H[,c(1, stringr::str_which(colnames(H), "^iBAQ_"))]
    }
    else if(data_type == "intensity"){
      idx <- stringr::str_which(colnames(H), "^iBAQ_|^Top3_")
      if(length(idx) > 0){
        H <- H[,-idx]
      }
    }
    else if(data_type == "all"){
      message("You chose to keep all numeric data, they may differ completly.")
    }
    else{
      stop("data_type can only be 'intensity', 'Top3', 'iBAQ' or 'all' .")
    }
  }

  if(transformation == "log2"){
    idx <- stringr::str_which(names(H), "^id$")
    H <- as.data.frame(cbind(H[,idx], log2(H[,-idx])))
    names(H)[1] <- "id"
    if(length((names(H)[-1])) == 1){
      H[,2] <- as.numeric(H[,2])
      names(H)[2] <- frac
    }
  }
  H <- tidyr::gather(H, fraction, intensity, -id)
  H <- H[order(H$id),]
  if(stringr::str_detect(transformation, "^z\\.score_")){
    if(length(unique(H$fraction)) == 1){
      message("No z.score transformation; need at least 3 fractions")
    }
    else{
      if(stringr::str_detect(transformation, "_protein$")){
        H <- H %>% dplyr::group_by(id) %>%
          dplyr::mutate(Mean = mean(intensity, na.rm = TRUE),
                        SD = sd(intensity, na.rm = TRUE))
        H$z.score <- (H$intensity - H$Mean)/H$SD

        H$intensity <- NULL
        colnames(H)[ncol(H)] <- "intensity"
      }
      else if(stringr::str_detect(transformation, "_fraction$")){
        H <- H %>% dplyr::group_by(fraction) %>%
          dplyr::mutate(Mean = mean(intensity, na.rm = TRUE),
                        SD = sd(intensity, na.rm = TRUE))
        H$z.score <- (H$intensity - H$Mean)/H$SD

        H$intensity <- NULL
        colnames(H)[ncol(H)] <- "intensity"
      }
      else{
        stop("If you want to perform a z.score transformation, please precise if you want to
           perform it on proteins or on fractions : 'z.score_protein' or 'z.score_fraction'.")
      }
    }
  }
  H_na <- length(unique(H$id))
  H$intensity[which(is.infinite(H$intensity))] <- NA
  H <- H %>% dplyr::group_by(id) %>%
    dplyr::filter(sum(is.na(intensity)) <= maxna)
  message(paste(H_na - length(unique(H$id)), "proteins has been removed because thay had more
                than", maxna, "missing values."))
  d_H <- H[,c("id", "fraction", "intensity")]
  d_H <- tidyr::spread(d_H, fraction, intensity)
  d_H <- dist(d_H[,-1])
  d_H[which(is.na(d_H))] <- 0
  prot_dend <- hclust(d_H)
  H$id <- factor(H$id, levels = unique(H$id)[prot_dend$order])
  H$id2 <- H$id
  H$id2 <- as.character(lapply(stringr::str_split(H$id2, ";"), function(x){
    if(length(x) > 1){
      if(length(x) > 3){
        l <- 3
      }
      else{
        l <- length(x)
      }
      paste(x[1:l], collapse = ";")
    }
    else{
      x
    }
  }))
  lv <- levels(H$id)
  lv <- as.character(lapply(stringr::str_split(lv, ";"), function(x){
    if(length(x) > 1){
      if(length(x) > 3){
        l <- 3
      }
      else{
        l <- length(x)
      }
      paste(x[1:l], collapse = ";")
    }
    else{
      x
    }
  }))
  lv <- unique(lv)
  H$id2 <- factor(H$id2, levels = lv)
  R <- range(H$intensity, na.rm = TRUE)
  R <- c(floor(R[1]), ceiling(R[2]))
  R <- c(R[1], mean(R), R[2])

  g <- ggplot2::ggplot(H, ggplot2::aes(fraction, id2, fill = intensity,
                                       text = paste("UniprotID:", id2,
                                                    "<br>Fraction:", fraction,
                                                    "<br>Intensity:", intensity)
                                       )
                       ) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = element_blank()) +
    ggplot2::scale_fill_gradientn(breaks = R,
                         colors = gradient_color,
                         limits = c(R[1], R[3]))

  if(print_val){
    g <- g +
      ggplot2::geom_text(ggplot2::aes(label = round(intensity,3)), color = "black", size = 3)
  }

  if(!static){
    g <- plotly::ggplotly(g, tooltip = c("text")) %>%
      plotly::layout(yaxis = list(range = c(1,30))) %>%
      plotly::config(displaylogo = FALSE,
                     displayModeBar = TRUE,
                     modeBarButtonsToRemove = list('sendDataToCloud',
                                                   'hoverClosestCartesian',
                                                   'hoverCompareCartesian',
                                                   'select2d',
                                                   'lasso2d',
                                                   'toggleSpikelines'
                     )
      )
  }

  return(g)
}
