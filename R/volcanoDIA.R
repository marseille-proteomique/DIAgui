#' volcanoDIA
#'
#' Plot volcano plot from your DIA data
#'
#' @param data Processed data from DIAnn (from iq processing or diann_matrix for example)
#' @param control The names of the controls in your data (needs to match column names in your data)
#' @param treated The names of the treated in your data (needs to match column names in your data)
#' @param id A character indicating the name of the coulumns that conatins the ids.
#'   If NULL, will check for the column named 'id' or the rownames.
#' @param transformation Which transformation do you want to apply (log2 or none)
#' @param tit Title of your plot
#' @param FDR The FDR for the volcano plot
#' @param FC_cut The fold-change cutoff
#' @param curvature The curvature used for the curve on the volcano plot
#' @param save_file Save results file
#'
#' @return A volcano plot
#'
#' @export

volcanoDIA <- function(data, control, treated, id = NULL,
                       transformation = c("none", "log2"), tit = "",
                       FDR = 0.01, FC_cut = 2.5, curvature = 0.1,
                       save_file = FALSE){
  if(!all(control %in% colnames(data))){
    not_indata <- control[!(control %in% colnames(data))]
    not_indata <- paste(not_indata, collapse = ", ")
    stop(paste(not_indata, ifelse(length(not_indata) > 1, "are", "is"),
               "not in your data ! Please check the control you typed"
               )
         )
  }
  if(!all(treated %in% colnames(data))){
    not_indata <- treated[!(treated %in% colnames(data))]
    not_indata <- paste(not_indata, collapse = ", ")
    stop(paste(not_indata, ifelse(length(not_indata) > 1, "are", "is"),
               "not in your data ! Please check the treated you typed"
               )
         )
  }

  check_num <- apply(data[,c(control, treated)], 2, class)
  if(!all(check_num == "numeric")){
    check_num <- names(check_num)[check_num != "numeric"]
    check_num <- paste(check_num, collapse = ", ")
    stop(paste(check_num, ifelse(length(check_num) > 1, "are", "is"),
               "not in your data !"
               )
         )
  }


  if(is.null(id)){
    if(!("id" %in% names(data))){
      if(sum(is.na(as.numeric(rownames(data)))) == 0){
        stop("Your data doesn't have any rownames (only numeric). Please provide the name of your column id via 'nm_id'.")
      }
      got_id <- apply(data, 2, function(x) sum(x == rownames(data), na.rm = TRUE) == length(x))
      if(sum(got_id) == 1){
        names(data)[got_id] <- "id"

      }
      else{
        data$id <- rownames(data)
      }
    }

    id <- "id"
  }
  else{
    if(class(id) != "character"){
      stop("You need to provide a character for the id column name")
    }
    else{
      if(!(id %in% colnames(data))){
        stop(paste(id, "is not in your data ! Please check your spelling"))
      }
    }
  }

  transformation <- match.arg(transformation)
  tit2 <- ""
  if(transformation == "log2"){
    if(stringr::str_length(tit) == 0){
      tit2 <- ", Log2 transformed"
    }
    data[,c(control, treated)] <- log2(data[,c(control, treated)])
  }
  if(stringr::str_length(tit) == 0){
    tit <- paste0("Volcano plot", tit2)
  }


  # define functions to get cutoff
  curve <- function(x, cut_neg, cut_pos, cut_p, curvature = curvature){
    y <- rep(NA, length(x))
    neg <- which(x <= cut_neg)
    pos <- which(x >= cut_pos)
    y[neg] <- curvature/abs(x[neg] - cut_neg) + -log10(cut_p)
    y[pos] <- curvature/abs(x[pos] - cut_pos) + -log10(cut_p)

    return(y)
  }
  find_cutoff <- function(x,y){
    id <- order(x)
    x <- x[id]
    y <- y[id]

    x <- x[which(x > y)]
    return(x[1])
  }

  volc <- mapply(function(ct, tr){
                    mean_diff <- mean(tr - ct, na.rm = TRUE)
                    pv <- tryCatch(t.test(tr, ct)$p.value,
                                   error = function(e) NA)

                    res <- data.frame(FC = mean_diff, pv = pv);
                    res
                  },
                 as.data.frame(t(data[, control])),
                 as.data.frame(t(data[, treated])),
                 SIMPLIFY = FALSE)

  volc <- as.data.frame(Reduce(rbind, volc))


  cutoff <- volc %>%
    dplyr::mutate(BH = (order(order(pv))/length(pv))*FDR) %>%
    dplyr::reframe(pval = find_cutoff(pv, BH),
                   FC_pos = FC_cut + median(FC[which(pv < quantile(pv, 0.5, na.rm = TRUE))], na.rm = TRUE),
                   FC_neg = -FC_cut - median(FC[which(pv < quantile(pv, 0.5, na.rm = TRUE))], na.rm = TRUE))

  volc <- volc %>%
    dplyr::mutate(criteria = pv <= cutoff$pval & (FC >= cutoff$FC_pos | FC <= cutoff$FC_neg),
                  curve = curve(FC, cutoff$FC_neg,
                                cutoff$FC_pos,
                                cutoff$pval,
                                curvature = curvature
                                ),
                  criteria_curve = -log10(pv) >= curve
                  )
  volc$criteria_curve <- tidyr::replace_na(volc$criteria_curve, FALSE)
  volc$id <- data[[id]]

  df_curve <- data.frame(FC = seq(min(volc$FC, na.rm = TRUE), max(volc$FC, na.rm = TRUE), 0.01))
  df_curve <- df_curve %>%
    dplyr::mutate(curve = curve(FC, cutoff$FC_neg,
                                cutoff$FC_pos, cutoff$pval,
                                curvature = curvature)
                  )


  g <- ggplot(volc, aes(FC, -log10(pv), color = criteria_curve)) +
    geom_point() +
    geom_line(data = df_curve, aes(x = FC, y = curve), linetype = "dashed", color = "black") +
    ylim(c(0, max(-log10(volc$pv), na.rm = TRUE))) +
    labs(title = tit,
         y = "-log10(p-value)",
         x = "fold-change") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggrepel::geom_label_repel(data = volc[volc$criteria_curve,],
                              aes(FC, -log10(pv), label = id), show.legend = FALSE)

  if(save_file){
    wd <- getwd()

    volc$criteria <- NULL
    volc$Curve <- NULL
    volc <- volc[,c("id", "FC", "pv", "criteria_curve")]
    colnames(volc) <- c(id, "fold-change", "p-value", "significant")
    openxlsx::write.xlsx(volc, paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M_"), "VolcanoDIA.xlsx"))

    cutoff_file <- paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "cutoff_VolcanoDIA.txt")
    readr::write_tsv(cutoff, cutoff_file)
    extra_info <- paste0("\nParameters: \nFC cutoff=", FC_cut, ", FDR=", FDR*100, "%, curvature=", curvature)
    write(extra_info, cutoff_file, sep = "\n", append = TRUE)
  }

  return(g)
}
