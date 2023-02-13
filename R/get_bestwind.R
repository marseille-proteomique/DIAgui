#' get_bestwind
#'
#' Get the best m/z windows for DIA according a number of window and the report-lib from DIA-nn.
#'
#' @param data The report-lib data from DIA-nn which contains the \emph{ProductMz} column.
#'             It can be a path to the file or the dataframe corresponding to the file
#'             (need to have \emph{FileName} and \emph{ProductMz} columns in that case).
#' @param bins The number of windows you want to have.
#' @param per_frac If \code{FALSE}, will select the best windows from all fraction without differentiation.
#'                 Else, it will select the best window for each fraction and then do the mean of those.
#'
#' @return A list containing a data frame with \emph{FileName}, \emph{ProductMz} and a \emph{bins} column;
#'         four plots corresponding to bar plots / histogramm from the distribution of the \emph{ProductMz}
#'         with and without the window selection.
#'
#' @details
#' In the end, the distribution of the \emph{ProductMz} per window will be uniform.
#'
#' @export

get_bestwind <- function(data, bins = 25, per_frac = FALSE){
  if(inherits(data, "character")){
    if(stringr::str_detect(data, "\\.tsv$"))
      data <-  diann_load(data)
    else{
      message("You need to import the report-lib file output from DIA-nn which is in tsv format.
              \nThe file you put is in another format; please provide the right file.")
      return()
    }
  }
  data <- as.data.frame(data)
  data <- data %>% dplyr::select(FileName, ProductMz)

  orig_hist <- ggplot2::ggplot(data, ggplot2::aes(x=ProductMz)) +
    ggplot2::geom_histogram(color="black", fill="white", bins = bins)
  orig_hist_frac <- ggplot2::ggplot(data, ggplot2::aes(x=ProductMz)) +
    ggplot2::geom_histogram(color="black", fill="white", bins = bins) +
    ggplot2::facet_wrap(~FileName)

  data <- data[order(data$ProductMz),]

  if(per_frac){
    best_wind <- list()
    frac = unique(data$FileName)
    ### first get best windows for each fraction
    for(i in frac){
      data_bis <- data %>% dplyr::select(FileName, ProductMz) %>%
        dplyr::filter(FileName == i)

      n = nrow(data_bis)
      n_average = n/bins
      n_average_down = floor(n/bins)
      n_average_up = ceiling(n/bins)

      a = round(-1*bins*(n_average - n_average_up))
      b = bins - a
      na = n_average_down
      nb = n_average_up
      if(a > b){
        b = a
        a = bins - b
        na = n_average_up
        nb = n_average_down
        # a will always be < b
      }

      med = median(1:bins)
      if(bins %% 2 == 0){ # if bins is even, end config with majority --> nb
        b = b-1
      }

      if(a == 0){ # if n/bins is an integer
        config = rep(nb, bins)
      }
      else if(b %% 2 == 0){ # b is even, a on the median
        config = c(nb, rep(na, (a-1)/2), rep(nb, (b-2)/2), na, rep(nb, (b-2)/2), rep(na, (a-1)/2), nb)
      }
      else if(b %% 2 == 1){ # a is even, b on the median
        config = c(nb, rep(na, (a-2)/2), rep(nb, (b-3)/2), na, nb, na, rep(nb, (b-3)/2), rep(na, (a-2)/2), nb)
      }

      if(bins %% 2 == 0){
        config = append(config, nb)
      }

      wind <- c()
      for(k in 1:bins){ # get the windows --> same amount for each
        if(k == 1){
          dep = 1
        }
        else{
          dep = sum(config[1:(k-1)])
        }
        val = data_bis$ProductMz[dep:sum(config[1:k])]
        wind <- append(wind, round(min(val)))
      }

      best_wind[[i]] <- c(wind, round(max(val)))*n # weight them by count in order to average the windows
    }

    n = nrow(data)
    # average them
    best_wind <- Reduce("+", best_wind)
    best_wind <- best_wind/n
    #get them
    best_wind <- c(floor(min(data$ProductMz)), round(best_wind[-c(1, bins+1)]), ceiling(max(data$ProductMz)))

    data$bins <- rep(NA, n)
    for(i in 1:bins){# write in the data
      down <- best_wind[i]
      up <- best_wind[i+1]
      idx <- which(data$ProductMz >= down & data$ProductMz < up)
      data$bins[idx] <- paste0(down, "-", up)
    }
  }
  else{
    ### get best windows from all data
    n = nrow(data)
    n_average = n/bins
    n_average_down = floor(n/bins)
    n_average_up = ceiling(n/bins)

    a = round(-1*bins*(n_average - n_average_up))
    b = bins - a
    na = n_average_down
    nb = n_average_up
    if(a > b){
      b = a
      a = bins - b
      na = n_average_up
      nb = n_average_down
      # a will always be < b
    }

    med = median(1:bins)
    if(bins %% 2 == 0){ # if bins is even, end config with majority --> nb
      b = b-1
    }

    if(a == 0){ # if n/bins is an integer
      config = rep(nb, bins)
    }
    else if(b %% 2 == 0){ # b is even, a on the median
      config = c(nb, rep(na, (a-1)/2), rep(nb, (b-2)/2), na, rep(nb, (b-2)/2), rep(na, (a-1)/2), nb)
    }
    else if(b %% 2 == 1){ # a is even, b on the median
      config = c(nb, rep(na, (a-2)/2), rep(nb, (b-3)/2), na, nb, na, rep(nb, (b-3)/2), rep(na, (a-2)/2), nb)
    }

    if(bins %% 2 == 0){
      config = append(config, nb)
    }

    data$bins <- rep(NA, n)
    for(i in 1:bins){ # get the windows --> same amount for each
      if(i == 1){
        dep = 1
      }
      else{
        dep = sum(config[1:(i-1)])
      }
      val = data$ProductMz[dep:sum(config[1:i])]
      if(i == 1){
        data$bins[dep:sum(config[1:i])] <- paste0(floor(min(val)), "-", round(max(val)))
      }
      else if(i == bins){
        data$bins[dep:sum(config[1:i])] <- paste0(round(min(val)), "-", ceiling(max(val)))
      }
      else{
        data$bins[dep:sum(config[1:i])] <- paste0(round(min(val)), "-", round(max(val)))
      }
    }
  }

  # print best windows
  un_wind <- unique(data$bins)
  message(paste("The best windows are", paste(un_wind, collapse = ", ")))

  # save plots
  new_hist <- ggplot2::ggplot(data, ggplot2::aes(x = factor(bins, levels = un_wind))) +
    ggplot2::geom_bar(stat="count", width=0.7, fill="steelblue") +
    ggplot2::labs(x = "bins")
  new_hist_frac <-ggplot2:: ggplot(data, ggplot2::aes(x = factor(bins, levels = un_wind))) +
    ggplot2::geom_bar(stat="count", width=0.7, fill="steelblue") +
    ggplot2::facet_wrap(~FileName) +
    ggplot2::labs(x = "bins")

  res <- list("data" = data, "orig_hist" = orig_hist, "orig_hist_perfrac" = orig_hist_frac,
              "new_hist" = new_hist, "new_hist_perfrac" = new_hist_frac)

  return(res)
}



