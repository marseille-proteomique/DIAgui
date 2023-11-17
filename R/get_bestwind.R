#' get_bestwind
#'
#' Get the best m/z windows for DIA according a number of window and the report-lib from DIA-nn.
#'
#' @param data The report-lib data from DIA-nn which contains the \emph{PrecursorMz} column.
#'             It can be a path to the file or the dataframe corresponding to the file
#'             (need to have \emph{FileName} and \emph{PrecursorMz} columns in that case).
#' @param n_window The number of windows you want to have.
#' @param per_frac If \code{FALSE}, will select the best windows from all fraction without differentiation.
#'                 Else, it will select the best window for each fraction and then do the mean of those.
#' @param overlap A fixed overlap between windows. Default is 0.
#' @param window_size A fix m/z window size. Default is NULL but if numeric, it will compute n
#'    windows of same size.
#'
#' @return A list containing a data frame with \emph{FileName}, \emph{PrecursorMz} and a \emph{bins} column;
#'         four plots corresponding to bar plots / histogramm from the distribution of the \emph{PrecursorMz}
#'         with and without the window selection.
#'
#' @details
#' In the end, the distribution of the \emph{PrecursorMz} per window will be uniform.
#'
#' @export

get_bestwind <- function(data, n_window = 25, per_frac = FALSE,
                         overlap = 0, window_size = NULL){
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
  data <- data %>% dplyr::select(FileName, PrecursorMz)

  orig_hist <- ggplot2::ggplot(data, ggplot2::aes(x=PrecursorMz)) +
    ggplot2::geom_histogram(color="black", fill="white", bins = n_window)
  orig_hist_frac <- ggplot2::ggplot(data, ggplot2::aes(x=PrecursorMz)) +
    ggplot2::geom_histogram(color="black", fill="white", bins = n_window) +
    ggplot2::facet_wrap(~FileName)

  data <- data[order(data$PrecursorMz),]

  if(is.numeric(window_size)){
    window_size <- abs(window_size)
    low <- floor(min(data$PrecursorMz, na.rm = TRUE))
    up <- ceiling(max(data$PrecursorMz, na.rm = TRUE))

    config <- unique(c(seq(low, up, window_size), up))
    data$bins <- NA
    for(i in 1:(length(config) - 1)){
      data$bins[which(data$PrecursorMz >= config[i] & data$PrecursorMz < config[i+1])] <- paste0(config[i], "-", config[i+1])
    }
  }
  else{
    if(per_frac){
      best_wind <- list()
      frac = unique(data$FileName)
      ### first get best windows for each fraction
      for(i in frac){
        data_bis <- data %>% dplyr::select(FileName, PrecursorMz) %>%
          dplyr::filter(FileName == i)

        n = nrow(data_bis)
        n_average = n/n_window
        n_average_down = floor(n/n_window)
        n_average_up = ceiling(n/n_window)

        a = round(-1*n_window*(n_average - n_average_up))
        b = n_window - a
        na = n_average_down
        nb = n_average_up
        if(a > b){
          b = a
          a = n_window - b
          na = n_average_up
          nb = n_average_down
          # a will always be < b
        }

        med = median(1:n_window)
        if(n_window %% 2 == 0){ # if n_window is even, end config with majority --> nb
          b = b-1
        }

        if(a == 0){ # if n/n_window is an integer
          config = rep(nb, n_window)
        }
        else if(b %% 2 == 0){ # b is even, a on the median
          config = c(nb, rep(na, (a-1)/2), rep(nb, (b-2)/2), na, rep(nb, (b-2)/2), rep(na, (a-1)/2), nb)
        }
        else if(b %% 2 == 1){ # a is even, b on the median
          config = c(nb, rep(na, (a-2)/2), rep(nb, (b-3)/2), na, nb, na, rep(nb, (b-3)/2), rep(na, (a-2)/2), nb)
        }

        if(n_window %% 2 == 0){
          config = append(config, nb)
        }

        config <- c(0, cumsum(config))
        un_wind <- list()
        for(k in 1:(length(config) - 1)){
          wind <- range(data_bis$PrecursorMz[(config[k]+1):config[k+1]], na.rm = TRUE)
          wind[1] <- floor(wind[1])
          wind[2] <- ceiling(wind[2])
          un_wind[[k]] <- wind
        }
        for(k in 1:(length(un_wind) - 1)){
          if(un_wind[[k]][2] == un_wind[[length(un_wind)]][2]){
            un_wind[[k]][2] <- un_wind[[length(un_wind)]][1]
          }
          if(un_wind[[k]][2] != un_wind[[k+1]][1]){
            un_wind[[k+1]][1] <- un_wind[[k]][2]
          }
        }

        best_wind[[i]] <- lapply(un_wind, "*", n) # weight them by count in order to average the windows
      }

      n = nrow(data)
      # average them
      best_wind <- lapply(best_wind, function(x) t(data.frame(x)))
      best_wind <- Reduce("+", best_wind)
      best_wind <- best_wind/n
      best_wind <- unname(apply(best_wind, 1, c, simplify = FALSE))

      data$bins <- rep(NA, n)
      for(i in best_wind){ # get the windows --> same amount for each
        data$bins[which(data$PrecursorMz >= i[1] & data$PrecursorMz < i[2])] <- paste(i, collapse = "-")
      }
    }
    else{
      ### get best windows from all data
      n = nrow(data)
      n_average = n/n_window
      n_average_down = floor(n/n_window)
      n_average_up = ceiling(n/n_window)

      a = round(-1*n_window*(n_average - n_average_up))
      b = n_window - a
      na = n_average_down
      nb = n_average_up
      if(a > b){
        b = a
        a = n_window - b
        na = n_average_up
        nb = n_average_down
        # a will always be < b
      }

      med = median(1:n_window)
      if(n_window %% 2 == 0){ # if n_window is even, end config with majority --> nb
        b = b-1
      }

      if(a == 0){ # if n/n_window is an integer
        config = rep(nb, n_window)
      }
      else if(b %% 2 == 0){ # b is even, a on the median
        config = c(nb, rep(na, (a-1)/2), rep(nb, (b-2)/2), na, rep(nb, (b-2)/2), rep(na, (a-1)/2), nb)
      }
      else if(b %% 2 == 1){ # a is even, b on the median
        config = c(nb, rep(na, (a-2)/2), rep(nb, (b-3)/2), na, nb, na, rep(nb, (b-3)/2), rep(na, (a-2)/2), nb)
      }

      if(n_window %% 2 == 0){
        config = append(config, nb)
      }

      config <- c(0, cumsum(config))

      un_wind <- list()
      for(i in 1:(length(config) - 1)){
        wind <- range(data$PrecursorMz[(config[i]+1):config[i+1]], na.rm = TRUE)
        wind[1] <- floor(wind[1])
        wind[2] <- ceiling(wind[2])
        un_wind[[i]] <- wind
      }
      for(i in 1:(length(un_wind) - 1)){
        if(un_wind[[i]][2] == un_wind[[length(un_wind)]][2]){
          un_wind[[i]][2] <- un_wind[[length(un_wind)]][1]
        }
        if(un_wind[[i]][2] != un_wind[[i+1]][1]){
          un_wind[[i+1]][1] <- un_wind[[i]][2]
        }
      }

      data$bins <- rep(NA, n)
      for(i in un_wind){ # get the windows --> same amount for each
        data$bins[which(data$PrecursorMz >= i[1] & data$PrecursorMz < i[2])] <- paste(i, collapse = "-")
      }
    }
  }

  # use overlap window
  un_wind <- unique(data$bins)
  if(is.numeric(overlap)){
    overlap <- abs(overlap)
    if(overlap > 0){
      mi <- min(as.numeric(unlist(strsplit(un_wind, "-"))))

      if(overlap > mi){
        message("The overlap you selected is greater than the minimal Precursor mz ! No overlap were applied.")
      }
      else{
        ma <- max(as.numeric(unlist(strsplit(un_wind, "-"))))
        un_wind <- sapply(un_wind,
                          function(x){
                            x <- as.numeric(strsplit(x, "-")[[1]])
                            if(x[1] > mi){
                              x[1] <- x[1] - overlap/2
                            }
                            if(x[2] < ma){
                              x[2] <- x[2] + overlap/2
                            }
                            x <- paste0(x, collapse = "-")
                          })

        bins_overlap <- list()
        for(i in un_wind){
          new_win <- as.numeric(strsplit(i, "-")[[1]])
          new_data <- data[which(data$PrecursorMz >= new_win[1] & data$PrecursorMz < new_win[2]),]
          new_data$bins <- i
          bins_overlap[[i]] <- new_data
        }
        data <- as.data.frame(Reduce(rbind, bins_overlap))
      }
    }
  }
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
