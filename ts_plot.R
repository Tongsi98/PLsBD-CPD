#' Plot change point detection results
#'
#' @param ts Time series input
#' @param step Starting point of detection
#' @param scores rPE scores
#' @param change_points Detected change point indices
#'
ts_plot <- function(ts, step, window_size, scores, change_points) {
    # dimension parameters
    ts_dim <- dim(ts)[1]
    n_time_points <- dim(ts)[2]

    # graphical parameters
    graphics::par(mfrow = c(ts_dim + 1, 1),
        oma = rep(1, 4),
        mar = rep(1.2, 4))

    # plot each time series
    for (s in 1:ts_dim) {
        # plot the time series
        plot(1:n_time_points, ts[s, ], type = 'l',
             xlab = "time", ylab = paste("y", s, sep = ''),
             xlim = c(1, n_time_points),
             ylim = c(min(ts[s, ]), max(ts[s, ])))

        # plot rectangles to highlight change points
        graphics::rect(xleft = change_points,
             ybottom = rep(min(ts[s, ]), length(change_points)),
             xright = change_points + 1,
             ytop = rep(max(ts[s, ]), length(change_points)),
             col = 'blue', border = NA)
    }
    
   # abline(v = c(101, 201, 301, 401, 501,601, 701, 801, 901), col = "blue")
   # dev.new() ## pop a new windows, can be removed if needed
    # plot scores
    plot(1:length(scores) + step + window_size, scores,      
         type = 'l',
         xlab = 'time', ylab = 'Score',
         xlim = c(1, n_time_points),
         col = 'black')
}

# Function to draw rectangles
draw_rectangles <- function(change_points, ts_dim, n_time_points) {
  for (cp in change_points) {
    rect(cp - 0.5, 0.5, cp + 0.5, ts_dim + 0.5,
         col = '#DD666666', border = NA)
  }
}

