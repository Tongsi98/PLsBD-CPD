#' Matrix of subsequences (sliding window)
#'
#' @param X Data matrix
#' @param window_size How many subsequences constitute a single window
#'
#' @return Matrix of subsequences
#'
sliding_window <- function(X, window_size) {
    # dimension parameters
    dim_ts <- dim(X)[1]
    N_points <- dim(X)[2]
    N_subsequences <- N_points - window_size + 1

    # create template matrix to fill in
    W <- matrix(0, nrow = dim_ts * window_size, ncol = N_subsequences)

    # fill in the result
    for (i in 1:N_subsequences) {
        W[ , i] <- as.vector(X[ , i:(i + window_size - 1)])
    }

    # return
    return(W)
}
