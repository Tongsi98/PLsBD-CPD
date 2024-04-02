#' Change point detection
#'
#' Detect change points in a time series using RelULSIF.
#'
#' @param ts Time series to detect change points in. Assumes this is a D by N matrix
#'           where D is the dimension of the time series and N is the number of
#'           time points. If a vector is provided, will assume the time series is
#'           one-dimensional.
#' @param window_size The length of the sub-sequences generated from the series. Default `5`.
#' @param step How many sub-sequences forward and backward to from a time point
#'             to compute a score from. Default is 10% of the length of the series
#'             if not specified.
#' @param alpha Relative parameter in [0, 1). Default `0.05`. Setting to `0` recovers
#'              ordinary unconstrained least squares importance fitting.
#' @param k Number of basis functions. Default is minimum of `100` and dimension of
#'          the time series.
#' @param n_folds Number of folds to use in determining optimal kernel bandwidth
#'                and lambda parameter in RULSIF.
#' @param thresh Scalar in (0, 1) indicating the percentile above which a score
#'               is considered a potential change-point. Lower values increase
#'               the sensitivity.
#' @param make_plot Logical. On the same figure, make a plot of each dimension
#'                  of the time series, the rPE scores, and highlight in the time
#'                  series plots the change points detected in red. Default `FALSE`.
#'
#'
#' @return List of 3:
#' - `step`: the step used
#' - `scores`: rPE scores
#' - `change_points`: time points that a change was detected at the given threshold
#' @export
#'
#' @examples
#' s <- c(rnorm(150, mean = 0), rnorm(150, mean = 5), rnorm(150, mean = 1))
#' s <- matrix(s, nrow = 1)
#' ts_detect(s)
ts_detect <- function(method, ts, window_size, step = NULL,
                      alpha, k, n_folds,
                      thresh, make_plot = FALSE) {

    # compatibility checks on input variables
    # ensure ts is a matrix
    if (is.vector(ts)) {
        ts <- matrix(ts, nrow = 1)
    } else if (!is.matrix(ts)) {
        stop("Parameter ts must be a matrix.")
    }

    # dimensional parameters
    dim_ts <- dim(ts)[1]
    N_time_points <- dim(ts)[2]

    # window_size
    if (window_size <= 0 || window_size %% 1 != 0) {
        stop("Parameter window_size must be a positive integer.")
    } else if (window_size > (N_time_points - 1)) {
        stop("Parameter window_size too large. Try reducing it.")
    }

    # step
    if (is.null(step)) {
        step <- floor(0.1 * N_time_points)
    } else if (step <= 0 || step %% 1 != 0) {
        stop("Parameter step must be a positive integer.")
    } else if (step > N_time_points - window_size) {
        stop("Parameter step is too large. Try reducing it.")
    }

    # alpha
    if (alpha < 0 || alpha >= 1 || length(alpha) != 1) {
        stop("Parameter alpha must be in [0, 1).")
    }

    # k
    if (k <= 0 || k %% 1 != 0) {
        stop("Parameter k must be a positive integer.")
    }

    # n_folds
    if (n_folds <= 0 || n_folds %% 1 != 0) {
        stop("Parameter n_folds must be a positive integer.")
    } else if (n_folds > N_time_points) {
        stop("Parameter n_folds exceeds number of points in the input time series.")
    }

    # tresh
    if (thresh <= 0 || thresh >= 1 || length(thresh) > 1) {
        stop("Parameter thresh should be a scalar in (0, 1).")
    }

    # make_plot
    if (!is.logical(make_plot) || length(make_plot) > 1) {
        stop("Parameter make_plot should be a single logical value.")
    }

    # constructing sliding window
    sw <- sliding_window(X = ts, window_size = window_size)
    n_samples <- dim(sw)[2]

    # compute change-point scores
    scores <- vector()
    t <- step + 1
    while(t + step - 1 <= n_samples) {
        y <- sw[ , (t - step):(step + t - 1), drop = FALSE]
        y_nu <- y[ , 1:step]
        y_de <- y[ , (step + 1):ncol(y)]

        out <- RelULSIF(method, y_nu, y_de, alpha = alpha, k = k, n_folds=n_folds)
        scores <- c(scores, out$rPE)
        t <- t + 1
    }

    # get change point indices
    change_points <- which(scores > (thresh * max(scores))) + step + window_size

    # make plot, if asked for
    if (make_plot) {
        ts_plot(ts, step, window_size, scores, change_points)
    }

    # return
    return(list(
        step = step,
        scores = scores,
        change_points = change_points
    ))
 
 }
