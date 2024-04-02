#' Compute median squared euclidean distance
#'
#' @param x Matrix
#'
#' @return Scalar: median squared euclidean distance
#'
comp_median <- function(x) {
    # dimension parameters
    d <- dim(x)[1]
    n <- dim(x)[2]

    # calculate median
    G <- matrix(rep(colSums(x ^ 2), times = n), nrow = n, byrow = TRUE)
    D <- G - (2 * t(x) %*% x) + t(G)
    D[upper.tri(D)] <- 0

    # return
    return(sqrt(0.5 * stats::median(D[D > 0])))
}
