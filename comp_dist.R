#' Squared euclidean distance for matrices
#'
#' @param x1 Matrix
#' @param x2 Matrix
#'
#' @return Squared euclidean distance between points in x1 and x2
#'
comp_dist <- function(x1, x2) {
    # dimensional parameters
    nx1 <- dim(x1)[2]
    nx2 <- dim(x2)[2]

    # compute distance matrix
    X1 <- matrix(rep(colSums(x1 ^ 2), times = nx2), nrow = nx2, byrow = TRUE)
    X2 <- matrix(rep(colSums(x2 ^ 2), times = nx1), nrow = nx1, byrow = TRUE)
    D <- t(X1) + X2 - (2 * t(x1) %*% x2)

    # return
    return(D)
}
