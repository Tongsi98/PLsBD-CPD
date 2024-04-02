#' Gaussian kernel function
#'
#' @param distance_matrix Matrix of squared euclidean matrices
#' @param sigma Kernel bandwidth parameter
#'
#' @return Gaussian kernel applied to the input matrix
#'
gaussian_kernel <- function(distance_matrix, sigma) {
    exp(-distance_matrix / (2 * sigma ^ 2))
}
