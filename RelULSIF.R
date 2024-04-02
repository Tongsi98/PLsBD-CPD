#' Estimate the ratio of probability densities that generated
#' data in two samples (Xnu and Xde).
#'
#' @param Xnu Samples from numerator probability density
#' @param Xde Samples from denomenator probability density
#' @param Xce Matrix of centers
#' @param sigma Scalar or vector; Gaussian kernel bandwidth(s). Positive.
#' @param lambda Scalar or vector; regularization parameter(s). Non-negative.
#' @param alpha Scalar; relative parameter in [0, 1)
#' @param k Positive integer; number of basis functions 
#' @param n_folds Integer; number of folds to use in cross fold validation 5
#'
#' @return List of 3
#' - opt_sigma: chosen sigma from CV
#' - opt_lambda: chosen lambda from CV
#' - rPE: relative Pearson divergence or Pearson-like scaled Bregman divergence
#' @export
#'
#' @examples
#' X <- matrix(
#'     c(rnorm(50), rnorm(50, mean = 5), rnorm(50, mean = -5),
#'     rnorm(100), rnorm(25, mean = 3), rnorm(25, mean = -1),
#'     rnorm(25), rnorm(75, mean = -2), rnorm(50, mean = 4)),
#'     nrow = 3, ncol = 150, byrow = TRUE
#' )
#' Xnu <- X[ , 1:50]
#' Xde <- X[ , 51:ncol(X)]
#' RelULSIF(Xnu, Xde)
RelULSIF <- function(method, Xnu, Xde, Xce = NULL, sigma = NULL, lambda = NULL,
                   alpha, k, n_folds) {
    # compatibility checks on data types
    if (is.vector(Xnu)) {
        Xnu <- as.matrix(Xnu)
    } else if (!is.matrix(Xnu)) {
        stop("Parameter Xnu should be a matrix.")
    }
    if (is.vector(Xde)) {
        Xde <- as.matrix(Xde)
    } else if (!is.matrix(Xnu)) {
        stop("Parameter Xde should be a matrix.")
    }

    # compatibility checks on alpha
    if (alpha < 0 || alpha >= 1 || length(alpha) != 1) {
        stop("Parameter alpha must be in [0, 1).")
    }

    # check on sigma
    if (!is.null(sigma)) {
        if (length(sigma) > 1) {
            if (any(sigma) <= 0) {
                stop("Values for sigma must be positive.")
            }
        } else {
            if (sigma <= 0) {
                stop("Parameter sigma must be positive.")
            }
        }
    }

    # check on lambda
    if (!is.null(lambda)) {
        if (length(lambda) > 1) {
            if (any(lambda) < 0) {
                stop("Values for sigma must be non-negative.")
            }
        } else {
            if (lambda <= 0) {
                stop("Parameter lambda must be non-negative.")
            }
        }
    }

    # check on k
    if (k <= 0 || length(k) > 1 || k %% 1 != 0) {
        stop("Parameter k must be a positive integer.")
    }

    # check on n_folds
    if (n_folds <= 0 || length(n_folds) > 1 || n_folds %% 1 != 0) {
        stop("Parameter n_folds must be a positive integer.")
    }

    # dimension compatibility checks
    dim_nu <- dim(Xnu)[1]
    dim_de <- dim(Xde)[1]
    n_nu <- dim(Xnu)[2]
    n_de <- dim(Xde)[2]
    if (dim_nu != dim_de) {
        stop("Number of rows must match between Xnu and Xde.")
    }

    # ensure that k does not exceed dimension of the data
    k <- min(k, dim_nu)

    # initialize centers for Gaussian kernel if NULL
    if (is.null(Xce)) {
        Xce <- Xnu[ , sample(n_nu, size = k, replace = TRUE), drop = FALSE]
    } else{
        # check compatibility
        if (dim(Xce)[1] != dim(Xnu)[1]) {
            stop("Dimensions of Xce and Xnu or Xde incompatible.")
        }
    }

    # distance matrices
    dist_nu <- comp_dist(Xnu, Xce)
    dist_de <- comp_dist(Xde, Xce)

    # check if search on sigma and lambda is needed
    if (length(sigma) == 1 && length(lambda) == 1) {
        # proceed with user specified sigma and lambda values
        opt_sigma <- sigma
        opt_lambda <- lambda
    } else {
        # first check if user specified sigma/lambda search ranges
        sigma_search_vec <- sigma
        lambda_search_vec <- lambda
        if (is.null(sigma)) {
            med <- comp_median(cbind(Xde, Xnu))
            sigma_search_vec <- med * seq(0.2, 1.6, 0.2)
        }
        if (is.null(lambda)) {
            lambda_search_vec <- 10 ^ seq(-3, 1, 1)
        }
        cv_out <- sigma_lambda_grid_search(n_nu, n_de, k, dist_nu, dist_de,
            sigma_search_vec, lambda_search_vec, alpha, n_folds)
        opt_sigma <- cv_out$opt_sigma
        opt_lambda <- cv_out$opt_lambda
    }

    # compute results
    ker_Xnu <- gaussian_kernel(t(dist_nu), opt_sigma)
    ker_Xde <- gaussian_kernel(t(dist_de), opt_sigma)
    
    
   ################################ squared loss 
    H <- ((1 - alpha) / n_de) * (ker_Xde %*% t(ker_Xde)) +
        (alpha / n_nu) * (ker_Xnu %*% t(ker_Xnu))
    h <- rowMeans(ker_Xnu)
    theta <- solve(H + diag(k) * opt_lambda) %*% h
    #################################
    
    ######below is absolute loss results    
    # hp <- rowMeans(ker_Xnu)
    # hq <- rowMeans(ker_Xde)    
    # theta <- 1/opt_lambda * (alpha * hp + (1-alpha)*hq)  
    #####################

    g_nu <- t(theta) %*% ker_Xnu
    g_de <- t(theta) %*% ker_Xde
    
    if (method=="PD"){
## Below is old Pearson divergence
    rPE <- mean(g_nu) - 0.5 * ((alpha * mean(g_nu ^ 2)) +
       (1 - alpha) * mean(g_de ^ 2)) - 0.5

#    rPE <- 0.5*mean(g_nu) - 0.5
    
  }

## Below is Pearson like scaled Bregman divergence
  if( method =="PL"){         
 rPE =   1/2 * mean(g_nu) - (2-alpha)/(2*(1-alpha)) * mean(g_de) +
             1/(2*(1-alpha))
 # rPE =   alpha * mean(g_nu ^ 2) + (1 - alpha) * mean(g_de ^ 2)   
 #           - 2* mean(g_de) - alpha/(1-alpha) * mean(g_de) +
 #             1/(1-alpha)

  }
 ##########################################################     
    
    # return
    return(list(
        opt_sigma = opt_sigma,
        opt_lambda = opt_lambda,
        rPE = rPE
    ))
}
