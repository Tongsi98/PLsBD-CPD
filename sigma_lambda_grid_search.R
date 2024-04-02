#' Grid search via cross validation on sigma and lambda parameters
#'
#' @param Nnu Number of points in numerator
#' @param Nde Number of points in denomenator
#' @param k Number of basis functions
#' @param dist_nu Distance matrix for numerator
#' @param dist_de Distance matrix for denomenator
#' @param sigma_list Sigma values to search over
#' @param lambda_list Lambda values to search over
#' @param alpha Relative parameter
#' @param n_folds Number of folds for CV
#'
#' @return Optimal sigma-lambda pair
#'
sigma_lambda_grid_search <- function(Nnu, Nde, k, dist_nu, dist_de,
    sigma_list, lambda_list, alpha, n_folds) {
    # parameters
    n_sigma <- length(sigma_list)
    n_lambda <- length(lambda_list)


    # create template score matrix
    scores <- matrix(0, nrow = n_sigma, ncol = n_lambda)

    # loop through each sigma
    for (s in 1:n_sigma) {
        # get kernel values at this sigma
        ker_nu <- gaussian_kernel(dist_nu, sigma_list[s])
        ker_de <- gaussian_kernel(dist_de, sigma_list[s])

        # loop through each lambda at this sigma
        for (l in 1:n_lambda) {
            # cross fold split
            folds_nu <- sample(cut(seq(1, Nnu), breaks = n_folds, labels = FALSE))
            folds_de <- sample(cut(seq(1, Nde), breaks = n_folds, labels = FALSE))

            # begin cross validation
            obj_score_sum <- 0
            for (f in 1:n_folds) {
                # get train and test folds
                nu_train <- t(ker_nu[folds_nu != f, , drop = FALSE])
                de_train <- t(ker_de[folds_de != f, , drop = FALSE])
                nu_test <- t(ker_nu[folds_nu == f, , drop = FALSE])
                de_test <- t(ker_de[folds_de == f, , drop = FALSE])

                # H matrix
                H_k <- ((1 - alpha)/dim(de_train)[2]) * (de_train %*% t(de_train)) +
                    (alpha / dim(nu_train)[2]) * (nu_train %*% t(nu_train))
                h_k <- rowMeans(nu_train) ## sample from P
                
                h_kq <- rowMeans(de_train)  ## sample from q

                # compute theta hat for L2
                theta_hat <- solve(H_k + diag(k) * lambda_list[l]) %*% h_k
                
                ####################
                # theta_hat for L1
                # theta_hatL1 <- -1/lambda_list[l]*(alpha * h_k + (1-alpha)* h_kq)
                

                # compute objective function value -- squared loss -L2
                J <- (alpha / 2) * mean((t(theta_hat) %*% nu_test) ^ 2) +
                    ((1 - alpha) / 2) * mean((t(theta_hat) %*% de_test) ^ 2) -
                    mean(t(theta_hat) %*% nu_test)
                    
                 # compute objective function value -- absolute loss -L1
                #J <- alpha * mean(t(theta_hatL1) %*% nu_test) + 
                #(1-alpha) * mean(t(theta_hatL1) %*% de_test) 
               
				####################################
				
				
                # add to overall sum
                obj_score_sum <- obj_score_sum + J
            }
            # record score for this sigma-lambda combo
            scores[s, l] <- obj_score_sum / k
        }
    }

    # identify optimal sigma and lambda from score matrix
    minimum_index <- arrayInd(which.min(scores), dim(scores))
    sigma_min_index <- minimum_index[1]
    lambda_min_index <- minimum_index[2]
    opt_sigma <- sigma_list[sigma_min_index]
    opt_lambda <- lambda_list[lambda_min_index]

    # return
    return(list(
        opt_sigma = opt_sigma,
        opt_lambda = opt_lambda
    ))
}
