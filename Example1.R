source("comp_dist.R")
source("comp_median.R")
source("gaussian_kernel.R")
source("sliding_window.R")
source("ts_detect.R")
source("ts_plot.R")
source("sigma_lambda_grid_search.R")
source("RelULSIF.R")



## Figure 1
series1 <- c(
  rnorm(100, mean = 0, sd = 1),
  # rnorm(25, mean = 8, sd = 1),
  rnorm(100, mean = 10, sd = 1),
  # rnorm(25, mean = 1, sd = 0.8),
  rnorm(100, mean = -5, sd = 1),
  # rnorm(100, mean = -5, sd = 0.2),
  # rnorm(50, mean = -2.5, sd = 0.4),
  rnorm(100, mean = 10, sd = 1)
)
# Experiment 1
d <- ts_detect(method="PL", series1, window_size = 5, alpha=0.5, k = 50, n_folds = 5,thresh = 0.9, make_plot = TRUE)
d$step
d$change_points
abline(v = c(101, 201, 301), col = "blue")


d <- ts_detect(method="PD", series1, window_size = 5, alpha=0.5, k = 50, n_folds = 5,thresh = 0.9, make_plot = TRUE)
d$change_points
abline(v = c(101, 201, 301), col = "blue")


d <- ts_detect(method="PD", series1, window_size = 10, alpha=0, k = 50, n_folds = 5,thresh = 0.9, make_plot = TRUE)
d$change_points
abline(v = c(101, 201, 301), col = "blue")


# Below are the code for some extra test not on the paper

#install.packages("foreach")
#install.packages("doParallel")
library(foreach)
library(doParallel)
# Set up parallel backend
# Adjust the number of cores as per your system
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)


niter = 50
#result = NULL
par(mfrow=c(3,1))
for (alpha in c(0.1)){
  result = NULL
  for (i in 1:niter){
    #d <- ts_detect(method="PL", series1, window_size = 5, alpha, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
    d <- ts_detect(method="PL", series1, window_size = 5, alpha, step = length(series1)/10, k = 50, n_folds = 5,thresh = 0.9, make_plot = FALSE)
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}



series2 <- c(
  rnorm(100, mean = 0, sd = 1),
  # rnorm(25, mean = 8, sd = 1),
  rnorm(100, mean = 10, sd = 2),
  # rnorm(25, mean = 1, sd = 0.8),
  rnorm(100, mean = -5, sd = 3),
  # rnorm(100, mean = -5, sd = 0.2),
  # rnorm(50, mean = -2.5, sd = 0.4),
  rnorm(100, mean = 10, sd = 1)
)

par(mfrow=c(2,2))
d <- ts_detect(method="PL", series2, window_size = 5, alpha=0.5, k = 50, n_folds = 5,thresh = 0.9, make_plot = TRUE)
d$change_points

d <- ts_detect(method="PD", series2, window_size = 5, alpha=0.5, k = 50, n_folds = 5,thresh = 0.9, make_plot = TRUE)
d$change_points





