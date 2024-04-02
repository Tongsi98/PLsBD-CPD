source("comp_dist.R")
source("comp_median.R")
source("gaussian_kernel.R")
source("sliding_window.R")
source("ts_detect.R")
source("ts_plot.R")
source("sigma_lambda_grid_search.R")
source("RelULSIF.R")


generate_data <- function(n) {
  y <- numeric(n)
  e <- rnorm(n, mean = 0, sd = 1)
  mu_N <- 0  # Initialize mu_N before the loop
  for (t in 3:n) {
    if (t %% 100 == 1) {
      N <- (t - 1) %/% 100 + 1
      #mu_N <- ifelse(N == 1, 0, mu_N + N/2)
      
      mu_N <- mu_N + 2
    }
    y[t] <- 0.6 * y[t-1] - 0.5 * y[t-2] + e[t] + mu_N
  }
  return(y)
}

## Figure 3
series <- generate_data(1000)  #step =length(series)/20
d <- ts_detect(method="PL",series, window_size = 5, step=50, alpha=0.5, k=50, n_folds = 5, 
               thresh = 0.9, make_plot = TRUE)
d$change_points
abline(v = c(101, 201, 301, 401, 501,601, 701, 801, 901), col = "blue")

d <- ts_detect(method="PD",series, window_size = 5, step=50, alpha=0.5, k=50, n_folds = 5, 
               thresh = 0.9, make_plot = TRUE)
d$change_points
abline(v = c(101, 201, 301, 401, 501,601, 701, 801, 901), col = "blue")

## Pearson only
d <- ts_detect(method="PD",series, window_size = 5, step=50, alpha=0, k=50, n_folds = 5, 
               thresh = 0.9, make_plot = TRUE)
d$change_points
abline(v = c(101, 201, 301, 401, 501,601, 701, 801, 901), col = "blue")



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

TPL = matrix(0, 20, 5)
j=1
for (alpha in seq(0.1, 0.9, 0.1)){
for (i in 1:20){
series <- generate_data(500)
d <- ts_detect(method="PL",series, window_size = 5, alpha, k=50, n_folds = 5, 
               thresh = 0.9, step = length(series)/10, make_plot = FALSE)
Predictedcpd = d$change_points
cpd1 = c(101, 201, 301, 401)
cpd = c(cpd1-1, cpd1, cpd1+1)
TPL[i,j]= length(intersect(cpd1, Predictedcpd))
}
  j= j+1
}
TPL/4
TPLratio = TPL/4
apply(TPLratio, 2,mean)



  niter = 50
#result = NULL
par(mfrow=c(3,1))
for (alpha in c(0.1)){
  result = NULL
  for (i in 1:niter){
    d <- ts_detect(method="PL", series, window_size = 5, step = length(series)/10, alpha, 
                   k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE)
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}

niter = 50
#result = NULL
par(mfrow=c(3,1))
for (alpha in c(0.1,0.5, 0.9)){
  result = NULL
  for (i in 1:niter){
    d <- ts_detect(method="PD", series, window_size = 5, step = length(series)/10, alpha, 
                   k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE)
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}

niter = 50
#result = NULL
par(mfrow=c(3,1))
for (alpha in c(0)){
  result = NULL
  for (i in 1:niter){
    d <- ts_detect(method="PL", series, window_size = 5, step = length(series)/10, alpha, 
                   k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE)
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}

TPD = matrix(0, 20, 6)
j=1
for (alpha in c(0, 0.1, 0.3, 0.5, 0.7, 0.9)){
  for (i in 1:20){
    series <- generate_data(500)
    d <- ts_detect(method="PD",series, window_size = 5, alpha, k=50, n_folds = 5, 
                   thresh = 0.9, step = length(series)/10, make_plot = FALSE)
    Predictedcpd = d$change_points
    cpd1 = c(101, 201, 301, 401)
    cpd = c(cpd1-1, cpd1, cpd1+1)
    TPD[i,j]= length(intersect(cpd, Predictedcpd))
  }
  j= j+1
}
TPDratio = TPD/4
apply(TPDratio, 2,mean)

