source("comp_dist.R")
source("comp_median.R")
source("gaussian_kernel.R")
source("sliding_window.R")
source("ts_detect.R")
source("ts_plot.R")
source("sigma_lambda_grid_search.R")
source("RelULSIF.R")


niter = 1000
resultM = matrix(0,niter,51)
result = NULL
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 40),]
  d <- ts_detect(method="PD", x1, window_size = 8, alpha =0.5, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig7a05.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 


niter = 1000
resultM = matrix(0,niter,51)
result = NULL
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 40),]
  d <- ts_detect(method="PD", x1, window_size = 8, alpha =0.1, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig7a01.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 

niter = 1000
resultM = matrix(0,niter,51)
result = NULL
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 40),]
  d <- ts_detect(method="PD", x1, window_size = 8, alpha =0.9, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig7a09.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 

niter = 1000
resultM = matrix(0,niter,51)
result = NULL
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 40),]
  d <- ts_detect(method="PD", x1, window_size = 8, alpha =0, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig8a05.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 

###### Two dimensional Multi-Normal, variance same
library(MASS)

mu1 = c(1,1)
sigma1 = matrix(c(1,0,0,1), 2,2)
mu2 = c(3,3)
sigma2 = matrix(c(1,0,0,1), 2,2)
series <- rbind(mvrnorm(200, mu1, sigma1), mvrnorm(200, mu2, sigma2))
#d <- ts_detect(method="PL", t(series), window_size = 10, step = 20, n_folds = 5, alpha=0.9, k=50, thresh = 0.9, make_plot = TRUE)

niter = 100
#result = NULL
par(mfrow=c(3,1))
for (alpha in c(0.1, 0.5, 0.9)){
  result = NULL
  for (i in 1:niter){
    d <- ts_detect(method="PL", t(series), window_size = 10, step = 20, alpha, k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE)
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}

###### Two dimensional Multi-Normal, variance different
library(MASS)
mu1 = c(1,1)
sigma1 = matrix(c(1,0.5,0.5,1), 2,2)
mu2 = c(15,13)
sigma2 = matrix(c(4,0.6,0.6,4), 2,2)

mu3 = c(30,40)
sigma3 = matrix(c(2,0.1,0.1,3), 2,2)

series <- rbind(mvrnorm(200, mu1, sigma1), mvrnorm(200, mu2, sigma2), mvrnorm(200, mu3, sigma3) )

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

# Example function to perform in parallel
my_function <- function(x) {
  return( ts_detect(method="PL", t(series), window_size = 10, alpha=x, step = 20, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE))
}
# Generate data
niter = 100
par(mfrow=c(3,1))
for (alpha in c(0.1, 0.5, 0.9)){
  # Perform computation in parallel
  result <- foreach(i=1:niter, .combine = c) %dopar% {
    my_function(alpha)$change_points
  }
  plot(table(unlist(result)))
}
# Stop parallel backend
stopCluster(cl)
# Result



niter = 100
#result = NULL
par(mfrow=c(3,1))
for (alpha in c(0.1, 0.5, 0.9)){
  result = NULL
  for (i in 1:niter){
    d <- ts_detect(method="PL", t(series), window_size = 10, step = 20, alpha, k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE)
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}

niter = 100
#result = NULL
par(mfrow=c(3,1))
for (alpha in c(0.1, 0.5, 0.9)){
  result = NULL
  for (i in 1:niter){
    d <- ts_detect(method="PD", t(series), window_size = 10, step = 20, alpha, k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE)
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}
