source("comp_dist.R")
source("comp_median.R")
source("gaussian_kernel.R")
source("sliding_window.R")
source("ts_detect.R")
source("ts_plot.R")
source("sigma_lambda_grid_search.R")
source("RelULSIF.R")



series = read.table("dros067.txt", sep="", header=T)
seriesT = t(series)
# baseline simulation 
d <- ts_detect(method="PL",seriesT, window_size = 8, step = 5, alpha=0.5, k=50, n_folds = 5, 
               thresh = 0.9, make_plot = TRUE)
length(d$scores)

## 11 genes
niter = 1000
resultM = matrix(0,niter,50)
x=seriesT
#  result = NULL
for (i in 1:niter){
  x1 = x[sample(1:11, 3),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha=0.5, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.8, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
pdf(file="Figure5a.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:63, out, type="l")
dev.off() 


## Pearson uLSIF
niter = 10
resultM = matrix(0,niter,51)
x=drosophila
#  result = NULL
for (i in 1:niter){
  x1 = x[sample(1:4028, 400),]
  d <- ts_detect(method="PD", x1, window_size = 8, alpha=0, step = 5, k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE)
  #result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
pdf(file="Test1.pdf")
# plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 

niter = 500
par(mfrow=c(3,1))

for (alpha in c(0.1, 0.5, 0.9)){
  result = NULL
  #resultM = matrix(0,niter,66)
  for (i in 1:niter){
    seriesT1 = seriesT[sample(1:11,4),]
    d <- ts_detect(method="PL", seriesT1, window_size = 8, alpha, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
    result[[i]] = d$change_points
    #resultM[,i] = d$scores
  }
  plot(table(unlist(result)))
  #out = apply(resultM, 1, mean)
  #plot(1:50, out, type="l")
}

niter = 100
par(mfrow=c(3,1))
for (alpha in c(0.1, 0.5, 0.9)){
  result = NULL
  for (i in 1:niter){
    d <- ts_detect(method="PD", seriesT, window_size = 2, alpha, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}


## 4028 genes, #32, 42, 60 change points
load("drosophila.rda")

niter = 100
par(mfrow=c(3,1))
for (alpha in c(0.1, 0.5, 0.9)){
  x=drosophila
  result = NULL
  for (i in 1:niter){
    x1 = x[sample(1:4028, 400),]
    # d <- ts_detect(method="PL", x1, window_size = 8, alpha, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
    d <- ts_detect(method="PD", x1, window_size = 8, alpha, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
    
    result[[i]] = d$change_points
  }
  plot(table(unlist(result)))
}

## Pearson uLSIF
niter = 500
resultM = matrix(0,niter,51)
x=drosophila
#  result = NULL
for (i in 1:niter){
  x1 = x[sample(1:4028, 400),]
  d <- ts_detect(method="PD", x1, window_size = 8, alpha=0, step = 5, k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE)
  #result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
# plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(7:60, out, type="l")






#31, 41, 59 change points

#ts_detect <- function(ts, window_size = 5, step = NULL,alpha=0.5, k = 50, n_folds = 5,thresh = 0.9, make_plot = TRUE)



#install.packages("foreach")
#install.packages("doParallel")
library(foreach)
library(doParallel)
# Set up parallel backend
# Adjust the number of cores as per your system
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

# function to perform in parallel
my_function <- function(x) {
  return( ts_detect(method="PL", seriesT, window_size = 8, alpha=x, step = 5, k = 50, n_folds = 5, thresh = 0.8, make_plot = FALSE))
}
# Generate data
niter = 1000
#for (alpha in c(0.5)){
  # Perform computation in parallel
  alpha = 0.5
  result1 = NULL
  result2 = NULL
  
  foreach(i=1:niter, .combine = c) %dopar% {
    result1 = my_function(alpha)$change_points
     result2 = my_function(alpha)$scores
  }
  result2 <- foreach(i=1:niter, .combine = c) %dopar% {
    my_function(alpha)$scores
  }
  pdf(file="Test10.pdf")
  par(mfrow=c(2,1))
  plot(table(unlist(result1)))
  resultM = matrix(result2,50, niter, byrow=T)
  out = apply(resultM, 1, mean)
  plot(14:63, out, type="l")
  dev.off() 
#}

# Stop parallel backend
stopCluster(cl)
# Result

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

# Example function to perform in parallel
my_function <- function(x) {
  return( ts_detect(method="PD", seriesT, window_size = 2, alpha=x, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE))
}
# Generate data
niter = 100
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

