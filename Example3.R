
source("comp_dist.R")
source("comp_median.R")
source("gaussian_kernel.R")
source("sliding_window.R")
source("ts_detect.R")
source("ts_plot.R")
source("sigma_lambda_grid_search.R")
source("RelULSIF.R")


## Figure 5
series = read.table("dros067.txt", sep="", header=T)
seriesT = t(series)
# baseline simulation 
#d <- ts_detect(method="PL",seriesT, window_size = 8, step = 5, alpha=0.5, k=50, n_folds = 5, 
#               thresh = 0.9, make_plot = TRUE)
#length(d$scores)


####### Windows size 8
## 11 genes
niter = 300. #300
result = NULL
resultM = matrix(0,niter,50) ## 50 if windows size is 8, it will be 14:63
x=seriesT

for (i in 1:niter){
  g = sample(3:5, 1)
  x1 = x[sample(1:11, g),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha=0.01, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.9, make_plot = FALSE)  # varying alpha 0.1, 0.5, 0.9
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
pdf(file="Fig5-001w8.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:63, out, type="l")  #1+step+size: length(score) + step+size
dev.off() 

for (i in 1:niter){
  g = sample(3:5, 1)
  x1 = x[sample(1:11, g),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha=0.1, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.9, make_plot = FALSE)  # varying alpha 0.1, 0.5, 0.9
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
pdf(file="Fig5-01w8.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:63, out, type="l")  #1+step+size: length(score) + step+size
dev.off() 

niter = 300
result = NULL
resultM = matrix(0,niter,50) ## 50 if windows size is 8, it will be 14:63
x=seriesT
for (i in 1:niter){
  g = sample(3:5, 1)
  x1 = x[sample(1:11, g),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha=0.5, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.9, make_plot = FALSE)  # varying alpha 0.1, 0.5, 0.9
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
pdf(file="Fig5-5w8.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:63, out, type="l")  #1+step+size: length(score) + step+size
dev.off() 

niter = 300
result = NULL
resultM = matrix(0,niter,50) ## 50 if windows size is 8, it will be 14:63  ## 53 if size 5, 11:63
x=seriesT
for (i in 1:niter){
  g = sample(3:5, 1)
  x1 = x[sample(1:11, g),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha=0.9, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.9, make_plot = FALSE)  # varying alpha 0.1, 0.5, 0.9
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
pdf(file="Fig5-9w8.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:63, out, type="l")  #1+step+size: length(score) + step+size
dev.off() 





###### Windows size = 5 Not on the paper, different window size


## 11 genes
niter = 500
result = NULL
resultM = matrix(0,niter,53) ## 50 if windows size=8, 14:63 
                              ## 54, win=4, 10:63
                             ## 53 if w=5, 11:63
                              ## 52 if window=6 12:63 
x=seriesT

for (i in 1:niter){
  g = sample(3:5, 1)
  x1 = x[sample(1:11, g),]
  d <- ts_detect(method="PL", x1, window_size = 5, alpha=0.01, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.9, make_plot = FALSE)  # varying alpha 0.1, 0.5, 0.9
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
pdf(file="Fig5-1.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(11:63, out, type="l")  #1+step+size: length(score) + step+size
dev.off() 

niter = 500
result = NULL
resultM = matrix(0,niter,53) ## 50 if windows size is 8, it will be 14:63
x=seriesT
for (i in 1:niter){
  g = sample(3:5, 1)
  x1 = x[sample(1:11, g),]
  d <- ts_detect(method="PL", x1, window_size = 5, alpha=0.5, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.9, make_plot = FALSE)  # varying alpha 0.1, 0.5, 0.9
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
pdf(file="Fig5-5.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(11:63, out, type="l")  #1+step+size: length(score) + step+size
dev.off() 

niter = 500
result = NULL
resultM = matrix(0,niter,53) ## 50 if windows size is 8, it will be 14:63
x=seriesT
for (i in 1:niter){
 # g = sample(3:5, 1)
  x1 = x[sample(1:11, 3),]
  d <- ts_detect(method="PD", x1, window_size = 5, alpha=0.9, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.9, make_plot = FALSE)  # varying alpha 0.1, 0.5, 0.9
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
#pdf(file="Fig5-9.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(11:63, out, type="l")  #1+step+size: length(score) + step+size
#dev.off() 





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
  
  
  

series = read.table("dros067.txt", sep="", header=T)
seriesT = t(series)
#d = ts_detect(seriesT, window_size = 2, alpha=0.5, k = 100, n_folds = 10,thresh = 0.9, make_plot = TRUE)

d <- ts_detect(method="PD",seriesT, window_size = 5, step = 5, alpha=0.5, k=50, n_folds = 5, 
               thresh = 0.9, make_plot = TRUE)


niter = 100
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

# Example function to perform in parallel
my_function <- function(x) {
  return( ts_detect(method="PL", seriesT, window_size = 2, alpha=x, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE))
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

