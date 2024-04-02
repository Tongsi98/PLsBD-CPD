source("comp_dist.R")
source("comp_median.R")
source("gaussian_kernel.R")
source("sliding_window.R")
source("ts_detect.R")
source("ts_plot.R")
source("sigma_lambda_grid_search.R")
source("RelULSIF.R")

## Figure 6
## 4028 genes  subset of 40 genes
 
niter = 1000
resultM = matrix(0,niter,51)
result = NULL
load("drosophila.rda")
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 40),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha =0.01, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig6a01D.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
sds = apply(resultM, 2, sd)
plot(14:64, out, type="l")
dev.off() 


niter = 1000
resultM = matrix(0,niter,51)
result = NULL
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 40),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha =0.5, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig6a05L.pdf")
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
  d <- ts_detect(method="PL", x1, window_size = 8, alpha =0.9, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig6a09L.pdf")
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



## Figure 6
series = read.table("dros067.txt", sep="", header=T)
seriesT = t(series)
# baseline simulation 
#d <- ts_detect(method="PL",seriesT, window_size = 8, step = 5, alpha=0.5, k=50, n_folds = 5, 
#               thresh = 0.9, make_plot = TRUE)
#length(d$scores)

## 11 genes
niter = 1000
result = NULL
resultM = matrix(0,niter,53)
x=seriesT

for (i in 1:niter){
  #g = sample(1:11, 1)
  x1 = x[sample(1:11, 5),]
  d <- ts_detect(method="PL", x1, window_size = 5, alpha=0.9, step = 5, k = 50, 
                 n_folds = 5, thresh = 0.9, make_plot = FALSE)  # varying alpha 0.1, 0.5, 0.9
  result[[i]] = d$change_points
  resultM[i,] = d$scores
  
}
#pdf(file="Fig1.pdf")
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
