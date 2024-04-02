source("comp_dist.R")
source("comp_median.R")
source("gaussian_kernel.R")
source("sliding_window.R")
source("ts_detect.R")
source("ts_plot.R")
source("sigma_lambda_grid_search.R")
source("RelULSIF.R")

## Figure 7
## 4028 genes  subset of 400 genes

niter = 1000
resultM = matrix(0,niter,51)
result = NULL
load("drosophila.rda")
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 400),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha =0.01, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig6a01-001-400L.pdf")
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
  x1 = x[sample(1:4028, 400),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha =0.5, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig6a05-400.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 

niter = 1000
resultM = matrix(0,niter,51)
result = NULL

load("drosophila.rda")
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 400),]
  d <- ts_detect(method="PL", x1, window_size = 8, alpha =0.9, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig6a09-400.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 

niter = 1000
resultM = matrix(0,niter,51)
result = NULL

load("drosophila.rda")
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 400),]
  d <- ts_detect(method="PD", x1, window_size = 8, alpha =0.1, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig7a01-400-PD.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 



niter = 1000
resultM = matrix(0,niter,51)
result = NULL

load("drosophila.rda")
x=drosophila

for (i in 1:niter){
  x1 = x[sample(1:4028, 400),]
  d <- ts_detect(method="PD", x1, window_size = 8, alpha =0, step = 5, k = 50, n_folds = 5, thresh = 0.9, make_plot = FALSE)
  result[[i]] = d$change_points
  resultM[i,] = d$scores
}

pdf(file="Fig7a00-400-PD.pdf")
par(mfrow=c(2,1))
plot(table(unlist(result)))
out = apply(resultM, 2, mean)
plot(14:64, out, type="l")
dev.off() 


## Below, we use INSPECT method, change threshold from 50, 100, 150
library(InspectChangepoint)
load("drosophila.rda")

par(mfrow=c(3,1))
niter = 1000
result = NULL

for (i in 1:niter){
  x=drosophila
  x = x[sample(1:4028, 400),]
  result[[i]] = inspect(x, threshold=150)$changepoints[,1]
  # printPercentage(i, niter)
}
plot(table(unlist(result)),  ylab="Frequency", xlab = "time points")
#, xlab = "time points", ylab="Frequency", main="(C): threshold = 175")



