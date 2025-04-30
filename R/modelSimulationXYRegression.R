#********************* SIMULATIONS FROM CLUSTERED MODELS **********************#
library(fda)
#library(scatterplot3d)
#library(shiny)
#library(shinydashboard)
library(MASS)
.fromClusteredSimulateX <- function(ncurves, p, mu, res_reg, K=2, prop=NULL,
                                    d=NULL, a=NULL, b=NULL, mu_div = 1.5){
  ####################################################################################
  #          #
  ####################################################################################
 
  
  # Class sizes
  N<-sum(ncurves)
  
  
  j<-sample(p, K)
  
  
  
  Q <- res_reg$Q1
  Wminv<-solve(res_reg$Wlist$W_m)
  #Wminv <- solve(sqrtm(inprod(res_reg$datax$basis, res_reg$datax$basis)))
  
  # Simulation
  S1<-vector(mode='list', length=K)
  
  for (i in 1:K)	{
    if(res_reg$model == "AKJBKQKDK") S1[[i]]<-tcrossprod(Wminv%*%Q[[i]]%*%sqrt(diag(c(a[i, 1:d[i]], rep(b[i], p-d[i])))))
    else if(res_reg$model == "AKJBQKDK") S1[[i]]<-tcrossprod(Wminv%*%Q[[i]]%*%sqrt(diag(c(a[i, 1:d[i]], rep(b, p-d[i])))))
    else if(res_reg$model == "AKBKQKDK") S1[[i]]<-tcrossprod(Wminv%*%Q[[i]]%*%sqrt(diag(c(rep(a[i], d[i]), rep(b[i], p-d[i])))))
    else if(res_reg$model == "AKBQKDK") S1[[i]]<-tcrossprod(Wminv%*%Q[[i]]%*%sqrt(diag(c(rep(a[i], d[i]), rep(b, p-d[i])))))
    else if(res_reg$model == "ABQKDK")S1[[i]]<-tcrossprod(Wminv%*%Q[[i]]%*%sqrt(diag(c(rep(a, d[i]), rep(b, p-d[i])))))
    else if(res_reg$model == "ABKQKDK")S1[[i]]<-tcrossprod(Wminv%*%Q[[i]]%*%sqrt(diag(c(rep(a, d[i]), rep(b[i], p-d[i])))))
  }
  
  
  cls<-X<-NULL
  clo <- rep(1, N)
  
  for (i in 1:K) for(j in 1:ncurves[i]){
    X<-rbind(X, mvrnorm(1, mu=mu[i, ]/mu_div,Sigma=S1[[i]]/2))
  }
  for (i in 1:K) cls<-c(cls, rep(i, ncurves[i]))
  
  ind <- sample(1:N, N)
  prms <- list(a=a, b=b, prop=prop, d=d, mu=mu, ncurves=ncurves)
  
  return(list(X = X,
              clx = cls,
              prms = prms))
  
}
# 
.fromClusteredSimulateY <- function(Xco, XW, fdX, res_reg, K=2, groupd=NULL){
  # note if nsplines = 70 a fourier basis will have 71 total functions (1 constant)
 
  Y <- NULL
  nsplines <- nrow(fdX$coefs)
  W <- inprod(fdX$basis, fdX$basis)
  for (i in 1:K){
    E <- rep(1, sum(groupd==i))
    
    coefk <- res_reg$gam[,,i][,1:6]
    # Cxk <- t(fdX$coefs[,groupd==i])
    Cxk <- Xco[groupd==i,]
    for(j in 1:sum(groupd==i)) {
      #Y<-rbind(Y, res_reg$gam[,,i][,7] +
      #          Cxk[j,]%*%W%*%t(coefk) +
      #         mvrnorm(n=1, rep(0, 6), res_reg$covy[,,i]))
      Y<-rbind(Y, res_reg$gam[,,i][,7] +
                 Cxk[j,]%*%XW%*%t(coefk) +
                 mvrnorm(n=1, rep(0, 6), res_reg$covy[,,i]))
      
    }
    
    
  }
  return(Y)
 
  
}

# ncurves is now given as a list accoring to the number K
genFromModelRegFD <- function(res_reg, ncurves = c(200, 300), mu_div = 1.5) {
  ####################################################################################
  mu <- res_reg$mu
  a <- res_reg$a
  b <- res_reg$b
  d <- res_reg$d
  
  # add 1 for constant
  coef <- .fromClusteredSimulateX(ncurves, 6, mu, res_reg, d=d, a=a, b=b, mu_div=mu_div)
  
  
  basis <- create.bspline.basis(rangeval = c(1,12), nbasis = 6)
  evalbas <- eval.basis(seq(1, 12, length.out = 24), basis)
  #print(dim(coef$X))
  #print(dim(coef$X %*% t(evalbas)))
  finaldata <- coef$X %*% t(evalbas)
  finaldata <- cbind(finaldata, coef$clx)
  bbasis <- create.bspline.basis(rangeval = c(1,12), nbasis = 6)
  
  fdx <- smooth.basis(argvals = seq(1, 12, length.out = 24),
                      y = t(finaldata[,1:24]),
                      fdParobj = bbasis)$fd
  
  coefY <- .fromClusteredSimulateY(coef$X,res_reg$Wlist$W, fdx, res_reg, groupd = coef$clx)
  
  
  basis <- create.bspline.basis(rangeval = c(13,24), nbasis = 6)
  evalbas <- eval.basis(seq(13, 24, length.out = 24), basis)
  finaldatay <- coefY %*% t(evalbas)
  bbasis <- create.bspline.basis(rangeval = c(13,24), nbasis = 6)
  
  fdy <- smooth.basis(argvals = seq(13, 24, length.out = 24),
                      y = t(finaldatay),
                      fdParobj = bbasis)$fd
  # fit the Y data
  return(list(xfd = fdx, yfd=fdy,groupd = coef$clx))
  
}
#Fourier basis
genFromModelRegFDF <- function(res_reg, ncurves = c(200, 300), mu_div = 1.5) {
  ####################################################################################
  mu <- res_reg$mu
  a <- res_reg$a
  b <- res_reg$b
  d <- res_reg$d
  
  # add 1 for constant
  coef <- .fromClusteredSimulateX(ncurves, 6, mu, res_reg, d=d, a=a, b=b, mu_div=mu_div)
  
  
  basis <- create.bspline.basis(rangeval = c(1,12), nbasis = 6)
  evalbas <- eval.basis(seq(1, 12, length.out = 24), basis)
  #print(dim(coef$X))
  #print(dim(coef$X %*% t(evalbas)))
  finaldata <- coef$X %*% t(evalbas)
  finaldata <- cbind(finaldata, coef$clx)
 # bbasis <- create.bspline.basis(rangeval = c(1,12), nbasis = 6)
  bbasis <-create.fourier.basis(rangeval = c(1,12), nbasis = 6)
  fdx <- smooth.basis(argvals = seq(1, 12, length.out = 24),
                      y = t(finaldata[,1:24]),
                      fdParobj = bbasis)$fd
  
  coefY <- .fromClusteredSimulateY(coef$X,res_reg$Wlist$W, fdx, res_reg, groupd = coef$clx)
  
  
  basis <- create.bspline.basis(rangeval = c(13,24), nbasis = 6)
  evalbas <- eval.basis(seq(13, 24, length.out = 24), basis)
  finaldatay <- coefY %*% t(evalbas)
  bbasis <-create.fourier.basis(rangeval = c(13,24), nbasis = 6)
#  bbasis <- create.bspline.basis(rangeval = c(13,24), nbasis = 6)
  
  fdy <- smooth.basis(argvals = seq(13, 24, length.out = 24),
                      y = t(finaldatay),
                      fdParobj = bbasis)$fd
  # fit the Y data
  return(list(xfd = fdx, yfd=fdy,groupd = coef$clx))
  
}
#plot the curves
plotRegFD <- function(gen_out) {
  groupd <- gen_out$groupd
  
  win.graph(width=9.5, height=6,pointsize=12)
  par(mfrow = c(2,1))
 
  plot.fd(gen_out$xfd, col = groupd, ylab = "X")
  plot.fd(gen_out$yfd, col = groupd, ylab = "Y")
}
