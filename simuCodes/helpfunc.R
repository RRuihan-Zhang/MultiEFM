
Diag<-function (vec){
  q <- length(vec)
  if (q > 1) {
    y <- diag(vec)
  }
  else {
    y <- matrix(vec, 1, 1)
  }
  return(y)
}


# Compared methods --------------------------------------------------------
mat2list <-function(z_int, nvec){
  
  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}

estimat.facs <- function(XList, hmu, hA, hB){
  FList <- list()
  HList <- list()
  p <- ncol(XList[[1]])
  S <- length(XList)
  for(s in 1:S){
    ns <- nrow(XList[[s]])
    Xms <- (XList[[s]]-matrix(hmu[,1], nrow=ns, ncol=p))
    FList[[s]] <- Xms %*% hA %*% qr.solve(t(hA)%*% hA)
    HList[[s]] <- (Xms - FList[[s]]%*% t(hA)) %*% hB[[s]] %*% qr.solve(t(hB[[s]])%*% hB[[s]])
  }
  return(list(F=FList, H=HList))
}

# MSFR.run <- function(XList, ZList, q, qs, maxIter=1e4, load.source=FALSE, dir.source=NULL){
#
#   # require(MSFA)
#   require(psych)
#   if(!load.source){
#     source(paste0(dir.source, "MSFR_main_R_MSFR_V1.R"))
#   }
#   #fa <- psych::fa
#   B_s <- ZList
#
#   X_s <- XList #
#   t1<- proc.time()
#   test <- start_msfa1(X_s, B_s, 5, k=q, j_s=qs, constraint = "block_lower2", method = "adhoc")
#   # EM_beta <- ecm_msfa(X_s, B_s, start=test, trace = FALSE, nIt=maxIter, constraint = "block_lower1")
#   EM_beta <- ecm_msfa1(X_s, B_s, start=test, trace = FALSE, nIt=maxIter)
#   t2<- proc.time()
#   EM_beta$Flist <- lapply(EM_beta$E_f, t)
#   EM_beta$Hlist <- lapply(EM_beta$E_l, t)
#   EM_beta$time.use <- t2[3] - t1[3]
#   return(EM_beta)
# }


#----------------------------------------------------------------------------------------
#  Spatial Multivariate Kendall' tau Matrix
#  Input:
#           x ------ n x p  data matrix (n: time dimension, p: cross-section)
#  Output:
#          TK ------ Sample Spatial Multivariate Kendall' tau Matrix
#----------------------------------------------------------------------------------------

SK<-function(X){
  p <- ncol(X)
  n <- nrow(X)
  TK <- matrix(0,p,p)
  for(i in 1:(n-2)){
    TT <- matrix(rep(X[i,],n-i),n-i,p,byrow = TRUE)-X[(i+1):n,]
    TT <- t(diag(1/diag(TT%*%t(TT)))%*%TT)%*%TT
    TK <- TK+TT
  }
  TT <- X[n-1,]-X[n,]
  TK <- TK+TT%*%t(TT)/sum(TT^2)
  TK <- 2/(n*(n-1))*TK
  return(TK)
}

#----------------------------------------------------------------------------------------
#  RTS Method
#  Input:
#           x ------ n x p  data matrix (n: time dimension, p: cross-section)
#           r ------ number of factors
#  Output:
#        Fhat ------ n x r  the estimated factor scores
#        Lhat ------ p x r  the estimated factor loadings
#----------------------------------------------------------------------------------------

RTS <- function(X,r){
  p <- ncol(X)
  n <- nrow(X)
  Khat <- SK(X)
  Lhat <- sqrt(p)*as.matrix(eigen(Khat)$vectors[,1:r])
  Fhat <- matrix(0,n,r)
  for (i in 1:n){
    Fhat[i,]=lm(X[i,]~Lhat-1)$coefficients
  }
  return(list(Fhat=Fhat,Lhat=Lhat))
}

### He Yong, 2022, JBES, Robust latent factor model without moment constraint.
RTS.run <- function(XList, q){
  nvec <- sapply(XList, nrow)
  X <- Reduce(rbind, XList)
  tic <- proc.time()
  fit <- RTS(X, r=q)
  toc <- proc.time()
  Flist <- mat2list(fit$Fhat, nvec)
  fit$Flist <- Flist
  fit$time.use <- toc[3] - tic[3]
  return(fit)
}

# Metrics -----------------------------------------------------------------
# mean.Fnorm <- function(x) sum(x^2)/ length(x)
normvec <- function(x) sqrt(sum(x^2)/ length(x))
colSD <- function(X) apply(X, 2, sd, na.rm=TRUE)
trace_statistic_fun <- function(H, H0){
  
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% qr.solve(t(H) %*% H) %*% t(H) %*% H0
  
  tr_fun(mat1) / tr_fun(t(H0) %*% H0)
  
}
trace_list_fun <- function(Hlist, H0list){
  trvec <- rep(NA, length(Hlist))
  for(i in seq_along(trvec)){
    trvec[i] <- trace_statistic_fun(Hlist[[i]], H0list[[i]])
  }
  return(mean(trvec))
}