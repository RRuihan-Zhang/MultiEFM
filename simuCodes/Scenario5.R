gendata_simu1_robust <- function (seed = 1, nvec = c(100,300), p = 50, q = 3,
                                  qs= rep(2, length(nvec)), fac.type = c("gaussian", 'mvt'),
                                  err.type=c("gaussian", 'mvt'),
                                  rho = c(1,1), sigma2_eps=0.1, nu=1){
  # seed = 1; nvec = c(100,300); p = 50; q = 3
  # qs= rep(2, length(nvec))
  # rho = c(1,1); sigma2_eps=0.1;
  
  
  if(length(nvec)<2) stop("nvec must have at least two elements!")
  fac.type <- match.arg(fac.type)
  err.type <- match.arg(err.type)
  S <- length(nvec)
  require(MASS)
  factor_term_A <- rho[1]
  factor_term_B <- rho[2]
  set.seed(1)
  Blist <- list()
  set.seed(1)
  bmu0 <- matrix(rnorm(p*S), p, S)
  Ztmp <- matrix(rnorm(p * (q+qs[1])), p, (q+qs[1]))
  A <- qr(Ztmp)
  A1 <- qr.Q(A) %*% Diag(seq(q+qs[1], 1, length=q+qs[1]))
  A1 <- A1 %*% Diag(sign(A1[1, ]))  ## Fixed B0 and mu0 for each repeat.
  A0 <- A1[,1:q] *factor_term_A;
  B1 <- A1[,(q+1):(q+qs[1])] *factor_term_B
  t(A0) %*% A0; t(B1) %*% B1
  Blist[[1]] <- B1
  for(r in 2:S){
    set.seed(r)
    Ztmp <- matrix(rnorm(p * (qs[r])), p, (qs[r]))
    ## Project on A0
    Ztmp <- Ztmp -A0 %*% qr.solve(t(A0) %*% A0) %*% t(A0) %*% Ztmp
    A <- qr(Ztmp)
    A1 <- qr.Q(A) %*% Diag(seq(qs[r], 1, length=qs[r]))  * factor_term_B
    t(A1) %*% A1
    Blist[[r]]   <- A1 %*% Diag(sign(A1[1, ]))
  }
  
  
  
  set.seed(seed)
  Xlist <- list()
  Flist <- list()
  Hlist <- list()
  for(s in 1:S){
    
    n <- nvec[s]
    if(err.type == 'gaussian'){
      epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
    }else if(err.type == 'mvt'){
      epsi <- mvtnorm::rmvt(n, sigma =sigma2_eps* diag(p), df=nu)
      #matrix(rt(n*p, df=nu), n, p) *
    }
    
    if(fac.type == 'gaussian'){
      message("Factor is from gaussian distribution.")
      FH <- mvrnorm(n, mu = rep(0, q+qs[s]), Diag(rep(1,q+qs[s])))
    }else if(fac.type == 'mvt'){
      message("Factor is from t distribution.")
      FH <-  mvtnorm::rmvt(n, sigma =Diag(rep(1,q+qs[s])), df=nu)
    }
    
    Flist[[s]] <- FH[,1:q]
    Hlist[[s]] <- FH[,(q+1):(q+qs[s])]
    AB1 <- cbind(A0, Blist[[s]])
    Xlist[[s]]  <- matrix(bmu0[,s], nrow=n, ncol=p, byrow = TRUE) + FH %*% t(AB1) + epsi
  }
  
  return(list(Xlist = Xlist, mu0=bmu0, A0 = A0, Blist0 = Blist,
              Flist = Flist, Hlist = Hlist,  q=q, qs=qs))
}


selectFac.MultiEFM <- function(XList, q_max=20, qs_max=20, n_threads = 4,
                               subsample_rate=1, verbose=TRUE, seed=1){
  
  ### Select the number of factors
  # n_threads <- 8;subsample_rate <- 1
  # q_max <- 10
  
  
  library(irlba)
  if(q_max<3) stop("selectFac.MultiEFM: q_max must be greater than 3!")
  S <- length(XList)
  XList <- lapply(XList, scale, scale=F)
  nvec <- sapply(XList, nrow)
  set.seed(seed)
  hK_list <- lapply(XList, SK_subsample_fast_cpp_parallel, rate=subsample_rate, n_threads=n_threads)
  hK_pool <- hK_list[[1]] * nvec[1]/sum(nvec)
  for(ii in 2:length(hK_list)){
    hK_pool <- hK_pool + hK_list[[ii]] * nvec[ii]/sum(nvec)
  }
  eigvals_A <- irlba(hK_pool, nv=q_max)$d[1:q_max]
  eig.ratio.A <- eigvals_A[1:(q_max-1)]/ eigvals_A[2:q_max]
  ### two signal
  top_ratio <- sort(eig.ratio.A, decreasing = TRUE)[1:2]
  hq1 <- which(eig.ratio.A==top_ratio[1])
  hq2 <- which(eig.ratio.A==top_ratio[2])
  hq <- min(hq1, hq2)
  qs_up <- abs(hq1-hq2)
  
  if(verbose){
    message("Estimated study-shared number of factors are: hq = ", hq)
  }
  
  
  
  reslist2 <- MultiEFM(XList, q=hq, qs_vec= rep(qs_max, S), loads_method = 'Residual', fac_est_method='Huber', n_threads = 8)
  eigvals_B <- reslist2$eigvals_B
  eig.ratio.B <- lapply(eigvals_B, function(x)  x[1:(qs_max-1)]/ x[2:qs_max])
  ## substract shared factor number
  hqs_vec <- rep(NA, S)
  for(s in 1:S){
    hqs_vec[s] <- (which.max(eig.ratio.B[[s]]))
    if(TRUE){
      message("Estimated study-specified number of factor is: hq",s, ' =', hqs_vec[s])
    }
  }
  
  
  return(list(hq = hq, hqs_vec = hqs_vec, eigvals_A=eigvals_A, eig_ratio_A=eig.ratio.A,
              eigvals_B=eigvals_B, eig_ratio_B=eig.ratio.B))
}

library(MultiEFM)
library(MultiRFM)
source("helpfunc.R")
nu.vec <- c(2, 3, 10)
# i <- commandArgs(TRUE) %>% as.integer()
j <- 1
nu <- nu.vec[j]
p <- 400
nvec <- c(250,300); d <- 3; q <- 3; qs <- c(2,1); S <- length(nvec)
sigma2_eps <- 1
N <- 100
method_run <- 'All'
methodNames <- c("MultiEFM_residual", "MultiRFM")
n_methods <- length(methodNames)
qqsMat <- array(dim= c(3, 2, N) ); row.names(qqsMat) <- c("q", "q1", "q2")
colnames(qqsMat) <- methodNames
q_max <- 10; q_s_max <- 6

for(i in 1:N){
  #i <- 2
  print("MultiRFM")
  message("i = ", i)
  datList <- gendata_simu1_robust(seed=i, nvec=nvec, p=p, q=q, qs=qs,
                                  rho=c(5,5), err.type='mvt',
                                  sigma2_eps = sigma2_eps, nu=nu)
  # str(datList)
  t(datList$A0) %*% datList$Blist0[[1]]
  XList <- datList$Xlist;
  S <- length(XList)
  
  qlist <- selectFac.MultiEFM(XList, q_max=15, qs_max = 4)
  qqsMat[,1,i] <- c(qlist$hq, qlist$hqs_vec)
  
  qlist_rfm <- MultiRFM::selectFac.MultiRFM(XList, q_max=15, qs_max=4)
  qqsMat[,2,i] <- c(qlist_rfm$q, qlist_rfm$qs)
}

apply(qqsMat, c(1,2), mean, na.rm=T)
apply(qqsMat, c(1,2), sd, na.rm=T)
