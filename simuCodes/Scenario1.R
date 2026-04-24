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



# Test peformance of different methods when fac.type = c("gaussian") ------------------------------------
library(MultiEFM)
library(MultiRFM)
dir.name <- "scenario1"
if(!dir.exists(dir.name)){
  dir.create(dir.name)
}
setwd(dir.name)
source("helpfunc.R")
nu.vec <- c(2, 3, 10)
# i <- commandArgs(TRUE) %>% as.integer()
j <- 1 # <--- Change this to 1, 2, or 3 to reproduce different nu settings in Table 1
nu <- nu.vec[j]
p <- 400
nvec <- c(250,300); d <- 3; q <- 3;qs <- c(2,2); S <- length(nvec)
sigma2_eps <- 1
N <- 100 
method_run <- 'All'
methodNames <- c("MultiEFM_residual", "MultiRFM", "MSFA-CAVI", "MSFA-SVI", "MSFA", "BMSFA", "MultiEFM_Ortho")
n_methods <- length(methodNames)
metricList <- list(A_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                   F_tr =matrix(NA,N, n_methods), H_tr =matrix(NA,N, n_methods),
                   timeMat = matrix(NA, N, n_methods))
for(ii in 1:length(metricList)) colnames(metricList[[ii]]) <- methodNames

for(i in 1:N){
  #i <- 1
  print("MultiRFM")
  message("i = ", i)
  datList <- gendata_simu1_robust(seed=i, nvec=nvec, p=p, q=q, qs=qs,
                                  rho=c(4,4), err.type='mvt',
                                  sigma2_eps = sigma2_eps, nu=nu)
  # str(datList)
  XList <- datList$Xlist;
  
  reslist1 <- MultiEFM(XList, q=q, qs_vec=qs, loads_method = 'Residual', fac_est_method='Huber', n_threads = 8)
  metricList$timeMat[i,1] <- reslist1$time_use
  metricList$A_tr[i,1] <- trace_statistic_fun(reslist1$A, datList$A0)
  metricList$B_tr[i,1] <- trace_list_fun(reslist1$B, datList$Blist0)
  metricList$F_tr[i,1] <- trace_list_fun(reslist1$F, datList$Flist)
  metricList$H_tr[i,1] <- trace_list_fun(reslist1$H, datList$Hlist)
  
  
  library(MultiRFM)
  res <- MultiRFM(XList, q=q, qs= qs)
  metricList$timeMat[i,2] <- res$time_use
  metricList$A_tr[i,2] <- trace_statistic_fun(res$A, datList$A0)
  metricList$B_tr[i,2] <- trace_list_fun(res$B, datList$Blist0)
  metricList$F_tr[i,2] <- trace_list_fun(res$F, datList$Flist)
  metricList$H_tr[i,2] <- trace_list_fun(res$H, datList$Hlist)
  sapply(metricList, colMeans, na.rm=T)
  sapply(metricList, colSD)
  
  
  
  X_s <- lapply(XList, scale, scale=FALSE)
  hmu <- sapply(XList, colMeans)
  library(VIMSFA)
  ### MSFA-CAVI
  print("MSFA-CAVI")
  tic <- proc.time()
  cavi_est <- cavi_msfa(X_s,  K=q, J_s=qs)
  toc <- proc.time()
  time_cavi <- toc[3] - tic[3]
  hLam <- Reduce(cbind, cavi_est$mean_psi_s)
  hF_cavi <- hH_cavi <- list()
  for(s in 1:S){
    # s <- 1
    hF_cavi[[s]] <- t(Reduce(cbind, cavi_est$mean_f[[s]]))
    hH_cavi[[s]] <- t(Reduce(cbind, cavi_est$mean_l[[s]]))
  }
  
  metricList$timeMat[i, 3] <- time_cavi
  metricList$A_tr[i, 3] <- trace_statistic_fun(cavi_est$mean_phi, datList$A0)
  metricList$B_tr[i, 3] <- trace_list_fun(cavi_est$mean_lambda_s, datList$Blist0)
  metricList$F_tr[i, 3] <- trace_list_fun(hF_cavi, datList$Flist)
  metricList$H_tr[i, 3] <- trace_list_fun(hH_cavi, datList$Hlist)
  
  ### MSFA-SVI
  print("MSFA-SVI")
  tic <- proc.time()
  svi_est <- svi_msfa(X_s, K=q, J_s=qs, verbose = 0)
  toc <- proc.time()
  time_svi <- toc[3] - tic[3]
  hLam <- Reduce(cbind, svi_est$mean_psi_s)
  hF_cavi <- hH_cavi <- list()
  for(s in 1:S){
    # s <- 1
    hF_cavi[[s]] <- t(Reduce(cbind, svi_est$mean_f[[s]]))
    hH_cavi[[s]] <- t(Reduce(cbind, svi_est$mean_l[[s]]))
  }
  metricList$timeMat[i, 4] <- time_svi
  metricList$A_tr[i, 4] <- trace_statistic_fun(svi_est$mean_phi, datList$A0)
  metricList$B_tr[i, 4] <- trace_list_fun(svi_est$mean_lambda_s, datList$Blist0)
  metricList$F_tr[i, 4] <- trace_list_fun(hF_cavi, datList$Flist)
  metricList$H_tr[i, 4] <- trace_list_fun(hH_cavi, datList$Hlist)
  
  
  library(MSFA)
  maxIter.msfa <- 1000
  ## MSFA
  X_s <- lapply(XList, scale, scale=FALSE)
  hmu <- sapply(XList, colMeans)
  require(MSFA)
  start_value <- start_msfa(X_s =X_s, k = q, j_s = qs)
  try({
    tic <- proc.time()
    mle <-  ecm_msfa(X_s, start_value, trace = TRUE, nIt = maxIter.msfa)
    toc <- proc.time()
    time.use <- toc[3] - tic[3]
    message("Time is :", time.use)
    str(mle)
    metricList$timeMat[i,5] <- time.use
    metricList$A_tr[i,5] <- trace_statistic_fun(mle$Phi, datList$A0)
    metricList$B_tr[i,5] <- trace_list_fun(mle$Lambda_s, datList$Blist0)
    res_fac_msfa <- estimat.facs(XList, hmu, mle$Phi, mle$Lambda_s)
    metricList$F_tr[i,5] <-trace_list_fun(res_fac_msfa$F, datList$Flist)
    metricList$H_tr[i,5] <-trace_list_fun(res_fac_msfa$H, datList$Hlist)
  }, silent = TRUE)
  
  
  
  # BMSFA
  set.seed(1971)
  tic <- proc.time()
  control.list <- sp_msfa_control(nrun=maxIter.msfa, burn=round(maxIter.msfa/2))
  out10_1010 <- sp_msfa(X_s,  k = q,  j_s = qs, trace = TRUE, control = control.list)
  toc <- proc.time()
  time.use_msfa.bayes <- toc[3] - tic[3]
  hA <- apply(out10_1010$Phi, c(1, 2), median)
  hBlist <- lapply(out10_1010$Lambda, function(x) apply(x, c(1,2), median))
  hLam <- sapply(out10_1010$psi, function(x) apply(x, c(1,2), median))
  metricList$timeMat[i,6] <- time.use_msfa.bayes
  metricList$A_tr[i,6] <- trace_statistic_fun(hA, datList$A0)
  metricList$B_tr[i,6] <- trace_list_fun(hBlist, datList$Blist0)
  res_fac_msfa.bayes <- estimat.facs(XList, hmu, hA, hBlist)
  metricList$F_tr[i,6] <-trace_list_fun(res_fac_msfa.bayes$F, datList$Flist)
  metricList$H_tr[i,6] <-trace_list_fun(res_fac_msfa.bayes$H, datList$Hlist)
  
  
  reslist1 <- MultiEFM(XList, q=q, qs_vec=qs, loads_method = 'Ortho', fac_est_method='Huber', n_threads = 4)
  metricList$timeMat[i,7] <- reslist1$time_use
  metricList$A_tr[i,7] <- trace_statistic_fun(reslist1$A, datList$A0)
  metricList$B_tr[i,7] <- trace_list_fun(reslist1$B, datList$Blist0)
  metricList$F_tr[i,7] <- trace_list_fun(reslist1$F, datList$Flist)
  metricList$H_tr[i,7] <- trace_list_fun(reslist1$H, datList$Hlist)
  
  
  
  sapply(metricList, colMeans, na.rm=T)
  sapply(metricList, colSD)
  
  save(metricList, file=paste0(dir.name,"_","nsum", sum(nvec),"nu",paste0(nu, collapse = "_"),
                               "q", q ,"p",p,"N",N, '_metricList.rds'))
  
  
}




# Test peformance of different methods when fac.type = c("mvt") ------------------------------------

library(MultiEFM)
nu.vec <- c(2, 3, 10)
# i <- commandArgs(TRUE) %>% as.integer()
j <- 1 # <--- Change this to 1, 2, or 3 to reproduce different nu settings in Table 1

nu <- nu.vec[j]
p <- 400
nvec <- c(250,300); d <- 3; q <- 3; qs <- c(2,2); S <- length(nvec)
sigma2_eps <- 1
N <- 1
method_run <- 'All'
methodNames <- c("MultiEFM", "MultiRFM", "MSFA-CAVI", "MSFA-SVI", "MSFA", "BMSFA","MultiEFM_Ortho")
n_methods <- length(methodNames)
metricList <- list(A_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                   F_tr =matrix(NA,N, n_methods), H_tr =matrix(NA,N, n_methods),
                   timeMat = matrix(NA, N, n_methods))
for(ii in 1:length(metricList)) colnames(metricList[[ii]]) <- methodNames

for(i in 1:N){
  #i <- 1
  print("MultiRFM")
  message("i = ", i)
  datList <- gendata_simu1_robust(seed=i, nvec=nvec, p=p, q=q, qs=qs,
                                 rho=c(5,5), err.type='mvt', fac.type = 'mvt',
                                 sigma2_eps = sigma2_eps, nu=nu)
  # str(datList)
  XList <- datList$Xlist;

  reslist1 <- MultiEFM(XList, q=q, qs_vec=qs, fac_est_method='Huber', n_threads = 6)
  metricList$timeMat[i,1] <- reslist1$time_use
  metricList$A_tr[i,1] <- trace_statistic_fun(reslist1$A, datList$A0)
  metricList$B_tr[i,1] <- trace_list_fun(reslist1$B, datList$Blist0)
  metricList$F_tr[i,1] <- trace_list_fun(reslist1$F, datList$Flist)
  metricList$H_tr[i,1] <- trace_list_fun(reslist1$H, datList$Hlist)

  library(MultiRFM)
  res <- MultiRFM(XList, q=q, qs= qs)
  metricList$timeMat[i,2] <- res$time_use
  metricList$A_tr[i,2] <- trace_statistic_fun(res$A, datList$A0)
  metricList$B_tr[i,2] <- trace_list_fun(res$B, datList$Blist0)
  metricList$F_tr[i,2] <- trace_list_fun(res$F, datList$Flist)
  metricList$H_tr[i,2] <- trace_list_fun(res$H, datList$Hlist)




  X_s <- lapply(XList, scale, scale=FALSE)
  hmu <- sapply(XList, colMeans)
  library(VIMSFA)
  ### MSFA-CAVI
  print("MSFA-CAVI")
  tic <- proc.time()
  cavi_est <- cavi_msfa(X_s,  K=q, J_s=qs)
  toc <- proc.time()
  time_cavi <- toc[3] - tic[3]
  hLam <- Reduce(cbind, cavi_est$mean_psi_s)
  hF_cavi <- hH_cavi <- list()
  for(s in 1:S){
    # s <- 1
    hF_cavi[[s]] <- t(Reduce(cbind, cavi_est$mean_f[[s]]))
    hH_cavi[[s]] <- t(Reduce(cbind, cavi_est$mean_l[[s]]))
  }

  metricList$timeMat[i, 3] <- time_cavi
  metricList$A_tr[i, 3] <- trace_statistic_fun(cavi_est$mean_phi, datList$A0)
  metricList$B_tr[i, 3] <- trace_list_fun(cavi_est$mean_lambda_s, datList$Blist0)
  metricList$F_tr[i, 3] <- trace_list_fun(hF_cavi, datList$Flist)
  metricList$H_tr[i, 3] <- trace_list_fun(hH_cavi, datList$Hlist)

  ### MSFA-SVI
  print("MSFA-SVI")
  tic <- proc.time()
  svi_est <- svi_msfa(X_s, K=q, J_s=qs, verbose = 0)
  toc <- proc.time()
  time_svi <- toc[3] - tic[3]
  hLam <- Reduce(cbind, svi_est$mean_psi_s)
  hF_cavi <- hH_cavi <- list()
  for(s in 1:S){
    # s <- 1
    hF_cavi[[s]] <- t(Reduce(cbind, svi_est$mean_f[[s]]))
    hH_cavi[[s]] <- t(Reduce(cbind, svi_est$mean_l[[s]]))
  }
  metricList$timeMat[i, 4] <- time_svi
  metricList$A_tr[i, 4] <- trace_statistic_fun(svi_est$mean_phi, datList$A0)
  metricList$B_tr[i, 4] <- trace_list_fun(svi_est$mean_lambda_s, datList$Blist0)
  metricList$F_tr[i, 4] <- trace_list_fun(hF_cavi, datList$Flist)
  metricList$H_tr[i, 4] <- trace_list_fun(hH_cavi, datList$Hlist)


  maxIter.msfa <- 1000
  ## MSFA
  X_s <- lapply(XList, scale, scale=FALSE)
  hmu <- sapply(XList, colMeans)
  require(MSFA)
  start_value <- start_msfa(X_s =X_s, k = q, j_s = qs)
  try({
    tic <- proc.time()
    mle <-  ecm_msfa(X_s, start_value, trace = TRUE, nIt = maxIter.msfa)
    toc <- proc.time()
    time.use <- toc[3] - tic[3]
    message("Time is :", time.use)
    str(mle)
    metricList$timeMat[i,5] <- time.use
    metricList$A_tr[i,5] <- trace_statistic_fun(mle$Phi, datList$A0)
    metricList$B_tr[i,5] <- trace_list_fun(mle$Lambda_s, datList$Blist0)
    res_fac_msfa <- estimat.facs(XList, hmu, mle$Phi, mle$Lambda_s)
    metricList$F_tr[i,5] <-trace_list_fun(res_fac_msfa$F, datList$Flist)
    metricList$H_tr[i,5] <-trace_list_fun(res_fac_msfa$H, datList$Hlist)
  }, silent = TRUE)



  # BMSFA
  set.seed(1971)
  tic <- proc.time()
  control.list <- sp_msfa_control(nrun=maxIter.msfa, burn=round(maxIter.msfa/2))
  out10_1010 <- sp_msfa(X_s,  k = q,  j_s = qs, trace = TRUE, control = control.list)
  toc <- proc.time()
  time.use_msfa.bayes <- toc[3] - tic[3]
  hA <- apply(out10_1010$Phi, c(1, 2), median)
  hBlist <- lapply(out10_1010$Lambda, function(x) apply(x, c(1,2), median))
  hLam <- sapply(out10_1010$psi, function(x) apply(x, c(1,2), median))
  metricList$timeMat[i,6] <- time.use_msfa.bayes
  metricList$A_tr[i,6] <- trace_statistic_fun(hA, datList$A0)
  metricList$B_tr[i,6] <- trace_list_fun(hBlist, datList$Blist0)
  res_fac_msfa.bayes <- estimat.facs(XList, hmu, hA, hBlist)
  metricList$F_tr[i,6] <-trace_list_fun(res_fac_msfa.bayes$F, datList$Flist)
  metricList$H_tr[i,6] <-trace_list_fun(res_fac_msfa.bayes$H, datList$Hlist)


  reslist1 <- MultiEFM(XList, q=q, qs_vec=qs, loads_method = 'Ortho', fac_est_method='Huber', n_threads = 4)
  metricList$timeMat[i,7] <- reslist1$time_use
  metricList$A_tr[i,7] <- trace_statistic_fun(reslist1$A, datList$A0)
  metricList$B_tr[i,7] <- trace_list_fun(reslist1$B, datList$Blist0)
  metricList$F_tr[i,7] <- trace_list_fun(reslist1$F, datList$Flist)
  metricList$H_tr[i,7] <- trace_list_fun(reslist1$H, datList$Hlist)

  sapply(metricList, colMeans, na.rm=T)
  sapply(metricList, colSD)


}


sapply(metricList, colMeans, na.rm=T)
sapply(metricList, function(x) apply(x, 2, sd, na.rm = TRUE))
