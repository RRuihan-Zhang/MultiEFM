# Maintenance Cheat Sheet (MultiEFM)
# ----------------------------------
# 1. Generate documentation and NAMESPACE
# devtools::document()

# 2. Build the full website (Docs folder)
# pkgdown::build_site()

# 3. Build specific tutorial articles
# pkgdown::build_article("simu_low_dim")

# 4. Final check for CRAN compatibility
# devtools::check(args = c('--no-manual', '--as-cran'))

# -------------------------------------------------------------------------
# Helper functions (Hidden from User Manual)
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd


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
mat2list <- function(z_int, nvec){

  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){

    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}
add_identifiability_fast <- function(H, B){

  n <- nrow(H)
  q <- ncol(H)
  svdH <- svd(H, nu = q, nv = q)
  tB <- B %*% Diag(svdH$d[1:q]) %*% t(svdH$v) /sqrt(n)
  e_max <- apply(tB, 2, function(v) v[which.max(abs(v))])
  signB1 <- sign(e_max )
  tH <-  sqrt(n) * svdH$u %*% Diag(signB1)
  tB <- tB %*% Diag(signB1)

  return(list(H = tH, B = tB))
}

# -------------------------------------------------------------------------
# Main Exported Functions
# -------------------------------------------------------------------------

#' Fit the high-dimensional multi-study elliptical factor model
#'
#' @description
#' Fit the high-dimensional multi-study elliptical factor model which learns latent features and accounts for the heterogeneity among sources.
#'
#' @param XList A length-S list, where each component represents a data matrix for a specific study.
#' @param q an integer, specify the number of study-shared factors.
#' @param qs_vec an integer vector with length S, specify the number of study-specified factors for each study.
#' @param loads_method a character string specifying the method for loading estimation. Options are "Residual" or "Ortho".
#' @param fac_est_method a character string specifying the method for factor estimation. Options are "Huber" or "LS" (Least Squares).
#' @param n_threads an integer, specify the number of threads for parallel computing; default as 4.
#' @param sample_pairs an integer, specify the number of sample pairs used for estimation; default as 0.
#' @param epsObj a positive value, tolerance for convergence, default as 1e-6.
#' @param maxIter the maximum iteration of the algorithm. The default is 20.
#' @param verbose a logical value, whether to output the information in iteration.
#' @param seed an integer, specify the random seed for reproducibility in initialization; default as 1.
#'
#' @return return a list including the following components:
#' \item{A}{the estimated loading matrix corresponding to study-shared factors;}
#' \item{F}{a list composed by the estimation of study-shared factor matrix for each study;}
#' \item{B}{a list composed by the loading matrices corresponding to the study-specified factors;}
#' \item{H}{a list composed by the estimation of study-specified factor matrix for each study.}
#'
#' @examples
#' \dontrun{
#' p <- 100
#' nvec <- c(150, 200)
#' qs <- c(2, 2)
#' datList <- gendata_simu_robust(seed=1, nvec=nvec, p=p, q=3, qs=qs, rho=c(5,5),
#'                                err.type='mvt', sigma2_eps = 1, nu=3)
#' XList <- datList$Xlist
#' res <- MultiEFM(XList, q=3, qs_vec=qs)
#' str(res)
#' }
#'
#' @export
#' @useDynLib MultiEFM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm sd
#' @importFrom MASS mvrnorm

MultiEFM <- function(XList, q, qs_vec, loads_method= c("Residual", "Ortho"), fac_est_method=c('Huber', 'LS'),
                     n_threads = 4,  sample_pairs=0, epsObj=1e-6, maxIter=20,
                     verbose=TRUE, seed=1){


  ## ---- XList ----
  if(!is.list(XList))
    stop("XList must be a list of numeric matrices.")

  S <- length(XList)
  if(S == 0)
    stop("XList must contain at least one matrix.")

  ## 检查每个矩阵
  for(i in seq_len(S)){
    Xi <- XList[[i]]

    if(!is.matrix(Xi))
      stop(sprintf("XList[[%d]] is not a matrix.", i))

    if(!is.numeric(Xi))
      stop(sprintf("XList[[%d]] must be numeric.", i))

    if(anyNA(Xi))
      stop(sprintf("XList[[%d]] contains NA. Please remove or impute.", i))

    if(nrow(Xi) < 2)
      stop(sprintf("XList[[%d]] must have at least 2 rows.", i))

    if(i > 1 && ncol(Xi) != ncol(XList[[1]]))
      stop("All matrices in XList must have the same number of columns.")
  }

  p <- ncol(XList[[1]])


  ## ---- q ----
  if(length(q) != 1 || !is.numeric(q) || q <= 0 || q != round(q))
    stop("q must be a positive integer.")

  if(q >= p)
    stop("q must be strictly smaller than the number of variables p.")


  ## ---- qs_vec ----
  if(length(qs_vec) != S)
    stop("qs_vec must have the same length as XList.")

  if(any(qs_vec <= 0) || any(qs_vec != round(qs_vec)))
    stop("All entries in qs_vec must be positive integers.")

  if(any(qs_vec + q >= p))
    stop("For each study, q + q_s must be strictly smaller than p.")


  ## ---- fac_est_method ----
  fac_est_method <- match.arg(fac_est_method)
  load_B_method <- match.arg(loads_method)

  nvec  <- sapply(XList, nrow)
  n_min <- min(nvec); n_max <- max(nvec)
  subsample_rate <- 1
  if(sample_pairs==0){
    sample_method <- 'subsample'
    subsample_rate <- 1
  }else{
    sample_method <- "sample_pairs"
    max_pairs <- n_min * (n_min - 1) / 2

    ## ---- sample_pairs ----
    if (length(sample_pairs) != 1 ||
        !is.numeric(sample_pairs) ||
        !is.finite(sample_pairs) ||
        sample_pairs != round(sample_pairs) ||
        sample_pairs < 0) {

      stop("sample_pairs must be a single finite integer >= 0.")
    }

    if (sample_pairs > max_pairs) {
      stop(paste0("sample_pairs must be <= ", max_pairs, "."))
    }
  }

  ## ---- n_threads ----
  if(length(n_threads) != 1 || n_threads <= 0 || n_threads != round(n_threads))
    stop("n_threads must be a positive integer.")





  ## ---- epsObj ----
  if(length(epsObj) != 1 || epsObj <= 0)
    stop("epsObj must be a positive number.")


  ## ---- maxIter ----
  if(length(maxIter) != 1 || maxIter <= 0 || maxIter != round(maxIter))
    stop("maxIter must be a positive integer.")


  ## ---- verbose ----
  if(!is.logical(verbose) || length(verbose) != 1)
    stop("verbose must be TRUE or FALSE.")


  ## ---- seed ----
  if(length(seed) != 1 || seed != round(seed))
    stop("seed must be an integer.")



  XList <- lapply(XList, scale, scale=F)
  factor_method <- switch(fac_est_method, Huber=2, LS=1)
  sample_method <- switch(sample_method, sample_pairs=2, subsample=1)
  set.seed(seed)
  tic <- proc.time()
  if(load_B_method=='Residual'){
    fit_km <- kenMultiRFM_cpp(XList, q, qs_vec=qs_vec, factor_method = factor_method, n_threads = n_threads,
                              sampling_method = sample_method, sample_pairs=sample_pairs,
                              subsample_rate=subsample_rate, epsObj=epsObj, maxIter=maxIter,
                              verbose=verbose)
  }else if(load_B_method=='Ortho'){
    fit_km <- kenMultiRFM_project_cpp(XList, q, qs_vec=qs_vec, factor_method = factor_method, n_threads = n_threads,
                                      subsample_rate=subsample_rate, epsObj=epsObj, maxIter=maxIter,
                                      verbose=verbose)
  }

  toc <- proc.time()
  fit_km$time_use <- toc[3] - tic[3]

  ### add identifiability conditions
  nvec <- sapply(XList, nrow)
  AFList <- add_identifiability_fast(do.call(rbind,fit_km$F), fit_km$A)
  fit_km$A <- AFList$B
  fit_km$F <- mat2list(AFList$H, nvec = nvec)
  rm(AFList)
  for(s in 1:S){
    BHList <- add_identifiability_fast(fit_km$H[[s]], fit_km$B[[s]])
    fit_km$B[[s]] <- BHList$B
    fit_km$H[[s]] <- BHList$H
  }
  rm(BHList)
  return(fit_km)
}

#' Select the number of factors for MultiEFM
#'
#' @description
#' Automatically estimate the number of study-shared factors and study-specific factors based on the eigenvalue ratio method.
#'
#' @param XList A length-S list, where each component represents a data matrix for a specific study.
#' @param q_max an integer, specify the maximum number of study-shared factors to consider; default as 20.
#' @param qs_max an integer, specify the maximum number of study-specified factors to consider; default as 20.
#' @param qd an optional integer, specify the true number of combined factors if known in advance; default is NULL.
#' @param n_threads an integer, specify the number of threads for parallel computing; default as 4.
#' @param sample_pairs an integer, specify the number of sample pairs used for estimation; default as 0.
#' @param verbose a logical value, whether to output the estimation information.
#' @param seed an integer, specify the random seed for reproducibility; default as 1.
#'
#' @return return a list including the following components:
#' \item{hq}{the estimated number of study-shared factors;}
#' \item{hqs_vec}{an integer vector representing the estimated number of study-specific factors for each study;}
#' \item{eigvals_A}{the eigenvalues used for determining the shared factor number;}
#' \item{eig_ratio_A}{the eigenvalue ratios used for determining the shared factor number;}
#' \item{eigvals_B}{a list of eigenvalues used for determining the specific factor numbers;}
#' \item{eig_ratio_B}{a list of eigenvalue ratios used for determining the specific factor numbers.}
#'
#' @importFrom irlba irlba
#' @export



selectFac.MultiEFM <- function(XList, q_max=20, qs_max=20, qd=NULL,   n_threads = 4,
                               sample_pairs=0, verbose=TRUE, seed=1){

  ### Select the number of factors
  # n_threads <- 8;subsample_rate <- 1
  # q_max <- 10


  #library(irlba)
  if(q_max<3) stop("selectFac.MultiEFM: q_max must be greater than 3!")
  S <- length(XList)
  XList <- lapply(XList, scale, scale=F)
  nvec <- sapply(XList, nrow)
  set.seed(seed)
  if(sample_pairs==0){
    hK_list <- lapply(XList, SK_subsample_fast_cpp_parallel, rate=1, n_threads=n_threads)
  }else{
    hK_list <- list()
    for(s in 1:S){

      if(sample_pairs == 1){
        num_sample_pairs <- as.integer(nvec[s] *log(nvec[s]))
        hK_list[[s]] <- SK_pair_MC_cpp_parallel(XList[[s]],  m=num_sample_pairs, n_threads=n_threads)
      }else{
        hK_list[[s]] <- SK_pair_MC_cpp_parallel(XList[[s]],  m=sample_pairs, n_threads=n_threads)
      }

    }
  }

  if(is.null(qd)){
    hK_pool <- hK_list[[1]] * nvec[1]/sum(nvec)
    if(S>=2){
      for(ii in 2:length(hK_list)){
        hK_pool <- hK_pool + hK_list[[ii]] * nvec[ii]/sum(nvec)
      }
    }
    eigvals_A <- irlba(hK_pool, nv=q_max)$d[1:q_max]
    eig.ratio.A <- eigvals_A[1:(q_max-1)]/ eigvals_A[2:q_max]
    ### two signal
    top_ratio <- sort(eig.ratio.A, decreasing = TRUE)[1:2]
    hq1 <- which(eig.ratio.A==top_ratio[1])
    hq2 <- which(eig.ratio.A==top_ratio[2])
    hqd <- max(hq1, hq2) ### Use aggregated information, is more robust.
  }else{
    eigvals_A <- NULL
    eig.ratio.A <- NULL
    hqd <- qd
  }

  eigList_B <- lapply(hK_list, function(K) irlba(K, nv=hqd))
  p <- ncol(XList[[1]])
  M <- matrix(0, p, p)
  for(s in 1:S){
    M <- M + eigList_B[[s]]$u[,1:hqd] %*% t(eigList_B[[s]]$u[,1:hqd])
  }
  M <- M%*%M - M
  eigvals_M <- irlba(M, nv=hqd)$d
  threshold <- (S-1)^2
  hq <- sum(eigvals_M>threshold)

  if(verbose){
    message("Estimated study-shared number of factors are: hq = ", hq)
  }

  reslist2 <- MultiEFM(XList, q=hq, qs_vec= rep(qs_max, S), loads_method = 'Residual',
                       fac_est_method='LS', sample_pairs = sample_pairs, n_threads = n_threads)
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

