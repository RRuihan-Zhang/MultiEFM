# Generate simulated data -------------------------------------------------
#' Generate Simulated Data for Robust Multi-Study Elliptical Factor Model
#'
#' @description
#' Generate simulated data for high-dimensional multi-study factor analysis under 
#' robust elliptical distributions (such as Multivariate t-distribution). The data follows 
#' a factor model with common factors (shared across studies) and study-specific factors, 
#' plus noise.
#'
#' @param seed Integer, default = 1. Random seed for reproducibility.
#' @param nvec Numeric vector (length >= 2). Sample sizes of each study.
#' @param p Integer, default = 50. Number of variables (features) in the data.
#' @param q Integer, default = 3. Number of common factors (shared across all studies).
#' @param qs Numeric vector. Number of study-specific factors for each study.
#' @param fac.type Character, default = "gaussian". Factor distribution type, one of "gaussian" or "mvt".
#' @param err.type Character, default = "gaussian". Error distribution type, one of "gaussian" or "mvt".
#' @param rho Numeric vector of length 2, default = `c(1,1)`. Scaling factors for common and specific factor loadings.
#' @param sigma2_eps Numeric, default = 0.1. Variance of the error term.
#' @param nu Integer, default = 1. Degrees of freedom for the multivariate t-distribution ("mvt").
#'
#' @return A list containing the simulated data and true parameter values:
#' \itemize{
#'   \item{\code{Xlist}: List of data matrices (ns × p) for each study.}
#'   \item{\code{mu0}: True mean vector matrix.}
#'   \item{\code{A0}: True common factor loadings matrix.}
#'   \item{\code{Blist0}: List of true study-specific factor loadings matrices.}
#'   \item{\code{Flist}: List of true common factor score matrices.}
#'   \item{\code{Hlist}: List of true study-specific factor score matrices.}
#'   \item{\code{q}: Number of common factors used.}
#'   \item{\code{qs}: Number of study-specific factors used.}
#' }
#' 
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
#' @importFrom stats rnorm sd
#' 
#' @export
gendata_simu_robust <- function (seed = 1, nvec = c(100,300), p = 50, q = 3,
                                 qs= rep(2, length(nvec)), fac.type = c("gaussian", 'mvt'),
                                 err.type=c("gaussian", 'mvt'),
                                 rho = c(1,1), sigma2_eps=0.1, nu=1){
  
  if(length(nvec)<2) stop("nvec must have at least two elements!")
  fac.type <- match.arg(fac.type)
  err.type <- match.arg(err.type)
  S <- length(nvec)
  
  # Removed require(MASS) to comply with CRAN standards. Handled by @importFrom.
  
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
    }
    
    if(fac.type == 'gaussian'){
      message("Factor is from gaussian distribution.")
      FH <- MASS::mvrnorm(n, mu = rep(0, q+qs[s]), Diag(rep(1,q+qs[s])))
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

# -------------------------------------------------------------------------
# Helper functions (Hidden from User Manual)
# -------------------------------------------------------------------------

#' @keywords internal
#' @noRd
normvec <- function(x) sqrt(sum(x^2)/ length(x))

#' @keywords internal
#' @noRd
colSD <- function(X) apply(X, 2, sd, na.rm=TRUE)

#' @keywords internal
#' @noRd
trace_statistic_fun <- function(H, H0){
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% qr.solve(t(H) %*% H) %*% t(H) %*% H0
  tr_fun(mat1) / tr_fun(t(H0) %*% H0)
}

#' @keywords internal
#' @noRd
trace_list_fun <- function(Hlist, H0list){
  trvec <- rep(NA, length(Hlist))
  for(i in seq_along(trvec)){
    trvec[i] <- trace_statistic_fun(Hlist[[i]], H0list[[i]])
  }
  return(mean(trvec))
}