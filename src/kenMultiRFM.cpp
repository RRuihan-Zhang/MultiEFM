// #define ARMA_64BIT_WORD 1
#include "RcppArmadillo.h"
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


/***  irlba wrapper  ***/
List irlbaCpp(const mat& X, const int& q){
  Environment irlba("package:irlba");
  Function f = irlba["irlba"];
  return f(Named("A") = X, Named("nv") = q);
}

// pair-sampling version
// [[Rcpp::export]]
arma::mat SK_pair_MC_cpp_parallel(const arma::mat& X,
                                  int m,
                                  int n_threads = 0) {

  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat TK(p, p, arma::fill::zeros);

#ifdef _OPENMP
  if(n_threads > 0) omp_set_num_threads(n_threads);

#pragma omp parallel
{
  arma::mat TK_local(p, p, arma::fill::zeros);

#pragma omp for schedule(static)
  for(int k = 0; k < m; k++) {

    // ---- sample i < j uniformly ----
    int i = std::rand() % n;
    int j = std::rand() % n;
    while(j == i)
      j = std::rand() % n;

    // ---- compute difference ----
    arma::rowvec d = X.row(i) - X.row(j);

    double norm = std::sqrt(arma::dot(d, d));
    if(norm > 1e-12) {
      d /= norm;
      TK_local += d.t() * d;
    }
  }

#pragma omp critical
  TK += TK_local;
}
#else
// fallback
for(int k = 0; k < m; k++) {

  int i = std::rand() % n;
  int j = std::rand() % n;
  while(j == i)
    j = std::rand() % n;

  arma::rowvec d = X.row(i) - X.row(j);
  double norm = std::sqrt(arma::dot(d, d));
  if(norm > 1e-12) {
    d /= norm;
    TK += d.t() * d;
  }
}
#endif

  TK /= m;

  return TK;
}



/***  Fast Kendall (subsampled)  ***/
// [[Rcpp::export]]
arma::mat SK_subsample_fast_cpp_parallel(const arma::mat& X, 
                                         double rate = 0.7,
                                         int n_threads = 0){
  
  int n = X.n_rows;
  int p = X.n_cols;
  int n_sub = std::round(n * rate);
  
  arma::mat Xsub;
  // random subsample
  if(n_sub==n){
    Xsub = X;
  }else{
    arma::uvec idx = arma::randperm(n, n_sub);
    Xsub = X.rows(idx);
  }
  
  arma::mat TK(p, p, arma::fill::zeros);
  
#ifdef _OPENMP
  if(n_threads > 0) omp_set_num_threads(n_threads);
#pragma omp parallel
{
  arma::mat TK_local(p, p, arma::fill::zeros);
  
#pragma omp for schedule(static) //dynamic
  for(int i = 0; i < n_sub - 1; i++){
    
    arma::mat D = Xsub.rows(i+1, n_sub-1);
    D.each_row() -= Xsub.row(i);
    
    arma::vec nr = arma::sum(arma::square(D), 1);
    nr = arma::sqrt(nr);
    D.each_col() /= nr;
    
    TK_local += D.t() * D;
  }
  
#pragma omp critical
  TK += TK_local;
}
#else
// fallback (no OpenMP)
for(int i = 0; i < n_sub - 1; i++){
  arma::mat D = Xsub.rows(i+1, n_sub-1);
  D.each_row() -= Xsub.row(i);
  arma::vec nr = arma::sum(arma::square(D), 1);
  nr = arma::sqrt(nr);
  D.each_col() /= nr;
  TK += D.t() * D;
}
#endif

TK *= (2.0 / (n_sub * (n_sub - 1.0)));

return TK;
}


arma::mat SK_subsample_fast_cpp(const arma::mat& X,
                                double rate = 0.7) {
  
  int n = X.n_rows;
  int p = X.n_cols;
  int n_sub = std::round(n * rate);
  
  arma::mat Xsub;
  // random subsample
  if(n_sub==n){
    Xsub = X;
  }else{
    arma::uvec idx = arma::randperm(n, n_sub);
    Xsub = X.rows(idx);
  }
  
  
  arma::mat TK(p, p, arma::fill::zeros);
  
  for(int i = 0; i < n_sub - 1; i++){
    
    arma::mat D = Xsub.rows(i+1, n_sub-1);
    D.each_row() -= Xsub.row(i);
    
    arma::vec nr = sqrt(sum(square(D), 1));
    
    D.each_col() /= nr;
    
    TK += D.t() * D;
  }
  
  TK *= (2.0 / (n_sub * (n_sub - 1.0)));
  
  return TK;
}


void update_FH_huber_cpp(
    const mat& Xs,
    const mat& A,
    const mat& B,
    mat& Fs,
    mat& Hs,
    const double c1 = 1.345,
    const double c2 = 0.6745,
    const int maxIter = 50,
    const double tol = 1e-5,
    const bool verbose = false
){
  
  int ns = Xs.n_rows;
  int p  = Xs.n_cols;
  
  mat W = ones<mat>(ns, p);     // weights
  mat R(ns,p);                  // residuals
  vec obj(maxIter, fill::zeros);
  
  for(int it = 0; it < maxIter; it++){
    
    //----------------------------------------
    // 1. compute residuals
    //----------------------------------------
    R = Xs - Fs * A.t() - Hs * B.t();
    
    //----------------------------------------
    // 2. compute Huber-type adaptive weights
    //----------------------------------------
    
    vec absR = vectorise(abs(R));
    double medAbsR = median(absR);
    
    double tau = c1 / (c2 * medAbsR + 1e-12);
    
    for(uword k = 0; k < R.n_elem; k++){
      double r = R(k);
      double ar = std::abs(r);
      
      if(ar <= tau) W(k) = 1.0;
      else W(k) = tau / ar;
    }
    
    //----------------------------------------
    // 3. update F_s (Hadamard weighted LS)
    //----------------------------------------
    mat XA = Xs - Hs * B.t();
    mat WA = W % XA;
    Fs = WA * A;
    
    //----------------------------------------
    // 4. update H_s (Hadamard weighted LS)
    //----------------------------------------
    mat XB = Xs - Fs * A.t();
    mat WB = W % XB;
    Hs = WB * B;
    
    //----------------------------------------
    // 5. compute objective (Huber loss)
    //----------------------------------------
    R = Xs - Fs * A.t() - Hs * B.t();
    double loss = 0.0;
    
    for(uword k = 0; k < R.n_elem; k++){
      double r = R(k);
      double ar = std::abs(r);
      if(ar <= tau) loss += 0.5 * r * r;
      else loss += tau * ar - 0.5 * tau * tau;
    }
    
    obj(it) = loss;
    
    if(verbose)
      Rprintf("iter %d  obj = %.6f\n", it+1, loss);
    
    if(it > 0 && std::abs(obj(it)-obj(it-1)) / obj(it-1) < tol)
      break;
  }
}


/***   Main Function   ***/

// [[Rcpp::export]]
Rcpp::List kenMultiRFM_cpp(const Rcpp::List& XsList,
                           const int& q,
                           const arma::vec& qs_vec,
                           const double& subsample_rate,
                           const int& factor_method, 
                           const int& sampling_method,
                           const int& sample_pairs, 
                           const int& n_threads, 
                           const double& epsObj,
                           const int& maxIter,
                           const bool& verbose){
  
  int n, S = XsList.size();
  
  field<mat> Xsf(S), hBf(S), hFf(S), hHf(S);
  
  int n_total = 0;
  
  for(int s = 0; s < S; ++s){
    mat Xtmp = as<mat>(XsList[s]);
    Xsf(s) = Xtmp;
    n_total += Xtmp.n_rows;
  }
  
  int p = Xsf(0).n_cols;
  
  
  /***********************
   *  Study-shared Loadings A
   ***********************/
  int num_sample_pairs;
  mat hK_pool(p,p,fill::zeros);
  mat K_tmp;
  for(int s = 0; s < S; ++s){
    n = Xsf(s).n_rows;
    if(sampling_method==1){
      K_tmp = SK_subsample_fast_cpp_parallel(Xsf(s), subsample_rate, n_threads);
    }else if(sampling_method ==2){
      if(sample_pairs==0){
        num_sample_pairs = n*(n-1)/2;
      }else if(sample_pairs==1){
        num_sample_pairs = int(n*log(n));
      }else{
        num_sample_pairs = sample_pairs;
      }
      K_tmp = SK_pair_MC_cpp_parallel(Xsf(s), num_sample_pairs, n_threads);
    }
    
    
    double w = double(Xsf(s).n_rows) / double(n_total);
    hK_pool += w * K_tmp;
  }
  
  List svdX = irlbaCpp(hK_pool, q);
  mat hA = as<mat>(svdX["v"]);     // p × q
  vec eigvals_A = as<vec>(svdX["d"]);
  
  
  
  /***********************
   *  Study-specific Loadings B_s
   ***********************/
  
  field<vec> eigvals_Bf(S);
  for(int s = 0; s < S; ++s){
    
    mat XA = Xsf(s) * hA;
    mat R = Xsf(s) - XA * hA.t();
    n = R.n_rows;
    
    if(sampling_method==1){
      K_tmp = SK_subsample_fast_cpp_parallel(R, subsample_rate, n_threads);
    }else if(sampling_method ==2){
      if(sample_pairs==0){
         num_sample_pairs = n*(n-1)/2;
      }else if(sample_pairs==1){
         num_sample_pairs = int(n*log(n));
      }else{
         num_sample_pairs = sample_pairs;
      }
      K_tmp = SK_pair_MC_cpp_parallel(R, num_sample_pairs, n_threads);
    }
    
    
    List svdR = irlbaCpp(K_tmp, qs_vec(s));
    
    hBf(s) = as<mat>(svdR["v"]);   // p × q_s
    // Rprintf("number of columns is: %d\n", hBf(s).n_cols);
    eigvals_Bf(s) = as<vec>(svdR["d"]); 
  }
  
  
  /***********************
   *  Initialize factors
   ***********************/
  for(int s = 0; s < S; ++s){
    hFf(s) = Xsf(s) * hA;
    hHf(s) = Xsf(s) * hBf(s);
  }
  
  
  /***********************
   *  Global Alternating Updates
   ***********************/
  vec Obj(maxIter);
  Obj(0) = -1e2;
  
  for(int iter = 1; iter < maxIter; ++iter){
    
    //--------------------------------
    // update factors for each study
    //--------------------------------
    if(factor_method == 1){
      
      // ----- Huber–IRLS -----
      for(int s = 0; s < S; ++s){
        
        update_FH_huber_cpp(
          Xsf(s),      // X_s
          hA,          // shared loading
          hBf(s),      // study-specific loading
          hFf(s),      // <-- UPDATED IN PLACE
          hHf(s),      // <-- UPDATED IN PLACE
          1.345,
          0.6745,
          50,
          1e-5,
          false        // verbose per study = false
        );
      }
      
    } else if(factor_method == 2){
      
      // ----- OLS 更新 -----
      for(int s = 0; s < S; ++s){
        
        mat& B = hBf(s);
        mat& X = Xsf(s);
        
        hFf(s) = (X - hHf(s) * B.t()) * hA;
        hHf(s) = (X - hFf(s) * hA.t()) * B;
      }
    }
    
    
    //--------------------------------
    // compute global objective
    //--------------------------------
    double obj_now = 0.0;
    
    for(int s = 0; s < S; ++s){
      mat R = Xsf(s) - hFf(s)*hA.t() - hHf(s)*hBf(s).t();
      obj_now += accu(square(R));
    }
    
    Obj(iter) = obj_now;
    
    
    //--------------------------------
    // print & stopping rule
    //--------------------------------
    if(verbose){
      Rprintf("Iter %d: Obj = %.6f, rel change = %.3e\n",
              iter,
              obj_now,
              std::abs((Obj(iter)-Obj(iter-1))/Obj(iter-1)));
    }
    
    if(std::abs((Obj(iter)-Obj(iter-1))/Obj(iter-1)) < epsObj){
      Obj = Obj.head(iter+1);
      break;
    }
  }
  
  /***********************
   *  Output
   ***********************/
  return List::create(
    Named("F")   = hFf,
    Named("H")   = hHf,
    Named("A")   = hA,
    Named("B")   = hBf,
    Named("eigvals_A")   = eigvals_A,
    Named("eigvals_B")   = eigvals_Bf,
    Named("Obj") = Obj
  );
}


// [[Rcpp::export]]
Rcpp::List kenMultiRFM_project_cpp(const Rcpp::List& XsList,
                           const int& q,
                           const arma::vec& qs_vec,
                           const double& subsample_rate,
                           const int& factor_method, 
                           const int& n_threads, 
                           const double& epsObj,
                           const int& maxIter,
                           const bool& verbose){
  
  int S = XsList.size();
  
  field<mat> Xsf(S), Ksf(S), hBf(S), hFf(S), hHf(S);
  field<vec> eigvals_Bf(S);
  int n_total = 0;
  
  for(int s = 0; s < S; ++s){
    mat Xtmp = as<mat>(XsList[s]);
    Xsf(s) = Xtmp;
    n_total += Xtmp.n_rows;
  }
  
  int p = Xsf(0).n_cols;
  
  
  /***********************
   *  Study-shared Loadings A
   ***********************/
  mat hK_pool(p,p,fill::zeros);
  for(int s = 0; s < S; ++s){
    
    Ksf(s) = SK_subsample_fast_cpp_parallel(Xsf(s), subsample_rate, n_threads);
    
    double w = double(Xsf(s).n_rows) / double(n_total);
    hK_pool += w * Ksf(s);
  }
  
  List svdX = irlbaCpp(hK_pool, q);
  mat hA = as<mat>(svdX["v"]);     // p × q
  vec eigvals_A = as<vec>(svdX["d"]);
  
  
  
  mat At = hA.t();
  
  /***********************
   *  Study-specific loadings B_s
   ***********************/
  for(int s = 0; s < S; ++s){
    
    mat& K = Ksf(s);
    
    mat AtK  = At * K;
    mat KA   = K * hA;
    mat AtKA = AtK * hA;
    
    mat K_res = K- hA * AtK- KA * At + hA * AtKA * At;
    
    List svdB = irlbaCpp(K_res, qs_vec(s));
    hBf(s) = as<mat>(svdB["v"]);
    eigvals_Bf(s) = as<vec>(svdB["d"]);
  }
  
  
  
  /***********************
   *  Initialize factors
   ***********************/
  for(int s = 0; s < S; ++s){
    hFf(s) = Xsf(s) * hA;
    hHf(s) = Xsf(s) * hBf(s);
  }
  
  
  /***********************
   *  Global Alternating Updates
   ***********************/
  vec Obj(maxIter);
  Obj(0) = -1e2;
  
  for(int iter = 1; iter < maxIter; ++iter){
    
    //--------------------------------
    // update factors for each study
    //--------------------------------
    if(factor_method == 1){
      
      // ----- Huber–IRLS -----
      for(int s = 0; s < S; ++s){
        
        update_FH_huber_cpp(
          Xsf(s),      // X_s
          hA,          // shared loading
          hBf(s),      // study-specific loading
          hFf(s),      // <-- UPDATED IN PLACE
          hHf(s),      // <-- UPDATED IN PLACE
          1.345,
          0.6745,
          50,
          1e-5,
          false        // verbose per study = false
        );
      }
      
    } else if(factor_method == 2){
      
      // ----- OLS 更新 -----
      for(int s = 0; s < S; ++s){
        
        mat& B = hBf(s);
        mat& X = Xsf(s);
        
        hFf(s) = (X - hHf(s) * B.t()) * hA;
        hHf(s) = (X - hFf(s) * hA.t()) * B;
      }
    }
    
    
    //--------------------------------
    // compute global objective
    //--------------------------------
    double obj_now = 0.0;
    
    for(int s = 0; s < S; ++s){
      mat R = Xsf(s) - hFf(s)*hA.t() - hHf(s)*hBf(s).t();
      obj_now += accu(square(R));
    }
    
    Obj(iter) = obj_now;
    
    
    //--------------------------------
    // print & stopping rule
    //--------------------------------
    if(verbose){
      Rprintf("Iter %d: Obj = %.6f, rel change = %.3e\n",
              iter,
              obj_now,
              std::abs((Obj(iter)-Obj(iter-1))/Obj(iter-1)));
    }
    
    if(std::abs((Obj(iter)-Obj(iter-1))/Obj(iter-1)) < epsObj){
      Obj = Obj.head(iter+1);
      break;
    }
  }
  
  /***********************
   *  Output
   ***********************/
  return List::create(
    Named("F")   = hFf,
    Named("H")   = hHf,
    Named("A")   = hA,
    Named("B")   = hBf,
    Named("eigvals_A")   = eigvals_A,
    Named("eigvals_B")   = eigvals_Bf,
    Named("Obj") = Obj
  );
}


