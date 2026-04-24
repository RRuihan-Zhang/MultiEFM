# High Dimensional Example of MultiEFM

This vignette introduces the usage of MultiEFM for the analysis of
high-dimensional multi-study multivariate data with heavy tail, by
comparison with other methods. In high-dimensional settings,
computational efficiency and robustness are particularly crucial.

The package can be loaded with the command, and define some metric
functions:

``` r
library(MultiEFM)
library(VIMSFA)
library(irlba)   
library(ggplot2)
library(cowplot) 

trace_statistic_fun <- function(H, H0){
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% qr.solve(t(H) %*% H) %*% t(H) %*% H0
  tr_fun(mat1) / tr_fun(t(H0) %*% H0)
}
trace_list_fun <- function(Hlist, H0list){
  trvec <- sapply(seq_along(Hlist), function(i) trace_statistic_fun(Hlist[[i]], H0list[[i]]))
  return(mean(trvec, na.rm = TRUE))
}
```

## Generate High-Dimensional Simulated Data

First, we generate the simulated data with heavy tail, where the error
term follows from a multivariate t-distribution with degree of freedom
3. To simulate a high-dimensional scenario, we significantly increase
the feature dimension to p = 500.

``` r
set.seed(1)
nu <- 3 # nu is set to
p <- 500 # High-dimensional setting
nvec <- c(150,200);  q <- 3; qs <- c(2,2); S <- length(nvec)
sigma2_eps <- 1
datList <- gendata_simu_robust(seed=1, nvec=nvec, p=p, q=q, qs=qs, rho=c(5,5), err.type='mvt', nu=nu)
XList <- datList$Xlist
```

## Fit the MultiEFM model

Fit the MultiEFM model using the function MultiEFM() in the R package
MultiEFM. Users can use ?MultiEFM to see the details about this
function. For two matrices $`\widehat D`$ and $`D`$, we use trace
statistic to measure their similarity. The trace statistic ranges from 0
to 1, with higher values indicating better performance.

``` r
methodNames <- c("MultiEFM", "MSFA-CAVI", "MSFA-SVI")
metricMat <- matrix(NA, nrow=length(methodNames), ncol=5)
colnames(metricMat) <- c('A_tr', 'B_tr',  'F_tr', 'H_tr', 'Time')
row.names(metricMat) <- methodNames

tic <- proc.time()
res <- MultiEFM(XList, q=q, qs_vec=qs)
toc <- proc.time()

metricMat["MultiEFM",'Time'] <- toc[3] - tic[3]
metricMat["MultiEFM",'A_tr'] <- trace_statistic_fun(res$A, datList$A0)
metricMat["MultiEFM",'B_tr'] <- trace_list_fun(res$B, datList$Blist0)
metricMat["MultiEFM",'F_tr'] <- trace_list_fun(res$F, datList$Flist)
metricMat["MultiEFM",'H_tr'] <- trace_list_fun(res$H, datList$Hlist)
```

## Compare with other methods

We compare MultiEFM with two prominent methods: MSFA-CAVI and MSFA-SVI

First, we implement MSFA-CAVI:

``` r
X_s <- lapply(XList, scale, scale=FALSE)
### MSFA-CAVI
print("MSFA-CAVI")
tic <- proc.time()
cavi_est <- cavi_msfa(X_s,  K=q, J_s=qs)
toc <- proc.time()
time_cavi <- toc[3] - tic[3]

hF_cavi <- hH_cavi <- list()
for(s in 1:S){
  hF_cavi[[s]] <- t(Reduce(cbind, cavi_est$mean_f[[s]]))
  hH_cavi[[s]] <- t(Reduce(cbind, cavi_est$mean_l[[s]]))
}

metricMat["MSFA-CAVI",'Time']  <- time_cavi
metricMat["MSFA-CAVI",'A_tr']  <- trace_statistic_fun(cavi_est$mean_phi, datList$A0)
metricMat["MSFA-CAVI",'B_tr']  <- trace_list_fun(cavi_est$mean_lambda_s, datList$Blist0)
metricMat["MSFA-CAVI",'F_tr'] <- trace_list_fun(hF_cavi, datList$Flist)
metricMat["MSFA-CAVI",'H_tr'] <- trace_list_fun(hH_cavi, datList$Hlist)
```

Next, we implement MSFA-SVI:

``` r
print("MSFA-SVI")
tic <- proc.time()
svi_est <- svi_msfa(X_s, K=q, J_s=qs, verbose = 0)
toc <- proc.time()
time_svi <- toc[3] - tic[3]

hF_svi <- hH_svi <- list()
for(s in 1:S){
  hF_svi[[s]] <- t(Reduce(cbind, svi_est$mean_f[[s]]))
  hH_svi[[s]] <- t(Reduce(cbind, svi_est$mean_l[[s]]))
}

metricMat["MSFA-SVI",'Time']  <- time_svi
metricMat["MSFA-SVI",'A_tr']  <- trace_statistic_fun(svi_est$mean_phi, datList$A0)
metricMat["MSFA-SVI",'B_tr']  <- trace_list_fun(svi_est$mean_lambda_s, datList$Blist0)
metricMat["MSFA-SVI",'F_tr'] <- trace_list_fun(hF_svi, datList$Flist)
metricMat["MSFA-SVI",'H_tr'] <- trace_list_fun(hH_svi, datList$Hlist)
```

## Visualize the comparison of performance

Next, we summarized the metrics for MultiEFM and other compared methods
in a data.frame object.

``` r
dat_metric <- data.frame(metricMat)
dat_metric$Method <- factor(row.names(dat_metric), levels=row.names(dat_metric))
```

Plot the results for MultiEFM and other methods, which suggests that
MultiEFM achieves better estimation accuracy for the study-shared
loading matrix A, study-specified loading matrix B and factor matrix H.
MultiEFM significantly outperforms the compared methods in terms of
estimation accuracy of B and H, as well as computational efficiency.

``` r
my_theme <- theme_bw(base_size = 14) + 
  theme(legend.position = "none",  
        plot.title = element_text(size = 14, hjust = 0.5), 
        axis.text.x = element_text(angle = 30, hjust = 1)) 

p1 <- ggplot(data=subset(dat_metric, !is.na(A_tr)), aes(x= Method, y=A_tr, fill=Method)) + geom_bar(stat="identity", width=0.6) + labs(title="Shared Loading (A)", x=NULL, y="Trace") + my_theme
p2 <- ggplot(data=subset(dat_metric, !is.na(F_tr)), aes(x= Method, y=F_tr, fill=Method)) + geom_bar(stat="identity", width=0.6) + labs(title="Shared Factor (F)", x=NULL, y="Trace") + my_theme
p3 <- ggplot(data=subset(dat_metric, !is.na(B_tr)), aes(x= Method, y=B_tr, fill=Method)) + geom_bar(stat="identity", width=0.6) + labs(title="Specific Loading (B)", x=NULL, y="Trace") + my_theme
p4 <- ggplot(data=subset(dat_metric, !is.na(H_tr)), aes(x= Method, y=H_tr, fill=Method)) + geom_bar(stat="identity", width=0.6) + labs(title="Specific Factor (H)", x=NULL, y="Trace") + my_theme
p5 <- ggplot(data=subset(dat_metric, !is.na(Time)), aes(x= Method, y=Time, fill=Method)) + geom_bar(stat="identity", width=0.6) + labs(title="Time (s)", x=NULL, y="Seconds") + my_theme

plot_grid(p1, p2, p3, p4, p5, nrow=2, ncol=3)
```

## Select the parameters

We applied the proposed TSP method to select the number of factors. The
results showed that the method has the potential to identify the true
values.

``` r
hq_res <- selectFac.MultiEFM(XList, q_max=10, qs_max=5, verbose = FALSE)
message("Estimated shared q = ", hq_res$hq, " VS true q = ", q)
```

Session Info

``` r
sessionInfo()
```
