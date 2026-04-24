# Generate Simulated Data for Robust Multi-Study Elliptical Factor Model

Generate simulated data for high-dimensional multi-study factor analysis
under robust elliptical distributions (such as Multivariate
t-distribution). The data follows a factor model with common factors
(shared across studies) and study-specific factors, plus noise.

## Usage

``` r
gendata_simu_robust(
  seed = 1,
  nvec = c(100, 300),
  p = 50,
  q = 3,
  qs = rep(2, length(nvec)),
  fac.type = c("gaussian", "mvt"),
  err.type = c("gaussian", "mvt"),
  rho = c(1, 1),
  sigma2_eps = 0.1,
  nu = 1
)
```

## Arguments

- seed:

  Integer, default = 1. Random seed for reproducibility.

- nvec:

  Numeric vector (length \>= 2). Sample sizes of each study.

- p:

  Integer, default = 50. Number of variables (features) in the data.

- q:

  Integer, default = 3. Number of common factors (shared across all
  studies).

- qs:

  Numeric vector. Number of study-specific factors for each study.

- fac.type:

  Character, default = "gaussian". Factor distribution type, one of
  "gaussian" or "mvt".

- err.type:

  Character, default = "gaussian". Error distribution type, one of
  "gaussian" or "mvt".

- rho:

  Numeric vector of length 2, default = \`c(1,1)\`. Scaling factors for
  common and specific factor loadings.

- sigma2_eps:

  Numeric, default = 0.1. Variance of the error term.

- nu:

  Integer, default = 1. Degrees of freedom for the multivariate
  t-distribution ("mvt").

## Value

A list containing the simulated data and true parameter values:

- `Xlist`: List of data matrices (ns × p) for each study.

- `mu0`: True mean vector matrix.

- `A0`: True common factor loadings matrix.

- `Blist0`: List of true study-specific factor loadings matrices.

- `Flist`: List of true common factor score matrices.

- `Hlist`: List of true study-specific factor score matrices.

- `q`: Number of common factors used.

- `qs`: Number of study-specific factors used.
