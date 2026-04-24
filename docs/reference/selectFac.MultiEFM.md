# Select the number of factors for MultiEFM

Automatically estimate the number of study-shared factors and
study-specific factors based on the eigenvalue ratio method.

## Usage

``` r
selectFac.MultiEFM(
  XList,
  q_max = 20,
  qs_max = 20,
  qd = NULL,
  n_threads = 4,
  sample_pairs = 0,
  verbose = TRUE,
  seed = 1
)
```

## Arguments

- XList:

  A length-S list, where each component represents a data matrix for a
  specific study.

- q_max:

  an integer, specify the maximum number of study-shared factors to
  consider; default as 20.

- qs_max:

  an integer, specify the maximum number of study-specified factors to
  consider; default as 20.

- qd:

  an optional integer, specify the true number of combined factors if
  known in advance; default is NULL.

- n_threads:

  an integer, specify the number of threads for parallel computing;
  default as 4.

- sample_pairs:

  an integer, specify the number of sample pairs used for estimation;
  default as 0.

- verbose:

  a logical value, whether to output the estimation information.

- seed:

  an integer, specify the random seed for reproducibility; default as 1.

## Value

return a list including the following components:

- hq:

  the estimated number of study-shared factors;

- hqs_vec:

  an integer vector representing the estimated number of study-specific
  factors for each study;

- eigvals_A:

  the eigenvalues used for determining the shared factor number;

- eig_ratio_A:

  the eigenvalue ratios used for determining the shared factor number;

- eigvals_B:

  a list of eigenvalues used for determining the specific factor
  numbers;

- eig_ratio_B:

  a list of eigenvalue ratios used for determining the specific factor
  numbers.
