# MultiEFM <a href='https://RRuihan-Zhang.github.io/MultiEFM/'><img src='docs/link.svg' align="right" height="30" /></a>

**MultiEFM** (Multi-study Elliptical Factor Model) is a robust statistical framework and R package designed for multi-study factor analysis. 

## 📖 Overview
Identifying study-shared and study-specific latent structures is central to multi-study factor analysis. However, existing methods heavily rely on Gaussian or finite second-moment assumptions, making them vulnerable to heavy-tailed distributions and outliers frequently encountered in modern applications like single-cell genomics.

**MultiEFM** addresses this limitation by operating under an elliptical factor model that does not require finite second moments. By leveraging geometric decoupling between radial and directional components, it constructs rank-based estimators via Kendall’s tau scatter matrices. 

## ✨ Key Features
* **Distribution-Free Robustness:** Achieves consistent estimation without moment conditions, making it highly resistant to heavy-tailed contamination.
* **Decoupled Estimation:** Recovers shared factor spaces through cross-study aggregation, while extracting study-specific factors via residual-based normalization.
* **Robust Factor Scores:** Utilizes a Huber-type estimator for reliable factor score computation.
* **Consistent Dimension Selection:** Introduces a frequency-based stability operator to accurately determine the number of latent factors.
* **High Performance:** Demonstrates superior statistical and computational efficiency compared to state-of-the-art methods.
* **Biological Interpretability:** Effectively separates biologically meaningful shared structures from condition- or platform-specific batch effects in single-cell RNA sequencing (scRNA-seq) datasets.

## 🚀 Installation

You can install the development version of MultiEFM from [GitHub](https://github.com/) with:

```r
# If you haven't installed devtools yet, please run: install.packages("devtools")
devtools::install_github("RRuihan-Zhang/MultiEFM")
