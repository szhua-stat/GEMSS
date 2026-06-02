# GEMSS: Generalization Error Minimization in SubSampling

**GEMSS** is an R package for optimal subdata selection in large-scale Gaussian Process (GP) regression. Accelerated by C++ (via RcppArmadillo), the algorithm sequentially identifies and selects the most informative data points to minimize generalization error, enabling efficient surrogate modeling for massive datasets.

## Installation

### Stable Release (CRAN)
You can install the released version of GEMSS from CRAN:

```R
install.packages("GEMSS")
```

### Development Version (GitHub)
You can install the development version from GitHub using the remotes package:
```R
# install.packages("remotes")
remotes::install_github("szhua-stat/GEMSS")
```

## Citation
If you use GEMSS in your research, please cite the following paper:

Chang, M. C., Hua, S. Z., & Wu, C. F. J. (2026). GEMSS-Driven Subsampling for Information Extraction and Redundancy Elimination. Technometrics, 1–20. https://doi.org/10.1080/00401706.2026.2670596
