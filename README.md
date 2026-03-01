# GEMSS: Generalization Error Minimization in SubSampling

**GEMSS** is an R package designed for subdata selection for Gaussian Process regression under large datasets. It uses C++ acceleration (via RcppArmadillo) to sequentially select the most informative points.

## Installation
You can install the development version from GitHub:
```R
# install.packages("remotes")
remotes::install_github("szhua-stat/GEMSS")
