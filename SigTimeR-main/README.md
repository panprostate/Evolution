# SigTimeR

An R package that collates methods for timing mutational signatures:
- Non-negative least squares (NNLS) as described in [Gerstung et al. (2020)](https://www.nature.com/articles/s41586-019-1907-7).
- Heirarchical Dirichlet Process (HDP), forked from [nicolaroberts/hdp](https://github.com/nicolaroberts/hdp).
- Dirichlet-Multinomial. Defined in this package.

## Installation
    
```R
# Install remotes to install packages from GitHub
install.packages("remotes")

# Install SigTimeR
remotes::install_github("NMNS93/SigTimeR")
```
## Quickstart

```R
library(SigTimeR)

# Load mutation count data for two epochs along with known active signatures
st_data = SigTimeR::STData(test_clonal, test_subclonal, test_sigs)

# Time signatures
st = SigTimeR::time(st_data, method = "nnls")
```
