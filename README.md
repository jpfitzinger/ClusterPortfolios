# ClusterPortfolios

`ClusterPortfolios` is an R package for constructing portfolios based on statistical clustering techniques. Clustering financial asset returns and allocating capital along cluster boundaries can increase robustness, decrease sampling sensitivity, improve diversification and enhance portfolio performance.

Methods implemented in this package:

	-	Hierarchical Risk Parity (Lopez de Prado, 2016)
	-	Constrained HRP (Pfitzinger & Katzke, 2019)
	
## Installation

```
devtools::install_github()
```

## Usage

```
library(ClusterPortfolios)
data("Industry_10")
rets <- Industry_10
sigma <- cov(rets)
HRP(sigma, UB = 0.15, tau = 0.5)
```

## References

Lopez de Prado, M. (2016).
Building Diversified Portfolios that Outperform Out-of-Sample.
 _SSRN Electronic Journal_.

Pfitzinger, J., Katzke, N. (2019).
A Constrained Hierarchical Risk Parity Algorithm with Cluster-Based Capital Allocation.
_Stellenbosch University, Department of Economics_. Working Paper 14/2019.