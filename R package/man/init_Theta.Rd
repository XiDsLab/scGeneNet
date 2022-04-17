\name{init_Theta}
\alias{init_Theta}
\title{
Initialization of precision matrices in MPLN.
}
\description{
Initial the precision matrix of each component in MPLN based on the pre-estimated model parameters and variational parameters by graphical lasso algorithm.
}
\usage{
init_Theta(scGeneNet_list, lambda_use, penalize_diagonal, Theta_threshold, core_num)
}
\arguments{
  \item{scGeneNet_list}{
  A list of scGeneNet object.
  }
  \item{lambda_use}{
  The regularization parameter for graphical lasso algorithm.
  }
  \item{penalize_diagonal}{
  Should diagonal of inverse covariance be penalized?
  }
  \item{Theta_threshold}{
  Threshold for convergence of graphical lasso algorithm.
  }
  \item{core_num}{
  The number of cores for parallel initialization.
  }
  
}
\value{
  A list of scGeneNet object with updated precision matrices.
}
\references{
Friedman, J., T. Hastie, and R. Tibshirani (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics 9(3), 432-441.


Meinshausen, N. and P. Buhlmann (2006). High-dimensional graphs and variable selection with the lasso. The annals of statistics 34(3), 1436-1462.

}
