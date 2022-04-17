\name{scGeneNet_main}
\alias{scGeneNet_main}
\title{
The main function of training VMPLN algorithm
}
\description{
After initialization, we use this function for optimize the parameters of VMPLN based on a block-wise descent algorithm.
}
\usage{
scGeneNet_main(scGeneNet_list,
           lambda_use,
           Global_Control = list(ELBO_threshold = 1e-4, minit = 5, maxit = 30,maxit_nonGRN = 3, MS_update_threshold = 1e-6),
           M_Control = list(ADMM_max_step = 1000, ADMM_threshold = 1e-3),
           S_Control = list(Newton_max_step = 1000, Newton_threshold = 1e-4),
           Theta_Control = list(penalize_diagonal = FALSE,Theta_threshold = 1e-4),
           U_fix = FALSE, 
           verbose = TRUE
           core_num = 1)
}
\arguments{
  \item{scGeneNet_list}{
A list of scGeneNet object with initialized parameters.
}
  \item{lambda_use}{The regularization parameter for precision matrices.}
  \item{Global_Control}{A list of hyper-parameters that controling the convergence of VMPLN's optimization.}
  \item{M_Control}{A list of hyper-parameters that controling the optimization of variational mean parameter M using ADMM algorithm.}
  \item{S_Control}{A list of hyper-parameters that controling the optimization of variational standard deviation parameter S using Newton-Rapson algorithm.}
  \item{Theta_Control}{A list of hyper-parameters that controling the optimization of precision matrices Theta using graphical-lasso algorithm.}
  \item{U_fix}{Should the updatation of vatiational posterior probability U be freezing? Default is FALSE.}
  \item{verbose}{Display the detail progress of VMPLN's optimization. Default is TRUE.}
  \item{core_num}{The number of cores for parallel initialization. Default is 1.}
}
\details{
The list of hyper-parameters "Global_Control" controls the convergence of VMPLN's optimization, with the following entries:
\itemize{
\item "ELBO_threshold": The threshold for optimization convergence. Default value is 1e-4. Optimization stop when relative changle of Evidence of Low BOund (ELBO) between two steps is less than "ELBO_threshold".
\item "minit": The minimal number of optimization step. Default value is 5.
\item "maxit": The maximal number of optimization step. Default value is 30.
\item "maxit_nonGRN": The maximal number of optimization step for the parameters related to non-GRN part. Default value is 3. The optimization for the parameters related to non-GRN part will be freezing when the optimization step is greater than "maxit_nonGRN".
\item "MS_update_threshold": The threshold for optimization of posterior parameter M and S. Default value is 1e-6. The optimization of each M and S will be stop when its corresponding relative change between two steps is less than "MS_update_threshold".
}
The list of hyper-parameters "M_Control" controls the optimization of variational mean parameter M using ADMM algorithm, with the following entries:
\itemize{
\item "ADMM_max_step": The maximal number of optimization step for sub-algorithm based on ADMM. Default value is 1000.
\item "ADMM_threshold": The threshold for optimization convergence for sub-algorithm based on ADMM. Default value is 1e-3. Optimization stop when relative changle of Evidence of Low BOund (ELBO) between two steps within sub-algorithm is less than "ADMM_threshold". 
}
The list of hyper-parameters "S_Control" controls the optimization of variational standard deviation parameter S using Newton-Rapson algorithm, with the following entries:
\itemize{
\item "Newton_max_step": The maximal number of optimization step for sub-algorithm based on Newton-Rapson. Default value is 1000.
\item "Newton_threshold": The threshold for optimization convergence for sub-algorithm based on Newton-Rapson. Default value is 1e-4. Optimization stop when relative changle of Evidence of Low BOund (ELBO) between two steps within sub-algorithm is less than "Newton_threshold". 
}
The list of hyper-parameters "Theta_Control" controls the optimization of precision matrices Theta using graphical-lasso algorithm, with the following entries:
\itemize{
\item "penalize_diagonal": Should diagonal of precision matrix be penalized? Dafault value is TRUE.
\item "Theta_threshold": The threshold for convergence of glasso algorithm. Default value is 1e-4.
}
}
\value{
A list of scGeneNet object with updated model parameters.
}
\references{

Jordan, M. I., Z. Ghahramani, T. S. Jaakkola, and L. K. Saul (1999). An introduction to variational methods for graphical models. Machine learning 37(2), 183-233.

Boyd, S., N. Parikh, and E. Chu (2011). Distributed optimization and statistical learning via the alternating direction method of multipliers. Now Publishers Inc.

Friedman, J., T. Hastie, and R. Tibshirani (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics 9(3), 432-441.

Meinshausen, N. and P. Buhlmann (2006). High-dimensional graphs and variable selection with the lasso. The annals of statistics 34(3), 1436-1462.

Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for clustering with the integrated completed likelihood. IEEE transactions on pattern analysis and machine intelligence 22(7), 719-725.
}
\examples{
##1. Package loading
##-------------------------------------------------------------
library(scGeneNet)
library(Seurat)
library(glasso)
##-------------------------------------------------------------

##2. load reference dataset
##-------------------------------------------------------------
data("example_data")
##-------------------------------------------------------------


##3. Initialization
##-------------------------------------------------------------
gene_GRN <- rownames(example_data)
scGeneNet_list<-scGeneNet_init(expression_profile = example_data,
                     celltypes_num = 3,
                     celltypes_ref = NULL,
                     ls_est = "TSS",
                     gene_GRN = gene_GRN,
                     HVG_model_num = 0,
                     zero_GRN = NULL,
                     preprocess_Control = list(HVG_num = length(gene_GRN),npc = 50,
                     run_umap = TRUE,label_umap = NULL,
                     cluster_method = "Kmeans",resolution = 0.8),
                     core_num = 1)
##-------------------------------------------------------------

##4. Run the main function for optimization of model
##-------------------------------------------------------------
lambda_use<-10
scGeneNet_list<-scGeneNet_main(scGeneNet_list = scGeneNet_list,lambda_use = lambda_use,
                     Theta_Control = list(penalize_diagonal = FALSE),
                     U_fix = FALSE,
                     verbose = TRUE,
                     core_num = 1)
##-------------------------------------------------------------

}

