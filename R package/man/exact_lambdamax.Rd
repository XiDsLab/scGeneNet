\name{exact_lambdamax}
\alias{exact_lambdamax}
\title{
Exact the suitable value of maximal value in the path of hyper-parameter lambda.
}
\description{
Exact a list of hyper-parameter lambda for each component, which induce completely sparse network.
}
\usage{
exact_lambdamax(scGeneNet_list,
                lambda_max_ref = NULL,
                Global_Control = list(ELBO_threshold = 1e-04, minit = 5, maxit = 30, maxit_nonGRN = 3, MS_update_threshold = 1e-06),
                M_Control = list(ADMM_max_step = 1000, ADMM_threshold = 0.001),
                S_Control = list(Newton_max_step = 1000, Newton_threshold = 1e-04),
                Theta_Control = list(penalize_diagonal = FALSE, Theta_threshold = 1e-04),
                U_fix = FALSE,
                core_num = 20)
}
\arguments{
  \item{scGeneNet_list}{
A list of scGeneNet object with initialized parameters.
}
  \item{lambda_max_ref}{The regularization parameter for precision matrices.}
  \item{Global_Control}{A list of hyper-parameters that controling the convergence of VMPLN's optimization.}
  \item{M_Control}{A list of hyper-parameters that controling the optimization of variational mean parameter M using ADMM algorithm.}
  \item{S_Control}{A list of hyper-parameters that controling the optimization of variational standard deviation parameter S using Newton-Rapson algorithm.}
  \item{Theta_Control}{A list of hyper-parameters that controling the optimization of precision matrices Theta using graphical-lasso algorithm.}
  \item{U_fix}{Should the updatation of vatiational posterior probability U be freezing? Default is FALSE.}
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

##4. Exact the suitable value of maximal value in the path of hyper-parameter lambda.
##-------------------------------------------------------------
lambdamax_vec<-exact_lambdamax(scGeneNet_list = scGeneNet_list,lambda_max_ref = 10,
                     Theta_Control = list(penalize_diagonal = FALSE),
                     U_fix = FALSE,
                     core_num = 1)
##-------------------------------------------------------------

}