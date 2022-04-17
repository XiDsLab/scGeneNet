\name{scGeneNet_init}
\alias{scGeneNet_init}
\title{
Initialization of scGeneNet algorithm.
}
\description{
Given the single cell RNA-seq transcriptomic profile, the number of cell types and the gene set for cell-type-specifc gene regularory network of interest, we first follow the standard Seurat pipeline for cells clustering by SNN or Kmeans. The gene set we need for VMPLN algorithm consists of two parts, including the gene set for cell-type-specifc gene regularory network of interest and a number of top highly variable gene detected by Seurat. After clustering, we initial the model parameters and variational parameters of VMPLN algorithm, and finally output a scGeneNet type list which consists the parameters we need.
}
\usage{
scGeneNet_init(expression_profile,
           celltypes_num,
           celltypes_ref = NULL,
           ls_est = "TSS",
           gene_GRN,
           HVG_model_num,
           zero_GRN = NULL,
           preprocess_Control = list(HVG_num = 2000,npc = 50,run_umap = TRUE,label_umap = NULL, cluster_method = "SNN"),
           core_num = 1)
}
\arguments{
  \item{expression_profile}{
  A (sparse) expression profile matrix with genes as rows and dimensions names.
  }
  \item{celltypes_num}{
  The number of cell types.
  }
  \item{celltypes_ref}{
  (Optional) A character vector that consists of the reference cell type label of cells. Default is NULL.
  }
  \item{ls_est}{
  The estimation methods of library size for normalizing the expression profile. Available methods are "TSS" (Total Sum Scaling, which is normalized by the sum of each cell divided by 10^4), "GMPR" (Geometric Mean of Pairwise Ratio, introduced in) and "none". Default is "TSS".
  }
  \item{gene_GRN}{
  A character vector as the gene set that consists of the genes within the gene regulatory network of interest.
  }
  \item{HVG_model_num}{
  The number of highly variable gene (HVG) for MPLN modeling.
  }
  \item{zero_GRN}{
  (Optional) A two-column matrix, and each row indicating the indices of elements in gene regulatory network of interest to be constrained to be zero. Default is NULL.
  }
  \item{preprocess_Control}{
  A list of hyper-parameters that controling the pre-process of VMPLN's initialization.
  }
  \item{core_num}{
  The number of cores for parallel optimization. Default is 1.
  }
}
\details{
The list of hyper-parameters "preprocess_Control" controls the pre-process of VMPLN's initialization, with the following entries:
\itemize{
\item "HVG_num": The number of highly variable gene (HVG) used in standard Seurat pipeline. Default is 2000.
\item "npc": The number of principal components used for dimensional reduction in standard Seurat pipeline.  Default is 50.
\item "run_umap": Should run the UMAP in standard Seurat pipeline?  Default is TRUE.
\item "label_umap": (Optional) The UMAP will be colored by label_umap. Default is NULL.
\item "cluster_method": The clustering method for cells based on the principal components representation. The available methods are "SNN" and "Kmeans" (default).
}
}
\value{
  A list of scGeneNet object with initialized parameters.
}
\references{

Friedman, J., T. Hastie, and R. Tibshirani (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics 9(3), 432-441.

Meinshausen, N. and P. Buhlmann (2006). High-dimensional graphs and variable selection with the lasso. The annals of statistics 34(3), 1436-1462.

Stuart, T., A. Butler, P. Hoffman, C. Hafemeister, E. Papalexi, W. M. Mauck III, Y. Hao, M. Stoeckius, P. Smibert, and R. Satija (2019). Comprehensive integration of single-cell data. Cell 177(7), 1888-1902.

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

}
