\name{scGeneNet-package}
\alias{scGeneNet-package}
\alias{scGeneNet}
\docType{package}
\title{
\packageTitle{scGeneNet}
}
\description{
\packageDescription{scGeneNet}
}
\details{
 
The DESCRIPTION file:
\packageDESCRIPTION{scGeneNet}
\packageIndices{scGeneNet}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{

\packageAuthor{scGeneNet}

Maintainer: \packageMaintainer{scGeneNet}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
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

##4. (Optional) Exact the path of hyper-parameter lambda
##-------------------------------------------------------------
if_exact_lambdapath<-FALSE ## TRUE for exact an example lambda path
if(if_exact_lambdapath){
lambdamax_vec<-exact_lambdamax(scGeneNet_list = scGeneNet_list,lambda_max_ref = 10,
                     Theta_Control = list(penalize_diagonal = FALSE),
                     U_fix = FALSE,
                     core_num = 1)
##Build an example lambda path
lambda_point<-c(1e-6,unique(lambdamax_vec))
lambda_path<-sort(c(lambda_point,lambda_point[-1] + lambda_point[-length(lambda_point)]))
}else{
lambda_path<-c(1e-6,5,10)
}
##-------------------------------------------------------------

##4. Run the main function for optimization of model
##-------------------------------------------------------------
scGeneNet_list_alongpath<-list()
for(l in 1:length(lambda_path)){
cat(paste(paste(rep("##",40),sep = "",collapse = ""),"\n",(paste(rep("##",40),sep = "",collapse = "")),"\n",sep = "",collapse = ""))
##
scGeneNet_list_alongpath[[l]]<-scGeneNet_main(scGeneNet_list = scGeneNet_list,lambda_use = lambda_path[l],
                                      Theta_Control = list(penalize_diagonal = FALSE),
                                      U_fix = FALSE,
                                      verbose = TRUE,
                                      core_num = 1)
}
names(scGeneNet_list_alongpath)<-paste("lambda: ",lambda_path,sep = "")
##-------------------------------------------------------------

}
