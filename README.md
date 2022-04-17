# scGeneNet: Variational inference for Mixture Poisson Log-Normal graphial model in Single-cell gene regulatory network analysis with mixed cell populations
> scGeneNet is a R package for cell-type-specifc gene regulatory inference for Single-cell RNA-seq data with mixed cell populations using a block-wise descent algorithm based on the variational inference .

## Installation
**scGeneNet** is available on [Github](https://github.com/XiDsLab/scGeneNet).

### R Package installation
- For the development version, use the github install
```{r package github, eval = FALSE}
remotes::install_github("https://github.com/XiDsLab/scGeneNet")
```

## Usage and main functions
The package comes with a simulation example single-cell RNA-seq data to present the functionality of main function.

### Preparsion
```{r load scGeneNet, eval = FALSE}

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

```

### Initialization
```{r, warning = FALSE}
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
```
### (Optional) Exact the path of hyper-parameter lambda
```{r, warning = FALSE}
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
```
###  Run the main function for optimization of model
```{r, warning = FALSE}
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
```

## Reference

Please cite our work using the following references:
