##1. Package required 
##--------------------------------------------------------------------------------------------
library(Matrix)
library(MASS)
library(glasso)
library(mclust)
library(mltools)
library(data.table)
library(doParallel)
library(PLNmodels)
library(scGeneNet)
library(cluster)
library(Seurat)
library(pryr)
library(philentropy)
library(umap)
library(igraph)
library(stats)
library(fpc)
library(ppcor)
library(GENIE3)
library(XMRF)
##--------------------------------------------------------------------------------------------

##2. Define some function
##--------------------------------------------------------------------------------------------
#tr
tr<-function(mat_1){
  return(sum(diag(as.matrix(mat_1))))
}

#l1_norm of matrix
l1_norm<-function(mat_1,mat_2){
  return(max(rowSums(abs(mat_1 - mat_2))))
}

#l2_norm of matrix
l2_norm<-function(mat_1,mat_2){
  diff_mat<-mat_1 - mat_2
  eigen_res<-eigen(diff_mat %*% t(diff_mat))
  return(sqrt(max(abs(eigen_res$values))))
}

#LF_norm of matrix
lF_norm<-function(mat_1,mat_2){
  return(sqrt(sum((mat_1 - mat_2)^2)))
}

#LF_norm of matrix without diagonal
lF_norm_offdiag<-function(mat_1,mat_2){
  mat_1_offdiag<-mat_1
  mat_2_offdiag<-mat_2
  diag(mat_1_offdiag)<-0
  diag(mat_2_offdiag)<-0
  return(sqrt(sum((mat_1_offdiag - mat_2_offdiag)^2)))
}

#log of det of matrix
log_det<-function(matrix_use){
  eigen_res<-eigen((matrix_use + t(matrix_use))/2)
  eig_value<-eigen_res$values
  eig_value<-ifelse(eig_value<=0,1e-8,eig_value)
  return(sum(log(eig_value)))
}

#BIC_glasso
BIC_glasso<-function(data_1,pre_mat_array){
  #data 1: n*p
  #pre_mat_array: L*p*p
  sample_size<-nrow(data_1)
  cov_mat<-cov(data_1)
  BIC_vec<-c()
  for(l in 1:dim(pre_mat_array)[1]){
    BIC_vec<-c(BIC_vec,sample_size*(-1*log_det(pre_mat_array[l,,])+tr(cov_mat%*%pre_mat_array[l,,]))+
                 log(sample_size)*sum(ifelse(pre_mat_array[l,,]!=0,1,0))/2)  
  }
  return(BIC_vec)
}

##prior tran
prior_tran<-function(p_GRN,group_num,gene_GRN_name,prior_nw_list = NULL, prior_TF = NULL,weight_prior_nw = 0.3,penalize_diagonal = TRUE){
  ################
  ##weighted network
  prior_nw_fill = list()
  if(!is.null(prior_nw_list)){
    # prior_nw_fill = prior_nw_list
    prior_nw_fill = lapply(X = 1:length(prior_nw_list),FUN = function(g){return(prior_nw_list[[g]][1:p_GRN,1:p_GRN])})
    for(g in 1:group_num){
      if(is.null(prior_nw_fill[[g]])){
        prior_nw_fill[[g]] <- matrix(0,nrow = p_GRN, ncol = p_GRN)
      }
    }
  }else{
    for(g in 1:group_num){
      prior_nw_fill[[g]] <- matrix(0,nrow = p_GRN, ncol = p_GRN)
    }
  }
  
  weight_mat_list<-list()
  for(g in 1:group_num){
    weight_mat0<-prior_nw_fill[[g]] * weight_prior_nw + (1 - prior_nw_fill[[g]])
    if(penalize_diagonal == TRUE){
      diag(weight_mat0)<-1
    }else{
      diag(weight_mat0)<-0
    }
    weight_mat_list[[g]]<-weight_mat0
    ##
    rm(weight_mat0);gc()
  }
  rm(prior_nw_fill)
  
  ##zero mat
  prior_adjoint_mat = matrix(1,ncol = p_GRN,nrow = p_GRN)
  if(!is.null(prior_TF)){
    index <- which(gene_GRN_name %in% prior_TF)
    index = setdiff(1:p_GRN, index) 
    prior_adjoint_mat[index,index] = 0
  }
  
  zero_mat<-NULL
  zero_num<-length(which(prior_adjoint_mat[upper.tri(prior_adjoint_mat)] == 0))
  if(zero_num > 0){
    zero_mat<-matrix(NA,nrow = 0,ncol = 2)
    for(i in 1:(p_GRN - 1)){
      for(j in (i+1):p_GRN){
        if(prior_adjoint_mat[i,j] == 0){
          zero_mat<-rbind(zero_mat,c(i,j))
        }
      }
    }
  }
  
  return(list(weight_mat_list = weight_mat_list,
              zero_mat = zero_mat))
}

##evaluator_1
evaluator_1<-function(input_array_hat,
                      type = 1,
                      Theta_array_true,
                      density_choose = NULL,
                      cluster_predict_true_vec,
                      uppertri_choose = NULL
){
  #####
  ##input_array_hat: group_num * p * p
  ##Theta_array_true: group_num * p * p
  ## performance criterion: signed AUROC AUPRC EP TPR TDR
  #####
  group_num<-length(cluster_predict_true_vec)
  res_mat<-matrix(NA,nrow = group_num,ncol = 12)
  if(is.null(uppertri_choose)){
    uppertri_choose <- 1:((dim(input_array_hat)[2]) * (dim(input_array_hat)[2] - 1)/2)
  }
  if(is.null(density_choose)){
    density_choose<-1
  }
  ##
  for(g in 1:group_num){
    if(type == 1){
      Theta_pre<-input_array_hat[g,,]
      corr_pre<-(-Theta_pre) / (matrix(sqrt(diag(Theta_pre)),ncol = 1) %*% matrix(sqrt(diag(Theta_pre)),nrow = 1))
      diag(corr_pre)<-1 
    }
    if(type == 2){
      corr_pre<-input_array_hat[g,,]
    }
    if(type == 3){
      weight_pre0<-input_array_hat[g,,]
      
      corr_pre<-ifelse(abs(weight_pre0)>abs(t(weight_pre0)),abs(weight_pre0),abs(t(weight_pre0)))
    }
    corr_pre_upper<-(corr_pre[upper.tri(corr_pre)])[uppertri_choose]
    ##
    corr_pre_upper[(order(abs(corr_pre_upper),decreasing = FALSE))[1:ceiling((1-density_choose) * length(corr_pre_upper))]]<-0
    ##
    Theta_true<-Theta_array_true[cluster_predict_true_vec[g],,]
    corr_true<-(-Theta_true) / (matrix(sqrt(diag(Theta_true)),ncol = 1) %*% matrix(sqrt(diag(Theta_true)),nrow = 1))
    diag(corr_true)<-1
    corr_true_upper<-(corr_true[upper.tri(corr_true)])[uppertri_choose]
    ##
    edge_num_pre<-length(which(corr_pre_upper != 0))
    edge_num_true<-length(which(corr_true_upper !=0))
    edge_density_true<-edge_num_true/length(corr_true_upper)
    edge_order_pre<-order(abs(corr_pre_upper),decreasing = TRUE)[1:edge_num_pre]
    edge_true<-which(corr_true_upper !=0)
    ##
    ##TPR PRE
    TPR_vec<-c()
    PRE_vec<-c()
    for(l in 1:length(edge_order_pre)){
      PRE_vec<-c(PRE_vec,length(which(corr_true_upper[edge_order_pre[1:l]] != 0))/l)
      TPR_vec<-c(TPR_vec,length(which(corr_true_upper[edge_order_pre[1:l]] != 0))/edge_num_true)
    }
    
    ## add (TPR PRE) by random predictor
    TPR_vec_add<-c()
    PRE_vec_add<-c()
    K_max<-length(corr_pre_upper) - edge_num_pre
    if(K_max>0){
      TPR_vec_add<-TPR_vec[length(TPR_vec)] + (1:K_max) * ((edge_num_true - length(which(corr_true_upper[edge_order_pre]!=0)))/((edge_num_true - length(which(corr_true_upper[edge_order_pre]!=0))) + (length(which(corr_true_upper ==0)) - length(which(corr_true_upper[edge_order_pre]==0)))))/edge_num_true
      
      PRE_vec_add<-(length(which(corr_true_upper[edge_order_pre]!=0)) + (((edge_num_true - length(which(corr_true_upper[edge_order_pre]!=0)))/((edge_num_true - length(which(corr_true_upper[edge_order_pre]!=0))) + (length(which(corr_true_upper ==0)) - length(which(corr_true_upper[edge_order_pre]==0)))))) * (1:K_max))/(length(edge_order_pre)+(1:K_max))
    }
    TPR_vec_whole<-c(TPR_vec,TPR_vec_add)
    PRE_vec_whole<-c(PRE_vec,PRE_vec_add)
    ##
    TPR_vec<-c(0,TPR_vec)
    PRE_vec<-c(1,PRE_vec)
    TPR_vec_whole<-c(0,TPR_vec_whole)
    PRE_vec_whole<-c(1,PRE_vec_whole)
    
    #AUPRC (ratio)
    AUPRC<-sum(diff(TPR_vec) * (PRE_vec[-1] + PRE_vec[-length(PRE_vec)])/2)
    AUPRC_whole<-sum(diff(TPR_vec_whole) * (PRE_vec_whole[-1] + PRE_vec_whole[-length(PRE_vec_whole)])/2)
    res_mat[g,1]<-AUPRC
    res_mat[g,2]<-AUPRC_whole
    AUPRC_random<-sum(diff(TPR_vec) * rep(edge_density_true,length(TPR_vec) - 1))
    AUPRC_whole_random<-sum(diff(TPR_vec_whole) * rep(edge_density_true,length(TPR_vec_whole) - 1))
    AUPRC_ratio<-AUPRC/AUPRC_random
    AUPRC_ratio_whole<-AUPRC_whole/AUPRC_whole_random
    res_mat[g,3]<-AUPRC_ratio
    res_mat[g,4]<-AUPRC_ratio_whole
    ##
    
    ##Early Precision (ratio) & Early TPR
    edge_num_choose<-min(edge_num_true,edge_num_pre)
    EP<-PRE_vec[edge_num_choose + 1]
    res_mat[g,5]<-EP
    EP_random<-edge_density_true
    EP_ratio<-EP/EP_random
    res_mat[g,6]<-EP_ratio
    ETPR<-TPR_vec[edge_num_choose + 1]
    res_mat[g,7]<-ETPR
    
    ##
    edge_num_choose<-edge_num_true
    EP<-PRE_vec_whole[edge_num_choose + 1]
    res_mat[g,8]<-EP
    EP_random<-edge_density_true
    EP_ratio<-EP/EP_random
    res_mat[g,9]<-EP_ratio
    ETPR<-TPR_vec_whole[edge_num_choose + 1]
    res_mat[g,10]<-ETPR
    ##
    ##TPR
    TPR<-TPR_vec[length(TPR_vec)]
    res_mat[g,11]<-TPR
    ##PRE
    PRE<-PRE_vec[length(PRE_vec)]
    res_mat[g,12]<-PRE
    
    
  }
  ##
  return(res_mat)
}

#evaluator_1_sign
evaluator_1_signed<-function(input_array_hat,
                      type = 1,
                      Theta_array_true,
                      cluster_predict_true_vec,
                      uppertri_choose = NULL
){
  #####
  ##input_array_hat: group_num * p * p
  ##Theta_array_true: group_num * p * p
  ## performance criterion: signed AUROC AUPRC EP TPR TDR
  #####
  group_num<-length(cluster_predict_true_vec)
  res_mat<-matrix(NA,nrow = group_num,ncol = 12)
  if(is.null(uppertri_choose)){
    uppertri_choose <- 1:((dim(input_array_hat)[2]) * (dim(input_array_hat)[2] - 1)/2)
  }
  ##
  if(type %in% c(1,2)){
    for(g in 1:group_num){
      if(type == 1){
        Theta_pre<-input_array_hat[g,,]
        corr_pre<-(-Theta_pre) / (matrix(sqrt(diag(Theta_pre)),ncol = 1) %*% matrix(sqrt(diag(Theta_pre)),nrow = 1))
        diag(corr_pre)<-1 
      }
      if(type == 2){
        corr_pre<-input_array_hat[g,,]
      }
      # if(type == 3){
      #   weight_pre0<-input_array_hat[g,,]
      #   
      #   corr_pre<-ifelse(abs(weight_pre0)>abs(t(weight_pre0)),abs(weight_pre0),abs(t(weight_pre0)))
      # }
      corr_pre_upper<-(corr_pre[upper.tri(corr_pre)])[uppertri_choose]
      ##
      ##
      Theta_true<-Theta_array_true[cluster_predict_true_vec[g],,]
      corr_true<-(-Theta_true) / (matrix(sqrt(diag(Theta_true)),ncol = 1) %*% matrix(sqrt(diag(Theta_true)),nrow = 1))
      diag(corr_true)<-1
      corr_true_upper<-(corr_true[upper.tri(corr_true)])[uppertri_choose]
      ##
      edge_num_pre<-length(which(corr_pre_upper != 0))
      edge_num_true<-length(which(corr_true_upper !=0))
      edge_density_true<-edge_num_true/length(corr_true_upper)
      edge_order_pre<-order(abs(corr_pre_upper),decreasing = TRUE)[1:edge_num_pre]
      edge_true<-which(corr_true_upper !=0)
      ##
      ##TPR PRE
      TPR_vec<-c()
      PRE_vec<-c()
      for(l in 1:length(edge_order_pre)){
        PRE_vec<-c(PRE_vec,length(which(sign(corr_true_upper[edge_order_pre[1:l]]) == sign(corr_pre_upper[edge_order_pre[1:l]])))/l)
        TPR_vec<-c(TPR_vec,length(which(sign(corr_true_upper[edge_order_pre[1:l]]) == sign(corr_pre_upper[edge_order_pre[1:l]])))/edge_num_true)
      }
      
      ## add (TPR PRE) by random predictor
      TPR_vec_add<-c()
      PRE_vec_add<-c()
      K_max<-length(corr_pre_upper) - edge_num_pre
      if(K_max>0){
        TPR_vec_add<-TPR_vec[length(TPR_vec)] + (1:K_max) * ((edge_num_true - length(which(corr_true_upper[edge_order_pre]!=0)))/((edge_num_true - length(which(corr_true_upper[edge_order_pre]!=0))) + (length(which(corr_true_upper ==0)) - length(which(corr_true_upper[edge_order_pre]==0)))))/edge_num_true
        
        PRE_vec_add<-(length(which(corr_true_upper[edge_order_pre]!=0)) + (((edge_num_true - length(which(corr_true_upper[edge_order_pre]!=0)))/((edge_num_true - length(which(corr_true_upper[edge_order_pre]!=0))) + (length(which(corr_true_upper ==0)) - length(which(corr_true_upper[edge_order_pre]==0)))))) * (1:K_max))/(length(edge_order_pre)+(1:K_max))
      }
      TPR_vec_whole<-c(TPR_vec,TPR_vec_add)
      PRE_vec_whole<-c(PRE_vec,PRE_vec_add)
      ##
      TPR_vec<-c(0,TPR_vec)
      PRE_vec<-c(1,PRE_vec)
      TPR_vec_whole<-c(0,TPR_vec_whole)
      PRE_vec_whole<-c(1,PRE_vec_whole)
      
      #AUPRC (ratio)
      AUPRC<-sum(diff(TPR_vec) * (PRE_vec[-1] + PRE_vec[-length(PRE_vec)])/2)
      AUPRC_whole<-sum(diff(TPR_vec_whole) * (PRE_vec_whole[-1] + PRE_vec_whole[-length(PRE_vec_whole)])/2)
      res_mat[g,1]<-AUPRC
      res_mat[g,2]<-AUPRC_whole
      AUPRC_random<-sum(diff(TPR_vec) * rep((edge_density_true * 0.5),length(TPR_vec) - 1))
      AUPRC_whole_random<-sum(diff(TPR_vec_whole) * rep((edge_density_true * 0.5),length(TPR_vec_whole) - 1))
      AUPRC_ratio<-AUPRC/AUPRC_random
      AUPRC_ratio_whole<-AUPRC_whole/AUPRC_whole_random
      res_mat[g,3]<-AUPRC_ratio
      res_mat[g,4]<-AUPRC_ratio_whole
      ##
      
      ##Early Precision (ratio) & Early TPR
      edge_num_choose<-min(edge_num_true,edge_num_pre)
      EP<-PRE_vec[edge_num_choose + 1]
      res_mat[g,5]<-EP
      EP_random<-edge_density_true
      EP_ratio<-EP/EP_random
      res_mat[g,6]<-EP_ratio
      ETPR<-TPR_vec[edge_num_choose + 1]
      res_mat[g,7]<-ETPR
      
      ##
      edge_num_choose<-edge_num_true
      EP<-PRE_vec_whole[edge_num_choose + 1]
      res_mat[g,8]<-EP
      EP_random<-edge_density_true
      EP_ratio<-EP/EP_random
      res_mat[g,9]<-EP_ratio
      ETPR<-TPR_vec_whole[edge_num_choose + 1]
      res_mat[g,10]<-ETPR
      ##
      ##TPR
      TPR<-TPR_vec[length(TPR_vec)]
      res_mat[g,11]<-TPR
      ##PRE
      PRE<-PRE_vec[length(PRE_vec)]
      res_mat[g,12]<-PRE
      
      
    } 
  }
  ##
  return(res_mat)
}

##evalulator_early_stability
evalulator_early_stability<-function(edge_list,
                                     cluster_predict_true_mat
){
  group_num<-ncol(cluster_predict_true_mat)
  reptime<-length(edge_list)
  Jaccard_mat<-matrix(NA,nrow = group_num,ncol = (reptime * (reptime - 1)/2))
  for(g in 1:group_num){
    ##
    edgeindex_group_list<-list()
    for(rep in 1:reptime){
      edgeindex_group_list[[rep]]<-edge_list[[rep]][[which(cluster_predict_true_mat[rep,] == g)]]
    }
    ##
    Jaccard_vec<-c()
    for(rep1 in 1:(reptime - 1)){
      for(rep2 in (rep1+1):reptime){
        dom_value<-length(union(edgeindex_group_list[[rep1]],edgeindex_group_list[[rep2]]))
        if(dom_value == 0){
          Jaccard_vec<-c(Jaccard_vec,NA)
        }else{
          Jaccard_vec<-c(Jaccard_vec,length(intersect(edgeindex_group_list[[rep1]],edgeindex_group_list[[rep2]]))/dom_value)
        }
      }
    }
    ##
    Jaccard_mat[g,]<-Jaccard_vec
  }
  return(Jaccard_mat)
}


##prior_generator
prior_generator<-function(adjoint_mat,eta=0.3){
  edge_index<-which(adjoint_mat[upper.tri(adjoint_mat)]==1)
  edge_index_sel<-sample(edge_index,floor(length(edge_index) * eta),replace = FALSE)
  prior_adjoint<-diag(nrow(adjoint_mat))
  prior_adjoint[upper.tri(prior_adjoint)][edge_index_sel]<-1
  prior_adjoint<-prior_adjoint + t(prior_adjoint) - diag(nrow(adjoint_mat))
  ##
  return(prior_adjoint)
}


##data_generator_new
data_generator<-function(n,dim_use,group_num,
                             network_type=c("ER-graph","AFF-graph","PA-graph","HUB-graph"),
                             densy_degree = 0.1,v = 0.3,p_v = 0.5,
                             num_de = 15,ub = 1,lb = -3.5, ub_non = 1, lb_non = -1,
                             l_mu = log(10),l_sd = 0.05,sigma_scale = 1,
                             de_ratio_bet = 0.5,alpha_com = 0.01,
                             if_fix = FALSE,
                             fix_element = NULL){
  if(if_fix == FALSE){
    #Generate Group-specific adjoint matrix & precision
    group_specific_adjoint<-array(0,dim=c(group_num,dim_use,dim_use))
    group_specific_precision<-array(0,dim=c(group_num,dim_use,dim_use))
    hub_index<-NULL
    if(network_type == "HUB-graph"){
      hub_index<-sample(1:dim_use, 2 * dim_use * densy_degree, replace = FALSE)
    }
    ##adjoint matrix
    ##adjoint
    if(network_type == "ER-graph"){
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(1:(dim_use * (dim_use - 1)/2),floor(dim_use * (dim_use - 1)/2 * densy_degree * alpha_com),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-which(common_adjoint[upper.tri(common_adjoint)]==0)
      for(r in 1:group_num){
        adjoin_mat_choose<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat_choose[upper.tri(adjoin_mat_choose)][sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1-alpha_com)),replace = FALSE)]<-1
        adjoin_mat_choose<-adjoin_mat_choose+t(adjoin_mat_choose)
        
        group_specific_adjoint[r,,]<- adjoin_mat_choose + common_adjoint
      }
    }
    
    if(network_type == "AFF-graph"){
      community_label<-as.vector(matrix(rep(1:4,each = dim_use/4),ncol = 4,byrow = FALSE))
      adjoin_mat_com<-diag(dim_use)
      for(k in 1:4){
        adjoin_mat_com[which(community_label == k),which(community_label == k)]<-1
      }
      edge_index<-which(adjoin_mat_com[upper.tri(adjoin_mat_com)]==1)
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(edge_index,floor(dim_use * (dim_use - 1)/2 * densy_degree * alpha_com),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-setdiff(edge_index,which(common_adjoint[upper.tri(common_adjoint)] == 1))
      for(r in 1:group_num){
        adjoin_mat_choose<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat_choose[upper.tri(adjoin_mat_choose)][sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1 - alpha_com)),replace = FALSE)]<-1
        adjoin_mat_choose<-adjoin_mat_choose+t(adjoin_mat_choose)
        group_specific_adjoint[r,,]<-adjoin_mat_choose + common_adjoint 
      }
    }
    if(network_type == "PA-graph"){
      #calculate the densy of single network
      g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
      adjoin_mat<-as_adjacency_matrix(g1)
      adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat)) + diag(dim_use)
      ##
      num_common<-floor(densy_degree * alpha_com/as.numeric(prop.table(table(adjoin_mat[upper.tri(adjoin_mat)]))[2]))
      num_noncommon<-floor(densy_degree/as.numeric(prop.table(table(adjoin_mat[upper.tri(adjoin_mat)]))[2])) - num_common
      
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      for(j in 1:num_common){
        g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
        adjoin_mat<-as_adjacency_matrix(g1)
        adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat))
        common_adjoint <- common_adjoint + adjoin_mat
      }
      common_adjoint<-ifelse(common_adjoint!=0,1,0) + diag(dim_use)
      
      for(r in 1:group_num){
        adjoin_mat_group<-matrix(0,nrow = dim_use,ncol = dim_use)
        for(j in 1:num_noncommon){
          g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
          adjoin_mat<-as_adjacency_matrix(g1)
          adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat))
          adjoin_mat_group <- adjoin_mat_group + adjoin_mat
        }
        group_specific_adjoint[r,,]<-ifelse(adjoin_mat_group + common_adjoint!=0,1,0)
      }
      
    }
    if(network_type == "HUB-graph"){
      adjoin_mat_choose<-matrix(1,nrow = dim_use,ncol = dim_use)
      adjoin_mat_choose[-hub_index,-hub_index]<-0
      diag(adjoin_mat_choose)<-1
      edge_choose<-which(adjoin_mat_choose[upper.tri(adjoin_mat_choose)]==1)
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(edge_choose,floor(dim_use * (dim_use - 1)/2 * densy_degree * (alpha_com)),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-setdiff(edge_choose,which(common_adjoint[upper.tri(common_adjoint)] == 1))
      
      for(r in 1:group_num){
        edge_choose1<-sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1-alpha_com)),replace = FALSE)
        adjoin_mat<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat[upper.tri(adjoin_mat)][edge_choose1]<-1
        adjoin_mat<-adjoin_mat + t(adjoin_mat)
        group_specific_adjoint[r,,]<-adjoin_mat + common_adjoint
      }
    }
    
    #precision
    ##1. common part
    # rbinom_mat<-diag(dim_use)
    # rbinom_mat[upper.tri(rbinom_mat)]<-(rbinom(dim_use * (dim_use-1)/2,1,1-p_v)-0.5)*2*v
    # rbinom_mat<-rbinom_mat + t(rbinom_mat) - diag(dim_use)
    # pre_common <- common_adjoint * rbinom_mat
    # common_index<-which(pre_common[upper.tri(pre_common)]!=0)
    for(r in 1:group_num){
      # diag_value<-1
      # group_specific_precision0<-matrix(0,nrow = dim_use,ncol = dim_use)
      # edge_allow<-setdiff(which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1),common_index)
      # group_specific_precision0[upper.tri(group_specific_precision0)][edge_allow]<-(rbinom(length(edge_allow),1,1 - p_v)-0.5)*2*v
      # group_specific_precision0<-group_specific_precision0 + t(group_specific_precision0)
      # group_specific_precision0<-group_specific_precision0 + pre_common
      # con_pos<-TRUE
      # while (con_pos) {
      #   group_specific_precision1<-group_specific_precision0
      #   diag(group_specific_precision1)<-diag_value
      #   
      #   eigen_value<-eigen(group_specific_precision1)$values
      #   if(min(eigen_value)>0){
      #     con_pos<-FALSE
      #     group_specific_precision[r,,]<-group_specific_precision1 * sigma_scale
      #   }else{
      #     diag_value<-diag_value * 1.2
      #   }
      # }
      #############################
      diag_value<-1
      group_specific_precision0<-diag(dim_use)
      group_specific_precision0[upper.tri(group_specific_precision0)][which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1)]<-(rbinom(length(which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1))
                                                                                                                                                              ,1,1 - p_v)-0.5)*2*v
      group_specific_precision0<-group_specific_precision0 + t(group_specific_precision0) - diag(dim_use)
      eigen_value<-eigen(group_specific_precision0)$values
      if(min(eigen_value)<0){
        group_specific_precision0<-group_specific_precision0 + diag(dim_use) * (abs(min(eigen_value))+1e-1)
      }else{
        group_specific_precision0<-group_specific_precision0 + diag(dim_use) * (1e-1)
      }
      
      group_specific_precision[r,,]<-group_specific_precision0 * sigma_scale
      # ##
      # con_pos<-TRUE
      # while (con_pos) {
      #   group_specific_precision1<-group_specific_precision0
      #   diag(group_specific_precision1)<-diag_value
      # 
      #   eigen_value<-eigen(group_specific_precision1)$values
      #   if(min(eigen_value)>0){
      #     con_pos<-FALSE
      #     group_specific_precision[r,,]<-group_specific_precision1 * sigma_scale
      #   }else{
      #     diag_value<-diag_value * 1.2
      #   }
      # }
      # ##
      
    }
    
    #Generate the sc-RNA data
    #1. mu_mat
    differential_gene_index<-sample(1:dim_use,num_de,replace = FALSE)
    nondifferential_gene_index<-setdiff(1:dim_use,differential_gene_index)
    mu_mat<-matrix(0,nrow = group_num,ncol = dim_use)
    mu_nondifferent<-runif(length(nondifferential_gene_index),min = lb_non,max = ub_non)
    for(r in 1:group_num){
      mu_mat[r,nondifferential_gene_index]<-mu_nondifferent
    }
    for(j in 1:length(differential_gene_index)){
      mu_mat[,differential_gene_index[j]]<-sample(c(lb,ub,(lb+ub)/2),group_num,replace = TRUE)
    }
    # group_num_true = group_num
    # # group_num = max(5,group_num_true)
    # rep_con<-TRUE
    # while (rep_con) {
    #   AAA<-matrix(NA,nrow = group_num,ncol = num_de)
    #   for(j in 1:num_de){
    #     aaa_vec<-sample(c(1:3),group_num,replace = TRUE)
    #     while(var(aaa_vec) == 0){
    #       aaa_vec[sample(1:group_num,1)]<-sample(c(1:3),1,replace = TRUE)
    #     }
    #     
    #     AAA[,j]<-aaa_vec
    #   }
    #   con_1<-TRUE
    #   count<-0
    #   while(con_1 & (count<20)){
    #     count<-count+1
    #     con_vec<-c()
    #     for(g in 1:group_num){
    #       if(min(rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))<(de_ratio_bet * num_de)){
    #         choose_index<-sample(which(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))[which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),]==0),num_de * de_ratio_bet - min(rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,])))))))
    #         for(i in 1:length(choose_index)){
    #           AAA[-g,][which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),choose_index[i]]<-sample(setdiff(c(1:3),AAA[-g,][which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),choose_index[i]]),1)
    #         }
    #         con_vec<-c(con_vec,1)
    #       }else{
    #         con_vec<-c(con_vec,0)
    #       }
    #     }
    #     if(length(which(con_vec==1))==0){
    #       con_1<-FALSE
    #     }
    #   }
    #   if(con_1 == TRUE){
    #     rep_con<-TRUE
    #   }else{
    #     rep_con<-FALSE
    #   }
    # }
    # AAA = AAA[1:group_num_true,]
    # group_num = group_num_true
    # value_choose<-seq(from = lb, to = ub,length.out = 3)
    # mu_matde<-ifelse(AAA==1,value_choose[1],ifelse(AAA==2,value_choose[2],
    #                                                value_choose[3]))
    # 
    # mu_mat[,differential_gene_index]<-mu_matde
  }else{
    group_specific_precision = fix_element$group_specific_precision
    common_adjoint = fix_element$common_adjoint
    group_specific_adjoint = fix_element$group_specific_adjoint
    mu_mat = fix_element$mu_mat
    hub_gene = fix_element$hub_gene
    ################
    n<-nrow(fix_element$obs_mat)
    group_num<-nrow(mu_mat)
    dim_use<-ncol(fix_element$obs_mat)
    hub_index<-NULL
    if(network_type == "HUB-graph"){
      hub_index<-sample(1:dim_use, 2 * dim_use * densy_degree, replace = FALSE)
    }
  }
  
  
  #5.2 library size
  ls_vec<-exp(rnorm(n,mean = log(l_mu),sd=l_sd))
  
  #5.3 pi
  pi_vec<-rep(1/group_num,group_num)
  # print("parameter design complete.")
  #
  #5.4 cluster_lable & obs_mat generation
  locator_vec<-c()
  obs_mat1<-matrix(0,nrow = n,ncol = dim_use)
  inv_group_pre<-array(NA,dim = dim(group_specific_precision))
  for(g in 1:group_num){
    inv_group_pre[g,,]<-solve(group_specific_precision[g,,])
  }
  ###########################
  locator_vec<-sample(1:group_num,n,prob = pi_vec,replace = TRUE)
  log_a_mat<-matrix(NA,nrow = n,ncol = dim_use)
  for(g in 1:group_num){
    log_a_mat[which(locator_vec == g),]<-mvrnorm(n = length(which(locator_vec == g)), (mu_mat[g,]), inv_group_pre[g,,])
  }
  log_a_mat<-log_a_mat + log(ls_vec)
  
  log_a_mat<-ifelse(log_a_mat>10,10,log_a_mat)
  a_mat<-exp(log_a_mat)
  ###########################
  for(i in 1:n){
    for(j in 1:dim_use){
      obs_mat1[i,j]<-rpois(1,a_mat[i,j])
    }
  }
  colnames(obs_mat1)<-paste("Gene",1:ncol(obs_mat1),sep = "")
  rownames(obs_mat1)<-paste("Cell",1:nrow(obs_mat1),sep = "")
  hub_gene<-colnames(obs_mat1)[hub_index]
  return(list(obs_mat = obs_mat1,locator_vec = locator_vec,
              group_specific_precision = group_specific_precision,
              common_adjoint = common_adjoint,
              group_specific_adjoint = group_specific_adjoint,
              ls_vec = ls_vec,mu_mat = mu_mat,hub_gene = hub_gene))
  
}

##
##data_generator_mulitnomial
data_generator_mulitnomial<-function(n,dim_use,group_num,
                                         network_type=c("ER-graph","AFF-graph","PA-graph","HUB-graph"),
                                         densy_degree = 0.1,v = 0.3,p_v = 0.5,
                                         num_de = 15,ub = 1,lb = -3.5, ub_non = 1, lb_non = -1,
                                         l_mu = log(10),l_sd = 0.05,sigma_scale = 1,
                                         de_ratio_bet = 0.5,alpha_com = 0.2,
                                         if_fix = FALSE,
                                         fix_element = NULL){
  if(if_fix == FALSE){
    #Generate Group-specific adjoint matrix & precision
    group_specific_adjoint<-array(0,dim=c(group_num,dim_use,dim_use))
    group_specific_precision<-array(0,dim=c(group_num,dim_use,dim_use))
    hub_index<-NULL
    if(network_type == "HUB-graph"){
      hub_index<-sample(1:dim_use, 2 * dim_use * densy_degree, replace = FALSE)
    }
    ##adjoint matrix
    ##adjoint
    if(network_type == "ER-graph"){
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(1:(dim_use * (dim_use - 1)/2),floor(dim_use * (dim_use - 1)/2 * densy_degree * alpha_com),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-which(common_adjoint[upper.tri(common_adjoint)]==0)
      for(r in 1:group_num){
        adjoin_mat_choose<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat_choose[upper.tri(adjoin_mat_choose)][sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1-alpha_com)),replace = FALSE)]<-1
        adjoin_mat_choose<-adjoin_mat_choose+t(adjoin_mat_choose)
        
        group_specific_adjoint[r,,]<- adjoin_mat_choose + common_adjoint
      }
    }
    
    if(network_type == "AFF-graph"){
      community_label<-as.vector(matrix(rep(1:4,each = dim_use/4),ncol = 4,byrow = FALSE))
      adjoin_mat_com<-diag(dim_use)
      for(k in 1:4){
        adjoin_mat_com[which(community_label == k),which(community_label == k)]<-1
      }
      edge_index<-which(adjoin_mat_com[upper.tri(adjoin_mat_com)]==1)
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(edge_index,floor(dim_use * (dim_use - 1)/2 * densy_degree * alpha_com),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-setdiff(edge_index,which(common_adjoint[upper.tri(common_adjoint)] == 1))
      for(r in 1:group_num){
        adjoin_mat_choose<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat_choose[upper.tri(adjoin_mat_choose)][sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1 - alpha_com)),replace = FALSE)]<-1
        adjoin_mat_choose<-adjoin_mat_choose+t(adjoin_mat_choose)
        group_specific_adjoint[r,,]<-adjoin_mat_choose + common_adjoint 
      }
    }
    if(network_type == "PA-graph"){
      #calculate the densy of single network
      g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
      adjoin_mat<-as_adjacency_matrix(g1)
      adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat)) + diag(dim_use)
      ##
      num_common<-floor(densy_degree * alpha_com/as.numeric(prop.table(table(adjoin_mat[upper.tri(adjoin_mat)]))[2]))
      num_noncommon<-floor(densy_degree/as.numeric(prop.table(table(adjoin_mat[upper.tri(adjoin_mat)]))[2])) - num_common
      
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      for(j in 1:num_common){
        g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
        adjoin_mat<-as_adjacency_matrix(g1)
        adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat))
        common_adjoint <- common_adjoint + adjoin_mat
      }
      common_adjoint<-ifelse(common_adjoint!=0,1,0) + diag(dim_use)
      
      for(r in 1:group_num){
        adjoin_mat_group<-matrix(0,nrow = dim_use,ncol = dim_use)
        for(j in 1:num_noncommon){
          g1 <- sample_pa_age(dim_use, pa.exp=0.5, aging.exp=-1, aging.bin=1000)
          adjoin_mat<-as_adjacency_matrix(g1)
          adjoin_mat<-as.matrix(adjoin_mat) + t(as.matrix(adjoin_mat))
          adjoin_mat_group <- adjoin_mat_group + adjoin_mat
        }
        group_specific_adjoint[r,,]<-ifelse(adjoin_mat_group + common_adjoint!=0,1,0)
      }
      
    }
    if(network_type == "HUB-graph"){
      adjoin_mat_choose<-matrix(1,nrow = dim_use,ncol = dim_use)
      adjoin_mat_choose[-hub_index,-hub_index]<-0
      diag(adjoin_mat_choose)<-1
      edge_choose<-which(adjoin_mat_choose[upper.tri(adjoin_mat_choose)]==1)
      ##common_adjoint
      common_adjoint<-matrix(0,nrow = dim_use,ncol = dim_use)
      common_adjoint[upper.tri(common_adjoint)][sample(edge_choose,floor(dim_use * (dim_use - 1)/2 * densy_degree * (alpha_com)),replace = FALSE)]<-1
      common_adjoint <- common_adjoint + t(common_adjoint) + diag(dim_use)
      avil_edge<-setdiff(edge_choose,which(common_adjoint[upper.tri(common_adjoint)] == 1))
      
      for(r in 1:group_num){
        edge_choose1<-sample(avil_edge,floor(dim_use * (dim_use - 1)/2 * densy_degree * (1-alpha_com)),replace = FALSE)
        adjoin_mat<-matrix(0,nrow = dim_use,ncol = dim_use)
        adjoin_mat[upper.tri(adjoin_mat)][edge_choose1]<-1
        adjoin_mat<-adjoin_mat + t(adjoin_mat)
        group_specific_adjoint[r,,]<-adjoin_mat + common_adjoint
      }
    }
    
    #precision
    ##1. common part
    # rbinom_mat<-diag(dim_use)
    # rbinom_mat[upper.tri(rbinom_mat)]<-(rbinom(dim_use * (dim_use-1)/2,1,1-p_v)-0.5)*2*v
    # rbinom_mat<-rbinom_mat + t(rbinom_mat) - diag(dim_use)
    # pre_common <- common_adjoint * rbinom_mat
    # common_index<-which(pre_common[upper.tri(pre_common)]!=0)
    for(r in 1:group_num){
      # diag_value<-1
      # group_specific_precision0<-matrix(0,nrow = dim_use,ncol = dim_use)
      # edge_allow<-setdiff(which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1),common_index)
      # group_specific_precision0[upper.tri(group_specific_precision0)][edge_allow]<-(rbinom(length(edge_allow),1,1 - p_v)-0.5)*2*v
      # group_specific_precision0<-group_specific_precision0 + t(group_specific_precision0)
      # group_specific_precision0<-group_specific_precision0 + pre_common
      # con_pos<-TRUE
      # while (con_pos) {
      #   group_specific_precision1<-group_specific_precision0
      #   diag(group_specific_precision1)<-diag_value
      #   
      #   eigen_value<-eigen(group_specific_precision1)$values
      #   if(min(eigen_value)>0){
      #     con_pos<-FALSE
      #     group_specific_precision[r,,]<-group_specific_precision1 * sigma_scale
      #   }else{
      #     diag_value<-diag_value * 1.2
      #   }
      # }
      #############################
      diag_value<-1
      group_specific_precision0<-diag(dim_use)
      group_specific_precision0[upper.tri(group_specific_precision0)][which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1)]<-(rbinom(length(which(group_specific_adjoint[r,,][upper.tri(group_specific_adjoint[r,,])]==1))
                                                                                                                                                              ,1,1 - p_v)-0.5)*2*v
      group_specific_precision0<-group_specific_precision0 + t(group_specific_precision0) - diag(dim_use)
      eigen_value<-eigen(group_specific_precision0)$values
      if(min(eigen_value)<0){
        group_specific_precision0<-group_specific_precision0 + diag(dim_use) * (abs(min(eigen_value))+1e-1)
      }else{
        group_specific_precision0<-group_specific_precision0 + diag(dim_use) * (1e-1)
      }
      
      group_specific_precision[r,,]<-group_specific_precision0 * sigma_scale
      # ##
      # con_pos<-TRUE
      # while (con_pos) {
      #   group_specific_precision1<-group_specific_precision0
      #   diag(group_specific_precision1)<-diag_value
      # 
      #   eigen_value<-eigen(group_specific_precision1)$values
      #   if(min(eigen_value)>0){
      #     con_pos<-FALSE
      #     group_specific_precision[r,,]<-group_specific_precision1 * sigma_scale
      #   }else{
      #     diag_value<-diag_value * 1.2
      #   }
      # }
      # ##
      
    }
    
    #Generate the sc-RNA data
    #1. mu_mat
    differential_gene_index<-sample(1:dim_use,num_de,replace = FALSE)
    nondifferential_gene_index<-setdiff(1:dim_use,differential_gene_index)
    mu_mat<-matrix(0,nrow = group_num,ncol = dim_use)
    mu_nondifferent<-runif(length(nondifferential_gene_index),min = lb_non,max = ub_non)
    for(r in 1:group_num){
      mu_mat[r,nondifferential_gene_index]<-mu_nondifferent
    }
    for(j in 1:length(differential_gene_index)){
      mu_mat[,differential_gene_index[j]]<-sample(c(lb,ub,(lb+ub)/2),group_num,replace = TRUE)
    }
    # group_num_true = group_num
    # # group_num = max(5,group_num_true)
    # rep_con<-TRUE
    # while (rep_con) {
    #   AAA<-matrix(NA,nrow = group_num,ncol = num_de)
    #   for(j in 1:num_de){
    #     aaa_vec<-sample(c(1:3),group_num,replace = TRUE)
    #     while(var(aaa_vec) == 0){
    #       aaa_vec[sample(1:group_num,1)]<-sample(c(1:3),1,replace = TRUE)
    #     }
    #     
    #     AAA[,j]<-aaa_vec
    #   }
    #   con_1<-TRUE
    #   count<-0
    #   while(con_1 & (count<20)){
    #     count<-count+1
    #     con_vec<-c()
    #     for(g in 1:group_num){
    #       if(min(rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))<(de_ratio_bet * num_de)){
    #         choose_index<-sample(which(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))[which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),]==0),num_de * de_ratio_bet - min(rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,])))))))
    #         for(i in 1:length(choose_index)){
    #           AAA[-g,][which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),choose_index[i]]<-sample(setdiff(c(1:3),AAA[-g,][which.min((rowSums(sign(abs(t(t(AAA[-g,]) - as.vector(AAA[g,]))))))),choose_index[i]]),1)
    #         }
    #         con_vec<-c(con_vec,1)
    #       }else{
    #         con_vec<-c(con_vec,0)
    #       }
    #     }
    #     if(length(which(con_vec==1))==0){
    #       con_1<-FALSE
    #     }
    #   }
    #   if(con_1 == TRUE){
    #     rep_con<-TRUE
    #   }else{
    #     rep_con<-FALSE
    #   }
    # }
    # AAA = AAA[1:group_num_true,]
    # group_num = group_num_true
    # value_choose<-seq(from = lb, to = ub,length.out = 3)
    # mu_matde<-ifelse(AAA==1,value_choose[1],ifelse(AAA==2,value_choose[2],
    #                                                value_choose[3]))
    # 
    # mu_mat[,differential_gene_index]<-mu_matde
  }else{
    group_specific_precision = fix_element$group_specific_precision
    common_adjoint = fix_element$common_adjoint
    group_specific_adjoint = fix_element$group_specific_adjoint
    mu_mat = fix_element$mu_mat
    hub_gene = fix_element$hub_gene
    ################
    n<-nrow(fix_element$obs_mat)
    group_num<-nrow(mu_mat)
    dim_use<-ncol(fix_element$obs_mat)
    hub_index<-NULL
    if(network_type == "HUB-graph"){
      hub_index<-sample(1:dim_use, 2 * dim_use * densy_degree, replace = FALSE)
    }
  }
  
  
  #5.2 library size
  ls_vec<-exp(rnorm(n,mean = log(l_mu),sd=l_sd))
  
  #5.3 pi
  pi_vec<-rep(1/group_num,group_num)
  # print("parameter design complete.")
  #
  #5.4 cluster_lable & obs_mat generation
  locator_vec<-c()
  obs_mat1<-matrix(0,nrow = n,ncol = dim_use)
  inv_group_pre<-array(NA,dim = dim(group_specific_precision))
  for(g in 1:group_num){
    inv_group_pre[g,,]<-solve(group_specific_precision[g,,])
  }
  ###########################
  locator_vec<-sample(1:group_num,n,prob = pi_vec,replace = TRUE)
  log_a_mat<-matrix(NA,nrow = n,ncol = dim_use)
  for(g in 1:group_num){
    log_a_mat[which(locator_vec == g),]<-mvrnorm(n = length(which(locator_vec == g)), (mu_mat[g,]), inv_group_pre[g,,])
  }
  log_a_mat<-log_a_mat + log(ls_vec)
  
  log_a_mat<-ifelse(log_a_mat>10,10,log_a_mat)
  a_mat<-exp(log_a_mat)
  ##
  N_vec<-ceiling(rowSums(a_mat))
  pi_count_mat<-a_mat/as.vector(rowSums(a_mat))
  ###########################
  for(i in 1:n){
    pseudo_sample<-sample(c(1:dim_use),size = N_vec[i],prob = as.vector(pi_count_mat[i,]),replace = TRUE)
    for(j in 1:dim_use){
      obs_mat1[i,j]<-length(which(pseudo_sample == j))
    }
    # ##
    # for(j in 1:dim_use){
    #   obs_mat1[i,j]<-rpois(1,a_mat[i,j])
    # }
  }
  colnames(obs_mat1)<-paste("Gene",1:ncol(obs_mat1),sep = "")
  rownames(obs_mat1)<-paste("Cell",1:nrow(obs_mat1),sep = "")
  hub_gene<-colnames(obs_mat1)[hub_index]
  return(list(obs_mat = obs_mat1,locator_vec = locator_vec,
              group_specific_precision = group_specific_precision,
              common_adjoint = common_adjoint,
              group_specific_adjoint = group_specific_adjoint,
              ls_vec = ls_vec,mu_mat = mu_mat,hub_gene = hub_gene))
  
}
##--------------------------------------------------------------------------------------------

##3. Define the evaluator function for each method
##--------------------------------------------------------------------------------------------
##Evaluator for scGeneNet
scGeneNet_evaluator<-function(penalize_diagonal = FALSE,
                          if_stability_eval = FALSE
){
  library_size_est = "TSS"
  ##
  repres_scGeneNet<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS','glasso',
                  'mclust','mltools','data.table','PLNmodels','dplyr',
                  "cluster","pryr","philentropy","umap","Seurat","igraph","Rcpp","scGeneNet","Matrix","stats","fpc","GMPR")
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    
    
    ##data & res_init
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    
    group_num<-ncol(res_init$U_mat)
    bic_object_list<-list()
    
    ## data initialization------------------
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    locator_base = apply(res_init$U_mat, 1, which.max)
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    ##
    purity_orgin<-sum(diag(table_reshape))/sum(table_reshape)
    ##
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    bic_object_list[["ARI_orgin"]]<-ARI_orgin
    bic_object_list[["purity_orgin"]]<-purity_orgin
    ##
    rm(cluster_predict_true_base);gc()
    rm(table_reshape);gc()
    
    ## method run-------------
    ######################################
    ## 1. prior information generation
    prior_nw_list<-list()
    for(g in 1:group_num){
      prior_nw_list[[g]] <- prior_generator(data$common_adjoint,eta = 0)
    }
    prior_list<-list(prior_nw_list = prior_nw_list,
                     prior_TF = NULL)
    rm(prior_nw_list);gc()
    ##
    prior_tran_res<-prior_tran(p_GRN = length(res_init$gene_GRN_index_use),
                               group_num = length(unique(res_init$celltypes_label)),
                               gene_GRN_name = res_init$gene_GRN_vec,
                               prior_nw_list = prior_list$prior_nw_list, prior_TF = prior_list$prior_TF,
                               weight_prior_nw = 0.3,penalize_diagonal = penalize_diagonal)
    weight_mat_list<-prior_tran_res$weight_mat_list
    zero_mat<-prior_tran_res$zero_mat
    rm(prior_tran_res)
    
    ## one step update for lambda = 1e-6 (U_fix = true)
    res_init<-scGeneNet_main(scGeneNet_list = res_init,lambda_use = 1e-6,
                           Theta_Control = list(penalize_diagonal = FALSE),
                           U_fix = TRUE,
                           verbose = FALSE,
                           core_num = 1)
    ##
    ######################################
    lambda_max_res <- exact_lambdamax(scGeneNet_list = res_init,
                                      Global_Control = list(ELBO_threshold = 1e-4, minit = 1, maxit = 50,
                                                            maxit_nonGRN = 10,MS_update_threshold = 1e-6),
                                      M_Control = list(ADMM_max_step = 1000, ADMM_threshold = 1e-4),
                                      S_Control = list(Newton_max_step = 1000, Newton_threshold = 1e-04),
                                      Theta_Control = list(penalize_diagonal = penalize_diagonal,
                                                           Theta_threshold = 1e-04),
                                      U_fix = FALSE,
                                      core_num = 1)
    ##
    lambda_min_res<-1e-6
    lambda_max_uni<-sort(unique(lambda_max_res))
    lambda_length0<-70
    lambda_length1<-40
    if(length(lambda_max_uni) == 1){
      
      lambda_vec<-exp(seq(from=log(lambda_min_res),to=log(lambda_max_uni[1]),length.out = lambda_length0))
      
      
    }else{
      lambda_max_uni<-c(lambda_min_res,lambda_max_uni)
      lambda_vec<-c()
      
      for(k in 1:(length(lambda_max_uni) - 1)){
        
        if(k == 1){
          lambda_vec<-c(lambda_vec,(seq(from=(lambda_max_uni[k]),to=(lambda_max_uni[k+1]),length.out = lambda_length0)))
        }else{
          lambda_vec<-c(lambda_vec,(seq(from=(lambda_max_uni[k]),to=(lambda_max_uni[k+1]),length.out = lambda_length1)))
          
        }
      }
    }
    print("lambda_vec has been determine.")
    
    ## 3. model run
    
    BIC_mat_shrink = matrix(NA,nrow = group_num,ncol = length(lambda_vec))
    BIC_mat_shrink1 = matrix(NA,nrow = group_num,ncol = length(lambda_vec))
    cluster_acc = NULL
    U_mat_array<-array(NA,dim = c(length(lambda_vec),nrow(obs_mat),group_num))
    Theta_mat_array_hat_all<-array(NA,dim = c(p,p,group_num,length(lambda_vec)))
    densy_mat<-matrix(NA,nrow = group_num,ncol = length(lambda_vec))
    soft_cs_vec<-c()
    purity_vec<-c()
    gc()
    ##
    
    
    time_cost<-c()
    res_scGeneNet_NN<-NULL
    for (l in 1:length(lambda_vec)) {
      print(paste("lambda(",l,") = ",lambda_vec[l],sep = ""))
      res_init0<-res_init
      ##
      rm(res_scGeneNet_NN);gc()
      ##
      time_MPLN_start<-Sys.time()
      res_scGeneNet_NN<-scGeneNet_main(scGeneNet_list = res_init0,lambda_use = lambda_vec[l],
                               Global_Control = list(ELBO_threshold = 1e-4, minit = 1, maxit = 50,maxit_nonGRN = 10, MS_update_threshold = 1e-6),
                               M_Control = list(ADMM_max_step = 1000, ADMM_threshold = 1e-4),
                               S_Control = list(Newton_max_step = 1000, Newton_threshold = 1e-4),
                               Theta_Control = list(penalize_diagonal = FALSE,Theta_threshold = 1e-4),
                               U_fix = FALSE,
                               verbose = TRUE,
                               core_num = 1)
    
      rm(res_init0);gc()
      time_MPLN_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_MPLN_end,time_MPLN_start,units = "hours")))
      
      BIC_mat_shrink[,l]<-res_scGeneNet_NN$scGeneNet_bic_VMICL
      BIC_mat_shrink1[,l]<-res_scGeneNet_NN$scGeneNet_bic_VICL
      ##
      
      for(g in 1:group_num){
        Theta_mat_array_hat_all[,,g,l]<-as.matrix(res_scGeneNet_NN$Theta_mat_list[[g]])
      }
      ##
      for(g in 1:group_num){
        densy_mat[g,l]<- length(which((Theta_mat_array_hat_all[,,g,l])[upper.tri(Theta_mat_array_hat_all[,,g,l])]!=0))/((nrow(Theta_mat_array_hat_all[,,g,l])*(nrow(Theta_mat_array_hat_all[,,g,l])-1))/2)
      }
      ##
      U_mat_array[l,,]<-res_scGeneNet_NN$U_mat
      ## Calculate the soft clustering strength
      U_mat_max<-matrix(0,nrow = nrow(res_scGeneNet_NN$U_mat),ncol = ncol(res_scGeneNet_NN$U_mat))
      for(i in 1:nrow(U_mat_max)){
        U_mat_max[i,which.max(res_scGeneNet_NN$U_mat[i,])]<-1
      }
      soft_cs<-mean(as.vector(abs(res_scGeneNet_NN$U_mat - U_mat_max)))
      soft_cs_vec<-c(soft_cs_vec,soft_cs)
      ##purity
      locator_cur<-apply(res_scGeneNet_NN$U_mat,MARGIN = 1,which.max)
      cluster_predict_true_cur<-as.vector(apply(table(locator_cur,cluster_true),MARGIN = 1,which.max))
      table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
      for(g1 in 1:group_num){
        for(g2 in 1:group_num){
          table_reshape[g1,g2]<-length(which(locator_cur == g1 & cluster_true == cluster_predict_true_cur[g2]))
        }
      }
      ##
      purity_cur<-sum(diag(table_reshape))/sum(table_reshape)
      purity_vec<-c(purity_vec,purity_cur)
      ##
      
      ##
      cluster_acc<-c(cluster_acc,adjustedRandIndex(apply(res_scGeneNet_NN$U_mat,MARGIN = 1,which.max),cluster_true))
      print(paste("ARI: ",round(cluster_acc[length(cluster_acc)],3),sep = ""))
      print(paste("purity: ",round(purity_cur,3),sep = ""))
    }
    ###################################
    ##save
    bic_object_list[["BIC_VMICL"]]<-BIC_mat_shrink
    bic_object_list[["BIC_VICL"]]<-BIC_mat_shrink1
    bic_object_list[["soft_cs_vec"]]<-soft_cs_vec
    bic_object_list[["purity_vec"]]<-purity_vec
    ##
    soft_cs_mean<-mean(soft_cs_vec)
    ##
    bic_object_list[["cluster_acc"]]<-cluster_acc
    ##
    rm(obs_mat,prior_list,res_init0);gc()
    time_use<-sum(time_cost)
    print("model has been run.")
    #############
    ##BIC_shrink_choose
    BIC_vec_shrink<-colSums(BIC_mat_shrink)
    BIC_sum_choose_shrink<-which.min(BIC_vec_shrink)
    BIC_sep_choose_shrink<-apply(BIC_mat_shrink,MARGIN = 1,which.min)
    ##
    BIC_vec_shrink1<-colSums(BIC_mat_shrink1)
    BIC_sum_choose_shrink1<-which.min(BIC_vec_shrink1)
    BIC_sep_choose_shrink1<-apply(BIC_mat_shrink1,MARGIN = 1,which.min)
    ##
    bic_object_list[["BIC_sep_choose_VMICL"]]<-BIC_sep_choose_shrink
    bic_object_list[["BIC_sum_choose_VMICL"]]<-BIC_sum_choose_shrink
    bic_object_list[["BIC_sep_choose_VICL"]]<-BIC_sep_choose_shrink1
    bic_object_list[["BIC_sum_choose_VICL"]]<-BIC_sum_choose_shrink1
    ###
    ##ARI & purity
    ARI_mean<-mean(cluster_acc)
    ARI_max<-max(cluster_acc)
    ARI_min<-min(cluster_acc)
    ARI_BIC_sum_shrink<-cluster_acc[BIC_sum_choose_shrink]
    
    purity_mean<-mean(purity_vec)
    purity_max<-max(purity_vec)
    purity_min<-min(purity_vec)
    purity_BIC_sum_shrink<-purity_vec[BIC_sum_choose_shrink]
    ###
    ####################
    Theta_array_hat_VICL<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      Theta_array_hat_VICL[g,,]<-Theta_mat_array_hat_all[,,g,which.min(BIC_mat_shrink1[g,])]
    }
    Theta_array_hat_VMICL<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      Theta_array_hat_VMICL[g,,]<-Theta_mat_array_hat_all[,,g,which.min(BIC_mat_shrink[g,])]
    }
    ##
    ##density choose
    densy_sep_choose<-c()
    for(g in 1:group_num){
      densy_sep_choose<-c(densy_sep_choose,(which(densy_mat[g,] > densy_degree))[which.min(abs(densy_mat[g,which(densy_mat[g,] > densy_degree)] - densy_degree))])
    }
    densy_sep_choose2<-c()
    for(g in 1:group_num){
      densy_sep_choose2<-c(densy_sep_choose2,(which(densy_mat[g,] > 2 * densy_degree))[which.min(abs(densy_mat[g,which(densy_mat[g,] > 2 * densy_degree)] - 2 * densy_degree))])
    }
    
    ##
    Theta_array_hat_density<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      ##select
      pcor_mat_trans<-Theta_mat_array_hat_all[,,g,densy_sep_choose[g]]/(matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose[g]])),ncol = 1) %*% matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose[g]])),nrow = 1))
      diag(pcor_mat_trans)<-1
      choose_mat<-diag(ncol(pcor_mat_trans))
      choose_mat[upper.tri(choose_mat)][order((abs(pcor_mat_trans))[upper.tri(abs(pcor_mat_trans))],decreasing = TRUE)[1:floor((ncol(pcor_mat_trans) * (ncol(pcor_mat_trans) - 1)/2) * densy_degree)]]<-1
      choose_mat<-choose_mat + t(choose_mat) - diag(ncol(pcor_mat_trans))
      Theta_array_hat_density[g,,]<-Theta_mat_array_hat_all[,,g,densy_sep_choose[g]] * choose_mat
    }
    Theta_array_hat_density2<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      ##select
      pcor_mat_trans<-Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]]/(matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]])),ncol = 1) %*% matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]])),nrow = 1))
      diag(pcor_mat_trans)<-1
      choose_mat<-diag(ncol(pcor_mat_trans))
      choose_mat[upper.tri(choose_mat)][order((abs(pcor_mat_trans))[upper.tri(abs(pcor_mat_trans))],decreasing = TRUE)[1:floor((ncol(pcor_mat_trans) * (ncol(pcor_mat_trans) - 1)/2) * (2 * densy_degree))]]<-1
      choose_mat<-choose_mat + t(choose_mat) - diag(ncol(pcor_mat_trans))
      Theta_array_hat_density2[g,,]<-Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]] * choose_mat
    }
    
    cluster_predict_true_vec<-as.vector(apply(table(apply(U_mat_array[1,,],MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max)) 
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    ##save
    save(Theta_array_hat_VICL,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/Theta_array_hat_VICL_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(Theta_array_hat_VMICL,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/Theta_array_hat_VMICL_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                 "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(Theta_array_hat_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/Theta_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                 "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(Theta_array_hat_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/Theta_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                            "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    ##
    evaluator_1_VICL<-evaluator_1(input_array_hat = Theta_array_hat_VICL,
                                  type = 1,
                Theta_array_true = Theta_mat_array_TRUE,
                density_choose = NULL,
                cluster_predict_true_vec = cluster_predict_true_vec,
                uppertri_choose = NULL)
    evaluator_1_VMICL<-evaluator_1(input_array_hat = Theta_array_hat_VMICL,
                                   type = 1,
                                  Theta_array_true = Theta_mat_array_TRUE,
                                  density_choose = NULL,
                                  cluster_predict_true_vec = cluster_predict_true_vec,
                                  uppertri_choose = NULL)
    evaluator_1_density<-evaluator_1(input_array_hat = Theta_array_hat_density,
                                     type = 1,
                                  Theta_array_true = Theta_mat_array_TRUE,
                                  density_choose = NULL,
                                  cluster_predict_true_vec = cluster_predict_true_vec,
                                  uppertri_choose = NULL)
    evaluator_1_density2<-evaluator_1(input_array_hat = Theta_array_hat_density2,
                                     type = 1,
                                     Theta_array_true = Theta_mat_array_TRUE,
                                     density_choose = NULL,
                                     cluster_predict_true_vec = cluster_predict_true_vec,
                                     uppertri_choose = NULL)
    
    evaluator_VICL_minsample<-evaluator_1_VICL[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_VMICL_minsample<-evaluator_1_VMICL[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density_minsample<-evaluator_1_density[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density2_minsample<-evaluator_1_density2[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    
    evaluator_VICL_minmisc<-evaluator_1_VICL[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_VMICL_minmisc<-evaluator_1_VMICL[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density_minmisc<-evaluator_1_density[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density2_minmisc<-evaluator_1_density2[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_VMICL<-rbind(evaluator_1_VMICL,evaluator_VMICL_minsample,evaluator_VMICL_minmisc)
    evaluator_2_VICL<-rbind(evaluator_1_VICL,evaluator_VICL_minsample,evaluator_VICL_minmisc)
    evaluator_2_density<-rbind(evaluator_1_density,evaluator_density_minsample,evaluator_density_minmisc)
    evaluator_2_density2<-rbind(evaluator_1_density2,evaluator_density2_minsample,evaluator_density2_minmisc)
    
    ##signed
    evaluator_1_VICL_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_VICL,
                                  type = 1,
                                  Theta_array_true = Theta_mat_array_TRUE,
                                  cluster_predict_true_vec = cluster_predict_true_vec,
                                  uppertri_choose = NULL)
    evaluator_1_VMICL_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_VMICL,
                                   type = 1,
                                   Theta_array_true = Theta_mat_array_TRUE,
                                   cluster_predict_true_vec = cluster_predict_true_vec,
                                   uppertri_choose = NULL)
    evaluator_1_density_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_density,
                                     type = 1,
                                     Theta_array_true = Theta_mat_array_TRUE,
                                     cluster_predict_true_vec = cluster_predict_true_vec,
                                     uppertri_choose = NULL)
    evaluator_1_density2_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_density2,
                                      type = 1,
                                      Theta_array_true = Theta_mat_array_TRUE,
                                      cluster_predict_true_vec = cluster_predict_true_vec,
                                      uppertri_choose = NULL)
    
    evaluator_VICL_minsample_signed<-evaluator_1_VICL_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_VMICL_minsample_signed<-evaluator_1_VMICL_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density_minsample_signed<-evaluator_1_density_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density2_minsample_signed<-evaluator_1_density2_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    
    evaluator_VICL_minmisc_signed<-evaluator_1_VICL_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_VMICL_minmisc_signed<-evaluator_1_VMICL_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density_minmisc_signed<-evaluator_1_density_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density2_minmisc_signed<-evaluator_1_density2_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_VMICL_signed<-rbind(evaluator_1_VMICL_signed,evaluator_VMICL_minsample_signed,evaluator_VMICL_minmisc_signed)
    evaluator_2_VICL_signed<-rbind(evaluator_1_VICL_signed,evaluator_VICL_minsample_signed,evaluator_VICL_minmisc_signed)
    evaluator_2_density_signed<-rbind(evaluator_1_density_signed,evaluator_density_minsample_signed,evaluator_density_minmisc_signed)
    evaluator_2_density2_signed<-rbind(evaluator_1_density2_signed,evaluator_density2_minsample_signed,evaluator_density2_minmisc_signed)
    
    

    #######################################################
    ##
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    ##
    res_vec<-c(
      ARI_orgin,
      ARI_mean,
      ARI_max,
      ARI_min,
      ARI_BIC_sum_shrink,
      ##1:5
      
      purity_orgin,
      purity_mean,
      purity_max,
      purity_min,
      purity_BIC_sum_shrink,
      ##6:10
      
      soft_cs_mean,time_use,dr,
      ##11:13
      
      ##non-signed
      ##VMICL part
      as.vector(evaluator_2_VMICL),
      ##14:73
      ##VICL part
      as.vector(evaluator_2_VICL),
      ##74:133
      ##density part
      as.vector(evaluator_2_density),
      ##134:193
      as.vector(evaluator_2_density2),
      ##194:253
      
      ##signed
      ##VMICL part
      as.vector(evaluator_2_VMICL_signed),
      ##254:313
      ##VICL part
      as.vector(evaluator_2_VICL_signed),
      ##314:373
      ##density part
      as.vector(evaluator_2_density_signed),
      ##374:433
      as.vector(evaluator_2_density2_signed),
      ##434:493
      ##
      edge_true_num
      ##494:496
    )
    ##
    return(res_vec)
  }
  
  ##stability
  early_stability_VMICL<-NULL
  early_stability_VICL<-NULL
  early_stability_density<-NULL
  early_stability_density2<-NULL
  early_stability_VMICL_signed<-NULL
  early_stability_VICL_signed<-NULL
  early_stability_density_signed<-NULL
  early_stability_density2_signed<-NULL
  ##
  if(if_stability_eval == TRUE){
    ##VMICL
    #load network and cluster_predict_true_vec
    edge_list<-list()
    edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_scGeneNet[rep,(ncol(repres_scGeneNet)-2):ncol(repres_scGeneNet)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/Theta_array_hat_VMICL_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_VMICL)[1]){
        corr_mat<-(-1) * Theta_array_hat_VMICL[g,,]/(matrix(sqrt(diag(Theta_array_hat_VMICL[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_VMICL[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_VMICL<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_VMICL_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##VICL
    #load network and cluster_predict_true_vec
    edge_list<-list()
    edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_scGeneNet[rep,(ncol(repres_scGeneNet)-2):ncol(repres_scGeneNet)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/Theta_array_hat_VICL_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_VICL)[1]){
        corr_mat<-(-1) * Theta_array_hat_VICL[g,,]/(matrix(sqrt(diag(Theta_array_hat_VICL[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_VICL[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_VICL<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_VICL_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density
    #load network and cluster_predict_true_vec
    edge_list<-list()
    edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_scGeneNet[rep,(ncol(repres_scGeneNet)-2):ncol(repres_scGeneNet)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/Theta_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_density)[1]){
        corr_mat<-(-1) * Theta_array_hat_density[g,,]/(matrix(sqrt(diag(Theta_array_hat_density[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_density[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density2
    #load network and cluster_predict_true_vec
    edge_list<-list()
    edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_scGeneNet[rep,(ncol(repres_scGeneNet)-2):ncol(repres_scGeneNet)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/Theta_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/scGeneNet/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_density2)[1]){
        corr_mat<-(-1) * Theta_array_hat_density2[g,,]/(matrix(sqrt(diag(Theta_array_hat_density2[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_density2[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density2<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density2_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
  }
  
  return(list(repres_scGeneNet = repres_scGeneNet,
              early_stability_VMICL = early_stability_VMICL,
              early_stability_VICL = early_stability_VICL,
              early_stability_density = early_stability_density,
              early_stability_density2 = early_stability_density2,
              early_stability_VMICL_signed = early_stability_VMICL_signed,
              early_stability_VICL_signed = early_stability_VICL_signed,
              early_stability_density_signed = early_stability_density_signed,
              early_stability_density2_signed = early_stability_density2_signed))
}

##Evaluator for VPLN
VPLN_evaluator<-function(penalize_diagonal = FALSE,self_lambda = TRUE,if_stability_eval = FALSE){
  library_size_est = "TSS"
  ##
  repres_VPLN<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS','glasso',
                  'mclust','mltools','data.table','PLNmodels','dplyr',
                  "cluster","pryr","philentropy","umap","Seurat","igraph","Rcpp","scGeneNet","Matrix","stats","GMPR")
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = "")) 
    ##
    bic_object_list<-list()
    ##
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dim_p = dim(obs_mat)[2]
    
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    ls_vec<-res_init$ls
    locator_base = apply(res_init$U_mat, 1, which.max) 
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    
    ###
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    
    diag(table_reshape)<-0
    rm(table_reshape);gc()
    ##################
    
    print("VPLN")
    
    ##################
    Theta_mat_array_hat_all<-array(NA,dim = c(p,p,group_num,100))
    ##
    densy_mat<-matrix(NA,nrow = group_num,ncol = 100)
    BIC_mat<-matrix(NA,nrow = group_num,ncol = 100)
    
    BIC_sep_choose<-c()
    
    original_list<-list(data_1=as.data.frame(obs_mat),
                        Covariate=matrix(0,nrow = nrow(obs_mat),ncol = 1))
    rownames(original_list$Covariate)<-rownames(obs_mat)
    # rm(obs_mat);gc()
    pre_data<-prepare_data(counts = original_list$data_1,
                           covariates = original_list$Covariate,
                           offset = "TSS")
    names(pre_data)[2]<-"covariates"
    ##change library size
    if(library_size_est == "none"){
      pre_data$Offset<-rep(1,length(res_init$ls))
    }
    if(library_size_est == "TSS"){
      pre_data$Offset<-rowSums(obs_mat)/10000
    }
    if(library_size_est == "GMPR"){
      pre_data$Offset<-as.vector(GMPR(obs_mat))
    }
    ##
    rm(original_list);gc()
    
    ##
    time_cost<-c()
    for (g2 in 1:group_num) {
      print(paste("group(",g2,") task of VPLN."))
      if(self_lambda == TRUE){
        ##lambda_max
        fits <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = pre_data,
                           subset = which(locator_base==g2),
                           control_main = list(trace=0,penalize_diagonal = penalize_diagonal,ftol_rel = 1e-4),
                           control_init = list(nPenalties = 1))
        lambda_vec_group<-(seq(from = (fits$penalties[1]), to = (1e-6),length.out = 100))
        rm(fits);gc()
        ##
        time_VPLN_start<-Sys.time()
        fits <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = pre_data,
                           subset = which(locator_base==g2),
                           penalties = lambda_vec_group,
                           control_main = list(trace=0,penalize_diagonal = penalize_diagonal,ftol_rel = 1e-4))
        time_VPLN_end<-Sys.time()
        time_cost<-c(time_cost,as.double(difftime(time_VPLN_end,time_VPLN_start,units = "hours")))
        ##
        bm<-fits$models 
      }else{
        time_VPLN_start<-Sys.time()
        fits <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = pre_data,
                           subset = which(locator_base==g2),
                           control_main = list(trace=0,penalize_diagonal = penalize_diagonal,ftol_rel = 1e-4),
                           control_init = list(nPenalties = 100,min.ratio = 1e-6))
        time_VPLN_end<-Sys.time()
        time_cost<-c(time_cost,as.double(difftime(time_VPLN_end,time_VPLN_start,units = "hours")))
        ##
        bm<-fits$models
      }
      BIC_sep_choose<-c(BIC_sep_choose,which.max(as.vector(fits$criteria[,4])))
      BIC_mat[g2,]<-(-1) * as.vector(fits$criteria[,4])
      rm(fits)
      
      
      
      for (l in 1:100) {
        # print(l)
        graph_VPLN<-bm[[l]]$model_par$Theta
        Theta_mat_array_hat_all[,,g2,l]<-graph_VPLN
        
        densy_mat[g2,l]<-length(which(graph_VPLN[upper.tri(graph_VPLN)]!=0))/(nrow(graph_VPLN) * (nrow(graph_VPLN)-1)/2)
        rm(graph_VPLN)
        
        gc()
        
      }
      rm(bm)
    }
    time_use<-sum(time_cost)
    
    ######################
    Theta_array_hat_BIC<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      Theta_array_hat_BIC[g,,]<-Theta_mat_array_hat_all[,,g,BIC_sep_choose[g]]
    }
    ##
    ##density choose
    densy_sep_choose<-c()
    for(g in 1:group_num){
      densy_sep_choose<-c(densy_sep_choose,(which(densy_mat[g,] > densy_degree))[which.min(abs(densy_mat[g,which(densy_mat[g,] > densy_degree)] - densy_degree))])
    }
    densy_sep_choose2<-c()
    for(g in 1:group_num){
      densy_sep_choose2<-c(densy_sep_choose2,(which(densy_mat[g,] > 2 * densy_degree))[which.min(abs(densy_mat[g,which(densy_mat[g,] > 2 * densy_degree)] - 2 * densy_degree))])
    }
    ##
    Theta_array_hat_density<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      ##select
      pcor_mat_trans<-Theta_mat_array_hat_all[,,g,densy_sep_choose[g]]/(matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose[g]])),ncol = 1) %*% matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose[g]])),nrow = 1))
      diag(pcor_mat_trans)<-1
      choose_mat<-diag(ncol(pcor_mat_trans))
      choose_mat[upper.tri(choose_mat)][order((abs(pcor_mat_trans))[upper.tri(abs(pcor_mat_trans))],decreasing = TRUE)[1:floor((ncol(pcor_mat_trans) * (ncol(pcor_mat_trans) - 1)/2) * densy_degree)]]<-1
      choose_mat<-choose_mat + t(choose_mat) - diag(ncol(pcor_mat_trans))
      Theta_array_hat_density[g,,]<-Theta_mat_array_hat_all[,,g,densy_sep_choose[g]] * choose_mat
      
    }
    Theta_array_hat_density2<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      ##select
      pcor_mat_trans<-Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]]/(matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]])),ncol = 1) %*% matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]])),nrow = 1))
      diag(pcor_mat_trans)<-1
      choose_mat<-diag(ncol(pcor_mat_trans))
      choose_mat[upper.tri(choose_mat)][order((abs(pcor_mat_trans))[upper.tri(abs(pcor_mat_trans))],decreasing = TRUE)[1:floor((ncol(pcor_mat_trans) * (ncol(pcor_mat_trans) - 1)/2) * (2 * densy_degree))]]<-1
      choose_mat<-choose_mat + t(choose_mat) - diag(ncol(pcor_mat_trans))
      Theta_array_hat_density2[g,,]<-Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]] * choose_mat
      
    }
    cluster_predict_true_vec<-as.vector(apply(table(apply(res_init$U_mat,MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max)) 
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    
    ##
    ##save
    save(Theta_array_hat_BIC,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/Theta_array_hat_BIC_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                 "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(Theta_array_hat_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/Theta_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(Theta_array_hat_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/Theta_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    ##non-signed
    evaluator_1_BIC<-evaluator_1(input_array_hat = Theta_array_hat_BIC,
                                 type = 1,
                                  Theta_array_true = Theta_mat_array_TRUE,
                                  cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_density<-evaluator_1(input_array_hat = Theta_array_hat_density,
                                     type = 1,
                                     Theta_array_true = Theta_mat_array_TRUE,
                                     cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_density2<-evaluator_1(input_array_hat = Theta_array_hat_density2,
                                      type = 1,
                                     Theta_array_true = Theta_mat_array_TRUE,
                                     cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_BIC_minsample<-evaluator_1_BIC[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density_minsample<-evaluator_1_density[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density2_minsample<-evaluator_1_density2[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    
    evaluator_BIC_minmisc<-evaluator_1_BIC[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density_minmisc<-evaluator_1_density[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density2_minmisc<-evaluator_1_density2[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_BIC<-rbind(evaluator_1_BIC,evaluator_BIC_minsample,evaluator_BIC_minmisc)
    evaluator_2_density<-rbind(evaluator_1_density,evaluator_density_minsample,evaluator_density_minmisc)
    evaluator_2_density2<-rbind(evaluator_1_density2,evaluator_density2_minsample,evaluator_density2_minmisc)
    
    
    ##signed
    evaluator_1_BIC_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_BIC,
                                 type = 1,
                                 Theta_array_true = Theta_mat_array_TRUE,
                                 cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_density_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_density,
                                     type = 1,
                                     Theta_array_true = Theta_mat_array_TRUE,
                                     cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_density2_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_density2,
                                      type = 1,
                                      Theta_array_true = Theta_mat_array_TRUE,
                                      cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_BIC_minsample_signed<-evaluator_1_BIC_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density_minsample_signed<-evaluator_1_density_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density2_minsample_signed<-evaluator_1_density2_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    
    evaluator_BIC_minmisc_signed<-evaluator_1_BIC_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density_minmisc_signed<-evaluator_1_density_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density2_minmisc_signed<-evaluator_1_density2_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_BIC_signed<-rbind(evaluator_1_BIC_signed,evaluator_BIC_minsample_signed,evaluator_BIC_minmisc_signed)
    evaluator_2_density_signed<-rbind(evaluator_1_density_signed,evaluator_density_minsample_signed,evaluator_density_minmisc_signed)
    evaluator_2_density2_signed<-rbind(evaluator_1_density2_signed,evaluator_density2_minsample_signed,evaluator_density2_minmisc_signed)
    
    ##
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    ##
    res_vec<-c(
      time_use,
      ##1
      
      ##non-signed
      ##BIC part
      as.vector(evaluator_2_BIC),
      ##2:61
      
      ##density part
      as.vector(evaluator_2_density),
      ##62:121
      
      as.vector(evaluator_2_density2),
      ##122:181
      
      ##signed
      ##BIC part
      as.vector(evaluator_2_BIC_signed),
      ##182:241
      
      ##density part
      as.vector(evaluator_2_density_signed),
      ##242:301
      
      as.vector(evaluator_2_density2_signed),
      ##302:361
      
      edge_true_num
      ##362:364
    )
    
   
    #####
    return(res_vec)
    
    
  }
  
  ##stability
  early_stability_BIC<-NULL
  early_stability_density<-NULL
  early_stability_density2<-NULL
  early_stability_BIC_signed<-NULL
  early_stability_density_signed<-NULL
  early_stability_density2_signed<-NULL
  
  if(if_stability_eval == TRUE){
    ##BIC
    #load network and cluster_predict_true_vec
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_VPLN[rep,(ncol(repres_VPLN)-2):ncol(repres_VPLN)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/Theta_array_hat_BIC_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_BIC)[1]){
        corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_BIC<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_BIC_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density
    #load network and cluster_predict_true_vec
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_VPLN[rep,(ncol(repres_VPLN)-2):ncol(repres_VPLN)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/Theta_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_density)[1]){
        corr_mat<-(-1) * Theta_array_hat_density[g,,]/(matrix(sqrt(diag(Theta_array_hat_density[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_density[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density2
    #load network and cluster_predict_true_vec
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_VPLN[rep,(ncol(repres_VPLN)-2):ncol(repres_VPLN)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/Theta_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/VPLN/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_density2)[1]){
        corr_mat<-(-1) * Theta_array_hat_density2[g,,]/(matrix(sqrt(diag(Theta_array_hat_density2[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_density2[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density2<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density2_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
  }
  
  return(list(repres_VPLN = repres_VPLN,
              early_stability_BIC = early_stability_BIC,
              early_stability_density = early_stability_density,
              early_stability_density2 = early_stability_density2,
              early_stability_BIC_signed = early_stability_BIC_signed,
              early_stability_density_signed = early_stability_density_signed,
              early_stability_density2_signed = early_stability_density2_signed))
  
}

##Evaluator for Glasso
Glasso_evaluator<-function(penalize_diagonal = FALSE,if_stability_eval = FALSE){
  library_size_est = "TSS"
  ##
  repres_Glasso<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS','glasso',
                  'mclust','mltools','data.table','PLNmodels','dplyr',
                  "cluster","pryr","philentropy","umap","Seurat","igraph","Rcpp","scGeneNet","Matrix","stats","GMPR")
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    ##
    bic_object_list<-list()
    ##
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    ls_vec<-res_init$ls
    ##
    if(library_size_est == "none"){
      ls_vec<-rep(1,length(res_init$ls))
    }
    if(library_size_est == "TSS"){
      ls_vec<-rowSums(obs_mat)/10000
    }
    if(library_size_est == "GMPR"){
      ls_vec<-as.vector(GMPR(obs_mat))
    }
    ##
    locator_base = apply(res_init$U_mat, 1, which.max) 
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    
    ###
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    ##
    diag(table_reshape)<-0
    rm(table_reshape);gc()
    ##################
    
    print("Glasso")
    
    ##################
    Theta_mat_array_hat_all<-array(NA,dim = c(p,p,group_num,100))
    ##
    densy_mat<-matrix(NA,nrow = group_num,ncol = 100)
    BIC_mat<-matrix(NA,nrow = group_num,ncol = 100)
    BIC_sep_choose<-c()
    
    ##
    glasso_data = log((1+obs_mat)/ls_vec)
    rm(obs_mat);gc()
    lambda_vec<-exp(seq(from=log((1e-6)),to=log((50)),length.out = 100))
    
    
    ##
    time_cost<-c()
    for (g2 in 1:group_num) {
      print(paste("group(",g2,") task of Glasso."))
      time_Glasso_start<-Sys.time()
      graph_lasso_array<-array(0,dim = c(length(lambda_vec),ncol(glasso_data),ncol(glasso_data)))
      cov_X = cov(glasso_data[locator_base==g2,])
      dim_p = dim(glasso_data)[2]
      cov_X = cov_X + 1e-3 * diag(dim_p)
      for (l in 1:length(lambda_vec)) {
        # print(l)
        # graph_glasso<-graph_glasso_list[,,l]
        graph_glasso<-glasso(s=cov_X,rho = lambda_vec[l],penalize.diagonal = penalize_diagonal,thr = 1e-3)$wi
        ##
        graph_lasso_array[l,,]<-graph_glasso
        Theta_mat_array_hat_all[,,g2,l]<-graph_glasso
        densy_mat[g2,l]<-length(which(graph_glasso[upper.tri(graph_glasso)]!=0))/(nrow(graph_glasso) * (nrow(graph_glasso)-1)/2)
        rm(graph_glasso);gc()
        
        gc()
      }
      rm(cov_X);gc()
      BIC_glasso_group<-BIC_glasso(glasso_data[locator_base==g2,],graph_lasso_array)
      ##
      BIC_mat[g2,]<-BIC_glasso_group
      ##
      rm(graph_lasso_array);gc()
      BIC_sep_choose<-c(BIC_sep_choose,which.min(BIC_glasso_group))
      
      time_Glasso_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_Glasso_end,time_Glasso_start,units = "hours")))
    }
    
    time_use<-sum(time_cost)
    
    ######################
    Theta_array_hat_BIC<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      Theta_array_hat_BIC[g,,]<-Theta_mat_array_hat_all[,,g,BIC_sep_choose[g]]
    }
    ##
    ##density choose
    densy_sep_choose<-c()
    for(g in 1:group_num){
      # densy_sep_choose<-c(densy_sep_choose,which.min(abs(densy_mat[g,] - densy_degree)))
      densy_sep_choose<-c(densy_sep_choose,(which(densy_mat[g,] > densy_degree))[which.min(abs(densy_mat[g,which(densy_mat[g,] > densy_degree)] - densy_degree))])
    }
    densy_sep_choose2<-c()
    for(g in 1:group_num){
      # densy_sep_choose2<-c(densy_sep_choose2,which.min(abs(densy_mat[g,] - 2* densy_degree)))
      densy_sep_choose2<-c(densy_sep_choose2,(which(densy_mat[g,] > 2 * densy_degree))[which.min(abs(densy_mat[g,which(densy_mat[g,] > 2 * densy_degree)] - 2 * densy_degree))])
    }
    ##
    Theta_array_hat_density<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      # Theta_array_hat_density[g,,]<-Theta_mat_array_hat_all[,,g,densy_sep_choose[g]]
      ##select
      pcor_mat_trans<-Theta_mat_array_hat_all[,,g,densy_sep_choose[g]]/(matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose[g]])),ncol = 1) %*% matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose[g]])),nrow = 1))
      diag(pcor_mat_trans)<-1
      choose_mat<-diag(ncol(pcor_mat_trans))
      choose_mat[upper.tri(choose_mat)][order((abs(pcor_mat_trans))[upper.tri(abs(pcor_mat_trans))],decreasing = TRUE)[1:floor((ncol(pcor_mat_trans) * (ncol(pcor_mat_trans) - 1)/2) * densy_degree)]]<-1
      choose_mat<-choose_mat + t(choose_mat) - diag(ncol(pcor_mat_trans))
      Theta_array_hat_density[g,,]<-Theta_mat_array_hat_all[,,g,densy_sep_choose[g]] * choose_mat
      
    }
    Theta_array_hat_density2<-array(NA,dim = c(group_num,p,p))
    for(g in 1:group_num){
      # Theta_array_hat_density2[g,,]<-Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]]
      ##select
      pcor_mat_trans<-Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]]/(matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]])),ncol = 1) %*% matrix(sqrt(diag(Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]])),nrow = 1))
      diag(pcor_mat_trans)<-1
      choose_mat<-diag(ncol(pcor_mat_trans))
      choose_mat[upper.tri(choose_mat)][order((abs(pcor_mat_trans))[upper.tri(abs(pcor_mat_trans))],decreasing = TRUE)[1:floor((ncol(pcor_mat_trans) * (ncol(pcor_mat_trans) - 1)/2) * (2 * densy_degree))]]<-1
      choose_mat<-choose_mat + t(choose_mat) - diag(ncol(pcor_mat_trans))
      Theta_array_hat_density2[g,,]<-Theta_mat_array_hat_all[,,g,densy_sep_choose2[g]] * choose_mat
      
    }
    cluster_predict_true_vec<-as.vector(apply(table(apply(res_init$U_mat,MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max)) 
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    
    ##
    save(Theta_array_hat_BIC,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/Theta_array_hat_BIC_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                 "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(Theta_array_hat_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/Theta_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(Theta_array_hat_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/Theta_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    ##
    evaluator_1_BIC<-evaluator_1(input_array_hat = Theta_array_hat_BIC,
                                 type = 1,
                                 Theta_array_true = Theta_mat_array_TRUE,
                                 cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_density<-evaluator_1(input_array_hat = Theta_array_hat_density,
                                     type = 1,
                                     Theta_array_true = Theta_mat_array_TRUE,
                                     cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_density2<-evaluator_1(input_array_hat = Theta_array_hat_density2,
                                     type = 1,
                                     Theta_array_true = Theta_mat_array_TRUE,
                                     cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_BIC_minsample<-evaluator_1_BIC[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density_minsample<-evaluator_1_density[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density2_minsample<-evaluator_1_density2[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    
    evaluator_BIC_minmisc<-evaluator_1_BIC[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density_minmisc<-evaluator_1_density[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density2_minmisc<-evaluator_1_density2[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_BIC<-rbind(evaluator_1_BIC,evaluator_BIC_minsample,evaluator_BIC_minmisc)
    evaluator_2_density<-rbind(evaluator_1_density,evaluator_density_minsample,evaluator_density_minmisc)
    evaluator_2_density2<-rbind(evaluator_1_density2,evaluator_density2_minsample,evaluator_density2_minmisc)
    
    
    ##signed
    evaluator_1_BIC_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_BIC,
                                 type = 1,
                                 Theta_array_true = Theta_mat_array_TRUE,
                                 cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_density_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_density,
                                     type = 1,
                                     Theta_array_true = Theta_mat_array_TRUE,
                                     cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_density2_signed<-evaluator_1_signed(input_array_hat = Theta_array_hat_density2,
                                      type = 1,
                                      Theta_array_true = Theta_mat_array_TRUE,
                                      cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_BIC_minsample_signed<-evaluator_1_BIC_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density_minsample_signed<-evaluator_1_density_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_density2_minsample_signed<-evaluator_1_density2_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    
    evaluator_BIC_minmisc_signed<-evaluator_1_BIC_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density_minmisc_signed<-evaluator_1_density_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_density2_minmisc_signed<-evaluator_1_density2_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_BIC_signed<-rbind(evaluator_1_BIC_signed,evaluator_BIC_minsample_signed,evaluator_BIC_minmisc_signed)
    evaluator_2_density_signed<-rbind(evaluator_1_density_signed,evaluator_density_minsample_signed,evaluator_density_minmisc_signed)
    evaluator_2_density2_signed<-rbind(evaluator_1_density2_signed,evaluator_density2_minsample_signed,evaluator_density2_minmisc_signed)
    #############################
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    ##
    res_vec<-c(
      time_use,
      ##1
      
      ##BIC part
      as.vector(evaluator_2_BIC),
      
      ##density part
      as.vector(evaluator_2_density),
      
      as.vector(evaluator_2_density2),
      
      ##sigend
      ##BIC part
      as.vector(evaluator_2_BIC_signed),
      
      ##density part
      as.vector(evaluator_2_density_signed),
      
      as.vector(evaluator_2_density2_signed),
      
      
      edge_true_num
    )
    #####
    return(res_vec)
    
    
  }
  
  ##stability
  early_stability_BIC<-NULL
  early_stability_density<-NULL
  early_stability_density2<-NULL
  early_stability_BIC_signed<-NULL
  early_stability_density_signed<-NULL
  early_stability_density2_signed<-NULL
  ##
  if(if_stability_eval == TRUE){
    ##BIC
    #load network and cluster_predict_true_vec
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_Glasso[rep,(ncol(repres_Glasso)-2):ncol(repres_Glasso)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/Theta_array_hat_BIC_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_BIC)[1]){
        corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_BIC<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_BIC_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density
    #load network and cluster_predict_true_vec
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_Glasso[rep,(ncol(repres_Glasso)-2):ncol(repres_Glasso)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/Theta_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_density)[1]){
        corr_mat<-(-1) * Theta_array_hat_density[g,,]/(matrix(sqrt(diag(Theta_array_hat_density[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_density[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density2
    #load network and cluster_predict_true_vec
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_Glasso[rep,(ncol(repres_Glasso)-2):ncol(repres_Glasso)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/Theta_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/Glasso/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(Theta_array_hat_density2)[1]){
        corr_mat<-(-1) * Theta_array_hat_density2[g,,]/(matrix(sqrt(diag(Theta_array_hat_density2[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_density2[g,,])),nrow = 1))
        diag(corr_mat)<-1
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density2<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density2_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
  }
  
  return(list(repres_Glasso = repres_Glasso,
              early_stability_BIC = early_stability_BIC,
              early_stability_density = early_stability_density,
              early_stability_density2 = early_stability_density2,
              early_stability_BIC_signed = early_stability_BIC_signed,
              early_stability_density_signed = early_stability_density_signed,
              early_stability_density2_signed = early_stability_density2_signed))
}

##Evaluator for ppcor
ppcor_evaluator<-function(penalize_diagonal = FALSE,method = "pearson",if_stability_eval = FALSE){
  library_size_est = "TSS"
  ##
  repres_ppcor<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS','glasso',
                  'mclust','mltools','data.table','PLNmodels','dplyr',
                  "cluster","pryr","philentropy","umap","Seurat","igraph","Rcpp","scGeneNet","Matrix","stats","GMPR","ppcor")
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    ##
    bic_object_list<-list()
    ##
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    ls_vec<-res_init$ls
    ##
    if(library_size_est == "none"){
      ls_vec<-rep(1,length(res_init$ls))
    }
    if(library_size_est == "TSS"){
      ls_vec<-rowSums(obs_mat)/10000
    }
    if(library_size_est == "GMPR"){
      ls_vec<-as.vector(GMPR(obs_mat))
    }
    ##
    locator_base = apply(res_init$U_mat, 1, which.max) 
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    
    ###
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    ##
    diag(table_reshape)<-0
    rm(table_reshape);gc()
    ##################
    
    print("ppcor")
    pvalue_threshold<-0.05
    ##################
    ##
    ppcor_data = log((1+obs_mat)/ls_vec)
    rm(obs_mat);gc()
    
    
    ##
    time_cost<-c()
    pcor_mat_array_hat<-array(NA,dim = c(group_num,p,p))
    for (g2 in 1:group_num) {
      print(paste("group(",g2,") task of ppcor."))
      time_ppcor_start<-Sys.time()
      pcor_res<-pcor(ppcor_data[locator_base==g2,],method = method)
      ##
      pcor_res_use<-pcor_res$estimate
      pcor_mat_array_hat[g2,,]<-pcor_res_use
      ##
      time_ppcor_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_ppcor_end,time_ppcor_start,units = "hours")))
    }
    time_use<-sum(time_cost)
    pcor_mat_array_hat_density<-pcor_mat_array_hat
    pcor_mat_array_hat_density2<-pcor_mat_array_hat
    for(g in 1:group_num){
      mat_use<-pcor_mat_array_hat_density[g,,]
      mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
      mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
      pcor_mat_array_hat_density[g,,]<-mat_use
    }
    for(g in 1:group_num){
      mat_use<-pcor_mat_array_hat_density2[g,,]
      mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
      mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
      pcor_mat_array_hat_density2[g,,]<-mat_use
    }
    cluster_predict_true_vec<-as.vector(apply(table(apply(res_init$U_mat,MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max)) 
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    
    ##
    save(pcor_mat_array_hat,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/pcor_mat_array_hat_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
    save(pcor_mat_array_hat_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/pcor_mat_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                               "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
    save(pcor_mat_array_hat_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/pcor_mat_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                               "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
    ##
    evaluator_1_pcor<-evaluator_1(input_array_hat = pcor_mat_array_hat,
                                  type = 2,
                                 Theta_array_true = Theta_mat_array_TRUE,
                                 cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_pcor_density<-evaluator_1(input_array_hat = pcor_mat_array_hat_density,
                                  type = 2,
                                  Theta_array_true = Theta_mat_array_TRUE,
                                  cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_pcor_density2<-evaluator_1(input_array_hat = pcor_mat_array_hat_density2,
                                          type = 2,
                                          Theta_array_true = Theta_mat_array_TRUE,
                                          cluster_predict_true_vec = cluster_predict_true_vec)
    
    evaluator_pcor_minsample<-evaluator_1_pcor[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_pcor_minsample_density<-evaluator_1_pcor_density[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_pcor_minsample_density2<-evaluator_1_pcor_density2[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_pcor_minmisc<-evaluator_1_pcor[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_pcor_minmisc_density<-evaluator_1_pcor_density[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_pcor_minmisc_density2<-evaluator_1_pcor_density2[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_pcor<-rbind(evaluator_1_pcor,evaluator_pcor_minsample,evaluator_pcor_minmisc)
    evaluator_2_pcor_density<-rbind(evaluator_1_pcor_density,evaluator_pcor_minsample_density,evaluator_pcor_minmisc_density)
    evaluator_2_pcor_density2<-rbind(evaluator_1_pcor_density2,evaluator_pcor_minsample_density2,evaluator_pcor_minmisc_density2)
    
    ##signed
    evaluator_1_pcor_signed<-evaluator_1_signed(input_array_hat = pcor_mat_array_hat,
                                  type = 2,
                                  Theta_array_true = Theta_mat_array_TRUE,
                                  cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_pcor_density_signed<-evaluator_1_signed(input_array_hat = pcor_mat_array_hat_density,
                                          type = 2,
                                          Theta_array_true = Theta_mat_array_TRUE,
                                          cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_pcor_density2_signed<-evaluator_1_signed(input_array_hat = pcor_mat_array_hat_density2,
                                           type = 2,
                                           Theta_array_true = Theta_mat_array_TRUE,
                                           cluster_predict_true_vec = cluster_predict_true_vec)
    
    evaluator_pcor_minsample_signed<-evaluator_1_pcor_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_pcor_minsample_density_signed<-evaluator_1_pcor_density_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_pcor_minsample_density2_signed<-evaluator_1_pcor_density2_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    
    evaluator_pcor_minmisc_signed<-evaluator_1_pcor_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_pcor_minmisc_density_signed<-evaluator_1_pcor_density_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_pcor_minmisc_density2_signed<-evaluator_1_pcor_density2_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_pcor_signed<-rbind(evaluator_1_pcor_signed,evaluator_pcor_minsample_signed,evaluator_pcor_minmisc_signed)
    evaluator_2_pcor_density_signed<-rbind(evaluator_1_pcor_density_signed,evaluator_pcor_minsample_density_signed,evaluator_pcor_minmisc_density_signed)
    evaluator_2_pcor_density2_signed<-rbind(evaluator_1_pcor_density2_signed,evaluator_pcor_minsample_density2_signed,evaluator_pcor_minmisc_density2_signed)
    
    ##
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    #############################
    res_vec<-c(
      time_use,
      ##1
      
      as.vector(evaluator_2_pcor),
      
      as.vector(evaluator_2_pcor_density),
      
      as.vector(evaluator_2_pcor_density2),
      
      ##signed
      as.vector(evaluator_2_pcor_signed),
      
      as.vector(evaluator_2_pcor_density_signed),
      
      as.vector(evaluator_2_pcor_density2_signed),
      
      edge_true_num

    )
    #####
    return(res_vec)
    
    
  }
  ##stability
  early_stability_BIC<-NULL
  early_stability_density<-NULL
  early_stability_density2<-NULL
  early_stability_BIC_signed<-NULL
  early_stability_density_signed<-NULL
  early_stability_density2_signed<-NULL
  ##
  if(if_stability_eval == TRUE){
    #load network and cluster_predict_true_vec
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_ppcor[rep,(ncol(repres_ppcor)-2):ncol(repres_ppcor)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/pcor_mat_array_hat_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(pcor_mat_array_hat)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        corr_mat<-pcor_mat_array_hat[g,,]
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_BIC<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_BIC_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##denisty1
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_ppcor[rep,(ncol(repres_ppcor)-2):ncol(repres_ppcor)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/pcor_mat_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(pcor_mat_array_hat_density)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        corr_mat<-pcor_mat_array_hat_density[g,,]
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##denisty2
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_ppcor[rep,(ncol(repres_ppcor)-2):ncol(repres_ppcor)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/pcor_mat_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/ppcor/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_method",method,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(pcor_mat_array_hat_density2)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        corr_mat<-pcor_mat_array_hat_density2[g,,]
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density2<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density2_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
  }
  
  return(list(repres_ppcor = repres_ppcor,
              early_stability_BIC = early_stability_BIC,
              early_stability_density = early_stability_density,
              early_stability_density2 = early_stability_density2,
              early_stability_BIC_signed = early_stability_BIC_signed,
              early_stability_density_signed = early_stability_density_signed,
              early_stability_density2_signed = early_stability_density2_signed))
  
}

##Evaluator for GENIE3 (count as input)
GENIE3_evaluator<-function(penalize_diagonal = FALSE,if_stability_eval = FALSE){
    ############
  library_size_est = "TSS"
  ##
    repres_GENIE3<-matrix(NA,nrow = 0,ncol = 184)
    for(rep in 1:reptime){
      ## load data-------------
      load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                 "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
      load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                 "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
      ##
      bic_object_list<-list()
      ##
      obs_mat = data[["obs_mat"]]
      colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
      rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
      Theta_mat_array_TRUE =data[["group_specific_precision"]]
      cluster_true = data[["locator_vec"]]
      # hub_gene<-data[["hub_gene"]]
      if(length(data[["hub_gene"]])>0){
        hub_gene<-data[["hub_gene"]]
      }else{
        hub_gene<-NULL
      }
      
      n<-dim(obs_mat)[1]
      p<-dim(obs_mat)[2]
      dr = length(which( obs_mat == 0 )) / (n*p)
      
      ls_vec<-res_init$ls
      ##
      if(library_size_est == "none"){
        ls_vec<-rep(1,length(res_init$ls))
      }
      if(library_size_est == "TSS"){
        ls_vec<-rowSums(obs_mat)/10000
      }
      if(library_size_est == "GMPR"){
        ls_vec<-as.vector(GMPR(obs_mat))
      }
      ##
      locator_base = apply(res_init$U_mat, 1, which.max) 
      ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
      cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
      
      ###
      table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
      for(g1 in 1:group_num){
        for(g2 in 1:group_num){
          table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
        }
      }
      bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
      bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
      ##
      diag(table_reshape)<-0
      rm(table_reshape);gc()
      ##################
      
      print("GENIE3")
      set.seed(123)
      # # regulators <- c(2, 4, 7)
      # weightMat<-GENIE3(as.matrix(obs_mat),
      #                   nCores = 1)
      
      
      ##
      time_cost<-c()
      gene_name_vec<-colnames(obs_mat)
      GENIE3_mat_array_hat<-array(NA,dim = c(group_num,p,p))
      for (g2 in 1:group_num) {
        print(paste("group(",g2,") task of GENIE3."))
        time_GENIE3_start<-Sys.time()
        ##
        # regulators <- c(2, 4, 7)
        weightMat0<-GENIE3(t(as.matrix(obs_mat[locator_base==g2,])),
                          nCores = num_core)
        ##reorder
        weightMat<-matrix(NA,nrow = nrow(weightMat0),ncol = ncol(weightMat0))
        for(i in 1:nrow(weightMat0)){
          for(j in 1:ncol(weightMat0)){
            weightMat[i,j]<-weightMat0[which(rownames(weightMat0) == gene_name_vec[i]),which(colnames(weightMat0) == gene_name_vec[j])]
          }
        }
        weightMat<-ifelse(is.nan(weightMat),0,weightMat)
        GENIE3_mat_array_hat[g2,,]<-weightMat
        ##
        time_GENIE3_end<-Sys.time()
        time_cost<-c(time_cost,as.double(difftime(time_GENIE3_end,time_GENIE3_start,units = "hours")))
      }
      time_use<-sum(time_cost)
      
      GENIE3_mat_array_hat_density<-GENIE3_mat_array_hat
      GENIE3_mat_array_hat_density2<-GENIE3_mat_array_hat
      for(g in 1:group_num){
        mat_use<-GENIE3_mat_array_hat[g,,]
        mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
        mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
        GENIE3_mat_array_hat_density[g,,]<-mat_use
      }
      for(g in 1:group_num){
        mat_use<-GENIE3_mat_array_hat_density2[g,,]
        mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
        mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
        GENIE3_mat_array_hat_density2[g,,]<-mat_use
      }
      
      cluster_predict_true_vec<-as.vector(apply(table(apply(res_init$U_mat,MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max)) 
      bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
      
      ##
      save(GENIE3_mat_array_hat,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3_mat_array_hat_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                 "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      save(GENIE3_mat_array_hat_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3_mat_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                   "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      save(GENIE3_mat_array_hat_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3_mat_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                   "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      evaluator_1_GENIE3<-evaluator_1(input_array_hat = GENIE3_mat_array_hat,
                                      type = 3,
                                      Theta_array_true = Theta_mat_array_TRUE,
                                      cluster_predict_true_vec = cluster_predict_true_vec)
      evaluator_1_GENIE3_density<-evaluator_1(input_array_hat = GENIE3_mat_array_hat_density,
                                      type = 3,
                                      Theta_array_true = Theta_mat_array_TRUE,
                                      cluster_predict_true_vec = cluster_predict_true_vec)
      evaluator_1_GENIE3_density2<-evaluator_1(input_array_hat = GENIE3_mat_array_hat_density2,
                                      type = 3,
                                      Theta_array_true = Theta_mat_array_TRUE,
                                      cluster_predict_true_vec = cluster_predict_true_vec)
      
      evaluator_GENIE3_minsample<-evaluator_1_GENIE3[which.min(bic_object_list[["cluster_predicted_sample"]]),]
      evaluator_GENIE3_minsample_density<-evaluator_1_GENIE3_density[which.min(bic_object_list[["cluster_predicted_sample"]]),]
      evaluator_GENIE3_minsample_density2<-evaluator_1_GENIE3_density2[which.min(bic_object_list[["cluster_predicted_sample"]]),]
      evaluator_GENIE3_minmisc<-evaluator_1_GENIE3[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
      evaluator_GENIE3_minmisc_density<-evaluator_1_GENIE3_density[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
      evaluator_GENIE3_minmisc_density2<-evaluator_1_GENIE3_density2[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
      
      evaluator_2_GENIE3<-rbind(evaluator_1_GENIE3,evaluator_GENIE3_minsample,evaluator_GENIE3_minmisc)
      evaluator_2_GENIE3_density<-rbind(evaluator_1_GENIE3_density,evaluator_GENIE3_minsample_density,evaluator_GENIE3_minmisc_density)
      evaluator_2_GENIE3_density2<-rbind(evaluator_1_GENIE3_density2,evaluator_GENIE3_minsample_density2,evaluator_GENIE3_minmisc_density2)
      
      edge_true_num<-c()
      for(g in 1:group_num){
        edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
      }
      #############################
      res_vec<-c(
        time_use,
        ##1
        
        as.vector(evaluator_2_GENIE3),
        
        as.vector(evaluator_2_GENIE3_density),
        
        as.vector(evaluator_2_GENIE3_density2),
        
        edge_true_num
        
      )
      ##
      repres_GENIE3<-rbind(repres_GENIE3,res_vec)
    }
    
    ##stability
    early_stability_BIC<-NULL
    early_stability_density<-NULL
    early_stability_density2<-NULL
    if(if_stability_eval == TRUE){
      #load network and cluster_predict_true_vec
      edge_list<-list()
      cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
      for(rep in 1:reptime){
        edge_true_num<-repres_GENIE3[rep,(ncol(repres_GENIE3)-2):ncol(repres_GENIE3)]
        load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3_mat_array_hat_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                         "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
        load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                         "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
        ##
        cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
        ##
        edge_list_rep<-list()
        for(g in 1:dim(GENIE3_mat_array_hat)[1]){
          # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
          # diag(corr_mat)<-1
          weight_mat_use<-ifelse(GENIE3_mat_array_hat[g,,]>t(GENIE3_mat_array_hat[g,,]),GENIE3_mat_array_hat[g,,],t(GENIE3_mat_array_hat[g,,]))
          corr_mat<-weight_mat_use
          edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
          edge_list_rep[[g]]<-edge_index_vec
        }
        edge_list[[rep]]<-edge_list_rep
        
      }
      ##
      early_stability_BIC<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
      
      ##density
      edge_list<-list()
      cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
      for(rep in 1:reptime){
        edge_true_num<-repres_GENIE3[rep,(ncol(repres_GENIE3)-2):ncol(repres_GENIE3)]
        load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3_mat_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                         "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
        load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                         "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
        ##
        cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
        ##
        edge_list_rep<-list()
        for(g in 1:dim(GENIE3_mat_array_hat_density)[1]){
          # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
          # diag(corr_mat)<-1
          weight_mat_use<-ifelse(GENIE3_mat_array_hat_density[g,,]>t(GENIE3_mat_array_hat_density[g,,]),GENIE3_mat_array_hat_density[g,,],t(GENIE3_mat_array_hat_density[g,,]))
          corr_mat<-weight_mat_use
          edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
          edge_list_rep[[g]]<-edge_index_vec
        }
        edge_list[[rep]]<-edge_list_rep
        
      }
      ##
      early_stability_density<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
      
      ##density2
      edge_list<-list()
      cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
      for(rep in 1:reptime){
        edge_true_num<-repres_GENIE3[rep,(ncol(repres_GENIE3)-2):ncol(repres_GENIE3)]
        load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3_mat_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                         "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
        load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                         "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
        ##
        cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
        ##
        edge_list_rep<-list()
        for(g in 1:dim(GENIE3_mat_array_hat_density2)[1]){
          # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
          # diag(corr_mat)<-1
          weight_mat_use<-ifelse(GENIE3_mat_array_hat_density2[g,,]>t(GENIE3_mat_array_hat_density2[g,,]),GENIE3_mat_array_hat_density2[g,,],t(GENIE3_mat_array_hat_density2[g,,]))
          corr_mat<-weight_mat_use
          edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
          edge_list_rep[[g]]<-edge_index_vec
        }
        edge_list[[rep]]<-edge_list_rep
        
      }
      ##
      early_stability_density2<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
      
    }
    
    return(list(repres_GENIE3 = repres_GENIE3,
                early_stability_BIC = early_stability_BIC,
                early_stability_density = early_stability_density,
                early_stability_density2 = early_stability_density2))
}

##Evaluator for GENIE3 (normalized matrix as input)
GENIE3_nor_evaluator<-function(penalize_diagonal = FALSE,if_stability_eval = FALSE){
  library_size_est = "TSS"
  ##
  repres_GENIE3<-matrix(NA,nrow = 0,ncol = 184)
  for(rep in 1:reptime){
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    ##
    bic_object_list<-list()
    ##
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    ls_vec<-res_init$ls
    ##
    if(library_size_est == "none"){
      ls_vec<-rep(1,length(res_init$ls))
    }
    if(library_size_est == "TSS"){
      ls_vec<-rowSums(obs_mat)/10000
    }
    if(library_size_est == "GMPR"){
      ls_vec<-as.vector(GMPR(obs_mat))
    }
    ##
    locator_base = apply(res_init$U_mat, 1, which.max) 
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    
    ###
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    ##
    diag(table_reshape)<-0
    rm(table_reshape);gc()
    ##################
    
    print("GENIE3")
    set.seed(123)
    # # regulators <- c(2, 4, 7)
    # weightMat<-GENIE3(as.matrix(obs_mat),
    #                   nCores = 1)
    
    
    ##
    time_cost<-c()
    gene_name_vec<-colnames(obs_mat)
    GENIE3_mat_array_hat<-array(NA,dim = c(group_num,p,p))
    for (g2 in 1:group_num) {
      print(paste("group(",g2,") task of GENIE3."))
      time_GENIE3_start<-Sys.time()
      ##
      # regulators <- c(2, 4, 7)
      weightMat0<-GENIE3(t(as.matrix(obs_mat[locator_base==g2,] / ls_vec[which(locator_base==g2)])),
                         nCores = num_core)
      ##reorder
      weightMat<-matrix(NA,nrow = nrow(weightMat0),ncol = ncol(weightMat0))
      for(i in 1:nrow(weightMat0)){
        for(j in 1:ncol(weightMat0)){
          weightMat[i,j]<-weightMat0[which(rownames(weightMat0) == gene_name_vec[i]),which(colnames(weightMat0) == gene_name_vec[j])]
        }
      }
      weightMat<-ifelse(is.nan(weightMat),0,weightMat)
      GENIE3_mat_array_hat[g2,,]<-weightMat
      ##
      time_GENIE3_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_GENIE3_end,time_GENIE3_start,units = "hours")))
    }
    time_use<-sum(time_cost)
    
    GENIE3_mat_array_hat_density<-GENIE3_mat_array_hat
    GENIE3_mat_array_hat_density2<-GENIE3_mat_array_hat
    for(g in 1:group_num){
      mat_use<-GENIE3_mat_array_hat[g,,]
      mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
      mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
      GENIE3_mat_array_hat_density[g,,]<-mat_use
    }
    for(g in 1:group_num){
      mat_use<-GENIE3_mat_array_hat_density2[g,,]
      mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
      mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
      GENIE3_mat_array_hat_density2[g,,]<-mat_use
    }
    
    cluster_predict_true_vec<-as.vector(apply(table(apply(res_init$U_mat,MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max)) 
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    
    ##
    save(GENIE3_mat_array_hat,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3nor_mat_array_hat_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                 "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(GENIE3_mat_array_hat_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3nor_mat_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                         "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(GENIE3_mat_array_hat_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3nor_mat_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                          "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    ##
    evaluator_1_GENIE3<-evaluator_1(input_array_hat = GENIE3_mat_array_hat,
                                    type = 3,
                                    Theta_array_true = Theta_mat_array_TRUE,
                                    cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_GENIE3_density<-evaluator_1(input_array_hat = GENIE3_mat_array_hat_density,
                                            type = 3,
                                            Theta_array_true = Theta_mat_array_TRUE,
                                            cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_GENIE3_density2<-evaluator_1(input_array_hat = GENIE3_mat_array_hat_density2,
                                             type = 3,
                                             Theta_array_true = Theta_mat_array_TRUE,
                                             cluster_predict_true_vec = cluster_predict_true_vec)
    
    evaluator_GENIE3_minsample<-evaluator_1_GENIE3[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_GENIE3_minsample_density<-evaluator_1_GENIE3_density[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_GENIE3_minsample_density2<-evaluator_1_GENIE3_density2[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_GENIE3_minmisc<-evaluator_1_GENIE3[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_GENIE3_minmisc_density<-evaluator_1_GENIE3_density[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_GENIE3_minmisc_density2<-evaluator_1_GENIE3_density2[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_GENIE3<-rbind(evaluator_1_GENIE3,evaluator_GENIE3_minsample,evaluator_GENIE3_minmisc)
    evaluator_2_GENIE3_density<-rbind(evaluator_1_GENIE3_density,evaluator_GENIE3_minsample_density,evaluator_GENIE3_minmisc_density)
    evaluator_2_GENIE3_density2<-rbind(evaluator_1_GENIE3_density2,evaluator_GENIE3_minsample_density2,evaluator_GENIE3_minmisc_density2)
    
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    #############################
    res_vec<-c(
      time_use,
      ##1
      
      as.vector(evaluator_2_GENIE3),
      
      as.vector(evaluator_2_GENIE3_density),
      
      as.vector(evaluator_2_GENIE3_density2),
      
      edge_true_num
      
    )
    ##
    repres_GENIE3<-rbind(repres_GENIE3,res_vec)
  }
  
  ##stability
  early_stability_BIC<-NULL
  early_stability_density<-NULL
  early_stability_density2<-NULL
  if(if_stability_eval == TRUE){
    #load network and cluster_predict_true_vec
    edge_list<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_GENIE3[rep,(ncol(repres_GENIE3)-2):ncol(repres_GENIE3)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3nor_mat_array_hat_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      for(g in 1:dim(GENIE3_mat_array_hat)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(GENIE3_mat_array_hat[g,,]>t(GENIE3_mat_array_hat[g,,]),GENIE3_mat_array_hat[g,,],t(GENIE3_mat_array_hat[g,,]))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
      }
      edge_list[[rep]]<-edge_list_rep
      
    }
    ##
    early_stability_BIC<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density
    edge_list<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_GENIE3[rep,(ncol(repres_GENIE3)-2):ncol(repres_GENIE3)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3nor_mat_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      for(g in 1:dim(GENIE3_mat_array_hat_density)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(GENIE3_mat_array_hat_density[g,,]>t(GENIE3_mat_array_hat_density[g,,]),GENIE3_mat_array_hat_density[g,,],t(GENIE3_mat_array_hat_density[g,,]))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
      }
      edge_list[[rep]]<-edge_list_rep
      
    }
    ##
    early_stability_density<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density2
    edge_list<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_GENIE3[rep,(ncol(repres_GENIE3)-2):ncol(repres_GENIE3)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/GENIE3nor_mat_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/GENIE3/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      for(g in 1:dim(GENIE3_mat_array_hat_density2)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(GENIE3_mat_array_hat_density2[g,,]>t(GENIE3_mat_array_hat_density2[g,,]),GENIE3_mat_array_hat_density2[g,,],t(GENIE3_mat_array_hat_density2[g,,]))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
      }
      edge_list[[rep]]<-edge_list_rep
      
    }
    ##
    early_stability_density2<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    
  }
  
  return(list(repres_GENIE3 = repres_GENIE3,
              early_stability_BIC = early_stability_BIC,
              early_stability_density = early_stability_density,
              early_stability_density2 = early_stability_density2))
}

##Evaluator for LPGM
LPGM_evaluator<-function(penalize_diagonal = FALSE,if_stability_eval = FALSE){
  library_size_est = "TSS"
  ##
  ############
  repres_LPGM<-matrix(NA,nrow = 0,ncol = 184)
  for(rep in 1:reptime){
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    ##
    bic_object_list<-list()
    ##
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    ls_vec<-res_init$ls
    ##
    if(library_size_est == "none"){
      ls_vec<-rep(1,length(res_init$ls))
    }
    if(library_size_est == "TSS"){
      ls_vec<-rowSums(obs_mat)/10000
    }
    if(library_size_est == "GMPR"){
      ls_vec<-as.vector(GMPR(obs_mat))
    }
    ##
    locator_base = apply(res_init$U_mat, 1, which.max) 
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    
    ###
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    ##
    diag(table_reshape)<-0
    rm(table_reshape);gc()
    ##################
    
    print("LPGM")
    set.seed(123)
    
    
    
    ##
    time_cost<-c()
    gene_name_vec<-colnames(obs_mat)
    LPGM_mat_array_hat<-array(NA,dim = c(group_num,p,p))
    LPGM_mat_array_hat_density<-array(NA,dim = c(group_num,p,p))
    LPGM_mat_array_hat_density2<-array(NA,dim = c(group_num,p,p))
    for (g2 in 1:group_num) {
      print(paste("group(",g2,") task of LPGM."))
      time_LPGM_start<-Sys.time()
      
      
      ##########
      lpgm.fit <- XMRF(X = t(as.matrix(obs_mat[locator_base==g2,])),
                       method="LPGM",
                       nlams=100,
                       beta = 0.01,
                       stability="STAR",
                       th=0.005,
                       parallel = TRUE,
                       N = 50,
                       lmin = 1e-5,
                       sth = 0.95,
                       nCpus = num_core)
      weightMat<-lpgm.fit$D[[lpgm.fit$opt.index]] * lpgm.fit$network[[lpgm.fit$opt.index]]
      weightMat<-ifelse(is.nan(weightMat),0,weightMat)
      #########
      LPGM_mat_array_hat[g2,,]<-weightMat
      
      ###
      density_vec<-c()
      for(k in 1:length(lpgm.fit$network)){
        density_vec<-c(density_vec,length(which((lpgm.fit$network[[k]])[upper.tri(lpgm.fit$network[[k]])]!=0))/length((lpgm.fit$network[[k]])[upper.tri(lpgm.fit$network[[k]])]))
      }
      
      #density1
      # density_choose<-which.min(abs(density_vec - 0.1))
      density_choose <- (which(density_vec > densy_degree))[which.min(abs(density_vec[which(density_vec > densy_degree)] - densy_degree))]
      if(length(density_choose) == 0){
        density_choose<-which.max(density_vec)
      }
      #
      weightMat<-lpgm.fit$D[[density_choose]] * lpgm.fit$network[[density_choose]]
      weightMat<-ifelse(is.nan(weightMat),0,weightMat)
      #########
      LPGM_mat_array_hat_density[g2,,]<-weightMat
      
      #density2
      # density_choose<-which.min(abs(density_vec - 0.2))
      density_choose <- (which(density_vec > 2*densy_degree))[which.min(abs(density_vec[which(density_vec > 2*densy_degree)] - 2*densy_degree))]
      if(length(density_choose) == 0){
        density_choose<-which.max(density_vec)
      }
      weightMat<-lpgm.fit$D[[density_choose]] * lpgm.fit$network[[density_choose]]
      weightMat<-ifelse(is.nan(weightMat),0,weightMat)
      #########
      LPGM_mat_array_hat_density2[g2,,]<-weightMat
      ##
      time_LPGM_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_LPGM_end,time_LPGM_start,units = "hours")))
    }
    time_use<-sum(time_cost)
    
    
    cluster_predict_true_vec<-as.vector(apply(table(apply(res_init$U_mat,MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max)) 
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    
    ##
    save(LPGM_mat_array_hat,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/LPGM_mat_array_hat_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                               "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(LPGM_mat_array_hat_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/LPGM_mat_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(LPGM_mat_array_hat_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/LPGM_mat_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                        "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
    ##
    evaluator_1_LPGM<-evaluator_1(input_array_hat = LPGM_mat_array_hat,
                                  type = 3,
                                  Theta_array_true = Theta_mat_array_TRUE,
                                  cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_LPGM_density<-evaluator_1(input_array_hat = LPGM_mat_array_hat_density,
                                          type = 3,
                                          Theta_array_true = Theta_mat_array_TRUE,
                                          cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_LPGM_density2<-evaluator_1(input_array_hat = LPGM_mat_array_hat_density2,
                                           type = 3,
                                           Theta_array_true = Theta_mat_array_TRUE,
                                           cluster_predict_true_vec = cluster_predict_true_vec)
    
    evaluator_LPGM_minsample<-evaluator_1_LPGM[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_LPGM_minsample_density<-evaluator_1_LPGM_density[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_LPGM_minsample_density2<-evaluator_1_LPGM_density2[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_LPGM_minmisc<-evaluator_1_LPGM[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_LPGM_minmisc_density<-evaluator_1_LPGM_density[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_LPGM_minmisc_density2<-evaluator_1_LPGM_density2[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_LPGM<-rbind(evaluator_1_LPGM,evaluator_LPGM_minsample,evaluator_LPGM_minmisc)
    evaluator_2_LPGM_density<-rbind(evaluator_1_LPGM_density,evaluator_LPGM_minsample_density,evaluator_LPGM_minmisc_density)
    evaluator_2_LPGM_density2<-rbind(evaluator_1_LPGM_density2,evaluator_LPGM_minsample_density2,evaluator_LPGM_minmisc_density2)
    #############################
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    ##
    res_vec<-c(
      time_use,
      ##1
      
      as.vector(evaluator_2_LPGM),
      
      as.vector(evaluator_2_LPGM_density),
      
      as.vector(evaluator_2_LPGM_density2),
      
      edge_true_num
      
    )
    ##
    repres_LPGM<-rbind(repres_LPGM,res_vec)
  }
  # #####
  # return(res_vec)
  
  
  # }
  
  early_stability_BIC<-NULL
  early_stability_density<-NULL
  early_stability_density2<-NULL
  if(if_stability_eval == TRUE){
    #load network and cluster_predict_true_vec
    edge_list<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_LPGM[rep,(ncol(repres_LPGM)-2):ncol(repres_LPGM)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/LPGM_mat_array_hat_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      for(g in 1:dim(LPGM_mat_array_hat)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(LPGM_mat_array_hat[g,,]>t(LPGM_mat_array_hat[g,,]),LPGM_mat_array_hat[g,,],t(LPGM_mat_array_hat[g,,]))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
      }
      edge_list[[rep]]<-edge_list_rep
      
    }
    ##
    early_stability_BIC<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    
    #density1
    edge_list<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_LPGM[rep,(ncol(repres_LPGM)-2):ncol(repres_LPGM)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/LPGM_mat_array_hat_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      for(g in 1:dim(LPGM_mat_array_hat_density)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(LPGM_mat_array_hat_density[g,,]>t(LPGM_mat_array_hat_density[g,,]),LPGM_mat_array_hat_density[g,,],t(LPGM_mat_array_hat_density[g,,]))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
      }
      edge_list[[rep]]<-edge_list_rep
      
    }
    ##
    early_stability_density<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    
    #density2
    edge_list<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_LPGM[rep,(ncol(repres_LPGM)-2):ncol(repres_LPGM)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/LPGM_mat_array_hat_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/LPGM/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_penalize_diagonal_",penalize_diagonal,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list()
      for(g in 1:dim(LPGM_mat_array_hat_density2)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(LPGM_mat_array_hat_density2[g,,]>t(LPGM_mat_array_hat_density2[g,,]),LPGM_mat_array_hat_density2[g,,],t(LPGM_mat_array_hat_density2[g,,]))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
      }
      edge_list[[rep]]<-edge_list_rep
      
    }
    ##
    early_stability_density2<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
  }
  
  return(list(repres_LPGM = repres_LPGM,
              early_stability_BIC = early_stability_BIC,
              early_stability_density = early_stability_density,
              early_stability_density2 = early_stability_density2))
  
}

##Evaluator for SCODE
run_SCODE<-function(obs_mat,
                    S_depth = NULL,
                    tfnum,
                    pnum,
                    maxite,
                    rep_time = 10){
  
  ##Step 1: run Monocle for pseudo time
  library(Seurat)
  library(monocle)
  cnum<-nrow(obs_mat)
  ##normal
  if(is.null(S_depth)){
    X<-t(log(obs_mat/as.vector(rowSums(obs_mat)) * 10000 + 1)) 
  }else{
    X<-t(log(obs_mat/S_depth + 1))
  }
  rownames(X)<-paste("Gene",1:nrow(X),sep = "")
  colnames(X)<-paste("Cell",1:ncol(X),sep = "")
  Seurat_object<-CreateSeuratObject(counts = X)
  
  ##
  pd = data.frame(gene_short_name = paste("Gene",1:nrow(X),sep = ""))
  rownames(pd)<-paste("Gene",1:nrow(X),sep = "")
  pd<-new('AnnotatedDataFrame', data = pd)
  fd = data.frame(gene_short_name = paste("Cell",1:ncol(X),sep = ""))
  rownames(fd)<-paste("Cell",1:ncol(X),sep = "")
  fd<-new('AnnotatedDataFrame', data = fd)
  HSMM <- newCellDataSet(as.matrix(X),
                         phenoData = fd, featureData = pd)
  ##
  HSMM <- estimateSizeFactors(HSMM)
  # HSMM <- estimateDispersions(HSMM) 
  # HSMM <- detectGenes(HSMM, min_expr = 3)
  # expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))
  ##
  HSMM <- reduceDimension(HSMM, max_components = 2,
                          method = 'DDRTree')
  HSMM <- orderCells(HSMM)
  pseudo_time<-HSMM$Pseudotime
  # plot_cell_trajectory(HSMM)
  pseudo_time_data<-cbind(rep(0,length(pseudo_time)),as.vector(pseudo_time))
  # write.table(pseudo_time_data,file = ftime,col.names = FALSE,row.names = FALSE,
  #             sep = "\t")
  # ##
  # disp_table <- dispersionTable(HSMM)
  # unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  # HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
  # ##
  # HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 10,reduction_method = 'tSNE', verbose = T)
  # HSMM <- clusterCells(HSMM, num_clusters = group_num)
  # ##
  # diff_test_res <- differentialGeneTest(cds = HSMM[expressed_genes,])
  
  ####################################################
  ##Step 2: use SCODE for Network inference
  print("Start run SCODE.")
  A_list<-list()
  for(k in 1:rep_time){
    print(paste("The current rep time is ",k,sep = ""))
    maxB <- 2.0
    minB <- -10.0
    
    # system(paste("mkdir", dir, sep=" "))
    
    # X <- as.matrix(read.table(fdata, sep="\t"))[1:tfnum,1:cnum]
    W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
    Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
    WZ <- matrix(nrow=tfnum, ncol=cnum)
    
    #read pseudo-time and normalize pseudo-time
    pseudotime<-pseudo_time
    pseudotime <- pseudotime/max(pseudotime)
    
    new_B <- rep(0, pnum)
    old_B <- rep(0, pnum)
    
    #initialization
    RSS <- Inf
    for(i in 1:pnum){
      new_B[i] <- runif(1, min=minB, max=maxB)
      old_B[i] <- new_B[i]
    }
    
    #function to sample Z
    sample_Z <- function(){
      for(i in 1:pnum){
        for(j in 1:cnum){
          Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
        }
      }
    }
    
    #optimize W and B iteratively
    for(ite in 1:maxite){
      #sampling B
      target <- floor(runif(1, min=1, max=pnum+1))
      new_B[target] <- runif(1, min=minB, max=maxB)
      
      #for last calculation
      if(ite == maxite){
        for(i in 1:pnum){
          new_B[i] <- old_B[i]
        }
      }
      
      #sample Z from new B
      sample_Z()
      
      #regression
      for(i in 1:tfnum){
        X.lm <- lm(X[i,] ~ t(Z)-1)
        for(j in 1:pnum){
          W[i,j] <- X.lm$coefficients[j]
        }
        WZ[i,] <- W[i,] %*% Z
      }
      
      #RSS
      tmp_RSS <- sum((X-WZ)**2)
      if(tmp_RSS < RSS){
        RSS <- tmp_RSS
        old_B[target] <- new_B[target]
      }
      else{
        new_B[target] <- old_B[target]
      }
    }
    
    # #output RSS
    # write.table(RSS, paste(dir,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")
    # 
    # #output W
    # write.table(W, paste(dir,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")
    
    #infer A
    B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
    for(i in 1:pnum){
      B[i,i] <- new_B[i]
    }
    invW <- MASS::ginv(W)
    A <- W %*% B %*% invW
    
    # #output A and B
    # write.table(A, paste(dir,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")
    # write.table(B, paste(dir,"/B.txt",sep=""), row.names=F, col.names=F, sep="\t")
    A_list[[k]]<-A
  }
  
  ##summarized
  A_mean<-Reduce("+",A_list)/rep_time
  return(A_mean)
}
##
SCODE_evaluator<-function(pnum = 4,maxite = 1000,if_stability_eval = FALSE){
  library_size_est = "TSS"
  ##
  repres_SCODE<-foreach (
    rep = 1:reptime,
    .combine = cfun,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS','mclust','Seurat','monocle')
  )%dopar%{
    ## load data-------------
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    load(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
               "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
    ##
    bic_object_list<-list()
    ##
    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }
    
    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dr = length(which( obs_mat == 0 )) / (n*p)
    
    ls_vec<-res_init$ls
    ##
    if(library_size_est == "none"){
      ls_vec<-rep(1,length(res_init$ls))
    }
    if(library_size_est == "TSS"){
      ls_vec<-rowSums(obs_mat)/10000
    }
    if(library_size_est == "GMPR"){
      ls_vec<-as.vector(GMPR(obs_mat))
    }
    ##
    locator_base = apply(res_init$U_mat, 1, which.max) 
    ##
    ARI_orgin = adjustedRandIndex(cluster_true,locator_base)
    cluster_predict_true_base<-as.vector(apply(table(locator_base,cluster_true),MARGIN = 1,which.max))
    
    ###
    table_reshape<-matrix(NA,nrow = group_num,ncol = group_num)
    for(g1 in 1:group_num){
      for(g2 in 1:group_num){
        table_reshape[g1,g2]<-length(which(locator_base == g1 & cluster_true == cluster_predict_true_base[g2]))
      }
    }
    bic_object_list[["cluster_predicted_sample"]]<-rowSums(table_reshape)
    bic_object_list[["cluster_predicted_TDR"]]<-diag(table_reshape)/rowSums(table_reshape)
    ##
    diag(table_reshape)<-0
    rm(table_reshape);gc()
    ##################
    
    print("SCODE")
    set.seed(123)
    
    
    ##
    time_cost<-c()
    gene_name_vec<-colnames(obs_mat)
    SCODE_mat_array_hat<-array(NA,dim = c(group_num,p,p))
    for (g2 in 1:group_num) {
      print(paste("group(",g2,") task of SCODE."))
      time_SCODE_start<-Sys.time()
      ##
      weightMat<-run_SCODE(obs_mat = obs_mat[locator_base==g2,],
                           S_depth = NULL,
                           tfnum = ncol(obs_mat),
                           pnum = pnum,
                           maxite = maxite,
                           rep_time = 5)
      weightMat<-ifelse(is.nan(weightMat),0,weightMat)
      SCODE_mat_array_hat[g2,,]<-weightMat
      ##
      time_SCODE_end<-Sys.time()
      time_cost<-c(time_cost,as.double(difftime(time_SCODE_end,time_SCODE_start,units = "hours")))
    }
    time_use<-sum(time_cost)
    ##
    cluster_predict_true_vec<-as.vector(apply(table(apply(res_init$U_mat,MARGIN = 1,which.max),cluster_true),MARGIN = 1,which.max)) 
    bic_object_list[["cluster_predict_true_vec"]]<-cluster_predict_true_vec
    ##
    ##
    SCODE_mat_array_hat_max<-SCODE_mat_array_hat
    for(g in 1:group_num){
      mat_input<-matrix(0,nrow = nrow(SCODE_mat_array_hat_max[g,,]),ncol = ncol(SCODE_mat_array_hat_max[g,,]))
      aaa0<-SCODE_mat_array_hat_max[g,,]
      mat_input[upper.tri(mat_input)]<-ifelse(abs(aaa0[upper.tri(aaa0)])>abs(t(aaa0)[upper.tri(t(aaa0))]),aaa0[upper.tri(aaa0)],t(aaa0)[upper.tri(t(aaa0))])
      mat_input<-mat_input + t(mat_input)
      SCODE_mat_array_hat_max[g,,]<-mat_input
    }
    ##
    SCODE_mat_array_hat_abs<-SCODE_mat_array_hat_max
    for(g in 1:group_num){
      SCODE_mat_array_hat_abs[g,,]<-abs(SCODE_mat_array_hat_abs[g,,])
    }
    SCODE_mat_array_hat_abs_density<-SCODE_mat_array_hat_abs
    SCODE_mat_array_hat_abs_density2<-SCODE_mat_array_hat_abs
    for(g in 1:group_num){
      mat_use<-SCODE_mat_array_hat_abs_density[g,,]
      mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
      mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
      SCODE_mat_array_hat_abs_density[g,,]<-mat_use
    }
    for(g in 1:group_num){
      mat_use<-SCODE_mat_array_hat_abs_density2[g,,]
      mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
      mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
      SCODE_mat_array_hat_abs_density2[g,,]<-mat_use
    }
    
    SCODE_mat_array_hat_max_density<-SCODE_mat_array_hat_max
    SCODE_mat_array_hat_max_density2<-SCODE_mat_array_hat_max
    for(g in 1:group_num){
      mat_use<-SCODE_mat_array_hat_max_density[g,,]
      mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
      mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.1))]]<-0
      SCODE_mat_array_hat_max_density[g,,]<-mat_use
    }
    for(g in 1:group_num){
      mat_use<-SCODE_mat_array_hat_max_density2[g,,]
      mat_use[upper.tri(mat_use)][order(abs(mat_use[upper.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
      mat_use[lower.tri(mat_use)][order(abs(mat_use[lower.tri(mat_use)]),decreasing = FALSE)[1:ceiling((ncol(mat_use) * (ncol(mat_use) - 1)/2) * (1-0.2))]]<-0
      SCODE_mat_array_hat_max_density2[g,,]<-mat_use
    }
    
    ##
    ##
    save(SCODE_mat_array_hat_max,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_max_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
    save(SCODE_mat_array_hat_max_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_max_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                            "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
    save(SCODE_mat_array_hat_max_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_max_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                             "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
    
    save(SCODE_mat_array_hat_abs,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_abs_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                               "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
    save(SCODE_mat_array_hat_abs_density,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_abs_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
    save(SCODE_mat_array_hat_abs_density2,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_abs_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                    "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
    save(cluster_predict_true_vec,file = paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                                                     "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
    ##
    
    
    evaluator_1_SCODE<-evaluator_1(input_array_hat = SCODE_mat_array_hat_abs,
                                   type = 3,
                                   Theta_array_true = Theta_mat_array_TRUE,
                                   cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_SCODE_density<-evaluator_1(input_array_hat = SCODE_mat_array_hat_abs_density,
                                   type = 3,
                                   Theta_array_true = Theta_mat_array_TRUE,
                                   cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_SCODE_density2<-evaluator_1(input_array_hat = SCODE_mat_array_hat_abs_density2,
                                   type = 3,
                                   Theta_array_true = Theta_mat_array_TRUE,
                                   cluster_predict_true_vec = cluster_predict_true_vec)
    
    evaluator_SCODE_minsample<-evaluator_1_SCODE[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_SCODE_minsample_density<-evaluator_1_SCODE_density[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_SCODE_minsample_density2<-evaluator_1_SCODE_density2[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_SCODE_minmisc<-evaluator_1_SCODE[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_SCODE_minmisc_density<-evaluator_1_SCODE_density[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_SCODE_minmisc_density2<-evaluator_1_SCODE_density2[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_SCODE<-rbind(evaluator_1_SCODE,evaluator_SCODE_minsample,evaluator_SCODE_minmisc)
    evaluator_2_SCODE_density<-rbind(evaluator_1_SCODE_density,evaluator_SCODE_minsample_density,evaluator_SCODE_minmisc_density)
    evaluator_2_SCODE_density2<-rbind(evaluator_1_SCODE_density2,evaluator_SCODE_minsample_density2,evaluator_SCODE_minmisc_density2)
    
    
    ##signed
    evaluator_1_SCODE_signed<-evaluator_1_signed(input_array_hat = SCODE_mat_array_hat_max,
                                   type = 2,
                                   Theta_array_true = Theta_mat_array_TRUE,
                                   cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_SCODE_density_signed<-evaluator_1_signed(input_array_hat = SCODE_mat_array_hat_max_density,
                                           type = 2,
                                           Theta_array_true = Theta_mat_array_TRUE,
                                           cluster_predict_true_vec = cluster_predict_true_vec)
    evaluator_1_SCODE_density2_signed<-evaluator_1_signed(input_array_hat = SCODE_mat_array_hat_max_density2,
                                            type = 2,
                                            Theta_array_true = Theta_mat_array_TRUE,
                                            cluster_predict_true_vec = cluster_predict_true_vec)
    
    evaluator_SCODE_minsample_signed<-evaluator_1_SCODE_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_SCODE_minsample_density_signed<-evaluator_1_SCODE_density_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    evaluator_SCODE_minsample_density2_signed<-evaluator_1_SCODE_density2_signed[which.min(bic_object_list[["cluster_predicted_sample"]]),]
    
    evaluator_SCODE_minmisc_signed<-evaluator_1_SCODE_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_SCODE_minmisc_density_signed<-evaluator_1_SCODE_density_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    evaluator_SCODE_minmisc_density2_signed<-evaluator_1_SCODE_density2_signed[which.min(bic_object_list[["cluster_predicted_TDR"]]),]
    
    evaluator_2_SCODE_signed<-rbind(evaluator_1_SCODE_signed,evaluator_SCODE_minsample_signed,evaluator_SCODE_minmisc_signed)
    evaluator_2_SCODE_density_signed<-rbind(evaluator_1_SCODE_density_signed,evaluator_SCODE_minsample_density_signed,evaluator_SCODE_minmisc_density_signed)
    evaluator_2_SCODE_density2_signed<-rbind(evaluator_1_SCODE_density2_signed,evaluator_SCODE_minsample_density2_signed,evaluator_SCODE_minmisc_density2_signed)
    
    ##
    edge_true_num<-c()
    for(g in 1:group_num){
      edge_true_num<-c(edge_true_num,length(which((Theta_mat_array_TRUE[g,,])[upper.tri((Theta_mat_array_TRUE[g,,]))] != 0)))
    }
    ##
    res_vec<-c(
      time_use,
      ##1
      
      as.vector(evaluator_2_SCODE),
      
      as.vector(evaluator_2_SCODE_density),
      
      as.vector(evaluator_2_SCODE_density2),
      
      ##signed
      as.vector(evaluator_2_SCODE_signed),
      
      as.vector(evaluator_2_SCODE_density_signed),
      
      as.vector(evaluator_2_SCODE_density2_signed),
      
      edge_true_num
      
    )
    ##
    return(res_vec)
  }
  ##
  early_stability_BIC<-NULL
  early_stability_density<-NULL
  early_stability_density2<-NULL
  early_stability_BIC_signed<-NULL
  early_stability_density_signed<-NULL
  early_stability_density2_signed<-NULL
  ##
  if(if_stability_eval == TRUE){
    #load network and cluster_predict_true_vec
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_SCODE[rep,(ncol(repres_SCODE)-2):ncol(repres_SCODE)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_abs_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(SCODE_mat_array_hat_abs)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(abs(SCODE_mat_array_hat_abs[g,,])>abs(t(SCODE_mat_array_hat_abs[g,,])),abs(SCODE_mat_array_hat_abs[g,,]),abs(t(SCODE_mat_array_hat_abs[g,,])))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:min(edge_true_num[g],length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0)))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_BIC<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_BIC_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density1
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_SCODE[rep,(ncol(repres_SCODE)-2):ncol(repres_SCODE)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_abs_density_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(SCODE_mat_array_hat_abs_density)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(abs(SCODE_mat_array_hat_abs_density[g,,])>abs(t(SCODE_mat_array_hat_abs_density[g,,])),abs(SCODE_mat_array_hat_abs_density[g,,]),abs(t(SCODE_mat_array_hat_abs_density[g,,])))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
    ##density2
    edge_list<-list();edge_list_signed<-list()
    cluster_predict_true_mat<-matrix(NA,nrow = 0,ncol = 3)
    for(rep in 1:reptime){
      edge_true_num<-repres_SCODE[rep,(ncol(repres_SCODE)-2):ncol(repres_SCODE)]
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/SCODE_mat_array_hat_abs_density2_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
      load(paste(paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_all/SCODE/cluster_predict_true_vec_1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                       "_nettype",network_type[net],"_ARI_",ARI_i,"_library_size_est_",library_size_est,"_rep",rep,".Rdata",sep = "")))
      ##
      cluster_predict_true_mat<-rbind(cluster_predict_true_mat,cluster_predict_true_vec)
      ##
      edge_list_rep<-list();edge_list_rep_signed<-list()
      for(g in 1:dim(SCODE_mat_array_hat_abs_density2)[1]){
        # corr_mat<-(-1) * Theta_array_hat_BIC[g,,]/(matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),ncol = 1) %*% matrix(sqrt(diag(Theta_array_hat_BIC[g,,])),nrow = 1))
        # diag(corr_mat)<-1
        weight_mat_use<-ifelse(abs(SCODE_mat_array_hat_abs_density2[g,,])>abs(t(SCODE_mat_array_hat_abs_density2[g,,])),abs(SCODE_mat_array_hat_abs_density2[g,,]),abs(t(SCODE_mat_array_hat_abs_density2[g,,])))
        corr_mat<-weight_mat_use
        edge_index_vec<-order(abs(corr_mat)[upper.tri(abs(corr_mat))],decreasing = TRUE)[1:length(which(abs(corr_mat)[upper.tri(abs(corr_mat))] !=0))]
        edge_list_rep[[g]]<-edge_index_vec
        edge_list_rep_signed[[g]]<-edge_index_vec * (sign(corr_mat[upper.tri(corr_mat)]))[edge_index_vec]
      }
      edge_list[[rep]]<-edge_list_rep
      edge_list_signed[[rep]]<-edge_list_rep_signed
      
    }
    ##
    early_stability_density2<-evalulator_early_stability(edge_list = edge_list,cluster_predict_true_mat = cluster_predict_true_mat)
    early_stability_density2_signed<-evalulator_early_stability(edge_list = edge_list_signed,cluster_predict_true_mat = cluster_predict_true_mat)
    
  }
  
  return(list(repres_SCODE = repres_SCODE,
              early_stability_BIC = early_stability_BIC,
              early_stability_density = early_stability_density,
              early_stability_density2 = early_stability_density2,
              early_stability_BIC_signed = early_stability_BIC_signed,
              early_stability_density_signed = early_stability_density_signed,
              early_stability_density2_signed = early_stability_density2_signed))
  
}
##--------------------------------------------------------------------------------------------


##4. Simulation setting
##--------------------------------------------------------------------------------------------

##4.1 overall setting
#####################################################################
##generate way
generator_type_vec<-c("MPLN","compositional")
##dimension
dim_use_vec<-c(100,300,500)

##sample size for each cell type
sample_size_mat<-t(sapply(X=1:length(dim_use_vec),FUN = function(i){
  return(c(1000,2000))
}))

mu_design = matrix(c(2.4,-0.1,0.9,-0.1,
                     1.4,-1.1,-0.1,-1.1,
                     0.7,-1.8,-0.8,-1.8),ncol = 4,byrow = T)
## represent different dropout levels: 10% 30% 50%

##number of group
group_num_vec<-c(3)

#network type
network_type=c("ER-graph","HUB-graph","AFF-graph","PA-graph")

##different mixing degree of component evaluated by ARI
ARI_deisgn<-matrix(c(0.9,1,0.75,0.85,0.65,0.75),ncol = 2,byrow = TRUE)

## The density of network and the element of precision matrix
densy_degree_vec<-c(0.1)
p_v_vec<-c(0.5)
v<-0.3

##Initialized clustering method
cluster_method_vec<-c("Kmeans","SNN")

#time of experiment
reptime<-100
#####################################################################

##4.2 choose the specific setting for testing
#####################################################################
generator_type<-"MPLN"
run_index<-c(1,1,1,1,1)
s_i<-run_index[1]
p_i<-run_index[2]
m_i<-run_index[3]
net<-run_index[4]
ARI_i<-run_index[5]
##
pv<-1
g_i<-1
clustermethod_i<-1
densy_degree_i<-1
if_fix<-FALSE
#####################################################################

##4.3 parallel setting (Optional)
#####################################################################
num_core<-20
cl <- makeCluster(num_core, outfile = paste("debug_",paste(run_index,collapse = ""),".txt",sep = ""))
registerDoParallel(cl)
#combine function
cfun<-function(a,b){
  rbind(a,b)
}
#####################################################################


##5. Run
##--------------------------------------------------------------------------------------------

##5.1 Set up the current parameter settings
#####################################################################
print(paste("s_i:",s_i,"p_i:",p_i," m_i:",m_i," net:",net," ARI_i:",ARI_i," clustermethod_i",clustermethod_i,sep = ""))
dim_use = dim_use_vec[p_i]
group_num = group_num_vec[g_i]
sample_size_vec<-rep(sample_size_mat[p_i,s_i],group_num)
n = sum(sample_size_vec)
p_v = p_v_vec[pv]
mu_design_vec = mu_design[m_i,]
nwtype = network_type[net]
cluster_method<-cluster_method_vec[clustermethod_i]
densy_degree<-densy_degree_vec[densy_degree_i]
ARI_low<-ARI_deisgn[ARI_i,1]
ARI_upper<-ARI_deisgn[ARI_i,2]
#####################################################################

##5.2 Generate the simulation data
#####################################################################
diff_gene_num<-10
for(rep in 1:reptime){
  ##Generate data ------------------
  print(rep)
  con_use<-TRUE
  iter0<-1
  while (con_use & iter0<100) {
    iter0<-iter0+1
    set.seed(123*rep + iter0)
    if((if_fix == TRUE) & (rep>1)){
      if(generator_type == "MPLN"){
        data = data_generator(n,dim_use,group_num,network_type = nwtype,
                              densy_degree = densy_degree,v = v,p_v = p_v,
                              num_de = diff_gene_num, ub =mu_design_vec[1],
                              lb = mu_design_vec[2],ub_non = mu_design_vec[3],
                              lb_non = mu_design_vec[4],
                              l_mu = log(10),sigma_scale = 1,alpha_com = 0.01,
                              if_fix = TRUE,fix_element = data_fix) 
      }else{
        data = data_generator_mulitnomial(n,dim_use,group_num,network_type = nwtype,
                                              densy_degree = densy_degree,v = v,p_v = p_v,
                                              num_de = diff_gene_num, ub =mu_design_vec[1],
                                              lb = mu_design_vec[2],ub_non = mu_design_vec[3],
                                              lb_non = mu_design_vec[4],
                                              l_mu = log(10),sigma_scale = 1,alpha_com = 0.01,
                                              if_fix = TRUE,fix_element = data_fix)
      }
    }else{
      if(generator_type == "MPLN"){
        data = data_generator(n,dim_use,group_num,network_type = nwtype,
                              densy_degree = densy_degree,v = v,p_v = p_v,
                              num_de = diff_gene_num, ub =mu_design_vec[1],
                              lb = mu_design_vec[2],ub_non = mu_design_vec[3],
                              lb_non = mu_design_vec[4],
                              l_mu = log(10),sigma_scale = 1,alpha_com = 0.01,
                              if_fix = FALSE) 
      }else{
        data = data_generator_mulitnomial(n,dim_use,group_num,network_type = nwtype,
                                              densy_degree = densy_degree,v = v,p_v = p_v,
                                              num_de = diff_gene_num, ub =mu_design_vec[1],
                                              lb = mu_design_vec[2],ub_non = mu_design_vec[3],
                                              lb_non = mu_design_vec[4],
                                              l_mu = log(10),sigma_scale = 1,alpha_com = 0.01,
                                              if_fix = FALSE)
      }
    }

    obs_mat = data[["obs_mat"]]
    colnames(obs_mat)<-paste("Gene",1:ncol(obs_mat),sep = "")
    rownames(obs_mat)<-paste("Cell",1:nrow(obs_mat),sep = "")
    Theta_mat_array_TRUE =data[["group_specific_precision"]]
    cluster_true = data[["locator_vec"]]
    # hub_gene<-data[["hub_gene"]]
    if(length(data[["hub_gene"]])>0){
      hub_gene<-data[["hub_gene"]]
    }else{
      hub_gene<-NULL
    }

    n<-dim(obs_mat)[1]
    p<-dim(obs_mat)[2]
    dr = length(which( obs_mat == 0 )) / (n*p)

    # 1.0 initial part
    time_init_start = Sys.time()
    gene_GRN_use<-rownames(as(t(obs_mat),"sparseMatrix"))
    suppressWarnings(res_init<-scGeneNet_init(expression_profile = as(t(obs_mat),"sparseMatrix"),
                                          celltypes_num = group_num,
                                          celltypes_ref = NULL,
                                          ls_est = "TSS",
                                          gene_GRN = gene_GRN_use,
                                          HVG_model_num = 0,
                                          zero_GRN = NULL,
                                          preprocess_Control = list(HVG_num = length(gene_GRN_use),npc = 50,
                                                                    run_umap = FALSE,label_umap = NULL,
                                                                    cluster_method = "Kmeans",resolution = 0.8),
                                          core_num = 1))


    time_init_end<-Sys.time()
    time_init<-as.double(difftime(time_init_end,time_init_start,units = "hours"))

    locator = res_init$celltypes_label
    ARI_orgin = adjustedRandIndex(cluster_true,locator)
    # ##
    cluster_predict_true_base<-as.vector(apply(table(locator,cluster_true),MARGIN = 1,which.max))
    
    #
    if(ARI_orgin<ARI_low){
      diff_gene_num<-diff_gene_num + 1
      diff_gene_num<-max(diff_gene_num,1)
    }else{
      if(ARI_orgin>ARI_upper){
        diff_gene_num<-diff_gene_num - 1
        diff_gene_num<-max(diff_gene_num,1)
      }else{
        #check whether have paired component predicted
        if(length(unique(cluster_predict_true_base)) == group_num){
          con_use<-FALSE
        }
      }
    }
    print(paste("diff_gene_num: ",diff_gene_num," con_use: ",con_use," ARI_orgin: ",ARI_orgin,sep = ""))

  }
  ##
  if((if_fix == TRUE) & (rep == 1)){
    data_fix<-data
  }
  ##add noise to diff_gene_num
  diff_gene_num<-diff_gene_num+sample(c(-1:1),1)
  diff_gene_num<-max(diff_gene_num,1)
  #save
  save(data,file = paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/data1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                         "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
  save(res_init,file = paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/data/res_init1015_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                             "_nettype",network_type[net],"_ARI:",ARI_i,"_clustermethod",clustermethod_i,"_rep",rep,".Rdata",sep = ""))
}
print("Data generation complete.")
#####################################################################

##5.3 evaluate the performance for each method
#####################################################################
repres_all<-list()
##5.3.1 scGeneNet (penalized diagonal)
print("scGeneNet start (penalized diagonal).")
repres_scGeneNet_pd<-scGeneNet_evaluator(penalize_diagonal = TRUE,if_stability_eval = FALSE)
repres_all[[1]]<-repres_scGeneNet_pd

##5.3.2 scGeneNet (non penalized diagonal)
print("scGeneNet start (not penalized diagonal).")
repres_scGeneNet_nonpd<-scGeneNet_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE)
repres_all[[2]]<-repres_scGeneNet_nonpd


##5.3.3 VPLN (penalized diagonal)
print("VPLN start (penalized diagonal).")
repres_VPLN<-VPLN_evaluator(penalize_diagonal = TRUE,
                            self_lambda = TRUE,if_stability_eval = FALSE)
repres_all[[3]]<-repres_VPLN

##5.3.4 VPLN (non penalized diagonal)
print("VPLN start (not penalized diagonal).")
repres_VPLN_n<-VPLN_evaluator(penalize_diagonal = FALSE,
                              self_lambda = TRUE,if_stability_eval = FALSE)
repres_all[[4]]<-repres_VPLN_n


##5.3.5 Glasso (penalized diagonal)
print("Glasso start (penalized diagonal).")
repres_Glasso<-Glasso_evaluator(penalize_diagonal = TRUE,if_stability_eval = FALSE)
repres_all[[5]]<-repres_Glasso

##5.3.6 Glasso (non penalized diagonal)
print("Glasso start (not penalized diagonal).")
repres_Glasso_n<-Glasso_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE)
repres_all[[6]]<-repres_Glasso_n

##5.3.7 GENIE3 (count as input)
print("GENIE3 start.")
repres_GENIE3<-GENIE3_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE)
repres_all[[7]]<-repres_GENIE3


##5.3.7 PPCOR
print("ppcor start.")
repres_ppcor<-ppcor_evaluator(penalize_diagonal = FALSE,method = "pearson",if_stability_eval = FALSE)
repres_all[[8]]<-repres_ppcor
repres_ppcor<-ppcor_evaluator(penalize_diagonal = FALSE,method = "spearman",if_stability_eval = FALSE)
repres_all[[9]]<-repres_ppcor
repres_ppcor<-ppcor_evaluator(penalize_diagonal = FALSE,method = "kendall",if_stability_eval = FALSE)
repres_all[[10]]<-repres_ppcor

##5.3.8 SCODE
print("SCODE start.")
repres_SCODE<-SCODE_evaluator(pnum = 4,maxite = 1000,if_stability_eval = FALSE)
repres_all[[11]]<-repres_SCODE


##5.3.9 LPGM
print("LPGM start.")
repres_LPGM<-LPGM_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE)
repres_all[[12]]<-repres_LPGM

##5.3.10 GENIE3 (normalized matrix as input)
print("GENIE3 nor start.")
repres_GENIE3<-GENIE3_nor_evaluator(penalize_diagonal = FALSE,if_stability_eval = FALSE)
repres_all[[13]]<-repres_GENIE3
#####################################################################

##5.4save the result
#####################################################################
save(repres_all,file = paste("/home/ruibinxi_pkuhpc/lustre1/tjj/pro2/simulation_1.27.1/result/result_whole/addnor/repres_all_sam",sample_size_mat[p_i,s_i],"_dim",dim_use_vec[p_i],"_drop",m_i,
                             "_nettype",network_type[net],"_ARI",ARI_i,"_library_size_est_",library_size_est_method,
                             ".Rdata",sep = ""))
#####################################################################
##--------------------------------------------------------------------------------------------