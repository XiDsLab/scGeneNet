# init_Theta
init_Theta<-function (scGeneNet_list,lambda_use,penalize_diagonal,Theta_threshold,core_num) 
{
  ######
  celltypes_num <- ncol(scGeneNet_list$U_mat)
  ###
  ##
  scGeneNet_list$Theta_mat_list<-Theta_fun(U_mat = scGeneNet_list$U_mat,
                                     Theta_mat_list = lapply(X = 1:celltypes_num,FUN = function(g){return(as(diag(ncol(scGeneNet_list$mu_mat)),"sparseMatrix"))}),
                                     m_mat_list = scGeneNet_list$m_mat_list,
                                     s_sq_mat_list = scGeneNet_list$s_sq_mat_list,
                                     p_GRN = length(scGeneNet_list$gene_GRN_index_use),
                                     lambda_use = lambda_use,
                                     zero_GRN = scGeneNet_list$zero_GRN,
                                     stop_threshold = Theta_threshold,
                                     var_index_Theta = rep(1,celltypes_num),
                                     penalize_diagonal = penalize_diagonal,
                                     zero_GRN_use = scGeneNet_list$zero_GRN_use,
                                     core_num = 1)
  ######
  ##
  log_det_vec<-NULL
  if((ncol(scGeneNet_list$obs_mat) - length(scGeneNet_list$gene_GRN_index_use)) == 0){
    log_det_vec <- log_det(matrix_list = scGeneNet_list$Theta_mat_list,celltypes_num = celltypes_num,core_num = min(core_num,celltypes_num))
  }else{
    log_det_vec <- log_det_block(matrix_list = scGeneNet_list$Theta_mat_list,celltypes_num = celltypes_num,p_GRN = length(scGeneNet_list$gene_GRN_index_use),core_num = min(core_num,celltypes_num))
  }
  ######################################################################################
  scGeneNet_list[["log_det_vec"]]<-log_det_vec
  
  ######################################################################################
  
  return(scGeneNet_list)
}

