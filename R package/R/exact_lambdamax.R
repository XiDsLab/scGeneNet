# exact_lambdamax
exact_lambdamax = function(scGeneNet_list,lambda_max_ref = NULL,
                      Global_Control = list(ELBO_threshold = 1e-4, minit = 5, maxit = 30,maxit_nonGRN = 3, MS_update_threshold = 1e-6),
                      M_Control = list(ADMM_max_step = 1000, ADMM_threshold = 1e-3),
                      S_Control = list(Newton_max_step = 1000, Newton_threshold = 1e-4),
                      Theta_Control = list(penalize_diagonal = FALSE,Theta_threshold = 1e-4),
                      U_fix = FALSE,
                      core_num = 20){
  ##Initial the hyper-parameter
  ##---------------------------------------------------------
  
  ##Global_Control
  #######################################
  ELBO_threshold<-ifelse(is.null(Global_Control$ELBO_threshold),1e-4,Global_Control$ELBO_threshold)
  minit<-ifelse(is.null(Global_Control$minit),5,Global_Control$minit)
  maxit<-ifelse(is.null(Global_Control$maxit),30,Global_Control$maxit)
  maxit_nonGRN<-ifelse(is.null(Global_Control$maxit_nonGRN),3,Global_Control$maxit_nonGRN)
  MS_update_threshold<-ifelse(is.null(Global_Control$MS_update_threshold),1e-6,Global_Control$MS_update_threshold)
  #######################################
  
  ##M_Control
  #######################################
  ADMM_max_step<-ifelse(is.null(M_Control$ADMM_max_step),1000,M_Control$ADMM_max_step)
  ADMM_threshold<-ifelse(is.null(M_Control$ADMM_threshold),1e-3,M_Control$ADMM_threshold)
  #######################################
  
  ##S_Control
  #######################################
  Newton_max_step<-ifelse(is.null(S_Control$Newton_max_step),1000,S_Control$Newton_max_step)
  Newton_threshold<-ifelse(is.null(S_Control$Newton_threshold),1e-4,S_Control$Newton_threshold)
  #######################################
  
  ##Theta_Control
  #######################################
  penalize_diagonal<-ifelse(is.null(Theta_Control$penalize_diagonal),FALSE,Theta_Control$penalize_diagonal)
  Theta_threshold<-ifelse(is.null(Theta_Control$Theta_threshold),1e-4,Theta_Control$Theta_threshold)
  #######################################
  
  #######################################
  
  n<-dim(scGeneNet_list$obs_mat)[1]
  celltypes_num<-ncol(scGeneNet_list$U_mat)
  p_GRN<-length(scGeneNet_list$gene_GRN)
  p_nonGRN<-dim(scGeneNet_list$obs_mat)[2] - length(scGeneNet_list$gene_GRN)
  
  ##
  
  if(is.null(lambda_max_ref)){
    base = (n / celltypes_num) / 3 
  }else{
    base = lambda_max_ref
  }
  ##Divide the group into two components
  lambda<-base
  print(paste("(",0,") lambda (base) = ",lambda,".",sep = ""))
  scGeneNet_list1<-init_Theta(scGeneNet_list,lambda_use = lambda,penalize_diagonal = FALSE,Theta_threshold = Theta_threshold, core_num = min(celltypes_num,core_num))
  
  ##
  scGeneNet_list1 <- scGeneNet_optim(celltypes_label_input = scGeneNet_list1$celltypes_label,
                             U_mat_input = scGeneNet_list1$U_mat,
                             mu_mat_input = scGeneNet_list1$mu_mat,
                             m_mat_list_input = scGeneNet_list1$m_mat_list,
                             s_sq_mat_list_input = scGeneNet_list1$s_sq_mat_list,
                             pi_vec_input = scGeneNet_list1$pi_vec,
                             obs_mat_input = scGeneNet_list1$obs_mat,
                             S_depth_input = scGeneNet_list1$S_depth,
                             gene_GRN_index_use_input = scGeneNet_list1$gene_GRN_index_use,
                             Theta_mat_list_input = scGeneNet_list1$Theta_mat_list,
                             zero_GRN_input = scGeneNet_list1$zero_GRN,
                             if_zero = scGeneNet_list1$zero_GRN_use,
                             log_det_vec_input = scGeneNet_list1$log_det_vec,
                             lambda_use = lambda,
                             penalize_diagonal = penalize_diagonal,
                             maxit_nonGRN = maxit_nonGRN,
                             MS_update_threshold = MS_update_threshold,
                             maxit = maxit,
                             minit = minit,
                             ELBO_threshold = ELBO_threshold,
                             verbose = FALSE,
                             ADMM_max_step = ADMM_max_step,
                             ADMM_threshold = ADMM_threshold,
                             Newton_max_step = Newton_max_step,
                             Newton_threshold = Newton_threshold,
                             Theta_threshold = Theta_threshold,
                             core_num = core_num,
                             U_fix = U_fix)
  ##
  
  Theta_nonzero_vec<-c()
  for(g in 1:ncol(scGeneNet_list1$U_mat)){
    Theta_mat_sub<-as.matrix(scGeneNet_list1$Theta_mat_list[[g]])[scGeneNet_list1$gene_GRN_index_use,scGeneNet_list1$gene_GRN_index_use]
    Theta_nonzero_vec<-c(Theta_nonzero_vec,length(which((Theta_mat_sub)[upper.tri(Theta_mat_sub)]!=0)))
    rm(Theta_mat_sub);gc()
  }
  group_less<-which(Theta_nonzero_vec==0)
  group_more<-which(Theta_nonzero_vec>0)
  ############
  lambda_max_vec_res<-rep(NA,(length(group_less) + length(group_more)))
  #####################
  ## For less group
  if(length(group_less)>0){
    print("Part of lambdas that are less than reference maximal lambda.")
    connect_vec<-rep(0,length(group_less))
    lambda_max_vec_less<-rep(NA,length(group_less))
    ###
    cnt = 0.8
    iter<-0
    con_0<-TRUE
    while(con_0){
      iter<-iter + 1
      lambda<-base * (cnt)^iter
      print(paste("(",iter,") lambda = ",lambda,".",sep = ""))
      ##
      scGeneNet_list1<-init_Theta(scGeneNet_list,lambda_use = lambda,penalize_diagonal = FALSE,Theta_threshold = Theta_threshold,
                              core_num = min(celltypes_num,core_num))
      
      ##
      scGeneNet_list1 <- scGeneNet_optim(celltypes_label_input = scGeneNet_list1$celltypes_label,
                                 U_mat_input = scGeneNet_list1$U_mat,
                                 mu_mat_input = scGeneNet_list1$mu_mat,
                                 m_mat_list_input = scGeneNet_list1$m_mat_list,
                                 s_sq_mat_list_input = scGeneNet_list1$s_sq_mat_list,
                                 pi_vec_input = scGeneNet_list1$pi_vec,
                                 obs_mat_input = scGeneNet_list1$obs_mat,
                                 S_depth_input = scGeneNet_list1$S_depth,
                                 gene_GRN_index_use_input = scGeneNet_list1$gene_GRN_index_use,
                                 Theta_mat_list_input = scGeneNet_list1$Theta_mat_list,
                                 zero_GRN_input = scGeneNet_list1$zero_GRN,
                                 if_zero = scGeneNet_list1$zero_GRN_use,
                                 log_det_vec_input = scGeneNet_list1$log_det_vec,
                                 lambda_use = lambda,
                                 penalize_diagonal = penalize_diagonal,
                                 maxit_nonGRN = maxit_nonGRN,
                                 MS_update_threshold = MS_update_threshold,
                                 maxit = maxit,
                                 minit = minit,
                                 ELBO_threshold = ELBO_threshold,
                                 verbose = FALSE,
                                 ADMM_max_step = ADMM_max_step,
                                 ADMM_threshold = ADMM_threshold,
                                 Newton_max_step = Newton_max_step,
                                 Newton_threshold = Newton_threshold,
                                 Theta_threshold = Theta_threshold,
                                 core_num = core_num,
                                 U_fix = U_fix)
      
      
      Theta_nonzero_vec<-c()
      for(g in 1:ncol(scGeneNet_list1$U_mat)){
        Theta_mat_sub<-as.matrix(scGeneNet_list1$Theta_mat_list[[g]])[scGeneNet_list1$gene_GRN_index_use,scGeneNet_list1$gene_GRN_index_use]
        Theta_nonzero_vec<-c(Theta_nonzero_vec,length(which((Theta_mat_sub)[upper.tri(Theta_mat_sub)]!=0)))
        rm(Theta_mat_sub);gc()
      }
      ##
      Theta_nonzero_vec_less<-Theta_nonzero_vec[group_less]
      for(g in 1:length(group_less)){
        if((connect_vec[g] == 0) && (Theta_nonzero_vec_less[g]>0)){
          lambda_max_vec_less[g]<-lambda/cnt
        }
      }
      ##
      connect_vec<-ifelse(Theta_nonzero_vec_less>0,1,0)
      ##
      if(sum(is.na(lambda_max_vec_less)) == 0){
        con_0<-FALSE
      }
      ##
    }
    ##
    lambda_max_vec_res[group_less]<-lambda_max_vec_less
  }
  ##
  ## For more group
  if(length(group_more)>0){
    print("Part of lambdas that are greater than reference maximal lambda.")
    connect_vec<-rep(1,length(group_more))
    lambda_max_vec_more<-rep(NA,length(group_more))
    ###
    cnt = 1.2
    iter<-0
    con_0<-TRUE
    while(con_0){
      iter<-iter + 1
      lambda<-base * (cnt)^iter
      print(paste("(",iter,") lambda = ",lambda,".",sep = ""))
      ##
      scGeneNet_list1<-init_Theta(scGeneNet_list,lambda_use = lambda,penalize_diagonal = FALSE,Theta_threshold = Theta_threshold,
                              core_num = min(celltypes_num,core_num))
      
      ##
      scGeneNet_list1 <- scGeneNet_optim(celltypes_label_input = scGeneNet_list1$celltypes_label,
                                 U_mat_input = scGeneNet_list1$U_mat,
                                 mu_mat_input = scGeneNet_list1$mu_mat,
                                 m_mat_list_input = scGeneNet_list1$m_mat_list,
                                 s_sq_mat_list_input = scGeneNet_list1$s_sq_mat_list,
                                 pi_vec_input = scGeneNet_list1$pi_vec,
                                 obs_mat_input = scGeneNet_list1$obs_mat,
                                 S_depth_input = scGeneNet_list1$S_depth,
                                 gene_GRN_index_use_input = scGeneNet_list1$gene_GRN_index_use,
                                 Theta_mat_list_input = scGeneNet_list1$Theta_mat_list,
                                 zero_GRN_input = scGeneNet_list1$zero_GRN,
                                 if_zero = scGeneNet_list1$zero_GRN_use,
                                 log_det_vec_input = scGeneNet_list1$log_det_vec,
                                 lambda_use = lambda,
                                 penalize_diagonal = penalize_diagonal,
                                 maxit_nonGRN = maxit_nonGRN,
                                 MS_update_threshold = MS_update_threshold,
                                 maxit = maxit,
                                 minit = minit,
                                 ELBO_threshold = ELBO_threshold,
                                 verbose = FALSE,
                                 ADMM_max_step = ADMM_max_step,
                                 ADMM_threshold = ADMM_threshold,
                                 Newton_max_step = Newton_max_step,
                                 Newton_threshold = Newton_threshold,
                                 Theta_threshold = Theta_threshold,
                                 core_num = core_num,
                                 U_fix = U_fix)
      
      
      Theta_nonzero_vec<-c()
      for(g in 1:ncol(scGeneNet_list1$U_mat)){
        Theta_mat_sub<-as.matrix(scGeneNet_list1$Theta_mat_list[[g]])[scGeneNet_list1$gene_GRN_index_use,scGeneNet_list1$gene_GRN_index_use]
        Theta_nonzero_vec<-c(Theta_nonzero_vec,length(which((Theta_mat_sub)[upper.tri(Theta_mat_sub)]!=0)))
        rm(Theta_mat_sub);gc()
      }
      ##
      Theta_nonzero_vec_more<-Theta_nonzero_vec[group_more]
      for(g in 1:length(group_more)){
        if((connect_vec[g] == 1) && (Theta_nonzero_vec_more[g]==0)){
          lambda_max_vec_more[g]<-lambda
        }
      }
      ##
      connect_vec<-ifelse(Theta_nonzero_vec_more>0,1,0)
      ##
      if(sum(is.na(lambda_max_vec_more)) == 0){
        con_0<-FALSE
      }
      ##
    }
    ##
    lambda_max_vec_res[group_more]<-lambda_max_vec_more
  }
  
  return(lambda_max_vec_res)
}

