##scGeneNet_main
scGeneNet_main<-function(scGeneNet_list,lambda_use,
                     Global_Control = list(ELBO_threshold = 1e-4, minit = 5, maxit = 30,maxit_nonGRN = 3, MS_update_threshold = 1e-6),
                     M_Control = list(ADMM_max_step = 1000, ADMM_threshold = 1e-3),
                     S_Control = list(Newton_max_step = 1000, Newton_threshold = 1e-4),
                     Theta_Control = list(penalize_diagonal = FALSE,Theta_threshold = 1e-4),
                     U_fix = FALSE,
                     verbose = TRUE,
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
  

  
  ##---------------------------------------------------------
  
  
  ##Re-initial the precision matrix with the current lambda_use and penalize_diagonal = FALSE
  ##---------------------------------------------------------
  scGeneNet_list<-init_Theta(scGeneNet_list,lambda_use = lambda_use, penalize_diagonal = FALSE, Theta_threshold = Theta_threshold,core_num = min(ncol(scGeneNet_list$U_mat),core_num))
  ##---------------------------------------------------------
  
  ##
  scGeneNet_list<-scGeneNet_optim(celltypes_label_input = scGeneNet_list$celltypes_label,
              U_mat_input = scGeneNet_list$U_mat,
              mu_mat_input = scGeneNet_list$mu_mat,
              m_mat_list_input = scGeneNet_list$m_mat_list,
              s_sq_mat_list_input = scGeneNet_list$s_sq_mat_list,
              pi_vec_input = scGeneNet_list$pi_vec,
              obs_mat_input = scGeneNet_list$obs_mat,
              S_depth_input = scGeneNet_list$S_depth,
              gene_GRN_index_use_input = scGeneNet_list$gene_GRN_index_use,
              Theta_mat_list_input = scGeneNet_list$Theta_mat_list,
              zero_GRN_input = scGeneNet_list$zero_GRN,
              if_zero = scGeneNet_list$zero_GRN_use,
              log_det_vec_input = scGeneNet_list$log_det_vec,
              lambda_use = lambda_use,
              penalize_diagonal = penalize_diagonal,
              maxit_nonGRN = maxit_nonGRN,
              MS_update_threshold = MS_update_threshold,
              maxit = maxit,
              minit = minit,
              ELBO_threshold = ELBO_threshold,
              verbose = verbose,
              ADMM_max_step = ADMM_max_step,
              ADMM_threshold = ADMM_threshold,
              Newton_max_step = Newton_max_step,
              Newton_threshold = Newton_threshold,
              Theta_threshold = Theta_threshold,
              core_num = core_num,
              U_fix = U_fix)
  
  ##
  return(scGeneNet_list)
  
}
