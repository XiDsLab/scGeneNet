##exact_zeroindex
exact_zeroindex<-function(support = NULL){
  zeroindex_res<-NULL
  if(!is.null(support)){
    zeroindex_res<-matrix(NA,nrow = 0,ncol = 2)
    ##
    zero_mat<-ifelse(support==0,1,0)
    for(i in 1:(nrow(zero_mat) - 1)){
      zero_vec_index<-which(zero_mat[i,] == 1)
      zero_vec_index<-setdiff(zero_vec_index,(1:i))
      if(length(zero_vec_index)>0){
        zeroindex_res<-rbind(zeroindex_res,cbind(rep(i,length(zero_vec_index)),zero_vec_index))
      }
    }
  }
  ##
  return(zeroindex_res)
}
