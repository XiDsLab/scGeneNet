// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <Eigen/Sparse>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXf;
using std::cout;
using std::endl;
static omp_lock_t lock;
//
//
//
const double epsilon1 = 1e-50;
const double epsilon2 = 1e-9;
const double epsilon3 = 1e-5;
//


Eigen::MatrixXd trun_fun(Eigen::MatrixXd mat_use,double trun_value){
  int n = mat_use.rows();
  int p = mat_use.cols();
  Eigen::MatrixXd mat_res = mat_use;
  for(int i =0;i<n;i++){
    for(int j =0;j<p;j++){
      if(mat_res(i,j)>trun_value){
        mat_res(i,j) = trun_value;
      }
    }
  }
  //
  return mat_res;
}


Eigen::VectorXd trun_vec_fun(Eigen::VectorXd vec_use,double trun_value){
  int l = vec_use.size();
  Eigen::VectorXd vec_res = vec_use;
  for(int i = 0;i<l;i++){
    if(vec_res(i)>trun_value){
      vec_res(i) = trun_value;
    }
  }
  //
  return vec_res;
}

Eigen::SparseMatrix<double> log_cul_mat(Eigen::SparseMatrix<int> mat){
  int num_row = mat.rows();
  int num_col = mat.cols();
  Eigen::SparseMatrix<double> mat_res(num_row,num_col);
  
  for(int k=0;k<mat.outerSize();++k){
    for(SparseMatrix<int>::InnerIterator it(mat,k);it;++it){
      Eigen::VectorXd vec_cur;
      int it_value = it.value();
      vec_cur.setZero(it_value);
      for(int j = 0;j<it_value;j++){
        vec_cur[j] = j+1;
      }
      
      mat_res.insert(it.row(),it.col()) = vec_cur.array().log().sum();
    }
  }
  //
  mat_res.makeCompressed();
  
  return mat_res;
}

Eigen::VectorXd rowSums_sparse(Eigen::SparseMatrix<double> mat){
  Eigen::MatrixXd mat_densy = Eigen::MatrixXd(mat);
  Eigen::VectorXd rowSums_res = mat_densy.rowwise().sum();
  return rowSums_res;
}

Eigen::VectorXd colSums_sparse(Eigen::SparseMatrix<double> mat){
  Eigen::MatrixXd mat_densy = Eigen::MatrixXd(mat);
  Eigen::VectorXd colSums_res = mat_densy.colwise().sum();
  return colSums_res;
}

Eigen::MatrixXd stack_fun(Rcpp::List list_use, int celltypes_num, int row, int col){
  Eigen::MatrixXd stack_mat(row,celltypes_num * col);
  // int count = 0;
  for(int g = 0;g<celltypes_num;g++){
    Eigen::MatrixXd mat_use = list_use[g];
    stack_mat.middleCols(g * col,col) = mat_use;
  }
  //
  return stack_mat;
}

Eigen::MatrixXd stack_sparse_fun(Rcpp::List list_use, int celltypes_num, int row, int col){
  Eigen::MatrixXd stack_mat(row,celltypes_num * col);
  // int count = 0;
  for(int g = 0;g<celltypes_num;g++){
    Eigen::SparseMatrix<double> mat_sparse_use = list_use[g];
    Eigen::MatrixXd mat_use = MatrixXd(mat_sparse_use);
    stack_mat.middleCols(g * col,col) = mat_use;
  }
  //
  return stack_mat;
}


Eigen::VectorXd ELBO_fast(Rcpp::List Theta_mat_list,Eigen::MatrixXd U_mat,Rcpp::List m_mat_list, Rcpp::List s_sq_mat_list, Eigen::SparseMatrix<int> obs_mat,
                          Eigen::MatrixXd mu_mat, Eigen::VectorXd S_depth, Eigen::VectorXd pi_vec,
                          int p_choose, Eigen::VectorXd log_cul_sum,
                          Eigen::VectorXd log_det_vec,int core_num) {
  int celltypes_num = U_mat.cols();
  int n = U_mat.rows();
  // int p = mu_mat.cols();
  int p = p_choose;
  //
  //reshape the List object
  Eigen::MatrixXd Theta_mat_stackall = stack_sparse_fun(Theta_mat_list,celltypes_num,mu_mat.cols(),mu_mat.cols());
  Eigen::MatrixXd m_mat_stackall = stack_fun(m_mat_list,celltypes_num,n,mu_mat.cols());
  Eigen::MatrixXd s_sq_mat_stackall = stack_fun(s_sq_mat_list,celltypes_num,n,mu_mat.cols());
  
  Eigen::VectorXd S_depth_log = S_depth.array().log();
  //
  Eigen::VectorXd ELBO_vec;
  //
  ELBO_vec.setZero(celltypes_num);
  // Eigen::VectorXd log_cul_mat_obs_sum = MatrixXd(log_cul_mat_obs_mat).rowwise().sum();
  Eigen::MatrixXd mu_mat_GRN = mu_mat.leftCols(p);
  Eigen::SparseMatrix<int> obs_mat_GRN = obs_mat.leftCols(p);
  Eigen::MatrixXi obs_mat_GRN_densy = MatrixXi(obs_mat_GRN);
  //
  omp_init_lock(&lock);
  omp_set_lock(&lock);
  //
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int g = 0;g<celltypes_num;g++){
    // Eigen::SparseMatrix<double> Theta_mat_sparse = Theta_mat_list[g];
    // Eigen::MatrixXd Theta_mat0 = MatrixXd(Theta_mat_sparse);
    Eigen::MatrixXd Theta_mat0 = Theta_mat_stackall.middleCols(g*mu_mat.cols(),mu_mat.cols());
    Eigen::MatrixXd Theta_mat = Theta_mat0.block(0,0,p,p);
    Eigen::VectorXd Theta_diag = Theta_mat.diagonal();
    //
    double ELBO_group = 0;
    // Eigen::MatrixXd mat10 = as<Eigen::MatrixXd>(m_mat_list[g]);
    Eigen::MatrixXd mat10 = m_mat_stackall.middleCols(g*mu_mat.cols(),mu_mat.cols()); 
    Eigen::MatrixXd mat1 = mat10.leftCols(p);
    for(int j = 0;j<p;j++){
      mat1.col(j) = mat1.col(j).array() + S_depth_log.array();
    }
    for(int i = 0;i<n;i++){
      mat1.row(i) = mat1.row(i).array() + mu_mat_GRN.row(g).array();
    }
    // Eigen::MatrixXd A_mat = (mat1.array() + ((as<Eigen::MatrixXd>(s_sq_mat_list[g])).array() * 0.5).array()).array().exp();
    
    // Eigen::MatrixXd B_mat = ((MatrixXd(obs_mat)).array() * mat1.array()).array() - (MatrixXd(log_cul_mat_obs_mat)).array() - A_mat.array() + ((0.5) * ((as<Eigen::MatrixXd>(s_sq_mat_list[g])).array().log()).array()).array();
    // Eigen::MatrixXd s_sq_mat = as<Eigen::MatrixXd>(s_sq_mat_list[g]);
    Eigen::MatrixXd s_sq_mat = s_sq_mat_stackall.middleCols(g*mu_mat.cols(),mu_mat.cols());
    Eigen::MatrixXd s_sq_mat_GRN = s_sq_mat.leftCols(p);
    Eigen::MatrixXd obs_mat_GRN_double = obs_mat_GRN_densy.cast<double>();
    Eigen::MatrixXd mat_for_trun1 = (mat1.array() + ((s_sq_mat_GRN).array() * 0.5).array());
    mat_for_trun1 = trun_fun(mat_for_trun1,50);
    Eigen::MatrixXd B_mat = ((obs_mat_GRN_double).array() * mat1.array()).array() - (mat_for_trun1.array().exp()).array() + ((0.5) * ((s_sq_mat_GRN).array().log()).array()).array();
    //
    Eigen::VectorXd B_vec = B_mat.rowwise().sum();
    //
    double entropy_group = 0;
    // Eigen::MatrixXd m_mat1 = as<Eigen::MatrixXd>(m_mat_list[g]);
    Eigen::MatrixXd m_mat1 = m_mat_stackall.middleCols(g*mu_mat.cols(),mu_mat.cols());
    Eigen::MatrixXd m_mat_GRN = m_mat1.leftCols(p);
    Eigen::MatrixXd m_mat_GRN1 = m_mat1.leftCols(p);
    Eigen::VectorXd s_sq_sum;
    s_sq_sum.setZero(p);
    // cout<<"p:"<<p<<endl;
    // cout<<s_sq_sum<<endl;
    // double log_cur_sum = 0;
    for(int i = 0;i<n;i++){
      // ELBO_group = ELBO_group + ((B_mat.row(i)).array() * U_mat(i,g)).sum();
      //
      if(U_mat(i,g) > epsilon1){
        entropy_group = entropy_group + U_mat(i,g) * log(U_mat(i,g));
      }
      //
      m_mat_GRN.row(i) = m_mat_GRN.row(i).array() * U_mat(i,g);
      //
      Eigen::VectorXd s_sq_vec = (s_sq_mat_GRN.row(i)).array() * U_mat(i,g);
      s_sq_sum = s_sq_sum.array() + s_sq_vec.array();
      //
      // log_cur_sum = log_cur_sum + log_cul_mat_obs_sum(i) * U_mat(i,g);
    }
    ELBO_group = ELBO_group + (B_vec.array() * (U_mat.col(g)).array()).sum();
    //
    // cout<<s_sq_sum<<endl;
    // //
    // ELBO_group = ELBO_group + ((p + log_det_vec(g)) * 0.5 + log(pi_vec(g))) * (U_mat.col(g).sum()) - entropy_group - log_cur_sum;
    ELBO_group = ELBO_group + ((p + log_det_vec(g)) * 0.5 + log(pi_vec(g))) * (U_mat.col(g).sum()) - entropy_group - (log_cul_sum.array() * (U_mat.col(g)).array()).sum();
    // //
    //
    // Theta_mat = (Theta_mat.array() + (Theta_mat.adjoint()).array())*0.5;
    double value1 = ((m_mat_GRN * Theta_mat).array() * (m_mat_GRN1).array()).sum() * 0.5;
    ELBO_group = ELBO_group - value1;
    // //
    double value2 =  0.5 * ((s_sq_sum.array() * Theta_diag.array()).sum());
    ELBO_group = ELBO_group - value2;
    ///////
    ELBO_vec(g) = ELBO_group;
  }
  omp_destroy_lock(&lock);
  //
  return(ELBO_vec);
}




Eigen::MatrixXd clustering_fun(Rcpp::List Theta_mat_list,Eigen::MatrixXd U_mat, 
                               Rcpp::List m_mat_list,Rcpp::List s_sq_mat_list, 
                               Eigen::SparseMatrix<int> obs_mat,
                               Eigen::MatrixXd mu_mat,
                               Eigen::VectorXd S_depth,
                               Eigen::VectorXd pi_vec,
                               Eigen::VectorXd log_det_vec){
  
  //////////////////////////////
  int celltypes_num = U_mat.cols();
  int n = U_mat.rows();
  int p = obs_mat.cols();
  Eigen::VectorXd S_depth_log = S_depth.array().log();
  Eigen::MatrixXi obs_mat_densy = MatrixXi(obs_mat);
  Eigen::MatrixXd obs_mat_densy_double = obs_mat_densy.cast <double> ();
  //
  log_det_vec = trun_vec_fun(log_det_vec,1e+8);
  //
  Eigen::MatrixXd U_mat_res(n,celltypes_num);
  Eigen::VectorXi celltypes_label_res(n);
  //
  Eigen::MatrixXd C_mat(n,celltypes_num);
  for(int g=0;g<celltypes_num;g++){
    //
    Eigen::SparseMatrix<double> Theta_mat0 = Theta_mat_list[g];
    Eigen::MatrixXd Theta_mat = MatrixXd(Theta_mat0);
    //
    Theta_mat = trun_fun(Theta_mat,1e+8);
    
    //
    
    Eigen::VectorXd Theta_diag = Theta_mat.diagonal();
    Eigen::MatrixXd m_mat = m_mat_list[g];
    Eigen::MatrixXd s_sq_mat = s_sq_mat_list[g];
    Eigen::VectorXd mu_vec = mu_mat.row(g);
    //
    Eigen::MatrixXd m_mat1 = m_mat;
    for(int j = 0;j<p;j++){
      Eigen::VectorXd m_vec = m_mat1.col(j);
      m_mat1.col(j) = m_vec.array() + S_depth_log.array();
    }
    for(int i = 0;i<n;i++){
      Eigen::VectorXd m_vec = m_mat1.row(i);
      m_mat1.row(i) = m_vec.array() + mu_vec.array();
    }
    Eigen::MatrixXd A_mat0 = (m_mat1.array() + (s_sq_mat.array() * 0.5).array());
    Eigen::MatrixXd A_mat = trun_fun(A_mat0,50);
    A_mat = A_mat.array().exp();
    //
    Eigen::VectorXd term1_vec = ((obs_mat_densy_double.array() * m_mat1.array()).array() - A_mat.array() + (0.5 * (s_sq_mat.array().log()).array()).array()).rowwise().sum();
    //
    Eigen::VectorXd term2_vec = (((m_mat * Theta_mat).array() * m_mat.array()).array() * (-0.5)).rowwise().sum();
    Eigen::VectorXd term3_vec = s_sq_mat * Theta_diag;
    term2_vec = term2_vec.array() - ((0.5) * term3_vec).array() + 0.5 * log_det_vec(g) + log(pi_vec(g));
    //
    C_mat.col(g) = term1_vec.array() + term2_vec.array();
    //
    
  }
  
  for(int i =0;i<n;i++){
    C_mat.row(i) = (C_mat.row(i)).array() - (C_mat.row(i)).maxCoeff();
    C_mat.row(i) = (C_mat.row(i)).array().exp();
    U_mat_res.row(i) = (C_mat.row(i)).array()/((C_mat.row(i)).sum());
  }
  
  return U_mat_res;
}



Eigen::VectorXi celltypes_label_fun(Eigen::MatrixXd U_mat){
  int n = U_mat.rows();
  int celltypes_num = U_mat.cols();
  Eigen::VectorXi celltypes_label_res(n);
  //
  for(int i =0;i<n;i++){
    Eigen::VectorXd U_vec = U_mat.row(i);
    Eigen::VectorXd::Index maxindex;
    double max = U_vec.maxCoeff(&maxindex);
    celltypes_label_res(i) = maxindex;
  }
  //
  return celltypes_label_res;
}


Eigen::VectorXi cluster_predict_true_fun(Eigen::VectorXi celltypes_label,Eigen::VectorXi celltypes_label_start,int celltypes_num){
  //
  int n = celltypes_label.size();
  //
  Eigen::MatrixXi table_mat;
  table_mat.setZero(celltypes_num,celltypes_num);
  for(int i =0;i<n;i++){
    //
    int row_index = celltypes_label(i);
    int col_index = celltypes_label_start(i);
    table_mat(row_index,col_index) = table_mat(row_index,col_index)+1;
  }
  //
  Eigen::VectorXi cluster_predict_true_res(celltypes_num);
  for(int g = 0;g<celltypes_num;g++){
    Eigen::VectorXi table_vec = table_mat.row(g);
    Eigen::VectorXd::Index maxindex;
    double max = table_vec.maxCoeff(&maxindex);
    cluster_predict_true_res(g) = maxindex;
  }
  //
  return cluster_predict_true_res;
}


int count_group_fun(Eigen::VectorXi celltypes_label){
  int count = 0;
  int group_max  = celltypes_label.size();
  for(int g = 0;g<group_max;g++){
    int count_0 = 0;
    for(int g1 = 0;g1<group_max;g1++){
      if(celltypes_label(g1) == g){
        count_0 = 1;
      }
    }
    count = count + count_0;
  }
  //
  return count;
}



double purity_fun(Eigen::VectorXi celltypes_label,Eigen::VectorXi celltypes_label_start,Eigen::VectorXi cluster_predict_true,int celltypes_num){
  //
  double count = 0;
  int n = celltypes_label.size();
  for(int g = 0;g<celltypes_num;g++){
    for(int i=0;i<n;i++){
      if((celltypes_label(i) == g) && (celltypes_label_start(i) == cluster_predict_true(g))){
        count = count + 1;
      }
    }
  }
  double purity_res = count/n;
  //
  return purity_res;
}


Eigen::VectorXd MS_update_ratio_fun(Eigen::MatrixXi MS_update_index){
  int celltypes_num = MS_update_index.cols();
  int n = MS_update_index.rows();
  Eigen::MatrixXi MS_update_index_densy0 = MS_update_index;
  Eigen::MatrixXd MS_update_index_densy = MS_update_index_densy0.cast<double>();
  Eigen::VectorXd MS_update_ratio_res = MS_update_index_densy.colwise().sum();
  MS_update_ratio_res = MS_update_ratio_res/n;
  // for(int g=0;g<celltypes_num;g++){
  //   MS_update_index_densy.col
  //   }
  //
  return MS_update_ratio_res;
}

// [[Rcpp::export]]
Eigen::VectorXd log_det(Rcpp::List matrix_list,int celltypes_num,
                        int core_num = 1){
  Eigen::VectorXd log_det_vec;
  log_det_vec.setZero(celltypes_num);
  //stack
  Eigen::SparseMatrix<double> mat_in_list0 = matrix_list[0];
  
  int row_num = mat_in_list0.rows();
  int col_num = mat_in_list0.cols();
  //
  Eigen::MatrixXd matrix_list_stackall = stack_sparse_fun(matrix_list,celltypes_num,row_num,col_num);
  //
  omp_init_lock(&lock);
  omp_set_lock(&lock);
  
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int i=0;i<celltypes_num;i++){
    // Eigen::SparseMatrix<double> matrix_use0 = matrix_list[i];
    // Eigen::MatrixXd matrix_use = MatrixXd(matrix_use0);
    Eigen::MatrixXd matrix_use = matrix_list_stackall.middleCols(i*col_num,col_num);
    //
    int p = matrix_use.cols();
    //
    EigenSolver<MatrixXd> eig(matrix_use);
    Eigen::VectorXd ev = eig.eigenvalues().real();
    double log_det_res = 0;
    for(int j=0;j<p;j++){
      if((ev(j)<0) || ((ev(j)==0))){
        log_det_res = log_det_res + log(epsilon2);
      }else{
        log_det_res = log_det_res + log(ev(j));
      }
    }
    //
    log_det_vec(i) = log_det_res;
    //
  }
  omp_destroy_lock(&lock);
  return log_det_vec;
}


// [[Rcpp::export]]
Eigen::VectorXd log_det_block(Rcpp::List matrix_list,int celltypes_num,
                              int p_GRN,
                              int core_num = 1){
  Eigen::VectorXd log_det_vec;
  log_det_vec.setZero(celltypes_num);
  //stack
  Eigen::SparseMatrix<double> mat_in_list0 = matrix_list[0];
  
  int row_num = mat_in_list0.rows();
  int col_num = mat_in_list0.cols();
  //
  Eigen::MatrixXd matrix_list_stackall = stack_sparse_fun(matrix_list,celltypes_num,row_num,col_num);
  //
  omp_init_lock(&lock);
  omp_set_lock(&lock);
  
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int i=0;i<celltypes_num;i++){
    // Eigen::SparseMatrix<double> matrix_use0 = matrix_list[i];
    // Eigen::MatrixXd matrix_use = MatrixXd(matrix_use0);
    Eigen::MatrixXd matrix_use = matrix_list_stackall.middleCols(i*col_num,col_num);
    //
    int p = matrix_use.cols();
    int p_nonGRN = p - p_GRN;
    //component1
    Eigen::MatrixXd matrix_use1 = matrix_use.block(0,0,p_GRN,p_GRN);
    Eigen::MatrixXd matrix_use2 = matrix_use.block(p_GRN,p_GRN,p_nonGRN,p_nonGRN);
    //
    EigenSolver<MatrixXd> eig1(matrix_use1);
    Eigen::VectorXd ev1 = eig1.eigenvalues().real();
    double log_det_res1 = 0;
    for(int j=0;j<p_GRN;j++){
      if((ev1(j)<0) || ((ev1(j)==0))){
        log_det_res1 = log_det_res1 + log(epsilon2);
      }else{
        log_det_res1 = log_det_res1 + log(ev1(j));
      }
    }
    //
    //
    EigenSolver<MatrixXd> eig2(matrix_use2);
    Eigen::VectorXd ev2 = eig2.eigenvalues().real();
    double log_det_res2 = 0;
    for(int j=0;j<p_nonGRN;j++){
      if((ev2(j)<0) || ((ev2(j)==0))){
        log_det_res2 = log_det_res2 + log(epsilon2);
      }else{
        log_det_res2 = log_det_res2 + log(ev2(j));
      }
    }
    
    //
    
    log_det_vec(i) = log_det_res1 + log_det_res2;
    //
  }
  omp_destroy_lock(&lock);
  return log_det_vec;
}


// [[Rcpp::export]]
Eigen::VectorXd log_det_GRN(Rcpp::List matrix_list,int celltypes_num,
                            int p_GRN,
                            int core_num = 1){
  Eigen::VectorXd log_det_vec;
  log_det_vec.setZero(celltypes_num);
  //stack
  Eigen::SparseMatrix<double> mat_in_list0 = matrix_list[0];
  
  int row_num = mat_in_list0.rows();
  int col_num = mat_in_list0.cols();
  //
  Eigen::MatrixXd matrix_list_stackall = stack_sparse_fun(matrix_list,celltypes_num,row_num,col_num);
  //
  omp_init_lock(&lock);
  omp_set_lock(&lock);
  
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int i=0;i<celltypes_num;i++){
    // Eigen::SparseMatrix<double> matrix_use0 = matrix_list[i];
    // Eigen::MatrixXd matrix_use = MatrixXd(matrix_use0);
    Eigen::MatrixXd matrix_use = matrix_list_stackall.middleCols(i*col_num,col_num);
    //
    int p = matrix_use.cols();
    int p_nonGRN = p - p_GRN;
    //component1
    Eigen::MatrixXd matrix_use1 = matrix_use.block(0,0,p_GRN,p_GRN);
    //
    EigenSolver<MatrixXd> eig1(matrix_use1);
    Eigen::VectorXd ev1 = eig1.eigenvalues().real();
    double log_det_res1 = 0;
    for(int j=0;j<p_GRN;j++){
      if((ev1(j)<0) || ((ev1(j)==0))){
        log_det_res1 = log_det_res1 + log(epsilon2);
      }else{
        log_det_res1 = log_det_res1 + log(ev1(j));
      }
    }
    //
    
    //
    
    log_det_vec(i) = log_det_res1;
    //
  }
  omp_destroy_lock(&lock);
  return log_det_vec;
}

// [[Rcpp::export]]
Eigen::MatrixXd ADMM1(Eigen::SparseMatrix<double> Theta_mat_sparse, Eigen::MatrixXd s_sq_mat_all, 
                      Eigen::SparseMatrix<int> obs_mat, Eigen::VectorXd mu_vec_all, 
                      Eigen::MatrixXd m_mat_start_all, Eigen::VectorXd S_depth, Eigen::VectorXi var_index_m, int GRN_num,
                      int core_num = 1, 
                      int ADMM_maxiter = 50,double th = 1e-3,
                      double rho_input = 0.1, 
                      bool is_rho_input = true,
                      
                      bool is_rho_constant_nonGRN = true,
                      bool is_update_nonGRNpart = true){
  //
  Eigen::MatrixXi obs_mat_all_int = MatrixXi(obs_mat);
  Eigen::MatrixXd obs_mat_all = obs_mat_all_int.cast<double>();
  Eigen::MatrixXd Theta_mat = MatrixXd(Theta_mat_sparse);
  //
  int n = obs_mat_all.rows();
  int p = m_mat_start_all.cols();
  int non_GRN_num = p - GRN_num;
  // int celltypes_num = U_mat.cols();
  double rho = 0.1 * double(p); 
  // Eigen::MatrixXd inv_share;
  Eigen::MatrixXd identity_mat;
  identity_mat.setIdentity(p,p);
  
  //GRN
  Eigen::MatrixXd m_mat_start = m_mat_start_all.leftCols(GRN_num);
  Eigen::MatrixXd s_sq_mat = s_sq_mat_all.leftCols(GRN_num);
  Eigen::MatrixXd obs_mat1 = obs_mat_all.leftCols(GRN_num);
  Eigen::VectorXd mu_vec = mu_vec_all.head(GRN_num);
  //
  
  // Eigen::MatrixXd Theta_struct;
  // inv_share.setZero(p,p);
  // choose rho
  Eigen::MatrixXd Theta = Theta_mat.block(0,0,GRN_num,GRN_num);
  // Q_sum.setZero(p,p);
  // for(int g = 0;g < celltypes_num; g++){
  //   Q_sum  += double((1.0/celltypes_num)) * as<MatrixXd>(Theta_mat_array[g]);
  //   Theta_struct[g] = as<MatrixXd>(Theta_mat_array[g]);
  // }
  if(is_rho_input){
    rho = rho_input;
  } else {
    EigenSolver<MatrixXd> eig(Theta);
    //
    Eigen::VectorXd ev = eig.eigenvalues().real();
    double l_ev = ev(0);
    double s_ev = ev(GRN_num-1);
    // double delta = double(S_depth.sum()/n);
    // delta *= double((s_sq_mat/2 + m_mat_start).array().exp().sum()/(n*p));
    double delta = 0.0;
    int count = 0;
    for(int i=0;i<n;i++){
      if(var_index_m(i) == 1){
        Eigen::VectorXd s_sq_vec = s_sq_mat.row(i)/2;
        Eigen::VectorXd m_cur = m_mat_start.row(i);
        Eigen::VectorXd sum_vec = s_sq_vec.array() + m_cur.array() + mu_vec.array();
        sum_vec = trun_vec_fun(sum_vec,50);
        delta = delta + S_depth(i) * ((sum_vec).array().exp().sum());
        count = count + 1;
      }
      // Eigen::VectorXd s_sq_vec = s_sq_mat.row(i_index)/2;
      // Eigen::VectorXd m_cur = m_mat_start.row(i_index);
      // delta = delta + S_depth(i_index) * ((s_sq_vec.array() + m_cur.array()
      //                                        + mu_vec.array()).array().exp().sum());
    }
    delta /= count;
    if(delta < s_ev) rho = sqrt(delta * s_ev+0.001);
    else if (delta > l_ev) rho = sqrt(delta * l_ev+0.001);
    else rho = delta;
    //
    if((rho<0.1)) rho = 0.1; 
  }
  
  //
  // EigenSolver<MatrixXd> eig(Theta);
  // Eigen::VectorXd ev = eig.eigenvalues().real();
  // double l_ev = ev(0);
  // double s_ev = ev(p-1);
  // // double delta = double(S_depth.sum()/n);
  // // delta *= double((s_sq_mat/2 + m_mat_start).array().exp().sum()/(n*p));
  // double delta = 0.0;
  // for(int i=0;i<var_index_m.size();i++){
  //   int i_index = var_index_m(i) - 1;
  //   Eigen::VectorXd s_sq_vec = s_sq_mat.row(i_index)/2;
  //   Eigen::VectorXd m_cur = m_mat_start.row(i_index);
  //   delta = delta + S_depth(i_index) * ((s_sq_vec.array() + m_cur.array()
  //                                          + mu_vec.array()).array().exp().sum());
  // }
  // delta /= var_index_m.size();
  // if(delta < s_ev) rho = sqrt(delta * s_ev+0.001);
  // else if (delta > l_ev) rho = sqrt(delta * l_ev+0.001);
  // else rho = delta;
  // //
  // if(rho_con & (rho<0.1)) rho = 0.1; 
  
  // share inverse
  
  //   omp_set_num_threads(core_num);
  // #pragma omp parallel for
  // for(int g = 0;g < celltypes_num; g++){
  //   Eigen::MatrixXd tmp = rho * identity_mat + Theta_struct[g];
  //   //VectorXi member = get_com(tmp);
  //   //inv_share[g] = block_inverse(tmp, member, p);
  //   inv_share[g] = tmp.inverse();
  // }
  //  block matrix inverse
  // Eigen::MatrixXd tmp = rho * identity_mat.block(0,0,GRN_num,GRN_num) + Theta.block(0,0,GRN_num,GRN_num);
  Eigen::MatrixXd tmp = (rho * (identity_mat.block(0,0,GRN_num,GRN_num)).array()).array() + Theta.array();
  // inv_share.block(0,0,GRN_num,GRN_num) = tmp.inverse();
  // Eigen::MatrixXd inv_A = inv_share.block(0,0,GRN_num,GRN_num);
  Eigen::MatrixXd inv_A = tmp.inverse();
  
  Eigen::MatrixXd m_mat_all = m_mat_start_all;
  Eigen::MatrixXd m_mat = m_mat_start;
  // m_mat.setZero(n,p);
  Eigen::MatrixXd n_mat = m_mat_start;
  // n_mat.setZero(n,p);
  int conver_num = 0;
  // ADMM
  omp_init_lock(&lock);
  omp_set_lock(&lock);
  
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int i=0;i<n;i++){
    // cout<<"i = "<<i<<endl;
    if(var_index_m(i) == 1){
      // //
      int ADMM_step = ADMM_maxiter;
      // int id_share = -1;
      // for(int g=0;g<celltypes_num;g++){
      //   if(abs(U_mat(i_index,g) -1) < 0.001){
      //     id_share = g;
      //   }
      // }
      // if(id_share == -1){
      //   Eigen::MatrixXd Q;
      //   Q.setZero(p,p);
      //   // P.setZero(p,p);
      //
      //   for(int g=0; g<celltypes_num;g++){
      //     Q += U_mat(i_index,g) * Theta_struct[g];
      //     // P += U_mat(i,g) * Theta_struct[g];
      //   }
      //
      //
      //   Q +=  rho * identity_mat;
      //   //VectorXi member = get_com(Q);
      //   //inv_A = block_inverse(Q, member, p);
      //
      //   inv_A = Q.inverse();
      // } else {
      //   inv_A = inv_share[id_share];
      // }
      // Eigen::VectorXd mu_sum;
      // mu_sum.setZero(p);
      // for(int g=0; g<celltypes_num;g++){
      // Eigen::VectorXd mu_vec = mu_mat.row(g);
      // mu_sum = mu_sum + U_mat(i_index,g) * (Theta_struct[g] * (mu_vec));
      // }
      //
      Eigen::VectorXd m_cur, n_cur, alpha;
      alpha.setZero(p);
      m_cur = m_mat_start.row(i);
      n_cur = m_mat_start.row(i);
      Eigen::VectorXd s_sq_vec = s_sq_mat.row(i)/2;
      Eigen::VectorXd obs_vec = obs_mat1.row(i);
      //
      double elbo_pre = 0;
      double lib = S_depth(i);
      Eigen::VectorXd sum_vec1 = mu_vec.array() + s_sq_vec.array()+ m_cur.array();
      sum_vec1 = trun_vec_fun(sum_vec1,50);
      elbo_pre = - (m_cur.array() * obs_vec.array()).sum() +
        lib * sum_vec1.array().exp().sum();
      elbo_pre = elbo_pre + 0.5 * (((Theta * (m_cur)).array() * (m_cur).array()).sum());
      // for(int i=GRN_num;i<p;i++){
      //   elbo_pre += 0.5 * Theta(i,i) * m_cur(i) * m_cur(i);
      // }
      
      while(ADMM_step>0){
        ADMM_step -= 1;
        int Newton_step = 10;
        Eigen::VectorXd x0 = n_cur;
        while(Newton_step > 0){
          Newton_step -= 1;
          Eigen::VectorXd sum_vec2 = (mu_vec.array()+ s_sq_vec.array()+ x0.array());
          sum_vec2 = trun_vec_fun(sum_vec2,50);
          x0 = (x0.array() - (lib * sum_vec2.array().exp() - alpha.array() - rho*(m_cur.array() - x0.array()) - obs_vec.array())/
            (lib * sum_vec2.array().exp() + rho)).matrix();
        }
        n_cur = x0;
        m_cur = inv_A * (- alpha + rho * n_cur).matrix();
        alpha = alpha + rho * 1.618 * (m_cur - n_cur);
        Eigen::VectorXd sum_vec3 = (mu_vec.array() + s_sq_vec.array()+ m_cur.array());
        sum_vec3 = trun_vec_fun(sum_vec3,50);
        double elbo_cur = 0;
        double lib = S_depth(i);
        elbo_cur = - (m_cur.array() * obs_vec.array()).sum() +
          lib * sum_vec3.array().exp().sum();
        elbo_cur = elbo_cur + 0.5 * (((Theta * (m_cur)).array() * (m_cur).array()).sum());
        // for(int i=GRN_num;i<p;i++){
        //   elbo_pre += 0.5 * Theta(i,i) * m_cur(i) * m_cur(i);
        // }
        if(abs((elbo_cur - elbo_pre)/elbo_pre) < th){
          // cout<<"ELBO has been converged."<<endl;
          conver_num = conver_num + 1;
        }
        
        if(abs((elbo_cur - elbo_pre)/elbo_pre) < th) break;
        else elbo_pre = elbo_cur;
      }
      m_mat.row(i) = m_cur;
      n_mat.row(i) = n_cur;
    }
  }
  omp_destroy_lock(&lock);
  m_mat_all.leftCols(GRN_num) = m_mat;
  // if(!is_update_nonGRNpart) {
  //   return(m_mat_all);
  // } 
  //
  if((non_GRN_num>0) && (is_update_nonGRNpart == true)){
    
    m_mat_start = m_mat_start_all.rightCols(non_GRN_num);
    s_sq_mat = s_sq_mat_all.rightCols(non_GRN_num);
    obs_mat1 = obs_mat_all.rightCols(non_GRN_num);
    mu_vec = mu_vec_all.tail(non_GRN_num);
    Theta = Theta_mat.block(GRN_num,GRN_num,non_GRN_num,non_GRN_num);
    // Q_sum.setZero(p,p);
    // for(int g = 0;g < celltypes_num; g++){
    //   Q_sum  += double((1.0/celltypes_num)) * as<MatrixXd>(Theta_mat_array[g]);
    //   Theta_struct[g] = as<MatrixXd>(Theta_mat_array[g]);
    // }
    if(is_rho_constant_nonGRN){
      rho = rho;
    } else {
      EigenSolver<MatrixXd> eig(Theta);
      //
      Eigen::VectorXd ev = eig.eigenvalues().real();
      double l_ev = ev(0);
      double s_ev = ev(non_GRN_num - 1);
      // double delta = double(S_depth.sum()/n);
      // delta *= double((s_sq_mat/2 + m_mat_start).array().exp().sum()/(n*p));
      double delta = 0.0;
      int count = 0;
      for(int i=0;i<n;i++){
        if(var_index_m(i) == 1){
          Eigen::VectorXd s_sq_vec = s_sq_mat.row(i)/2;
          Eigen::VectorXd m_cur = m_mat_start.row(i);
          Eigen::VectorXd sum_vec4 = (s_sq_vec.array() + m_cur.array()
                                        + mu_vec.array());
          sum_vec4 = trun_vec_fun(sum_vec4,50);
          delta = delta + S_depth(i) * (sum_vec4.array().exp().sum());
          count = count + 1;
        }
        // Eigen::VectorXd s_sq_vec = s_sq_mat.row(i_index)/2;
        // Eigen::VectorXd m_cur = m_mat_start.row(i_index);
        // delta = delta + S_depth(i_index) * ((s_sq_vec.array() + m_cur.array()
        //                                        + mu_vec.array()).array().exp().sum());
      }
      delta /= count;
      if(delta < s_ev) rho = sqrt(delta * s_ev+0.001);
      else if (delta > l_ev) rho = sqrt(delta * l_ev+0.001);
      else rho = delta;
      //
      if((rho<0.1)) rho = 0.1; 
    }
    // tmp = rho * identity_mat.block(GRN_num,GRN_num,non_GRN_num,non_GRN_num) + Theta.block(GRN_num,GRN_num,non_GRN_num,non_GRN_num);
    Eigen::MatrixXd tmp = (rho * (identity_mat.block(GRN_num,GRN_num,non_GRN_num,non_GRN_num)).array()).array() + Theta.array();
    // inv_share.block(GRN_num,GRN_num,non_GRN_num,non_GRN_num) = tmp.inverse();
    // inv_A = inv_share.block(GRN_num,GRN_num,non_GRN_num,non_GRN_num);
    inv_A = tmp.inverse();    
    m_mat = m_mat_start;
    // m_mat.setZero(n,p);
    n_mat = m_mat_start;
    // n_mat.setZero(n,p);
    conver_num = 0;
    // ADMM
    omp_init_lock(&lock);
    omp_set_lock(&lock);
    
    omp_set_num_threads(core_num);
#pragma omp parallel for
    for(int i=0;i<n;i++){
      // cout<<"i = "<<i<<endl;
      if(var_index_m(i) == 1){
        // //
        int ADMM_step = ADMM_maxiter;
        // int id_share = -1;
        // for(int g=0;g<celltypes_num;g++){
        //   if(abs(U_mat(i_index,g) -1) < 0.001){
        //     id_share = g;
        //   }
        // }
        // if(id_share == -1){
        //   Eigen::MatrixXd Q;
        //   Q.setZero(p,p);
        //   // P.setZero(p,p);
        //
        //   for(int g=0; g<celltypes_num;g++){
        //     Q += U_mat(i_index,g) * Theta_struct[g];
        //     // P += U_mat(i,g) * Theta_struct[g];
        //   }
        //
        //
        //   Q +=  rho * identity_mat;
        //   //VectorXi member = get_com(Q);
        //   //inv_A = block_inverse(Q, member, p);
        //
        //   inv_A = Q.inverse();
        // } else {
        //   inv_A = inv_share[id_share];
        // }
        // Eigen::VectorXd mu_sum;
        // mu_sum.setZero(p);
        // for(int g=0; g<celltypes_num;g++){
        // Eigen::VectorXd mu_vec = mu_mat.row(g);
        // mu_sum = mu_sum + U_mat(i_index,g) * (Theta_struct[g] * (mu_vec));
        // }
        //
        Eigen::VectorXd m_cur, n_cur, alpha;
        alpha.setZero(p);
        m_cur = m_mat_start.row(i);
        n_cur = m_mat_start.row(i);
        Eigen::VectorXd s_sq_vec = s_sq_mat.row(i)/2;
        Eigen::VectorXd obs_vec = obs_mat1.row(i);
        //
        double elbo_pre = 0;
        double lib = S_depth(i);
        Eigen::VectorXd sum_vec5 = (mu_vec.array() + s_sq_vec.array()+ m_cur.array());
        sum_vec5 = trun_vec_fun(sum_vec5,50);
        elbo_pre = - (m_cur.array() * obs_vec.array()).sum() +
          lib * sum_vec5.array().exp().sum();
        elbo_pre = elbo_pre + 0.5  * (((Theta * (m_cur)).array() * (m_cur).array()).sum());
        // for(int i=GRN_num;i<p;i++){
        //   elbo_pre += 0.5 * Theta(i,i) * m_cur(i) * m_cur(i);
        // }
        
        while(ADMM_step>0){
          ADMM_step -= 1;
          int Newton_step = 10;
          Eigen::VectorXd x0 = n_cur;
          while(Newton_step > 0){
            Newton_step -= 1;
            Eigen::VectorXd sum_vec6 = (mu_vec.array()+ s_sq_vec.array()+ x0.array());
            sum_vec6 = trun_vec_fun(sum_vec6,50);
            x0 = (x0.array() - (lib * sum_vec6.array().exp() - alpha.array() - rho*(m_cur.array() - x0.array()) - obs_vec.array())/
              (lib * sum_vec6.array().exp() + rho)).matrix();
          }
          n_cur = x0;
          m_cur = inv_A * (- alpha + rho * n_cur).matrix();
          alpha = alpha + rho * 1.618 * (m_cur - n_cur);
          double elbo_cur = 0;
          double lib = S_depth(i);
          Eigen::VectorXd sum_vec7 = (mu_vec.array() + s_sq_vec.array()+ m_cur.array());
          sum_vec7 = trun_vec_fun(sum_vec7,50);
          elbo_cur = - (m_cur.array() * obs_vec.array()).sum() +
            lib * sum_vec7.array().exp().sum();
          elbo_cur = elbo_cur + 0.5 * (((Theta * (m_cur)).array() * (m_cur).array()).sum());
          // for(int i=GRN_num;i<p;i++){
          //   elbo_pre += 0.5 * Theta(i,i) * m_cur(i) * m_cur(i);
          // }
          if(abs((elbo_cur - elbo_pre)/elbo_pre) < th){
            // cout<<"ELBO has been converged."<<endl;
            conver_num = conver_num + 1;
          }
          
          if(abs((elbo_cur - elbo_pre)/elbo_pre) < th) break;
          else elbo_pre = elbo_cur;
        }
        m_mat.row(i) = m_cur;
        n_mat.row(i) = n_cur;
      }
    }
    omp_destroy_lock(&lock);
    //
    m_mat_all.rightCols(non_GRN_num) = m_mat;
  }
  
  //   // cout<<"the number of convergence is "<<conver_num<<endl;
  //   // return m_mat;
  //   // return n_mat;
  return m_mat_all;
}

// [[Rcpp::export]]
Eigen::MatrixXd Newtown(Eigen::VectorXd diag_Theta_all, Eigen::MatrixXd s_sq_mat_start_all,
                        Eigen::MatrixXd m_mat_all, Eigen::VectorXd S_depth_use, Eigen::VectorXi var_index_s,
                        int p_GRN,
                        int core_num = 1,
                        int newtown_maxiter = 50, double th = 1e-3, double low_bound = 1e-8,
                        bool is_update_nonGRNpart = true){
  int n = m_mat_all.rows();
  int p = m_mat_all.cols();
  int non_GRN = p - p_GRN;
  // Eigen::VectorXi var_index_s_densy = Eigen::VectorXi(var_index_s);
  Eigen::VectorXi var_index_s_densy = var_index_s;
  Eigen::MatrixXd s_sq_mat_all = s_sq_mat_start_all;
  //GRN
  Eigen::MatrixXd s_sq_mat;
  s_sq_mat = s_sq_mat_start_all.leftCols(p_GRN);
  // s_sq_mat.setZero(n,p_GRN);
  Eigen::MatrixXd s_sq_mat_start = s_sq_mat_start_all.leftCols(p_GRN);
  Eigen::MatrixXd m_mat = m_mat_all.leftCols(p_GRN);
  Eigen::VectorXd diag_Theta = diag_Theta_all.head(p_GRN);
  //
  omp_init_lock(&lock);
  omp_set_lock(&lock);
  //
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int i=0;i<n;i++){
    if(var_index_s_densy(i) == 1){
      int Newstep = newtown_maxiter;
      Eigen::VectorXd m_vec = m_mat.row(i);
      Eigen::VectorXd c = (S_depth_use(i) * m_vec.array().exp()).matrix();
      Eigen::VectorXd x0 = s_sq_mat_start.row(i);
      double elbo_pre = 0;
      Eigen::VectorXd vec_use1= ((x0).array()/2);
      vec_use1 = trun_vec_fun(vec_use1,50);
      elbo_pre = (c.array() * vec_use1.array().exp()).sum() + 0.5 * ((x0).array() * 
        diag_Theta.array()).sum() - 0.5 * x0.array().log().sum();
      while (Newstep>0) {
        Newstep -= 1;
        Eigen::VectorXd vec_use2= ((x0).array()/2);
        vec_use2 = trun_vec_fun(vec_use2,50);
        Eigen::ArrayXd s_exp = vec_use2.array().exp();
        Eigen::ArrayXd firsr_d = 1/x0.array();
        Eigen::ArrayXd secd = (1/x0.array())*(1/x0.array());
        
        x0 = (x0.array() - (((c/2).array()) * s_exp + 0.5 * diag_Theta.array() -
          0.5 * firsr_d) / (((c/4).array()) * s_exp + 0.5 * secd)).matrix();
        for(int j=0;j<p_GRN;j++)
          if(x0(j) <= 0)  x0(j) = low_bound;
          Eigen::VectorXd vec_use3= ((x0).array()/2);
          vec_use3 = trun_vec_fun(vec_use3,50);
          double elbo_cur = 0;
          elbo_cur = (c.array() * vec_use3.array().exp()).sum() + 0.5 * ((x0).array() * 
            diag_Theta.array()).sum() - 0.5 * x0.array().log().sum();
          if((abs(elbo_cur - elbo_pre)/(abs(elbo_pre) + (1e-2))) < th) break;
          else elbo_pre = elbo_cur;
      }
      s_sq_mat.row(i) = x0;
      //
    }
    //
  }
  omp_destroy_lock(&lock);
  //
  s_sq_mat_all.leftCols(p_GRN) = s_sq_mat;
  //
  // if(!is_update_nonGRNpart) return(s_sq_mat_all);
  //   
  ///  non GRN
  if((non_GRN>0) && (is_update_nonGRNpart == true)){
    // s_sq_mat.setZero(n,non_GRN);
    s_sq_mat = s_sq_mat_start_all.rightCols(non_GRN);
    s_sq_mat_start = s_sq_mat_start_all.rightCols(non_GRN);
    m_mat = m_mat_all.rightCols(non_GRN);
    diag_Theta = diag_Theta_all.tail(non_GRN);
    //
    omp_init_lock(&lock);
    omp_set_lock(&lock);
    //
    omp_set_num_threads(core_num);
#pragma omp parallel for
    for(int i=0;i<n;i++){
      if(var_index_s_densy(i) == 1){
        int Newstep = newtown_maxiter;
        Eigen::VectorXd m_vec = m_mat.row(i);
        Eigen::VectorXd c = (S_depth_use(i) * m_vec.array().exp()).matrix();
        Eigen::VectorXd x0 = s_sq_mat_start.row(i);
        double elbo_pre = 0;
        Eigen::VectorXd vec_use4= ((x0).array()/2);
        vec_use4 = trun_vec_fun(vec_use4,50);
        elbo_pre = (c.array() * vec_use4.array().exp()).sum() + 0.5 * ((x0).array() *
          diag_Theta.array()).sum() - 0.5 * x0.array().log().sum();
        while (Newstep>0) {
          Newstep -= 1;
          Eigen::VectorXd vec_use5= ((x0).array()/2);
          vec_use5 = trun_vec_fun(vec_use5,50);
          Eigen::ArrayXd s_exp = vec_use5.array().exp();
          Eigen::ArrayXd firsr_d = 1/x0.array();
          Eigen::ArrayXd secd = (1/x0.array())*(1/x0.array());
          
          x0 = (x0.array() - (((c/2).array()) * s_exp + 0.5 * diag_Theta.array() -
            0.5 * firsr_d) / (((c/4).array()) * s_exp + 0.5 * secd)).matrix();
          for(int j=0;j<non_GRN;j++)
            if(x0(j) <= 0)  x0(j) = low_bound;
            double elbo_cur = 0;
            Eigen::VectorXd vec_use6= ((x0).array()/2);
            vec_use6 = trun_vec_fun(vec_use6,50);
            elbo_cur = (c.array() * vec_use6.array().exp()).sum() + 0.5 * ((x0).array() *
              diag_Theta.array()).sum() - 0.5 * x0.array().log().sum();
            if((abs(elbo_cur - elbo_pre)/(abs(elbo_pre) + (1e-2))) < th) break;
            else elbo_pre = elbo_cur;
        }
        s_sq_mat.row(i) = x0;
        //
      }
      //
    }
    omp_destroy_lock(&lock);
    //
    s_sq_mat_all.rightCols(non_GRN) = s_sq_mat;
    //
  }
  //   
  //   
  return s_sq_mat_all;
}



Eigen::VectorXd upper_sparse(Eigen::SparseMatrix<double> mat_use){
  //
  int p = mat_use.cols();
  Eigen::VectorXd upper_vec;
  upper_vec.setZero((p*(p - 1))/2);
  //
  int count = 0;
  for(int i = 1;i<p;i++){
    for(int j = 0;j<i;j++){
      upper_vec(count) = mat_use.coeff(j,i);
      count = count + 1;
    }
  }
  // // //
  // Eigen::SparseVector<double> upper_vec_sparse((p*(p - 1))/2);
  // upper_vec_sparse  = upper_vec.sparseView();
  // upper_vec_sparse.makeCompressed();
  // return upper_vec_sparse;
  return upper_vec;
}


// [[Rcpp::export]]
Rcpp::List Theta_fun(Eigen::MatrixXd U_mat,
                     Rcpp::List Theta_mat_list,
                     Rcpp::List m_mat_list,Rcpp::List s_sq_mat_list,
                     int p_GRN,
                     double lambda_use,
                     Eigen::MatrixXd zero_GRN,
                     double stop_threshold,
                     Eigen::VectorXi var_index_Theta,
                     bool penalize_diagonal = false,
                     bool zero_GRN_use = false,
                     int core_num = 1){
  
  //////////////////////////////
  int celltypes_num = U_mat.cols();
  int n = U_mat.rows();
  Eigen::MatrixXd m_mat0 = m_mat_list[0];
  int p = (m_mat0).cols();
  int p_nonGRN = p - p_GRN;
  Rcpp::List Theta_mat_res_list(celltypes_num);
  //
  //stack
  Eigen::MatrixXd Theta_mat_stackall = stack_sparse_fun(Theta_mat_list,celltypes_num,p,p);
  Eigen::MatrixXd m_mat_stackall = stack_fun(m_mat_list,celltypes_num,n,p);
  Eigen::MatrixXd s_sq_mat_stackall = stack_fun(s_sq_mat_list,celltypes_num,n,p);
  //
  Environment pkg = Environment::namespace_env("glasso");
  Function glasso = pkg["glasso"];
  //
  /////////////////
  omp_init_lock(&lock);
  omp_set_lock(&lock);
  
  omp_set_num_threads(core_num);
#pragma omp parallel for
  for(int g=0;g<celltypes_num;g++){
    //
    Eigen::MatrixXd Theta_mat_res;
    Theta_mat_res.setZero(p,p);
    // Eigen::MatrixXd m_mat = m_mat_list[g];
    Eigen::MatrixXd m_mat = m_mat_stackall.middleCols(g * p,p);
    // Eigen::MatrixXd s_sq_mat = s_sq_mat_list[g];
    Eigen::MatrixXd s_sq_mat = s_sq_mat_stackall.middleCols(g * p,p);
    Eigen::VectorXd U_vec = U_mat.col(g);
    //
    Eigen::MatrixXd m_mat1 = m_mat;
    Eigen::VectorXd s_sq_sum;
    s_sq_sum.setZero(p);
    double sample_sum = 0;
    //
    for(int i=0;i<n;i++){
      m_mat1.row(i) = (m_mat1.row(i)).array() * U_mat(i,g);
      Eigen::VectorXd s_sq_vec = (s_sq_mat.row(i)).array() * U_mat(i,g);
      s_sq_sum = s_sq_sum.array() + s_sq_vec.array();
      sample_sum = sample_sum + U_mat(i,g);
    }
    Eigen::MatrixXd sigma_mat = ((m_mat1.adjoint()) * m_mat);
    sigma_mat.diagonal() = (sigma_mat.diagonal()).array() + (s_sq_sum).array();
    sigma_mat = sigma_mat.array() / sample_sum;
    sigma_mat.diagonal() = (sigma_mat.diagonal()).array() + epsilon3;
    // //
    Eigen::MatrixXd sigma_hat_GRN = sigma_mat.block(0,0,p_GRN,p_GRN);
    Eigen::MatrixXd Theta_mat_glasso;
    if(var_index_Theta(g) == 1){
      Rcpp::List res_list;
      if(zero_GRN_use){
        res_list = glasso(Named("s",Rcpp::wrap(sigma_hat_GRN)),
                          Named("rho",Rcpp::wrap((lambda_use * 2/(sample_sum)))),
                          Named("penalize.diagonal",Rcpp::wrap(penalize_diagonal)),
                          Named("zero",Rcpp::wrap(zero_GRN)),
                          Named("thr",stop_threshold));
      }else{
        res_list = glasso(Named("s",Rcpp::wrap(sigma_hat_GRN)),
                          Named("rho",Rcpp::wrap((lambda_use * 2/(sample_sum)))),
                          Named("penalize.diagonal",Rcpp::wrap(penalize_diagonal)),
                          Named("thr",stop_threshold));
      }
      // NumericMatrix Theta_mat_hat = as<NumericMatrix> (res_list["wi"]);
      Theta_mat_glasso = as<Eigen::MatrixXd>(res_list["wi"]);
    }else{
      // Eigen::SparseMatrix<double> Theta_mat_use0 = Theta_mat_list[g];
      // Eigen::MatrixXd Theta_mat_use = MatrixXd(Theta_mat_use0);
      Eigen::MatrixXd Theta_mat_use = Theta_mat_stackall.middleCols(g * p,p);
      // Eigen::MatrixXd Theta_mat_use = Theta_mat_list[g];
      Eigen::MatrixXd Theta_mat_GRN = Theta_mat_use.block(0,0,p_GRN,p_GRN);
      Theta_mat_glasso = Theta_mat_GRN;
    }
    
    Theta_mat_res.block(0,0,p_GRN,p_GRN) = Theta_mat_glasso;
    if(p_nonGRN > 0){
      Eigen::MatrixXd sigma_mat_nonGRN = sigma_mat.block(p_GRN,p_GRN,p_nonGRN,p_nonGRN);
      Eigen::MatrixXd Theta_mat_nonGRN = sigma_mat_nonGRN.inverse();
      Theta_mat_res.block(p_GRN,p_GRN,p_nonGRN,p_nonGRN) = Theta_mat_nonGRN;
    }
    
    //
    
    // Eigen::SparseMatrix<double> Theta_mat_res_sparse = Theta_mat_res.sparseView();
    // Theta_mat_res_list[g] = Theta_mat_res_sparse;
    Theta_mat_stackall.middleCols(g * p,p) = Theta_mat_res;
  }
  omp_destroy_lock(&lock);
  //
  //reshape to list
  for(int g = 0;g<celltypes_num;g++){
    Eigen::MatrixXd Theta_mat_res_tran = Theta_mat_stackall.middleCols(g * p,p);
    Eigen::SparseMatrix<double> Theta_mat_res_sparse = Theta_mat_res_tran.sparseView();
    Theta_mat_res_list[g] = Theta_mat_res_sparse;
  }
  //
  //
  return Theta_mat_res_list;
}


Eigen::VectorXd BIC_fun(Rcpp::List Theta_mat_list,
                        Eigen::MatrixXd U_mat,
                        Rcpp::List m_mat_list,Rcpp::List s_sq_mat_list,
                        Eigen::SparseMatrix<int> obs_mat,
                        Eigen::MatrixXd mu_mat,
                        Eigen::VectorXd S_depth,
                        Eigen::VectorXd pi_vec,
                        Eigen::VectorXd log_cul_sum,
                        int p_GRN,
                        Eigen::VectorXd log_det_GRN_vec,
                        int core_num){
  //
  Eigen::VectorXd gene_GRN_index_use(p_GRN);
  for(int j = 0;j<p_GRN;j++){
    gene_GRN_index_use(j) = j + 1;
  }
  //
  int celltypes_num = U_mat.cols();
  Eigen::VectorXd bic_value_vec(celltypes_num);
  //
  int core_num_use = core_num;
  if(celltypes_num < core_num){
    core_num_use = celltypes_num;
  }
  //
  Eigen::VectorXd loglikelihood_vec = ELBO_fast(Theta_mat_list,
                                                U_mat,
                                                m_mat_list,
                                                s_sq_mat_list,
                                                obs_mat,
                                                mu_mat,
                                                S_depth,
                                                pi_vec,
                                                p_GRN,
                                                log_cul_sum,
                                                log_det_GRN_vec,
                                                core_num_use);
  //
  for(int g=0;g<celltypes_num;g++){
    //
    Eigen::SparseMatrix<double> Theta_mat_sparse = Theta_mat_list[g];
    Eigen::MatrixXd Theta_mat = MatrixXd(Theta_mat_sparse);
    Eigen::MatrixXd Theta_mat_con = (Theta_mat.array().sign()).array().abs();
    double num_nonzero = (Theta_mat_con.sum() - p_GRN) * 0.5;
    //
    double sample_group = (U_mat.col(g)).sum();
    // int num_nonzero = 0;
    // for(int i =0;j<p_GRN;j++){
    //   for(int j = 0;j<p_GRN;j++){
    //     if()
    //     }
    //   }
    bic_value_vec(g) = (-2) * loglikelihood_vec(g) + log(sample_group) * num_nonzero;
  }
  //
  // return Rcpp::List::create(Named("bic_value_vec") = bic_value_vec,
  //                           Named("loglikelihood_vec") = loglikelihood_vec);
  return bic_value_vec;
}






/////////////////////////////////

// [[Rcpp::export]]
Rcpp::List scGeneNet_optim(Eigen::VectorXi celltypes_label_input,
                           Eigen::MatrixXd U_mat_input,
                           Eigen::MatrixXd mu_mat_input,
                           Rcpp::List m_mat_list_input,
                           Rcpp::List s_sq_mat_list_input,
                           Eigen::VectorXd pi_vec_input,
                           Eigen::SparseMatrix<int> obs_mat_input,
                           Eigen::VectorXd S_depth_input,
                           Eigen::VectorXi gene_GRN_index_use_input,
                           Rcpp::List Theta_mat_list_input,
                           Eigen::MatrixXd zero_GRN_input,
                           bool if_zero,
                           Eigen::VectorXd log_det_vec_input,
                           double lambda_use,
                           bool penalize_diagonal,
                           int maxit_nonGRN,
                           double MS_update_threshold,
                           int maxit,
                           int minit,
                           double ELBO_threshold,
                           bool verbose,
                           int ADMM_max_step,
                           double ADMM_threshold,
                           int Newton_max_step,
                           double Newton_threshold,
                           double Theta_threshold,
                           int core_num,
                           bool U_fix){
  //////////////////////
  //
  Eigen::VectorXi celltypes_label = (celltypes_label_input);
  Eigen::MatrixXd U_mat = (U_mat_input);
  Eigen::MatrixXd mu_mat = (mu_mat_input);
  Rcpp::List m_mat_list = clone(m_mat_list_input);
  Rcpp::List s_sq_mat_list = clone(s_sq_mat_list_input);
  Eigen::VectorXd pi_vec = (pi_vec_input);
  Eigen::SparseMatrix<int> obs_mat = (obs_mat_input);
  Eigen::VectorXd S_depth = (S_depth_input);
  Eigen::VectorXi gene_GRN_index_use = (gene_GRN_index_use_input);
  Rcpp::List Theta_mat_list = clone(Theta_mat_list_input);
  Eigen::MatrixXd zero_GRN = (zero_GRN_input);
  Eigen::VectorXd log_det_vec = (log_det_vec_input);
  //////////////////////////
  int celltypes_num = U_mat.cols();
  int n = U_mat.rows();
  int p = mu_mat.cols();
  int p_GRN = gene_GRN_index_use.size();
  int p_nonGRN = p - p_GRN;
  int core_min_use = core_num;
  if(celltypes_num < core_num){
    core_min_use = celltypes_num;
  }
  
  //
  Eigen::VectorXi celltypes_label_trans;
  Eigen::MatrixXd U_mat_trans;
  Eigen::MatrixXd mu_mat_trans;
  Rcpp::List m_mat_list_trans;
  Rcpp::List s_sq_mat_list_trans;
  Eigen::VectorXd pi_vec_trans;
  Eigen::SparseMatrix<int> obs_mat_trans;
  Eigen::VectorXd S_depth_trans;
  Eigen::VectorXi gene_GRN_index_use_trans;
  Rcpp::List Theta_mat_list_trans;
  Eigen::MatrixXd zero_GRN_trans;
  bool zero_GRN_use_trans;
  Eigen::VectorXd log_det_vec_trans;
  //
  // Eigen::VectorXd diag_Theta_all_save;
  // Eigen::MatrixXd s_sq_mat_start_save;
  // Eigen::MatrixXd m_mat_start_save;
  // Eigen::VectorXd S_depth_save;
  // Eigen::VectorXi MS_update_vec_save;
  // int p_GRN_save;
  // int core_num_save;
  // int Newton_max_step_save;
  // double Newton_threshold_save;
  // double low_bound_save;
  ///////
  ///////
  // Rcpp::List Theta_mat_list_pd(celltypes_num);
  //
  Eigen::MatrixXi MS_update_index = Eigen::MatrixXi::Ones(n,celltypes_num);
  //
  // if(sample_filter == true){
  //   Eigen::MatrixXd U_mat1 = clustering_fun(Theta_mat_list, U_mat, m_mat_list, s_sq_mat_list, obs_mat,
  //                                           mu_mat, S_depth, pi_vec, log_det_vec); 
  //   for(int i =0;i<n;i++){
  //     for(int g=0;g<celltypes_num;g++){
  //       if(U_mat1(i,g)>(MS_update_threshold * MS_update_threshold_ratio)){
  //         MS_update_index(i,g) = 1;
  //       }else{
  //         MS_update_index(i,g) = 0;
  //       }
  //     }
  //   }
  // }
  // Eigen::MatrixXi MS_update_index;
  // MS_update_index.resize(n,celltypes_num);
  // for(int i =0;i<n;i++){
  //   for(int g=0;g<celltypes_num;g++){
  //     if(U_mat(i,g)>MS_update_threshold){
  //       MS_update_index(i,g) = 1;
  //     }else{
  //       MS_update_index(i,g) = 0;
  //     }
  //   }
  // }
  // MS_update_index.makeCompressed();
  //
  int block_iter_max = 10;
  double block_threshold = ELBO_threshold;
  double low_bound = 1e-8;
  //Initial
  int iter = 0;
  bool con_run = true;
  //
  Eigen::SparseMatrix<double>log_cul_mat_obs_mat = log_cul_mat(obs_mat);
  
  Eigen::VectorXd log_cul_mat_obs_sum_all = rowSums_sparse(log_cul_mat_obs_mat);
  Eigen::VectorXd log_cul_mat_obs_sum_GRN;
  if(p_nonGRN == 0){
    log_cul_mat_obs_sum_GRN = log_cul_mat_obs_sum_all;
  }else{
    log_cul_mat_obs_sum_GRN = rowSums_sparse(log_cul_mat_obs_mat.leftCols(p_GRN));
  }
  //
  Eigen::VectorXd log_det_vec_GRN = log_det_GRN(Theta_mat_list,celltypes_num,p_GRN,core_min_use);
  Eigen::VectorXd log_det_vec_GRN_cur = log_det_vec_GRN;
  //
  
  Eigen::VectorXi var_index_vec(celltypes_num);
  for(int g=0;g<celltypes_num;g++){
    var_index_vec(g) = 1;
  }
  //
  Eigen::VectorXd ELBO_vec_start = ELBO_fast(Theta_mat_list,U_mat,m_mat_list,s_sq_mat_list, obs_mat,
                                             mu_mat,S_depth, pi_vec,p,
                                             log_cul_mat_obs_sum_all,
                                             log_det_vec,core_min_use);
  Eigen::VectorXd ELBO_vec_cur = ELBO_vec_start;
  //
  double diff1_sign_start = 1e-4;
  bool non_pendiag_finish = false;
  bool non_pendiag_finish_cur = false;
  bool penalize_diagonal_use = false;
  double purity_threshold = 0.99;
  //
  double purity_start = 0;
  double purity_cur = 0;
  //
  //
  while(con_run){
    //
    iter = iter + 1;
    //
    if(iter>1){
      ELBO_vec_start = ELBO_vec_cur;
    }
    //
    if(non_pendiag_finish == true){
      if(penalize_diagonal == true){
        if(non_pendiag_finish_cur == false){
          //
          penalize_diagonal_use = true;
          //
          for(int g=0;g<celltypes_num;g++){
            var_index_vec(g) = 1;
          }
          // if(verbose == true){
          //   cout<<"Start run penalized diagonal version."<<endl;
          // }
        }
      }
    }
    //
    // bool non_pendiag_finish_cur = non_pendiag_finish;
    non_pendiag_finish_cur = non_pendiag_finish;
    // cout<<"non_pendiag_finish:"<<non_pendiag_finish<<endl;
    //
    if(verbose == true){
      if(iter == 1){
        cout<<"lambda: "<<lambda_use<<endl; 
        cout<<"----------------------------------------------------------------"<<endl;
        cout<<"iter: "<<iter<<endl;
      }else{
        cout<<"iter: "<<iter<<endl;
      }
    }
    //
    //////////////////////////////////////////////////////////////////
    //Step 1: update result of clustering
    // cout<<"Part 1."<<endl;
    //
    if(iter>1){
      purity_start = purity_cur;
    }
    //
    bool con_pair = true;
    Eigen::VectorXi cluster_predict_true(celltypes_num);
    for(int g=0;g<celltypes_num;g++){
      cluster_predict_true(g) = g;
    }
    Eigen::VectorXi celltypes_label_start = celltypes_label;
    //
    // cout<<"Part 1.1"<<endl;
    if(iter>1){
      if((purity_cur > purity_threshold) || (iter>maxit_nonGRN)){
        non_pendiag_finish = true;
        // cout<<"point1."<<endl;
      }
      
      //
      bool con_clustering;
      // if(clustering_con == true){
      //   con_clustering = (purity_cur<purity_threshold) || (purity_cur == purity_threshold);
      // }else{
      //   con_clustering = ((purity_cur<purity_threshold) || (purity_cur == purity_threshold)) && (iter<(maxit_nonGRN + 1)) && (non_pendiag_finish == false);
      // }
      // cout<<"purity_cur:"<<purity_cur<<endl;
      // cout<<"purity_threshold:"<<purity_threshold<<endl;
      // cout<<"iter:"<<iter<<endl;
      // cout<<"non_pendiag_finish:"<<non_pendiag_finish<<endl;
      //
      con_clustering = ((purity_cur<purity_threshold) || (purity_cur == purity_threshold)) && (iter<(maxit_nonGRN + 1)) && (non_pendiag_finish == false);
      //
      if(con_clustering && (U_fix == false)){
        // if(verbose == true){
        // cout<<"Run clustering."<<endl;
        // }
        // for(int g = 0;g<celltypes_num;g++){
        //   cout<<"log_det_vec("<<g<<"):"<<log_det_vec(g)<<endl;
        //   }
        Eigen::MatrixXd U_mat0 = clustering_fun(Theta_mat_list, U_mat, m_mat_list, s_sq_mat_list, obs_mat,
                                                mu_mat, S_depth, pi_vec, log_det_vec);
        // cout<<"Part 1.2"<<endl;
        Eigen::VectorXi celltypes_label0 = celltypes_label_fun(U_mat0);
        // cout<<"Part 1.3"<<endl;
        //
        Eigen::VectorXi cluster_predict_true = cluster_predict_true_fun(celltypes_label0,celltypes_label_start,celltypes_num);
        // cout<<"Part 1.4"<<endl;
        int count_group = count_group_fun(cluster_predict_true);
        if(!(count_group == celltypes_num)){
          con_pair = false;
        }
        purity_cur = purity_fun(celltypes_label0,celltypes_label_start,cluster_predict_true,celltypes_num);
        // cout<<"Part 1.5"<<endl;
        //
        if(con_pair){
          // cout<<"Part 1.6.1"<<endl;
          //permutation
          int count_per = 0;
          for(int g = 0;g<cluster_predict_true.size();g++){
            if(!(cluster_predict_true(g) == g)){
              count_per = count_per + 1;
            }
          }
          //
          if(count_per>0){
            for(int g = 0;g<celltypes_num;g++){
              U_mat.col(cluster_predict_true(g)) = U_mat0.col(g);
            }
            celltypes_label = celltypes_label_fun(U_mat);
            pi_vec = (U_mat.colwise().sum()).array() / n;
            
          }else{
            U_mat = U_mat0;
            celltypes_label = celltypes_label0;
            pi_vec = (U_mat.colwise().sum()).array() / n;
          }
          //
          // cout<<"Part 1.6.2"<<endl;
        }else{
          non_pendiag_finish = true;
          // cout<<"point2."<<endl;
          purity_cur = 1;
        }
        //
        for(int i =0;i<n;i++){
          for(int g=0;g<celltypes_num;g++){
            if(U_mat(i,g)>MS_update_threshold){
              MS_update_index(i,g) = 1;
            }else{
              MS_update_index(i,g) = 0;
            }
          }
        }
        // MS_update_index.makeCompressed();
        //
      }else{
        con_pair = true;
        purity_cur = 1;
      }
    }
    // cout<<"Part 1.7"<<endl;
    //
    Eigen::VectorXd MS_update_ratio  = MS_update_ratio_fun(MS_update_index);
    // if(verbose == true){
    //   for(int g=0;g<celltypes_num;g++){
    //     cout<<"MS_update_ratio("<<g<<"): "<<MS_update_ratio(g)<<endl;
    //   }
    // }
    //
    //
    //Step 2: update mu M S 
    // cout<<"Part 2."<<endl;
    int block_iter = 1;
    Eigen::VectorXd ELBO_rel_change(celltypes_num);
    Eigen::VectorXi var_index_group(celltypes_num);
    for(int g=0;g<celltypes_num;g++){
      var_index_group(g) = 1;
      ELBO_rel_change(g) = 1e+8;
    }
    //
    Eigen::VectorXd ELBO_rel_start;
    if(iter == 1){
      ELBO_rel_start = ELBO_vec_start;
    }else{
      if(non_pendiag_finish == true){
        ELBO_rel_start = ELBO_fast(Theta_mat_list,U_mat, m_mat_list, s_sq_mat_list, obs_mat,
                                   mu_mat, S_depth, pi_vec,p_GRN, log_cul_mat_obs_sum_GRN,
                                   log_det_vec_GRN,core_min_use);
        
      }else{
        ELBO_rel_start = ELBO_fast(Theta_mat_list,U_mat, m_mat_list, s_sq_mat_list, obs_mat,
                                   mu_mat, S_depth, pi_vec, p, log_cul_mat_obs_sum_all,
                                   log_det_vec,core_min_use);
      }
    }
    Eigen::VectorXd ELBO_rel_cur = ELBO_rel_start;
    //
    while(((block_iter < block_iter_max) || (block_iter == block_iter_max)) && (ELBO_rel_change.maxCoeff() > block_threshold) && (var_index_group.sum()>0)){
      // cout<<"block_iter:"<<block_iter<<endl;
      if(block_iter>1){
        ELBO_rel_start = ELBO_rel_cur;
      }
      //
      //update mu
      for(int g = 0;g<celltypes_num;g++){
        ////
        if((var_index_group(g) == 1) && (var_index_vec(g) == 1)){
          if(non_pendiag_finish == false){
            Eigen::MatrixXd m_mat = as<Eigen::MatrixXd> (m_mat_list[g]);
            Eigen::VectorXd U_vec = U_mat.col(g);
            //
            Eigen::MatrixXd m_mat1 = m_mat;
            for(int i =0 ;i<n;i++){
              m_mat1.row(i) = (m_mat1.row(i)).array() * U_vec(i);
            }
            Eigen::VectorXd mu_vec = m_mat1.colwise().sum();
            Eigen::VectorXd mu_vec1 = mu_mat.row(g);
            mu_vec = mu_vec.array()/(U_vec.sum());
            mu_vec = mu_vec.array() + mu_vec1.array();
            mu_mat.row(g) = mu_vec;
            //
          }else{
            Eigen::MatrixXd m_mat0 = as<Eigen::MatrixXd> (m_mat_list[g]);
            Eigen::MatrixXd m_mat = m_mat0.leftCols(p_GRN);
            //
            Eigen::VectorXd U_vec = U_mat.col(g);
            //
            Eigen::MatrixXd m_mat1 = m_mat;
            for(int i =0 ;i<n;i++){
              m_mat1.row(i) = (m_mat1.row(i)).array() * U_vec(i);
            }
            Eigen::VectorXd mu_vec = m_mat1.colwise().sum();
            Eigen::VectorXd mu_vec10 = mu_mat.row(g);
            Eigen::VectorXd mu_vec1 = mu_vec10.head(p_GRN);
            mu_vec = mu_vec.array()/(U_vec.sum());
            mu_vec = mu_vec.array() + mu_vec1.array();
            //
            for(int j = 0;j<p_GRN;j++){
              mu_mat(g,j) = mu_vec(j);
            }
            // mu_mat.block(g,0,1,p_GRN) = mu_vec;
          }
          //
          // Eigen::MatrixXd m_mat = as<Eigen::MatrixXd> (m_mat_list[g]);
          // Eigen::VectorXd U_vec = U_mat.col(g);
          // //
          // Eigen::MatrixXd m_mat1 = m_mat;
          // for(int i =0 ;i<n;i++){
          //   m_mat1.row(i) = (m_mat1.row(i)).array() * U_vec(i);
          // }
          // Eigen::VectorXd mu_vec = m_mat1.colwise().sum();
          // Eigen::VectorXd mu_vec1 = mu_mat.row(g);
          // mu_vec = mu_vec.array()/(U_vec.sum());
          // mu_vec = mu_vec.array() + mu_vec1.array();
          // mu_mat.row(g) = mu_vec;
          
          //
        }
        ///////
      }
      
      //////
      for(int g=0;g<celltypes_num;g++){
        Eigen::MatrixXd s_sq_mat_start = as<Eigen::MatrixXd>(s_sq_mat_list[g]);
        s_sq_mat_start = trun_fun(s_sq_mat_start,50);
        s_sq_mat_list[g] = s_sq_mat_start;
      }
      //////
      //update M
      for(int g = 0;g<celltypes_num;g++){
        if((var_index_group(g) == 1) && (var_index_vec(g) == 1)){
          Eigen::VectorXd mu_vec = mu_mat.row(g);
          Eigen::VectorXi MS_update_vec = MS_update_index.col(g);
          //
          if(p_nonGRN>0){
            if(non_pendiag_finish == false){
              Eigen::MatrixXd m_mat_start = as<Eigen::MatrixXd>(m_mat_list[g]);
              Eigen::MatrixXd m_mat_res = ADMM1(Theta_mat_list[g],
                                                s_sq_mat_list[g],
                                                             obs_mat,mu_vec,
                                                             m_mat_start,
                                                             S_depth,
                                                             MS_update_vec,
                                                             p_GRN,
                                                             core_num,
                                                             ADMM_max_step,
                                                             ADMM_threshold,
                                                             0.1,
                                                             false,
                                                             true,
                                                             true);
              m_mat_list[g] = (m_mat_res);
            }else{
              Eigen::MatrixXd m_mat_start = as<Eigen::MatrixXd>(m_mat_list[g]);
              Eigen::MatrixXd m_mat_res = ADMM1(Theta_mat_list[g],
                                                s_sq_mat_list[g],
                                                             obs_mat,mu_vec,
                                                             m_mat_start,
                                                             S_depth,
                                                             MS_update_vec,
                                                             p_GRN,
                                                             core_num,
                                                             ADMM_max_step,
                                                             ADMM_threshold,
                                                             0.1,
                                                             false,
                                                             true,
                                                             false);
              m_mat_list[g] = (m_mat_res);
            }
          }else{
            Eigen::MatrixXd m_mat_start = as<Eigen::MatrixXd>(m_mat_list[g]);
            Eigen::MatrixXd m_mat_res = ADMM1(Theta_mat_list[g],
                                              s_sq_mat_list[g],
                                                           obs_mat,mu_vec,
                                                           m_mat_start,
                                                           S_depth,
                                                           MS_update_vec,
                                                           p_GRN,
                                                           core_num,
                                                           ADMM_max_step,
                                                           ADMM_threshold,
                                                           0.1,
                                                           false,
                                                           true,
                                                           true);
            m_mat_list[g] = (m_mat_res);
          }
          //
        } 
      }
      
      
      //update S
      for(int g = 0;g<celltypes_num;g++){
        if((var_index_group(g) == 1) && (var_index_vec(g) == 1)){
          //
          Eigen::MatrixXd m_mat_start = as<Eigen::MatrixXd>(m_mat_list[g]);
          // m_mat_start = trun_fun(m_mat_start,50);
          Eigen::VectorXd mu_vec = mu_mat.row(g);
          for(int i = 0;i<n;i++){
            if(MS_update_index(i,g) == 1){
              Eigen::VectorXd m_vec = m_mat_start.row(i);
              m_mat_start.row(i) = m_vec.array() + mu_vec.array();
            }
            // Eigen::VectorXd m_vec = m_mat_start.row(i);
            // m_mat_start.row(i) = m_vec.array() + mu_vec.array();
          }
          Eigen::MatrixXd s_sq_mat_start = as<Eigen::MatrixXd>(s_sq_mat_list[g]);
          // s_sq_mat_start = trun_fun(s_sq_mat_start,50);
          Eigen::SparseMatrix<double> Theta_mat_sparse = Theta_mat_list[g];
          Eigen::MatrixXd Theta_mat0 = MatrixXd(Theta_mat_sparse);
          Eigen::VectorXd diag_Theta_all = Theta_mat0.diagonal();
          //
          Eigen::VectorXi MS_update_vec = MS_update_index.col(g);
          //
          // if((block_iter == 1) && (g == 0) && (iter == 1)){
          //   diag_Theta_all_save = diag_Theta_all;
          //   s_sq_mat_start_save = s_sq_mat_start;
          //   m_mat_start_save = m_mat_start;
          //   S_depth_save = S_depth;
          //   MS_update_vec_save = MS_update_vec;
          //   p_GRN_save = p_GRN;
          //   core_num_save = core_num;
          //   Newton_max_step_save = Newton_max_step;
          //   Newton_threshold_save = Newton_threshold;
          //   low_bound_save = low_bound;
          //   }
          //
          if(p_nonGRN>0){
            if(non_pendiag_finish == false){
              Eigen::MatrixXd s_sq_mat_res = Newtown(diag_Theta_all, s_sq_mat_start,
                                                     m_mat_start, S_depth,MS_update_vec,
                                                     p_GRN,
                                                     core_num,
                                                     Newton_max_step, Newton_threshold, low_bound,
                                                     true);
              s_sq_mat_list[g] = (s_sq_mat_res);
            }else{
              Eigen::MatrixXd s_sq_mat_res = Newtown(diag_Theta_all, s_sq_mat_start,
                                                     m_mat_start, S_depth,MS_update_vec,
                                                     p_GRN,
                                                     core_num,
                                                     Newton_max_step, Newton_threshold, low_bound,
                                                     false);
              s_sq_mat_list[g] = (s_sq_mat_res);
            }
          }else{
            Eigen::MatrixXd s_sq_mat_res = Newtown(diag_Theta_all, s_sq_mat_start,
                                                   m_mat_start, S_depth, MS_update_vec,
                                                   p_GRN,
                                                   core_num,
                                                   Newton_max_step, Newton_threshold, low_bound,
                                                   true);
            s_sq_mat_list[g] = (s_sq_mat_res);
          }
          //
        }
      }
      
      
      
      
      
      block_iter = block_iter + 1;
      //
      if(non_pendiag_finish == true){
        ELBO_rel_cur = ELBO_fast(Theta_mat_list,U_mat, m_mat_list, s_sq_mat_list, obs_mat,
                                 mu_mat, S_depth, pi_vec, p_GRN, log_cul_mat_obs_sum_GRN,
                                 log_det_vec_GRN,core_min_use);
      }else{
        ELBO_rel_cur = ELBO_fast(Theta_mat_list,U_mat, m_mat_list, s_sq_mat_list, obs_mat,
                                 mu_mat, S_depth, pi_vec, p, log_cul_mat_obs_sum_all,
                                 log_det_vec,core_min_use);
      }
      //
      for(int g=0;g<celltypes_num;g++){
        ELBO_rel_change(g) = abs(ELBO_rel_cur(g) - ELBO_rel_start(g))/(abs(ELBO_rel_start(g)) + 1e-2);
      }
      //
      for(int g=0;g<celltypes_num;g++){
        if(ELBO_rel_change(g)>block_threshold){
          var_index_group(g) = 1;
        }else{
          var_index_group(g) = 0;
        }
      }
      //
    }
    ///////////////////////////
    // cout<<"Part 3."<<endl;
    //Step 3: update Theta
    Eigen::MatrixXd Theta_mat_sign_start(celltypes_num,((p_GRN * (p_GRN - 1))/2));
    for(int g = 0;g<celltypes_num;g++){
      Eigen::SparseMatrix<double> Theta_mat0 = Theta_mat_list[g];
      Eigen::SparseMatrix<double> Theta_mat1 = Theta_mat0.block(0,0,p_GRN,p_GRN);
      Eigen::VectorXd upper_sparse_vec = upper_sparse(Theta_mat1);
      Eigen::VectorXd upper_sparse_vec_sign = upper_sparse_vec.array().sign();
      Theta_mat_sign_start.row(g) = upper_sparse_vec_sign;
    }
    //
    if(var_index_vec.sum()>0){
      // bool penalize_diagonal_use;
      // if(non_pendiag_finish == false){
      //   penalize_diagonal_use = false;
      // }else{
      //   // if(non_pendiag_finish_cur == false){
      //   //   penalize_diagonal_use == false;
      //   // }else{
      //   //   penalize_diagonal_use = penalize_diagonal; 
      //   // }
      //     penalize_diagonal_use = penalize_diagonal;
      // }
      // cout<<"penalize_diagonal_use:"<<penalize_diagonal_use<<endl;
      Theta_mat_list = Theta_fun(U_mat,
                                 Theta_mat_list,
                                 m_mat_list,
                                 s_sq_mat_list,
                                 p_GRN,
                                 lambda_use,
                                 zero_GRN,
                                 Theta_threshold,
                                 var_index_vec,
                                 penalize_diagonal_use,
                                 if_zero,
                                 1);
    }
    //
    log_det_vec_GRN_cur = log_det_vec_GRN;
    log_det_vec_GRN = log_det_GRN(Theta_mat_list,celltypes_num,p_GRN,core_min_use);
    //
    if(var_index_vec.sum()>0){
      if(p_nonGRN == 0){
        log_det_vec = log_det(Theta_mat_list,celltypes_num,core_min_use);
      }else{
        if(non_pendiag_finish == false){
          log_det_vec = log_det_block(Theta_mat_list,celltypes_num,p_GRN,core_min_use);
        }else{
          log_det_vec = log_det_vec.array() - log_det_vec_GRN_cur.array() + log_det_vec_GRN.array();
        }
      }
    }
    //
    Eigen::MatrixXd Theta_mat_sign(celltypes_num,((p_GRN * (p_GRN - 1))/2));
    for(int g = 0;g<celltypes_num;g++){
      Eigen::SparseMatrix<double> Theta_mat0 = Theta_mat_list[g];
      Eigen::SparseMatrix<double> Theta_mat1 = Theta_mat0.block(0,0,p_GRN,p_GRN);
      Eigen::VectorXd upper_sparse_vec = upper_sparse(Theta_mat1);
      Eigen::VectorXd upper_sparse_vec_sign = upper_sparse_vec.array().sign();
      Theta_mat_sign.row(g) = upper_sparse_vec_sign;
    }
    
    Eigen::VectorXi diff1_vec_sign(celltypes_num);
    for(int g = 0;g<celltypes_num;g++){
      int count = 0;
      for(int j = 0;j<Theta_mat_sign.cols();j++){
        double value1 = Theta_mat_sign(g,j);
        double value2 = Theta_mat_sign_start(g,j);
        if(!(value1 == value2)){
          count = count + 1;
        }
      }
      diff1_vec_sign(g) = count;
    }
    int diff1_sign;
    if(var_index_vec.sum()>0){
      Eigen::VectorXi diff1_vec_sel = diff1_vec_sign.array() * var_index_vec.array();
      diff1_sign = diff1_vec_sel.maxCoeff();
    }else{
      diff1_sign = 0;
    }
    //
    // 
    //convergence criterion
    if(non_pendiag_finish == true){
      ELBO_vec_cur = ELBO_fast(Theta_mat_list, U_mat, m_mat_list, s_sq_mat_list, obs_mat,
                               mu_mat, S_depth, pi_vec, p_GRN, log_cul_mat_obs_sum_GRN,
                               log_det_vec_GRN,core_min_use);
    }else{
      ELBO_vec_cur = ELBO_fast(Theta_mat_list, U_mat, m_mat_list, s_sq_mat_list, obs_mat,
                               mu_mat, S_depth, pi_vec, p, log_cul_mat_obs_sum_all,
                               log_det_vec,core_min_use);
    }
    //
    Eigen::VectorXd ELBO_vec_change =((ELBO_vec_start.array() - ELBO_vec_cur.array()).array().abs()).array()/ ((ELBO_vec_start.array().abs()).array() + 1e-2).array();
    double diff_ELBO;
    if((var_index_vec.sum())>0){
      Eigen::VectorXd ELBO_vec_change_sel;
      ELBO_vec_change_sel.setZero(celltypes_num);
      for(int g=0;g<celltypes_num;g++){
        if(var_index_vec(g) == 1){
          ELBO_vec_change_sel(g) = ELBO_vec_change(g);
        }
      }
      
      diff_ELBO = ELBO_vec_change_sel.maxCoeff();
    }else{
      diff_ELBO = 0;
    }
    //
    for(int g =0;g<celltypes_num;g++){
      if(ELBO_vec_change(g) > ELBO_threshold){
        var_index_vec(g) = 1;
      }else{
        var_index_vec(g) = 0;
      }
    }
    Eigen::VectorXi Theta_nonzero_vec(celltypes_num);
    for(int g =0;g<celltypes_num;g++){
      Eigen::SparseMatrix<double> Theta_mat = Theta_mat_list[g];
      Eigen::SparseMatrix<double> Theta_mat_GRN = Theta_mat.block(0,0,p_GRN,p_GRN);
      Theta_nonzero_vec(g) = Theta_mat_GRN.nonZeros();
    }
    Theta_nonzero_vec = (Theta_nonzero_vec.array() - p_GRN).array() /2;
    //
    // non_pendiag_finish_cur = non_pendiag_finish;
    //
    if((iter > minit) || (iter == minit)){
      if(penalize_diagonal == false){
        if(((diff_ELBO < ELBO_threshold) && con_pair && ((purity_cur>purity_threshold) || (purity_cur==purity_threshold))) || ((diff1_sign == 0) && ((purity_cur>purity_threshold) || (purity_cur==purity_threshold)))){
          con_run=false;
        }
      }else{
        if(non_pendiag_finish == false){
          if(((diff_ELBO < ELBO_threshold) && con_pair && ((purity_cur>purity_threshold) || (purity_cur==purity_threshold))) || ((diff1_sign == 0) && ((purity_cur>purity_threshold) || (purity_cur==purity_threshold)))){
            non_pendiag_finish=true;
            // cout<<"point3."<<endl;
          }
        }else{
          if(non_pendiag_finish_cur == true){
            if((diff_ELBO < ELBO_threshold) || ((diff1_sign == 0))){
              con_run = false;
            }
          }
        }
      }
    }
    //
    if(iter > maxit){
      if((penalize_diagonal == true) && (non_pendiag_finish == false)){
        non_pendiag_finish = true;
        // cout<<"point4."<<endl;
      }else{
        con_run = false;
      }
    }
    //
    // cout<<"penalize_diagonal:"<<penalize_diagonal<<endl;
    if((penalize_diagonal == true) && (non_pendiag_finish_cur == false) && (non_pendiag_finish == true)){
      // cout<<"clone"<<endl;
      celltypes_label_trans = celltypes_label;
      U_mat_trans = U_mat;
      mu_mat_trans = mu_mat;
      m_mat_list_trans = clone(m_mat_list);
      s_sq_mat_list_trans = clone(s_sq_mat_list);
      pi_vec_trans = pi_vec;
      obs_mat_trans = obs_mat;
      S_depth_trans = S_depth;
      gene_GRN_index_use_trans = gene_GRN_index_use;
      Theta_mat_list_trans = clone(Theta_mat_list);
      zero_GRN_trans = zero_GRN;
      zero_GRN_use_trans = if_zero;
      log_det_vec_trans = log_det_vec;
    }
    //
    if(penalize_diagonal == false){
      celltypes_label_trans = celltypes_label;
      U_mat_trans = U_mat;
      mu_mat_trans = mu_mat;
      m_mat_list_trans = clone(m_mat_list);
      s_sq_mat_list_trans = clone(s_sq_mat_list);
      pi_vec_trans = pi_vec;
      obs_mat_trans = obs_mat;
      S_depth_trans = S_depth;
      gene_GRN_index_use_trans = gene_GRN_index_use;
      Theta_mat_list_trans = clone(Theta_mat_list);
      zero_GRN_trans = zero_GRN;
      zero_GRN_use_trans = if_zero;
      log_det_vec_trans = log_det_vec;
    }
    //
    diff1_sign_start = diff1_sign;
    if(verbose == true){
      cout<<"The convergence citerion:"<<endl;
      
      cout<<"(a) Purity: "<<purity_cur<<endl;
      cout<<"(b) The maximal number of changed sign in precision matices among all component:"<<diff1_sign<<endl;
      cout<<"(c) The relative change of ELBO: "<<diff_ELBO<<endl;
      cout<<"----------------------------------------------------------------"<<endl;
      
      // cout<<"The convergence citerion: Purity: "<<purity_cur<<", con_pair:"<<con_pair<<", Theta_sign:"<<diff1_sign<<", diff_ELBO: "<<diff_ELBO<<endl;
      // for(int g = 0;g<celltypes_num;g++){
      //   cout<<"Theta_nonzero("<<g<<"):"<<Theta_nonzero_vec(g)<<endl;
      // }
      // for(int g = 0;g<celltypes_num;g++){
      //   cout<<"diff1_vec_sign("<<g<<"):"<<diff1_vec_sign(g)<<endl;
      // }
      // for(int g = 0;g<celltypes_num;g++){
      //   cout<<"ELBO_vec_change("<<g<<"):"<<ELBO_vec_change(g)<<endl;
      // }
      
    }
  }
  
  //
  Eigen::VectorXd scGeneNet_bic_VMICL;
  // Eigen::VectorXd scGeneNet_loglikelihood_VMICL;
  Eigen::VectorXd scGeneNet_bic_VICL;
  // Eigen::VectorXd scGeneNet_loglikelihood_VICL;
  scGeneNet_bic_VMICL.setZero(celltypes_num);
  // scGeneNet_loglikelihood_VMICL.setZero(celltypes_num);
  scGeneNet_bic_VICL.setZero(celltypes_num);
  // scGeneNet_loglikelihood_VICL.setZero(celltypes_num);
  
  //
  Eigen::VectorXd log_det_GRN_vec;
  if(p_nonGRN == 0){
    log_det_GRN_vec = log_det_vec;
  }else{
    log_det_GRN_vec = log_det_GRN(Theta_mat_list,celltypes_num,p_GRN,core_min_use);
  }
  //
  scGeneNet_bic_VMICL = BIC_fun(Theta_mat_list,
                                U_mat,
                                m_mat_list,
                                s_sq_mat_list,
                                obs_mat,
                                mu_mat,
                                S_depth,
                                pi_vec,
                                log_cul_mat_obs_sum_GRN,
                                p_GRN,
                                log_det_GRN_vec,
                                core_min_use);
  //
  Eigen::MatrixXd U_max;
  U_max.setZero(n,celltypes_num);
  for(int i = 0;i<n;i++){
    U_max(i,celltypes_label(i)) = 1;
  }
  scGeneNet_bic_VICL = BIC_fun(Theta_mat_list,
                               U_max,
                               m_mat_list,
                               s_sq_mat_list,
                               obs_mat,
                               mu_mat,
                               S_depth,
                               pi_vec,
                               log_cul_mat_obs_sum_GRN,
                               p_GRN,
                               log_det_GRN_vec,
                               core_min_use);
  
  //
  Eigen::MatrixXi obs_mat_trans_densy = MatrixXi(obs_mat);
  Eigen::MatrixXd obs_mat_trans_densy_double = obs_mat_trans_densy.cast<double>();
  Eigen::SparseMatrix<double> obs_mat_trans_double = obs_mat_trans_densy_double.sparseView();
  //
  if(penalize_diagonal == false){
    return Rcpp::List::create(Named("celltypes_label") = celltypes_label_trans,
                              Named("U_mat") = U_mat_trans,
                              Named("mu_mat") = mu_mat_trans,
                              Named("m_mat_list") = m_mat_list_trans,
                              Named("s_sq_mat_list") = s_sq_mat_list_trans,
                              Named("pi_vec") = pi_vec_trans,
                              Named("obs_mat") = obs_mat_trans_double,
                              Named("S_depth") = S_depth_trans,
                              Named("gene_GRN_index_use") = gene_GRN_index_use_trans,
                              Named("Theta_mat_list") = Theta_mat_list_trans,
                              Named("zero_GRN") = zero_GRN_trans,
                              Named("zero_GRN_use") = zero_GRN_use_trans,
                              Named("log_det_vec") = log_det_vec_trans,
                              Named("scGeneNet_bic_VMICL") = scGeneNet_bic_VMICL,
                              Named("scGeneNet_bic_VICL") = scGeneNet_bic_VICL);
  }else{
    return Rcpp::List::create(Named("celltypes_label") = celltypes_label_trans,
                              Named("U_mat") = U_mat_trans,
                              Named("mu_mat") = mu_mat_trans,
                              Named("m_mat_list") = m_mat_list_trans,
                              Named("s_sq_mat_list") = s_sq_mat_list_trans,
                              Named("pi_vec") = pi_vec_trans,
                              Named("obs_mat") = obs_mat_trans_double,
                              Named("S_depth") = S_depth_trans,
                              Named("gene_GRN_index_use") = gene_GRN_index_use_trans,
                              Named("Theta_mat_list") = Theta_mat_list_trans,
                              Named("zero_GRN") = zero_GRN_trans,
                              Named("zero_GRN_use") = zero_GRN_use_trans,
                              Named("log_det_vec") = log_det_vec_trans,
                              Named("scGeneNet_bic_VMICL") = scGeneNet_bic_VMICL,
                              Named("scGeneNet_bic_VICL") = scGeneNet_bic_VICL,
                              Named("Theta_mat_list_pd") = Theta_mat_list);
  }
  
}
