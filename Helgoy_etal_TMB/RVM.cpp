// Type II MLE of sparse learning model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(t);  
  DATA_MATRIX(x);
  DATA_INTEGER(clas_ind);  // 1: classification, 0: regression
  DATA_IVECTOR(A);         // Indices of active set 
  DATA_MATRIX(x_pred);     //  x values for prediction of f(x) and associated uncertainty.
                           
  PARAMETER_VECTOR(w);
  PARAMETER_VECTOR(log_alpha);
  PARAMETER(log_gamma);
  PARAMETER(log_sigma2);

  // Transform parameters
  vector<Type> alpha = exp(log_alpha);
  Type gamma = exp(log_gamma);
  Type sigma2 = exp(log_sigma2);

  Type nll = 0;  // Will hold (negative)logarithm of eqn (10) in manuscript

  // Gaussian prior on w
  nll -= dnorm(w,Type(0),sqrt(1.0/alpha),true).sum(); // Eqn (5)
  
  // Calculate Phi*w for Gaussian kernel implemented via dnorm()
  int N = x.rows();               
  int D = x.cols(); 
  Type c_N01 = D*0.5*log(2*PI);  
  vector<Type> Phi_w(N);
  for(int i=0; i<N; i++){ 
     Phi_w(i) = w(0);                 // Intercept
     for(int j=0;j<A.size(); j++){ 
        int J = A[j];
        vector<Type> x_i = x.row(i);
        vector<Type> x_j = x.row(J-1);
        Phi_w(i) += exp(2*gamma*(dnorm(x_i,x_j,Type(1),true).sum()+c_N01))*w(J);  
      }
  }
       
  // Likelihood
  if(clas_ind==1){    // Classification
    vector<Type> p = invlogit(Phi_w);
    nll -= dbinom(t,Type(1),p,true).sum();  // Eqn (17)
  }
  else{               // Regression
    nll -= dnorm(t,Phi_w, sqrt(sigma2),true).sum();   // Eqn (15)
  }

  // Calculate prediction f(x); structurally the same that for Phi*w
  int N_pred = x_pred.rows();     
  vector<Type> f(N_pred);
  for(int i=0; i<N_pred; i++){ 
     f(i) = w(0);                 
     for(int j=0;j<A.size(); j++){ 
        int J = A[j];
        vector<Type> x_i = x_pred.row(i);
        vector<Type> x_j = x_pred.row(J-1);
        f(i) += exp(2*gamma*(dnorm(x_i,x_j,Type(1),true).sum()+c_N01))*w(J); 
      }
  }
  
  REPORT(w);    // Report back to R witout standard deviation
  ADREPORT(f);  // Report back to R with standard deviation
  
  return nll;   // Return negative joint log likelihood to R session
}
