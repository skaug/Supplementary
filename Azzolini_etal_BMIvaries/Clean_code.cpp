// A mixture-model for twin-data where the only restriction is equal marginals
// There is a gender effect in the mean vector

#include <TMB.hpp>


// Utility function which maps unconstrained tdelta into constrained p via logistic transformation
template<class Type>
vector<Type> delta_w2n(int m, vector<Type> tdelta){  
  vector<Type> p(m);
  p(0) = Type(1);       // set first element to one
  if(m>1){
    p.tail(m - 1) = exp(tdelta); // Fill in the last m-1 elements with working parameters
    p = p/p.sum(); // normalize
  }
  return p;
}

// Function which calculates the correlation curve at point "value"
template<class Type>
Type corr_curve(int m, Type value, vector<Type> mu, vector<Type> sigma, vector<Type> rho, vector<Type> p, Type marginal_sd){
  //Components of the correlation
  vector<Type> CV_vec(m);   // elements of conditional variance 
  vector<Type> CM_vec(m);   // elements of conditional mean 
  vector<Type> Beta_A_vec(m); //elements the of derivative of the conditional mean 
  vector<Type> Beta_B_vec(m); //elements of the derivative of the conditional mean
  vector<Type> Beta_C_vec(m); //elements of the derivative of the conditional mean

  // create conditional probabilities P(p | Y_2 = y) = p*
  vector<Type> p_star(m);
  vector<Type> pk_norm(m);

  for(int k = 0; k < m; k++){
    pk_norm(k) = p(k)*dnorm(value, mu(k), sigma(k), false);
  }

  for(int k = 0; k < m; k++){
    p_star(k) = pk_norm(k)/pk_norm.sum();  // Berentsen et al, equation 22
  }

  // create the conditional mean of mu - rho(y - mu)  with respect to p* ---  Berentsen et al, equation 23  
  for(int k = 0; k < m; k++){
    CM_vec(k) = p_star(k)*(mu(k) + rho(k)*(value - mu(k)));  
  }

  Type CM = CM_vec.sum();

  // Create conditional variance --- Berentsen et al, equation 24
  for(int k = 0; k < m; k++){
    CV_vec(k) = sigma(k)*sigma(k)*(1 - rho(k)*rho(k))*p_star(k) + p_star(k)*(mu(k) + rho(k)*(value - mu(k)) - CM)*(mu(k) + rho(k)*(value - mu(k)) - CM);
  } 
  
 
  Type CV = CV_vec.sum();
  
 
  // Create derivative of conditional mean ---  Berentsen et al, equation 25
  for(int k = 0; k < m; k++){
    Beta_A_vec(k) = p_star(k)*(rho(k) - (mu(k) + rho(k)*(value - mu(k)))*(value - mu(k))/(sigma(k)*sigma(k)));
    Beta_B_vec(k) = p_star(k)*(mu(k) + rho(k)*(value - mu(k)));
    Beta_C_vec(k) = p_star(k)*(value - mu(k))/(sigma(k)*sigma(k));
  }
  
  Type Beta = Beta_A_vec.sum() + Beta_B_vec.sum()*Beta_C_vec.sum();


  // correlation curve at point "value" ---  Berentsen et al, equation 9
  Type sigma_Beta = marginal_sd*Beta;
  Type sigma_Beta_squared = sigma_Beta*sigma_Beta;
  Type cor_value = sigma_Beta/sqrt(sigma_Beta_squared + CV); 
  return(cor_value);
}

// negative log likelihood 

template<class Type>
Type objective_function<Type>::operator() ()
{
  // IMPORT DATA
  DATA_MATRIX(y_mz);     // Manifested variables MZ-twins (twin 1, 1st column, twin 2, 2nd column)
  DATA_MATRIX(y_dz);     // Manifested variables DZ-twins (twin 1, 1st column, twin 2, 2nd column)
  DATA_INTEGER(m);       // Number of mixing-components
  DATA_VECTOR(y_points); // Points to calculate heratibility curve in
  DATA_IVECTOR(genderMZ);    //Tells if the genderMZ is 1 or 2
  DATA_IVECTOR(genderDZ);    //Tells if the genderDZ is 1 or 2
  
  // READ THE PARAMETERS  
  PARAMETER_VECTOR(alpha);	                                                 // Vector of marginal means mu = (mu1, mu2, ..., mum)
  PARAMETER_VECTOR(log_sigma);  vector<Type> sigma = log_sigma.exp();    // Vector of marginal sds sigma = (sigma1, sigma2, ..., sigmam) 
  PARAMETER_VECTOR(rhoMZ);                                               // Correlations in each kernel MZ (rhoMZ1, rhoMZ2, ..., rhoMZm)
  PARAMETER_VECTOR(rhoDZ);                                               // Correlation in each kernel DZ (rhoDZ1, rhoDZ2, ..., rhoDZm)    
  PARAMETER_VECTOR(tdelta);   vector<Type> p = delta_w2n(m, tdelta); // Mixing probability
  PARAMETER(Beta_gender);                                                // covariate effect of gender on the mean vector
 
  
  // creation of the vector mu from alpha
  vector<Type> mu(m);
  mu(0)=exp(alpha(0));
  for (int i=1; i<m; i++){
    
    mu(i)=mu(i-1) + exp(alpha(i));
  }
  
  // reporting the parameters back to R
  ADREPORT(mu);
  ADREPORT(sigma);
  ADREPORT(rhoMZ);
  ADREPORT(rhoDZ);
  ADREPORT(p);
  ADREPORT(Beta_gender);
  
  
  // Misc quantities 
  int n_mz = y_mz.rows();               // No. of MZ-twins
  int n_dz = y_dz.rows();               // No. of DZ-twins
  vector<Type> sigma2 = sigma*sigma;    // variances
  
  // mean vector for male data
  vector<Type> mu_M(m);
  for (int i=0; i<m; i++){
    
    mu_M(i)=mu(i) + 0.5*Beta_gender;
  }
  
  // mean vector for female data
  vector<Type> mu_F(m);
  for (int i=0; i<m; i++){
    
    mu_F(i)=mu(i) - 0.5*Beta_gender;
  }
  
  
  //Create uncentered mixing densities and mean vectors
  vector<vector<Type> > mu_vec_M(m);
  vector<vector<Type> > mu_vec_F(m);
  matrix<Type> covMZ(2,2);
  matrix<Type> covDZ(2,2);
  matrix<Type> I(2,2);
  matrix<Type> U(2,2);
  U = U.setOnes();
  I = I.setIdentity();
  matrix<Type> V = U - I; // off-diagonal matrix with 1's
  
  using namespace density;
  vector<MVNORM_t<Type> > mvdnorm_MZ(m);
  vector<MVNORM_t<Type> > mvdnorm_DZ(m);
  
  
  // Set up the components of the Gaussian mixture
  for(int i = 0; i < m; i++){
    //covariance matrix of component i 
    covMZ =  sigma2(i)*I + sigma2(i)*rhoMZ(i)*V;
    covDZ =  sigma2(i)*I + sigma2(i)*rhoDZ(i)*V;
    
    
    //mean vector of component i for males 
    vector<Type> mu_vec_M_i(2);
    mu_vec_M_i(0) = mu_M(i);
    mu_vec_M_i(1) = mu_M(i);
    mu_vec_M(i) = mu_vec_M_i;
    
    //mean vector of component i for females
    vector<Type> mu_vec_F_i(2);
    mu_vec_F_i(0) = mu_F(i);
    mu_vec_F_i(1) = mu_F(i);
    mu_vec_F(i) = mu_vec_F_i;
    
    
    // uncentered density of component i 
    mvdnorm_MZ(i) = MVNORM_t<Type>(covMZ);
    mvdnorm_DZ(i) = MVNORM_t<Type>(covDZ);
  }
  
  
 // ------ Evaluate negative log-likelihood -log(L) ----------------
 
 // Define negative log likelihood
 Type nll=0;			
 
 // Contribution MZ-twins
 for (int i = 0; i < n_mz; i++)
 {
   
   // Twin pair no.i
   vector<Type> Y(2);
   Y(0) = y_mz(i,0);
   Y(1) = y_mz(i,1);
   
   
   // Evaluate mixture for obs no. i
   // In the loop, we take into consideration the gender of the twin (1= male, 2=female)
   Type prob_i;
   
   if (genderMZ(i) ==1) {
     // loop over mixture components
     for(int j = 0; j < m; j++){
       
       prob_i += p(j)*exp(-mvdnorm_MZ(j)(Y - mu_vec_M(j)));
     }
   }
   
   else if (genderMZ(i) ==2) {
     
     for(int j = 0; j < m; j++){
       
       prob_i += p(j)*exp(-mvdnorm_MZ(j)(Y - mu_vec_F(j)));
     }
     
     
   }
   
   // Contribution
   nll -= log(prob_i);
   
   
 }
 
 
 // Contribution DZ-twins
 for (int i = 0; i < n_dz; i++)
 {
   
   // Twin pair no.i
   vector<Type> Y(2);
   Y(0) = y_dz(i,0);
   Y(1) = y_dz(i,1);
   
   
   
   // Evaluate mixture for obs no. i
   // In the loop, we take into consideration the gender of the twin
   Type prob_i;
   
   
   
   if (genderDZ(i) ==1) {
     

     for(int j = 0; j < m; j++){
       
       prob_i += p(j)*exp(-mvdnorm_DZ(j)(Y - mu_vec_M(j)));
     }
   }
   
   else if (genderDZ(i) ==2) {
     
     for(int j = 0; j < m; j++){
       
       prob_i += p(j)*exp(-mvdnorm_DZ(j)(Y - mu_vec_F(j)));
     }
     
     
   }
   
   // Contribution
   nll -= log(prob_i);
   
 } // end nll loop
 
 // ------ Evaluate heritability curves --------------------- 
 

   
  int n_points = y_points.size();
  //Curves for male data
  vector<Type> cor_curve_MZ_M(n_points); 
  vector<Type> cor_curve_DZ_M(n_points);
  vector<Type> her_curve_ACE_M(n_points);
  vector<Type> her_curve_ADE_M(n_points);
  vector<Type> env_curve_M(n_points);
  vector<Type> c_env_curve_M(n_points);
  vector<Type> dom_curve_M(n_points);

  
  //Curves for female data
  vector<Type> cor_curve_MZ_F(n_points);
  vector<Type> cor_curve_DZ_F(n_points);
  vector<Type> her_curve_ACE_F(n_points);
  vector<Type> her_curve_ADE_F(n_points);
  vector<Type> env_curve_F(n_points);
  vector<Type> c_env_curve_F(n_points);
  vector<Type> dom_curve_F(n_points);

  
  //Marginal variance for male data
  Type mu_bar_M = (p*mu_M).sum();                                           
  vector<Type> mu_minus_mu_bar_M = mu_M - mu_bar_M;
  vector<Type> mu_minus_mu_bar_squared_M = mu_minus_mu_bar_M*mu_minus_mu_bar_M;
  Type marginal_variance_M = (p*sigma2).sum() + (p*mu_minus_mu_bar_squared_M).sum();
  Type marginal_sd_M = sqrt(marginal_variance_M);
  

  //Marginal variance for female data
  Type mu_bar_F = (p*mu_F).sum();                                           
  vector<Type> mu_minus_mu_bar_F = mu - mu_bar_F;
  vector<Type> mu_minus_mu_bar_squared_F = mu_minus_mu_bar_F*mu_minus_mu_bar_F;
  Type marginal_variance_F = (p*sigma2).sum() + (p*mu_minus_mu_bar_squared_F).sum();
  Type marginal_sd_F = sqrt(marginal_variance_F);
  
  
  // Loop for evaluating the curves at every point in the vector y_points
    for(int i = 0; i < n_points; i++){
   
  //curves for male data
  
    Type cor_curve_MZ_i_M= corr_curve(m, y_points(i), mu_M, sigma, rhoMZ, p, marginal_sd_M);
    cor_curve_MZ_M(i)=cor_curve_MZ_i_M;
    Type cor_curve_DZ_i_M= corr_curve(m, y_points(i), mu_M, sigma, rhoDZ, p, marginal_sd_M);
    cor_curve_DZ_M(i)=cor_curve_DZ_i_M;
    her_curve_ACE_M(i) = 2*(cor_curve_MZ_i_M - cor_curve_DZ_i_M); // Berentsen et al, eq. 14
    c_env_curve_M(i)= 2*cor_curve_DZ_i_M - cor_curve_MZ_i_M;      // Berentsen et al, eq. 14
    her_curve_ADE_M(i) = 4*cor_curve_DZ_i_M - cor_curve_MZ_i_M;   // Berentsen et al, eq. 11
    dom_curve_M(i)= 2*cor_curve_MZ_i_M - 4*cor_curve_DZ_i_M;      // Berentsen et al, eq. 12
    env_curve_M(i) =1 - her_curve_ACE_M(i) - c_env_curve_M(i);    // Berentsen et al, eq. 13
   
   
  //curves for female data
  

    Type cor_curve_MZ_i_F= corr_curve(m, y_points(i), mu_F, sigma, rhoMZ, p, marginal_sd_F);
    cor_curve_MZ_F(i)=cor_curve_MZ_i_F;
    Type cor_curve_DZ_i_F= corr_curve(m, y_points(i), mu_F, sigma, rhoDZ, p, marginal_sd_F);
    cor_curve_DZ_F(i)=cor_curve_DZ_i_F;
    her_curve_ACE_F(i) = 2*(cor_curve_MZ_i_F - cor_curve_DZ_i_F); // Berentsen et al, eq. 14
    c_env_curve_F(i)= 2*cor_curve_DZ_i_F - cor_curve_MZ_i_F;      // Berentsen et al, eq. 14
    her_curve_ADE_F(i) = 4*cor_curve_DZ_i_F - cor_curve_MZ_i_F;   // Berentsen et al, eq. 11
    dom_curve_F(i)= 2*cor_curve_MZ_i_F - 4*cor_curve_DZ_i_F;      // Berentsen et al, eq. 12
    env_curve_F(i) =1 - her_curve_ACE_F(i) - c_env_curve_F(i);    // Berentsen et al, eq. 13
   } // end of iteration over y_points
   
   // reporting the curves back to R
   ADREPORT(cor_curve_MZ_M);
   ADREPORT(cor_curve_DZ_M);
   ADREPORT(cor_curve_MZ_F);
   ADREPORT(cor_curve_DZ_F);
   ADREPORT(her_curve_ACE_M);
   ADREPORT(her_curve_ADE_M);
   ADREPORT(her_curve_ACE_F);
   ADREPORT(her_curve_ADE_F);
   ADREPORT(env_curve_M);
   ADREPORT(env_curve_F);
   ADREPORT(c_env_curve_M);
   ADREPORT(c_env_curve_F);
   ADREPORT(dom_curve_M);
   ADREPORT(dom_curve_F);

   return nll;

}

