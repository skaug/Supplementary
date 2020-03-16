// Analysis (partial) of BMI of twins in Section 4.1 of Berentsen et al. (2020)

#include <TMB.hpp>

// Utility function which maps unconstrained tdelta into constrained delta via logistic transformation
template<class Type>
vector<Type> delta_w2n(int m, vector<Type> tdelta){  
  vector<Type> delta(m);
  delta(0) = Type(1);       // set first element to one
  delta.tail(m - 1) = exp(tdelta); // Fill in the last m-1 elements with working parameters
  delta = delta/delta.sum(); // normalize
  return delta;
}

// Objective function (negative log-likehood of 3-dimensional Gaussian mixture)
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data from R
  DATA_MATRIX(y_mz);     // Trait values for MZ-twins (columns: twin 1,  twin 2)
  DATA_MATRIX(y_dz);     // Trait values for DZ-twins (columns: twin 1,  twin 2)
  DATA_INTEGER(m);       // Number of mixture components
  DATA_VECTOR(y_points); // Points to calculate heratibility curve in
    
  // Parameter values from R
  PARAMETER_VECTOR(alpha);	                                             // Reparameterization of mu (see below)
  PARAMETER_VECTOR(log_sigma);  vector<Type> sigma = log_sigma.exp();    // Vector of marginal sds sigma = (sigma1, sigma2, ..., sigmam) 
  PARAMETER_VECTOR(rhoMZ);                                               // Correlations in each kernel MZ (rhoMZ1, rhoMZ2, ..., rhoMZm)
  PARAMETER_VECTOR(rhoDZ);                                               // Correlation in each kernel DZ (rhoDZ1, rhoDZ2, ..., rhoDZm)    
  PARAMETER_VECTOR(tdelta);   vector<Type> p = delta_w2n(m, tdelta); // Mixture probabilities
  
  // In this code (unlike the paper) the mu's are constrained to be increasing in value
  // This is achived by using an underlying (unconstrained) parameter vector alpha.
  vector<Type> mu(m);
  mu(0)=exp(alpha(0));
  for (int i=1; i<m; i++){
    mu(i)=mu(i-1) + exp(alpha(i));
  }
  
  // Report back to R
  ADREPORT(mu);
  ADREPORT(sigma);
  ADREPORT(rhoMZ);
  ADREPORT(rhoDZ);
  ADREPORT(p);
    
  // Misc quantities 
  int n_mz = y_mz.rows();               // No. of MZ-twins
  int n_dz = y_dz.rows();               // No. of DZ-twins
  vector<Type> sigma2 = sigma*sigma;    // Variances
   
  // Create uncentered mixture densities and mean vectors
  vector<vector<Type> > mu_vec(m);
  matrix<Type> covMZ(2,2);
  matrix<Type> covDZ(2,2);
  matrix<Type> I(2,2);
  matrix<Type> U(2,2);
  U = U.setOnes();
  I = I.setIdentity();
  matrix<Type> V = U - I; // off-diagonal matrix with 1's
  
  using namespace density;
  vector<MVNORM_t<Type> > mvdnorm_MZ(m);  // List of bivariate gaussian densities to be populated below
  vector<MVNORM_t<Type> > mvdnorm_DZ(m);
  
  for(int i = 0; i < m; i++){             // Loop over mixture components
    covMZ =  sigma2(i)*I + sigma2(i)*rhoMZ(i)*V;
    covDZ =  sigma2(i)*I + sigma2(i)*rhoDZ(i)*V;
    
    mvdnorm_MZ(i) = MVNORM_t<Type>(covMZ);  // Zero-mean bivariate gaussian mixture component with covariance covMZ
    mvdnorm_DZ(i) = MVNORM_t<Type>(covDZ);

    //mean vector of i'th mixture component
    vector<Type> mu_vec_i(2);
    mu_vec_i(0) = mu(i);
    mu_vec_i(1) = mu(i);          // The two twins have the same mean
    mu_vec(i) = mu_vec_i;
  }

  
// ------ Evaluate negative log-likelihood -log(L) ----------------

  Type nll=0;
  
  // Contribution from MZ-twins
  for (int i = 0; i < n_mz; i++)
  {
    // Twin pair no.i
    vector<Type> Y(2);
    Y(0) = y_mz(i,0);
    Y(1) = y_mz(i,1);

    // Evaluate mixture density (L_i) for i'th twin pair 
    Type L_i = Type(0);
    for(int j = 0; j < m; j++){
      L_i += p(j)*exp(-mvdnorm_MZ(j)(Y - mu_vec(j)));
    }
    
    nll -= log(L_i);  // Add contribution to objective function
  }
  
  // Contribution DZ-twins
  for (int i = 0; i < n_dz; i++)
  {
    // Twin pair no.i
    vector<Type> Y(2);
    Y(0) = y_dz(i,0);
    Y(1) = y_dz(i,1);
    
    // Evaluate mixture density (L_i) for i'th twin pair 
    Type L_i = Type(0);    
    for(int j = 0; j < m; j++){   
      L_i += p(j)*exp(-mvdnorm_DZ(j)(Y - mu_vec(j)));
    }
    
    // Contribution
    nll -= log(L_i);  // Add contribution to objective function
  }

  
// ------ Evaluate heratibility curves ---------------------

  // Define objects
  int n_points = y_points.size();
  vector<Type> cor_curve_MZ(n_points);    // rho(y)
  vector<Type> cor_curve_DZ(n_points);
  vector<Type> her_curve(n_points);       // a^2(y)
  vector<Type> dom_curve(n_points);       // c^2(y)
  vector<Type> env_curve(n_points);       // e^2(y)  xxx check
  vector<Type> CVMZ(n_points);            // sigma^2(y)
  vector<Type> CVDZ(n_points);
  vector<Type> CMMZ(n_points);            // sigma^2(y)
  vector<Type> CMDZ(n_points); 
  vector<Type> BetaMZ(n_points);          // beta(y)
  vector<Type> BetaDZ(n_points);
  
  //Marginal variance
  Type mu_bar = (p*mu).sum();                                           
  vector<Type> mu_minus_mu_bar = mu - mu_bar;
  vector<Type> mu_minus_mu_bar_squared = mu_minus_mu_bar*mu_minus_mu_bar;
  Type marginal_variance = (p*sigma2).sum() + (p*mu_minus_mu_bar_squared).sum();
  Type marginal_sd = sqrt(marginal_variance);
  
  // Conditional variances and derivative of conditional expectations
  
  for(int i = 0; i < n_points; i++){
      
    vector<Type> CVMZ_vec(m);   // elements of conditional variance MZ
    vector<Type> CVDZ_vec(m);   // elements of conditional variance DC
    vector<Type> CMMZ_vec(m);   // elements of conditional mean MZ
    vector<Type> CMDZ_vec(m);   // elements of conditional mean DZ
    vector<Type> BetaMZ_A_vec(m); //
    vector<Type> BetaMZ_B_vec(m); //
    vector<Type> BetaMZ_C_vec(m); //
    vector<Type> BetaDZ_A_vec(m); //
    vector<Type> BetaDZ_B_vec(m); //
    vector<Type> BetaDZ_C_vec(m); //
       
    // create conditional probabilities P(p | Y_2 = y) = p*
    vector<Type> p_star(m);
    vector<Type> pk_norm(m);
    
    for(int k = 0; k < m; k++){
      pk_norm(k) = p(k)*dnorm(y_points(i), mu(k), sigma(k), false);
    }
    
    for(int k = 0; k < m; k++){
      p_star(k) = pk_norm(k)/pk_norm.sum();
    }
    REPORT(p_star);
    
    // create the conditional mean of mu - rhoMZ(y - mu) and mu - rhoDZ(y - mu) with respect to p*
    
    for(int k = 0; k < m; k++){
      CMMZ_vec(k) = p_star(k)*(mu(k) + rhoMZ(k)*(y_points(i) - mu(k)));
      CMDZ_vec(k) = p_star(k)*(mu(k) + rhoDZ(k)*(y_points(i) - mu(k)));
    }
    
    Type CMMZ_i = CMMZ_vec.sum();
    Type CMDZ_i = CMDZ_vec.sum();
    CMMZ(i) = CMMZ_i;
    CMDZ(i) = CMDZ_i; 
    
    // create elements in the Conditional Variance sums
    for(int k = 0; k < m; k++){
      CVMZ_vec(k) = sigma2(k)*(1 - rhoMZ(k)*rhoMZ(k))*p_star(k) + p_star(k)*(mu(k) + rhoMZ(k)*(y_points(i) - mu(k)) - CMMZ_i)*(mu(k) + rhoMZ(k)*(y_points(i) - mu(k)) - CMMZ_i);
      CVDZ_vec(k) = sigma2(k)*(1 - rhoDZ(k)*rhoDZ(k))*p_star(k) + p_star(k)*(mu(k) + rhoDZ(k)*(y_points(i) - mu(k)) - CMDZ_i)*(mu(k) + rhoDZ(k)*(y_points(i) - mu(k)) - CMDZ_i);
    }
    
    // Create conditional variances
    Type CVMZ_i = CVMZ_vec.sum();
    Type CVDZ_i = CVDZ_vec.sum();
    CVMZ(i) = CVMZ_i;
    CVDZ(i) = CVDZ_i;
    
    
    // Create derivatives of conditional means
    for(int k = 0; k < m; k++){
      BetaMZ_A_vec(k) = p_star(k)*(rhoMZ(k) - (mu(k) + rhoMZ(k)*(y_points(i) - mu(k)))*(y_points(i) - mu(k))/sigma2(k));
      BetaMZ_B_vec(k) = p_star(k)*(mu(k) + rhoMZ(k)*(y_points(i) - mu(k)));
      BetaMZ_C_vec(k) = p_star(k)*(y_points(i) - mu(k))/sigma2(k);
      BetaDZ_A_vec(k) = p_star(k)*(rhoDZ(k) - (mu(k) + rhoDZ(k)*(y_points(i) - mu(k)))*(y_points(i) - mu(k))/sigma2(k));
      BetaDZ_B_vec(k) = p_star(k)*(mu(k) + rhoDZ(k)*(y_points(i) - mu(k)));
      BetaDZ_C_vec(k) = p_star(k)*(y_points(i) - mu(k))/sigma2(k);
    }
    
    Type BetaMZ_i = BetaMZ_A_vec.sum() + BetaMZ_B_vec.sum()*BetaMZ_C_vec.sum();
    Type BetaDZ_i = BetaDZ_A_vec.sum() + BetaDZ_B_vec.sum()*BetaDZ_C_vec.sum();
    BetaMZ(i) = BetaMZ_i;
    BetaDZ(i) = BetaDZ_i;
    
    // correlation curves
    Type sigma_BetaMZ = marginal_sd*BetaMZ_i;
    Type sigma_BetaMZ_squared = sigma_BetaMZ*sigma_BetaMZ;
    Type cor_curve_MZ_i = sigma_BetaMZ/sqrt(sigma_BetaMZ_squared + CVMZ_i);
    cor_curve_MZ(i) = cor_curve_MZ_i;
    
    Type sigma_BetaDZ = marginal_sd*BetaDZ_i;
    Type sigma_BetaDZ_squared = sigma_BetaDZ*sigma_BetaDZ;
    Type cor_curve_DZ_i = sigma_BetaDZ/sqrt(sigma_BetaDZ_squared + CVDZ_i);
    cor_curve_DZ(i) = cor_curve_DZ_i;
    
    // heritability curve
    her_curve(i) = 4*cor_curve_DZ_i - cor_curve_MZ_i;
    
    //dominant genetic curve
    dom_curve(i)= 2*cor_curve_MZ_i - 4*cor_curve_DZ_i;
    
    //environmental curve
    env_curve(i) =1 - her_curve(i) - dom_curve(i);
    
  } // end of iteration over i
  
  ADREPORT(cor_curve_MZ);
  ADREPORT(cor_curve_DZ);
  ADREPORT(her_curve);
  ADREPORT(dom_curve);
  ADREPORT(env_curve);
  ADREPORT(CVMZ);
  ADREPORT(CVDZ);
  ADREPORT(BetaMZ);
  ADREPORT(BetaDZ);
  ADREPORT(CMMZ);
  ADREPORT(CMDZ);
  
  REPORT(CVMZ);
  REPORT(CVDZ);
  REPORT(CMMZ);
  REPORT(CMDZ);
  REPORT(BetaMZ);
  REPORT(BetaDZ);
  REPORT(marginal_sd);
  
  return nll;
  
}