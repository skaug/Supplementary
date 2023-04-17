# RTMB code for "A TMB approach to study spatial variation in weather-generated claims in insurance" by Thorsen et al. (2023) 
# We use a simplified simulated data set as the example

library(RTMB)
library(INLA)
library(rgdal)
library(gdata)

set.seed(1234)

# SPDE approximation of the precision matrix (inverse Matern covariance matrix)
Q_spde <- function(spde, kappa) {  # "spde" is a list of sparse matrices imported from INLA (see below)
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2    ## M0=G0, M1=G1, M2=G2
}

### 1: Generate and plot INLA mesh for location of Norwegian municipalities
data <- read.xls("Kommune_punkt_2014_korrigert.xls", sheet="Komm_punkt_2014", header=T) 
longlat <- cbind(data[,5:6])
colnames(longlat) <- c("Longitude", "Latitude")
coords <- as.matrix(longlat[,1:2])
coords.utm <- project(coords, "+proj=utm +zone=33 ellps=WGS84")/1000
colnames(coords.utm) <- c("Easting", "Northing")
mesh_nor1 <- inla.mesh.2d(loc=coords.utm, max.edge=c(130,200))
m = length(mesh_nor1$idx$loc) 
plot(mesh_nor1, main="")
points(coords.utm, pch=21, bg=1, col="white", cex=1.8)
spde <- inla.spde2.matern(mesh_nor1, alpha=2)$param.inla[c("M0", "M1", "M2")] # Input to Q_spde()
n_s=nrow(spde$M0)

### 2: Simulate data
n = 10000     # Number of observations
beta0 = 0.1   # Intercept
beta1 = 0.2   # Slope of "rainfall"
p = 0.3       # Zero inflation probability
kappa=exp(-3.204) # Taken from Table 2 in Thorsen et al (2023)
sigma2 = 0.398    # Taken from Table 2 in Thorsen et al (2023)
tau=1/sqrt(4*pi*kappa^2*sigma2)  # The SPDE is internally parameterized in terms of "tau" rather than "sigma2"
COV = solve(Q_spde(spde,kappa=kappa)) # Covariance matrix from precission matrix
u = tau^(-2)*chol(COV)%*%rnorm(n_s)      # Simulate latent Gaussian field
dat = NULL
dat$rainfall = rnorm(n,0,0.2)
dat$meshidxloc = sample(x=mesh_nor1$idx$loc,size=n,replace=T)             # Sample location randomly from Municipalities
dat$y = rpois(n,lambda=exp(beta0+beta1*dat$rainfall+u[dat$meshidxloc]))   # Poisson response
dat$y[sample(c(T,F),n,replace=T,prob=c(p,1-p))] = 0                       # Add zero inflation

### Section 3: Specify and fit model in RTMB

# Because dzipois is not yet implemented in RTMB
log_dzipois = function(x,lambda,zip){
  if(x==0){
    logres <- log(zip + (1-zip)*dpois(x, lambda,log=FALSE))
  } 
  else{
    logres <- log(1-zip) + dpois(x, lambda, log=TRUE)
  }
  logres
}

# TMB model. A specification g(u) in Appendix B of the paper.
f <- function(parms) {
  
  # Makes both parameters and data local to the function, so that we can use
  getAll(parms, dat)  
  n <- length(y)
  tau <- exp(log_tau)
  kappa <- exp(log_kappa)

  ## GMRF prior
  Q <- tau^2*Q_spde(spde, kappa)
  nll <- -dgmrf(u, 0, Q, log=TRUE)  # Eqn. (6) in paper

  # Zero inflated Poisson likelihood
  eta = beta0 + beta1*rainfall + u[meshidxloc] # Eqn. (5) in paper
  for(i in 1:n)
    nll = nll - log_dzipois(y[i], exp(eta[i]), p) # Eqn. (4) in paper
    
  # Correlation range
  nu <- 1.0                # nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren
  rho <- sqrt(8*nu)/kappa  # # Eqn. (A.2) in paper
  ADREPORT(rho)
  
  return(nll)
}

# Fit the model using true parameters values as starting values
par_init <- list(beta0=beta0, beta1=beta1, log_tau=log(tau), log_kappa=log(kappa), p=p, u=rep(0.0,n_s))
obj <- MakeADFun(f, par_init, random="u")   # Build TMB model from R code. Argument "random" invokes the Laplace approximation for u
opt <- nlminb(obj$par, obj$fn, obj$gr)      # Fit the model using stand R optimizer
sdreport(obj) # Calculate standard deviations ("rho" can also be obtained)
