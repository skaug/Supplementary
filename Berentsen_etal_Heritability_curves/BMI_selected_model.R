library(TMB)      # For estimation
library(ggplot2)  # For figures 
library(OpenMx)   # For the dataset "twinData"
source("utils.R") # Utility functions (transformation, plotting etc.)

#### The Dataset -----
data(twinData)                 # Load Data
selVars   <- c('bmi1','bmi2')  # Select Variables for Analysis: bmi for twin 1 and 2

# Select Data for the analysis. We choose only zyg=1 and zyg=3, that is MZ female twins and DZ female twins
mzData    <- subset(twinData, zyg==1, selVars)
dzData    <- subset(twinData, zyg==3, selVars)
y_mz <- as.matrix(mzData[complete.cases(mzData), ])
y_dz <- as.matrix(dzData[complete.cases(dzData), ])

ade.global <- ade_global(y_mz, y_dz)  ##computation of the global a^2, d^2 and e^2 values


# Compile and load cpp file 
compile("BMItwins.cpp")
dyn.load(dynlib("BMItwins"))


###### Estimation (for m=2) ------

# Initial values
m <- 2                        # Number of mixture components
alpha <- c(3.05, rep(0.006,m-1))      # mu = (exp(alpha1), exp(alpha1)+exp(alpha2))
sigma <- rep(2,m)             # Sigma = (sigma1, sigma2)
rhoMZ <- rep(0.8, m)          # Correlations MZ twins
rhoDZ <- rep(0.4, m)          # Correlations DZ twins
p <- rep(1,m)/m               # Mixture probabilities
tdelta <- delta.n2w(m, p)     # Working mixture probabilities

y_points <- seq(19.8, 24, .1) # Points to evaluate correlation and heritability curves in plots


# Estimation using TMB
parameters <- list(alpha = alpha, log_sigma  = log(sigma), rhoMZ = rhoMZ, rhoDZ = rhoDZ, tdelta = tdelta)
dat <- list(y_mz = y_mz, y_dz = y_dz, m = m, y_points = y_points)
model <- MakeADFun(data = dat, parameters = parameters, DLL="BMItwins", silent = TRUE)
L = c(-5,-5,-2,-2,rep(-.995,4),-5) # Lower bound for nlminb()
U = -L
fit <- nlminb(model$par, model$fn, model$gr,model$he,lower=L,upper=U)
adrep <- summary(sdreport(model, par.fixed = fit$par), "report") # Evaluate standard devitiations

#### Table 2 (reformated)
print(round(adrep[1:10,],2)) # Second column are standard deviations not shown in Table 2


# AIC and BIC values
n <- nrow(y_mz) + nrow(y_dz)
AIC_m2 <- 2*length(fit$par) + 2*fit$objective # AIC
BIC_m2<- log(n)*length(fit$par) + 2*fit$objective# BIC


##### Plots of the curves: You need to uncomment the ggsave() lines to actually save the pdf's

# Correlation curves
p.cors<-ggplot.corMZDZ(MakeADFun.obj=model, adreport.obj = adrep, rMZ = ade.global$r.mz, rDZ = ade.global$r.dz, alpha =  0.95, x.lab = "BMI")
#ggsave("rho_y_twins.pdf", plot=p.cors, width = 13, height=7, units= "cm")

# Heritability, dominant component, non-shared environment
p.curves3<-ggplot.3curves(MakeADFun.obj=model, adreport.obj=adrep, g.1=ade.global$a2, g.2=ade.global$d2, g.3=ade.global$e2, member1="her_curve", 
                             member2="dom_curve", member3="env_curve", type1="heritability", type2="dominant component",
                             type3="non-shared environment", alpha = 0.95, x.lab = "BMI")
#ggsave("curves_plot_twins.pdf", p.curves3, width=13, height=13, units="cm")
