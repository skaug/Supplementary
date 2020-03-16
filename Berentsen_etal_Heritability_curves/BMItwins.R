options(warn=-1) #Supress warnings in pdf
library(TMB)      # For estimation
library(ggplot2)  # For figures 
library(egg)      # For figures
library(OpenMx)   # For twindata
source("utils.R") # Utility functions (transformation, plotting etc.)

# Load Data
data(twinData)

# Select Variables for Analysis
selVars   <- c('bmi1','bmi2')

# Select Data for Analysis
mzData    <- subset(twinData, zyg==1, selVars)
dzData    <- subset(twinData, zyg==3, selVars)

y_mz <- as.matrix(mzData[complete.cases(mzData), ])
y_dz <- as.matrix(dzData[complete.cases(dzData), ])

n<-dim(y_mz)[1]+dim(y_dz)[1]

# Compile and load cpp file with the first method
compile("BMItwins.cpp")
dyn.load(dynlib("BMItwins"))


#### Model with two m=2 mixture-components ###############

# Initial values
m <- 4                        # Number of components
alpha <- c(3.05, rep(0.006,m-1))      # mu = (mu1, mu2)
sigma <- rep(2,m)              # sigma = (sigma1, sigma2)
rhoMZ <- rep(0.8, m)          # correlations MZ
rhoDZ <- rep(0.4, m)          # correlations DZ
p <- rep(1,m)/m              # mixture probabilities
tdelta <- delta.n2w(m, p) # working mixing probabilities
y_points <- seq(19.8, 24, .1)      # points to evaluate herability curve



# Estimation
map = list(log_sigma=as.factor(rep(1,m)))
parameters <- list(alpha = alpha, log_sigma  = log(sigma), rhoMZ = rhoMZ, rhoDZ = rhoDZ, tdelta = tdelta)
dat <- list(y_mz = y_mz, y_dz = y_dz, m = m, y_points = y_points)
model <- MakeADFun(data = dat, parameters = parameters, map=map, DLL="BMItwins", silent = T)
fit <- nlminb(model$par, model$fn, model$gr,
              control = list(iter.max = 7000, eval.max = 3000))

adrep <- summary(sdreport(model, par.fixed = fit$par), "report")

# Parameters
print("Parameter estimates for m = 2:")
round(adrep[1:12,], 4)

# AIC
print("AIC for m = 2:")
AICm2 <-2*length(fit$par) + 2*fit$objective
AICm2

#BIC
print("BIC for m=2")
BICm2<- log(n)*length(fit$par) + 2*fit$objective
BICm2


print("model with 2 mixture components preferable according to BIC")
##################### Plots ##############################

##evaluate the global coefficients
ade.global <- ade_global(y_mz, y_dz)

# plot of the correlation curves
pcor <- ggplot.corMZDZ(MakeADFun.obj=model, adreport.obj=adrep,
                       rMZ = ade.global$r.mz, rDZ = ade.global$r.dz, x.lab = "BMI")
pcor




#Plot of the heritability, dominant component and environment curves
p.her<-ggplot.curves(MakeADFun.obj = model, adreport.obj = adrep, member="her_curve",
                     global.value  = ade.global$a2,  alpha = 0.95, x.lab = "BMI") +
  theme(plot.margin = unit(c(1,0,1,0), "lines"))+theme_Publication()+ ylab("heritability")
p.dom<-ggplot.curves(MakeADFun.obj = model, adreport.obj = adrep, member="dom_curve",
                     global.value  = ade.global$d2,  alpha = 0.95, x.lab = "BMI") +
  theme(plot.margin = unit(c(1,0,1,0), "lines"))+theme_Publication()+ ylab("dominant component")
p.env<-ggplot.curves(MakeADFun.obj = model, adreport.obj = adrep, member="env_curve",
                     global.value = ade.global$e2,  alpha = 0.95, x.lab = "BMI")+
  theme(plot.margin = unit(c(1,0,0,0), "lines"))+theme_Publication()+ ylab("non-shared environment")

p.her.dom.env <- ggarrange( p.her, p.dom, p.env, ncol = 1, nrow = 3 )
p.her.dom.env




