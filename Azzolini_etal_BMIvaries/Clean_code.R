library(tidyverse)
library(TMB)      # For estimation
source("utils.R") # Various useful functions

load("Sim_dataset.RData")
summary(Dat)

# Vector used to estimate the mean quantiles
y<-c(Dat[,1], Dat[,2])

# data matrix used in the estimation process
MZ_dat<-Dat%>%filter(Zygosity=="MZ")
DZ_dat<-Dat%>%filter(Zygosity=="DZ")
v.mz<-as.matrix(MZ_dat[,1:2])
v.dz<-as.matrix(DZ_dat[,1:2])

# gender vector for the estimation process
genderMZ<-MZ_dat$Gender
genderDZ<-DZ_dat$Gender
levels(genderMZ)<-c(1,2)
levels(genderDZ)<-c(1,2)
genderMZ<-as.numeric(genderMZ)
genderDZ<-as.numeric(genderDZ)

#visual representation of the simulated points

p<-ggplot(Dat)+geom_point(aes(x=First, y=Second), alpha=0.6, size=.1)+
  facet_grid(Zygosity~Gender)+theme_Publication()+
  xlim(13,47)+ylim(13, 47)
p

#ggsave(filename = "simulated_data.pdf", plot=p, height=15, width = 15, units="cm")


#Compile the c++ code
compile("Clean_code.cpp")
dyn.load(dynlib("Clean_code"))



# Best model ------------------------------------------
m=3
# Initial values
quant<-quantile(y, c(1/4, 2/4, 3/4))
alpha<-c()                               # for the mean vector 
alpha[1]<-log(quant[[1]])                #exp(alpha[1])=mu[1]
alpha[2]<-log(quant[[2]]-quant[[1]])     #mu[2]=mu[1]+exp(alpha[2])
alpha[3]<-log(quant[[3]]-quant[[2]])     #mu[3]=mu[2]+exp(alpha[3])
sigma <- rep(2,3)                        # sigma = (sigma1, sigma2)
rhoMZ <- rep(0.7, 3)                     # MZ correlations 
rhoDZ <- rep(0.3, 3)                     # DZ correlations 
p_0 <- rep(1,3)/3                          # mixture probabilities
tdelta <- delta.n2w(3, p_0)                # working mixing probabilities
Beta_gender<-1.6                         # gender effect on the mean vector
y_points <- seq(18, 30, length.out=500)  # points to evaluate heritability curve

# Estimation
parameters <- list(alpha = alpha, log_sigma  = log(sigma), rhoMZ = rhoMZ, rhoDZ = rhoDZ, tdelta = tdelta, Beta_gender=Beta_gender)
dat <- list(y_mz = v.mz, y_dz = v.dz, m = m, y_points = y_points, genderMZ=genderMZ, genderDZ=genderDZ)
model <- MakeADFun(data = dat, parameters = parameters, DLL="Clean_code")         

fit <- nlminb(model$par, model$fn, model$gr,model$he,
              control = list(iter.max = 7000, eval.max = 3000)) 
fit
adrep1<-summary(sdreport(model, par.fixed = fit$par), "report")			

#parameter estimates
print(round(adrep1[1:16,],2)) 


#AIC value
2*length(fit$par) + 2*fit$objective 

#BIC value
n<-dim(Dat)[1]
log(n)*length(fit$par) + 2*fit$objective


