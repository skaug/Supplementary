library(TMB)
library(ISLR)
library(kerndwd)
library(graphics)
library(ggplot2)

# Compile and link the C++ into the R session
compile("RVM.cpp") 
dyn.load(dynlib("RVM"))

# Function that fit the RVM and calculate the prediction f(x)
RVM <- function (delta, nu, B_size, train_x, test_x, train_t, clas_ind) {
  
    ML = numeric()                 
    N = length(train_t)
    active = rep(FALSE,(N + 1))
    
    #  Step 2. in ALgorithm 1: set initial values   
    log_sigma2 = 0
    log_gamma = 0
    log_alpha = rep(5,(N + 1))
    w = rep(0, (N + 1)) 
    
    
    for (iter in 1:1e4) {              
      # MAIN LOOP
      
        #  Step 3. (and 7) in Algorithm 1:                                     
                                                                                          
        dat = list(t = train_t, x = train_x, clas_ind = clas_ind, A = 1:N, x_pred = test_x)                   
        map = list(log_alpha = factor(rep(NA,(N +1 ))), log_gamma = factor(NA), log_sigma2 = factor(NA))
        
        obj_w <- MakeADFun(data = dat,
                           parameters = list(w=w, log_alpha=log_alpha, log_gamma=log_gamma, log_sigma2=log_sigma2),
                           silent = TRUE,
                           map = map,
                           DLL = "RVM")

        
        fit_w <- nlminb(obj_w$par, obj_w$fn, obj_w$gr, control=list(iter.max = 1000,eval.max = 1000,rel.tol = 1e-10))  
        fit_tmp=fit_w$par
        
        if (iter ==1){
          
        # Step 4. in Algorithm 1:   
          num_start = B_size           
          ix = order(abs(fit_tmp))
          active[ix[(length(ix) - num_start):length(ix)]] = TRUE
          w[active] = obj_w$env$report()$w[active]
          log_alpha[active] = 0  
          
        }
        
        # Step 9. in Algorithm 1: 
        else{
          active[abs(fit_tmp) > delta  & ! active] = TRUE
        }
        
        
        n = sum(active)
        map_na = rep(NA,(N +1))
        map_na[active] = 1:n
        
        if(clas_ind==1){  # 1: classification 
            
            map = list(log_alpha = factor(map_na), w = factor(n+map_na),log_sigma2 = factor(NA))
            lower = c(rep(-5,sum(active)+1))
            upper = c(rep(nu,sum(active)+1)) 
            
          }
        else{  # 0: regression
    
            map = list(log_alpha = factor(map_na), w = factor(n+map_na))
            lower = c(rep(-5,sum(active)+2))
            upper = c(rep(nu,sum(active)+2))
          } 
           
     
        # Step 5 in Algorithm 1: 
        dat = list(t = train_t, x = train_x, clas_ind = clas_ind, A = which(active[-1]), x_pred = test_x)
        
        obj_2 <- MakeADFun(data = dat,
                           parameters = list(w=w,log_alpha=log_alpha,log_gamma=log_gamma,log_sigma2=log_sigma2),
                           random = "w", 
                           silent = TRUE,
                           map = map,
                           DLL = "RVM")
        
        fit <- nlminb(obj_2$par, obj_2$fn, obj_2$gr,lower=lower,upper=upper, control=list(iter.max=1000,eval.max=1000,rel.tol=1e-10))  
        
        
        ML[iter] = fit$objective
        w[active] = obj_2$env$report()$w[active]
        log_alpha[active] = fit$par[1:n]
        log_gamma = fit$par[(n+1)]
        log_sigma2 = fit$par[(n+2)]
          
         #Step 6 in Algorithm 1: 
         active[log_alpha ==  nu & active] = FALSE 
          
         #  SD evaluation
         rep_SD = sdreport(obj_2,ignore.parm.uncertainty = FALSE)
         sum_rep = summary(rep_SD)
         f = sum_rep[row.names(sum_rep)=="f",]
          
         #             Checking if terminates
         if ((iter > 2) && (abs(ML[iter] - ML[iter-1]) < 1e-6)) {   
            break;
          }
        
      }  
    
#   END OF MAIN LOOP
    
    return(list( "f" = f[,1],"f_sd" = f[,2], "RVindex"= which(active), "sigma2"= exp(log_sigma2) ))
   
}


Dat = Auto
y = Dat$mpg
x = Dat[5]

#    Scaling
y = y[order(as.numeric(scale(x)))]
x = sort(as.numeric(scale(x[,1])))


train_x = as.matrix(x[seq(1:300)])
test_x = as.matrix(x[seq(1:392)])
train_y = y[seq(1:300)]
test_y = y[seq(1:392)]


#  parameters from algorithm 1
delta = 0.001
nu = 5
B_size = 60
clas_in = 0  


RVM_auto = RVM(delta, nu, B_size, train_x, test_x, train_y, clas_in)


## C.I.
f_U = RVM_auto$f + 1.96* RVM_auto$f_sd
f_L = RVM_auto$f - 1.96* RVM_auto$f_sd
f = RVM_auto$f


##   PRED I.
f_Up = RVM_auto$f + 1.96*sqrt( RVM_auto$f_sd^2 + RVM_auto$sigma2)
f_Lp = RVM_auto$f - 1.96*sqrt( RVM_auto$f_sd^2 + RVM_auto$sigma2)





RV = RVM_auto$RVindex  
RV = RV[!RV==1] # remove intercept

xi = train_x[RV]
yi = train_y[RV]

dat = data.frame(test_x,test_y, f_U,f_L,f)
d1 = data.frame(xi,yi)
d2 = data.frame(train_x,train_y)
d3 = data.frame(test_x[301:392], test_y[301:392])


plot1 = ggplot(dat) + 
  ylab("mpg") +
  xlab("Mass")+  
  theme(text=element_text(size=25), ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))+
  geom_ribbon(aes(x = test_x, ymax = f_U, ymin = f_L), alpha = 0.9, fill = "skyblue")+
  geom_ribbon(aes(x = test_x, ymax = f_Up, ymin = f_Lp), alpha = 0.3, fill = "skyblue")+
  geom_line(aes(test_x, f), group = 1, col = "black") +
  geom_point(aes(train_x,train_y),col = "black", shape = 19, size = 1.1, data = d2) +  
  geom_point(aes(xi,yi), col = "#FFB90F", shape = 19, size = 2, data = d1) +    
  geom_point(aes(test_x[301:392], test_y[301:392]), col = "red", shape = 19, size = 1.1, data = d3) 
auto_plot = plot1 + theme(panel.background = element_rect(fill = "white", colour = "grey50"))
print(auto_plot)










