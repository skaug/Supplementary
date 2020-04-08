###Publication theme
theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = rel(0.8)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text()
    ))
  
}

# Function to transform natural parameters to working

delta.n2w <- function(m, delta){

  foo <- log(delta/delta[1])
  tdelta <- as.vector(tail(foo, m - 1))
  return(tdelta) 
}

delta.w2n <- function(m, tdelta){
  
  # set first element to one and fill in the last m - 1 elements with working parameters and take exp
  delta <- c(1, exp(tdelta))
  
  # normalize
  delta = delta/sum(delta)
  
  return(delta)
}


##Formulas for the global coefficients  rhoMZ, rhoDZ, a^2, d^2, and e^2
ade_global<-function(y.mz, y.dz){
  r.mz <-cor(y.mz)[1,2]  
  r.dz <- cor(y.dz)[1,2]
  ybar <- mean(rbind(y.mz,y.dz))
  yvar <- var(rbind(y.mz,y.dz))[1,2]
  
  
  # Global use of ade formula
  W <- c(r.mz, r.dz, 1)
  Q <- rbind(c(1, 1, 0), c(.5, 0.25, 0), c(1, 1, 1))
  abc <- solve(Q, W)
  a2 <- abc[1] # heratibility
  d2 <- abc[2] # dominant component
  e2 <- abc[3] # non-shared environment
  
  # variance of a2 by delta-method
  sd.r.mz <- sqrt((1 - r.mz^2)/length(y.mz))
  sd.r.dz <- sqrt((1 - r.dz^2)/length(y.dz))
  
  res <- list(r.mz = r.mz, r.dz = r.dz, ybar = ybar, yvar = yvar, a2 = a2, d2 = d2, e2 = e2, sd.r.mz = sd.r.mz, sd.r.dz = sd.r.dz)
  return(res)
}


###function to plot MZ and DZ correlation curves together

ggplot.corMZDZ <- function(MakeADFun.obj, adreport.obj, rDZ, rMZ, alpha = 0.95, x.lab = "response"){
  
  parMZ <- adreport.obj[rownames(adreport.obj) == "cor_curve_MZ", ]
  parDZ <- adreport.obj[rownames(adreport.obj) == "cor_curve_DZ", ]
  correlationMZ <- parMZ[,1]
  correlationDZ <- parDZ[,1]
  lowerMZ <- correlationMZ - qnorm(alpha)*parMZ[,2]
  upperMZ <- correlationMZ + qnorm(alpha)*parMZ[,2]
  lowerDZ <- correlationDZ - qnorm(alpha)*parDZ[,2]
  upperDZ <- correlationDZ + qnorm(alpha)*parDZ[,2]
  y <- MakeADFun.obj$env$data$y_points
  dat <- c(MakeADFun.obj$env$data$y_mz, MakeADFun.obj$env$data$y_dz)
  # Create df for ggplot
  df <- data.frame(y = c(y, y), correlation = c(correlationMZ, correlationDZ),
                   lower = c(lowerMZ, lowerDZ), upper = c(upperMZ, upperDZ), 
                   relationship = factor(rep(c("MZ", "DZ"), each = length(y))))
  
  df2 <- data.frame(global = c(rMZ, rDZ), relationship = factor(c("MZ", "DZ")))
  
  p <- ggplot(data = df, aes(x = y, y = correlation)) +
    geom_line(aes(color = relationship)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, group=relationship), alpha = 0.3) +
    geom_hline(data = df2, mapping = aes(yintercept = global, color = relationship), linetype = 2) +
    geom_vline(xintercept = quantile(dat, probs = c(0.05, 0.95), na.rm =TRUE), col=3) +
    xlab(x.lab)
  p <- p +  theme_Publication()
  p
}
  
###function to plot heritability, dominant component, and non-shared environment together
ggplot.3curves<-function(MakeADFun.obj, adreport.obj, g.1, g.2, g.3, member1, member2, member3, type1, type2, type3, alpha = 0.95, x.lab = "response"){
  y_points <- MakeADFun.obj$env$data$y_points
  q=length(y_points)
  g.values<-c(rep(g.1, q), rep(g.2, q), rep(g.3, q))
  BMI=rep(y_points, 3)
  type=as.factor(c(rep(type1, q), rep(type2, q), rep(type3, q)))
  par1 <- adrep[rownames(adrep) == member1, ]
  curve1<- par1[,1]
  lower1 <- curve1 - qnorm(alpha)*par1[,2]
  upper1 <- curve1 + qnorm(alpha)*par1[,2]
  par2 <- adrep[rownames(adrep) == member2, ]
  curve2<- par2[,1]
  lower2 <- curve2 - qnorm(alpha)*par2[,2]
  upper2 <- curve2 + qnorm(alpha)*par2[,2]
  par3 <- adrep[rownames(adrep) == member3, ]
  curve3<- par3[,1]
  lower3 <- curve3 - qnorm(alpha)*par3[,2]
  upper3 <- curve3 + qnorm(alpha)*par3[,2]
  value<-c(curve1, curve2, curve3)
  
  lower<-c(lower1, lower2, lower3)
  upper<-c(upper1, upper2, upper3)
  dat <- c(MakeADFun.obj$env$data$y_mz, MakeADFun.obj$env$data$y_dz)
  
  
  df<-data.frame(value=value, BMI=BMI, lower=lower, upper=upper, type=type, g.values=g.values)
  p <- ggplot(data = df, aes(x = BMI, y = value)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_hline(aes(yintercept = g.values), linetype = 2, col = "red") +
    geom_vline(xintercept = quantile(dat, probs = c(0.05, 0.95), na.rm =TRUE), col=3) +
    facet_wrap(~type, nrow = 3, scales = "free_y")
  p <- p +  theme_Publication()
  p
  
  
}
