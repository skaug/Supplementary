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





##Formulas for the global coefficients
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
  d2 <- abc[2] # dominance
  e2 <- abc[3] # residual environmental
  
  # variance of a2 by delta-method
  sd.r.mz <- sqrt((1 - r.mz^2)/length(y.mz))
  sd.r.dz <- sqrt((1 - r.dz^2)/length(y.dz))
  
  res <- list(r.mz = r.mz, r.dz = r.dz, ybar = ybar, yvar = yvar, a2 = a2, d2 = d2, e2 = e2, sd.r.mz = sd.r.mz, sd.r.dz = sd.r.dz)
  return(res)
}



# function for plotting heritability, dominant component and environment curves
ggplot.curves <- function(MakeADFun.obj, adreport.obj, global.value, member,  alpha = 0.95, x.lab = "response"){
  
  par <- adreport.obj[rownames(adreport.obj) == member, ]
  curve<- par[,1]
  lower <- curve - qnorm(alpha)*par[,2]
  upper <- curve + qnorm(alpha)*par[,2]
  y <- MakeADFun.obj$env$data$y_points
  dat <- c(MakeADFun.obj$env$data$y_mz, MakeADFun.obj$env$data$y_dz) 
  
  # Create df for ggplot
  df <- data.frame(y = y, curve=curve, lower = lower, upper = upper)
  
  # Construct plot
  p <- ggplot(data = df, aes(x = y, y = curve)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_hline(aes(yintercept = global.value), linetype = 2, col = "red") +
    geom_vline(xintercept = quantile(dat, probs = c(0.05, 0.95), na.rm =T), col=3) +
    xlab(x.lab)
  p <- p +  theme_Publication()
  p
}



###function to plot MZ and DZ correlation curves

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
  
  # Create df for ggplot
  dfMZ <- data.frame(y = y, correlation = correlationMZ, lower = lowerMZ, upper = upperMZ,
                     zygosity = factor(rep("MZ", length(y))))
  dfDZ <- data.frame(y = y, correlation = correlationDZ, lower = lowerDZ, upper = upperDZ,
                     zygosity = factor(rep("DZ", length(y))))
  df <- rbind(dfMZ, dfDZ)
  # Construct plot
  p <- ggplot(data = df, aes(x = y, y = correlation)) +
    geom_line(aes(color = zygosity)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, group = zygosity), alpha = 0.3) +
    geom_hline(aes(yintercept = rDZ), linetype = 2, color = "red") +
    geom_hline(aes(yintercept = rMZ), linetype = 2, color = "red") +
    xlab(x.lab) +
    scale_color_manual(values=c("green", "blue")) +
    theme_Publication()
  p
}
  

ggplot.cor.mc.fc <- function(MakeADFun.obj, adreport.obj, cor.global.mc, cor.global.fc,  alpha = 0.95, x.lab = "response"){
  
  par.mc <- adreport.obj[rownames(adreport.obj) == "cor_curve_MC", ]
  par.fc <- adreport.obj[rownames(adreport.obj) == "cor_curve_FC", ]
  correlation.mc <- par.mc[,1]
  correlation.fc <- par.fc[,1]
  lower.mc <- correlation.mc - qnorm(alpha)*par.mc[,2]
  upper.mc <- correlation.mc + qnorm(alpha)*par.mc[,2]
  lower.fc <- correlation.fc - qnorm(alpha)*par.fc[,2]
  upper.fc <- correlation.fc + qnorm(alpha)*par.fc[,2]
  
  y <- MakeADFun.obj$env$data$y_points
  dat <- c(MakeADFun.obj$env$data$y[,1], MakeADFun.obj$env$data$y[,2])
  
  
  # Create df for ggplot
  df <- data.frame(y = c(y, y), correlation = c(correlation.mc, correlation.fc),
                   lower = c(lower.mc, lower.fc), upper = c(upper.mc, upper.fc), 
                   relationship = factor(rep(c("MC", "FC"), each = length(y))))
  
  df2 <- data.frame(global = c(cor.global.mc, cor.global.fc), relationship = factor(c("MC", "FC")))
  
  # Construct plot
  p <- ggplot(data = df, aes(x = y, y = correlation)) +
    geom_line(aes(color = relationship)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, group=relationship), alpha = 0.3) +
    geom_hline(data = df2, mapping = aes(yintercept = global, color = relationship), linetype = 2) +
    geom_vline(xintercept = quantile(dat, probs = c(0.05, 0.95), na.rm =TRUE), col=3) +
    xlab(x.lab)
  p <- p +  theme_Publication()
  p
}

