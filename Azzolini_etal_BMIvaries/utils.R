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
  if(m>1){
    tdelta <- as.vector(tail(foo, m - 1))
  } else{
    tdelta <- 0 # When m==1 we still need to pass a vector to the C++ code
  }
  return(tdelta) 
}

delta.w2n <- function(m, tdelta){
  
  # set first element to one and fill in the last m - 1 elements with working parameters and take exp
  delta <- c(1, exp(tdelta))
  
  # normalize
  delta = delta/sum(delta)
  
  return(delta)
}


