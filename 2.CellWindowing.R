###
#Calculate the intersection of the density curves of vectors a and b
###
#'@param a vectors a
#'@param b vectors b
#'@param filename the file name of the density plot
#'result the density plot and the intersection 
densityintersection <- function(a, b, filename){
  xlim = c(min(c(a,b)), max(c(a,b))) 
  df <- merge(
    as.data.frame(density(a, from = xlim[1], to = xlim[2])[c("x", "y")]), 
    as.data.frame(density(b, from = xlim[1], to = xlim[2])[c("x", "y")]),
    by = "x", suffixes = c(".a", ".b")
  )
  df$comp <- as.numeric(df$y.a > df$y.b) 
  df$cross <- c(NA, diff(df$comp))       
  intersection.point <- df[which(df$cross != 0), "x"]
  
  library(ggplot2)
  plotdata <- rbind.data.frame(data.frame(lines=rep("a",length(a)), time=a), 
                               data.frame(lines=rep("b",length(b)),time=b))
  p <- ggplot(data=plotdata, aes(x=time))+
    geom_density(aes(color=lines))+
    geom_vline(xintercept = intersection.point, color = "red")
  ggsave(p, filename = filename)
  
  return(intersection.point)
}

