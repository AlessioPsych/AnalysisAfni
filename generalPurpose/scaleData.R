scaleData <- function(x, newMax, newMin) {

oldMax <- max(x)
oldMin <- min(x)
  
delta <- (newMax-newMin)/(oldMax-oldMin);
xScaled <- delta*(x-oldMin) + newMin;

return(xScaled)

}