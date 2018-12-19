vonMisesDist <- function(x,mu,k,scaleOutput,scalePar) {
	if (scaleOutput==0) {
		return( exp( k*cos(x-mu) ) / ( 2*pi*besselI( k, 0, expon.scaled = FALSE) ) )
	}
	if (scaleOutput==1) {
		originalFun <- exp( k*cos(x-mu) ) / ( 2*pi*besselI( k, 0, expon.scaled = FALSE) ) 
		oldMax <- max( originalFun )
		oldMin <- min( originalFun )
		newMax <- 1
		newMin <- 0	
		delta <- (newMax - newMin) / ( oldMax - oldMin );
		xScaled <- delta*(originalFun-oldMin) + newMin;
		return( xScaled )
	}
	if (scaleOutput==2) {
		originalFun <- exp( k*cos(x-mu) ) / ( 2*pi*besselI( k, 0, expon.scaled = FALSE) ) 
		oldMax <- max( originalFun )
		oldMin <- min( originalFun )
		newMax <- 1
		newMin <- 0	
		delta <- (newMax - newMin) / ( oldMax - oldMin );
		xScaled <- delta*(originalFun-oldMin) + newMin;
		xScaled <- xScaled - scalePar
		return( xScaled )
	}
}
