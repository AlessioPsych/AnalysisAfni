powerSpectra <- function( inputTs, inputTr ) {
  
  library(pracma)
  tsStepX <- seq(0,length(inputTs)*inputTr,by=inputTr)
  fftOut <- fft( inputTs )
  nt <- length(tsStepX);
  dt <- tsStepX[2]-tsStepX[1];
  dc <- fftOut[1]/nt;
  amp <- 2*abs(fftOut)/nt;
  ph <- -180*angle(fftOut)/pi;
  id <- 2:(ceil(nt/2)+1)
  freq <- (1:length(id)) / (nt*dt);
  amp <- amp[id]
  ph <- ph[id]
  return( data.frame(amp,ph,freq) )
  #plot(amp~freq, type='l')
  
}
