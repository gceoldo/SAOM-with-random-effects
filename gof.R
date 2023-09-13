n.iterations.used = 2500L

outDegDistr = function(outDegrees, maxDeg) {
  distr = rep(NA, times = maxDeg+1)
  for(i in 0:maxDeg) distr[i+1] = sum(outDegrees == i)
  return(distr/sum(distr))
} 

pvalueGOF = function(phase3, maxDeg, which_out_degrees, correction_factor = 0) {
  simOutDegrees = phase3$sf2[1:n.iterations.used,1,which_out_degrees]
  obsOutDegrees = phase3$targets2[which_out_degrees,1]
  
  simOutDegDistr = t(apply(simOutDegrees, 1, function(x) outDegDistr(x, maxDeg)))
  obsOutDegDistr = outDegDistr(obsOutDegrees, maxDeg)
  
  meaOutDegDistr = colMeans(simOutDegDistr)
  varOutDegDistr = var(simOutDegDistr) + diag(correction_factor, ncol(simOutDegDistr))
  #varOutDegDistr = (1 - correction_factor) * varOutDegDistr + correction_factor * diag(diag(varOutDegDistr))
  
  deviance = function(x) as.vector(t(x - meaOutDegDistr) %*% solve(varOutDegDistr) %*% (x - meaOutDegDistr))
  simDeviance = apply(simOutDegDistr, 1, deviance)
  obsDeviance = deviance(obsOutDegDistr)
  return(mean(simDeviance > obsDeviance))
  
  #deviances = apply(simOutDegDistr, 1, function(x) as.vector((x - obsOutDegDistr) %*% solve(varOutDegDistr) %*% (x - obsOutDegDistr)))
  #return(mean(deviances))
}

maxDeg = 20
correction_factor = 0.2

pvalueGOF(ph3.stand, maxDeg, 9:47, correction_factor) 
pvalueGOF(p3c.stand, maxDeg, 9:47, correction_factor) 

pvalueGOF(p3c.fullt, maxDeg, 9:47, correction_factor) 

pvalueGOF(ph3.notrt, maxDeg, 9:47, correction_factor) 
pvalueGOF(p3c.notrt, maxDeg, 9:47, correction_factor) 

pvalueGOF(ph3.nosts, maxDeg, 9:47, correction_factor) 
pvalueGOF(p3c.nosts, maxDeg, 9:47, correction_factor) 

