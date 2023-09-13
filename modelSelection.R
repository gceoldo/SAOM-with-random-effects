n.iterations.used = 2500L
simulated.statistics <- function(modelp3, theta.locations) cbind(modelp3$sf2[1:n.iterations.used,1,theta.locations], compute.simulated.variances(modelp3, theta.locations))
observed.statistics  <- function(modelp3, theta.locations) c(modelp3$targets2[theta.locations],var(modelp3$targets2[-(1:8)]))
simulated.deviations <- function(modelp3, theta.locations) t(apply(simulated.statistics(modelp3, theta.locations), 1, function(x) x - observed.statistics(modelp3, theta.locations)))


degrees.freedom <- function(n) 9
penaltyAIC <- function(df) 2L
penaltyBIC <- function(df) log(df)

p0 = 8L
q0 = 1L

psc <- function(modelp3, p, has.ranef, penalty = penaltyAIC, theta.locations = 1:8, n = 39L) {
  q = as.integer(has.ranef)
  second.term = (p0 - p + q0 - q) * penalty(degrees.freedom(n))
  
  inverted.variance = solve(var(simulated.statistics(modelp3, theta.locations)))
  first.term  = mean(apply(simulated.deviations(modelp3, theta.locations), 1, function(x) as.vector(t(x) %*% inverted.variance %*% x))) * degrees.freedom(n)
  
  return(first.term - second.term)
}

psc(ph3.stand, 7L, TRUE)
psc(ph3.notrt, 6L, TRUE)
psc(ph3.nosts, 4L, TRUE)
psc(p3c.stand, 7L, FALSE)
psc(p3c.fullt, 8L, FALSE)
psc(p3c.notrt, 6L, FALSE)
psc(p3c.nosts, 4L, FALSE)

psc(ph3.stand, 7L, TRUE, penaltyBIC)
psc(ph3.notrt, 6L, TRUE, penaltyBIC)
psc(ph3.nosts, 4L, TRUE, penaltyBIC)
psc(p3c.stand, 7L, FALSE, penaltyBIC)
psc(p3c.fullt, 8L, FALSE, penaltyBIC)
psc(p3c.notrt, 6L, FALSE, penaltyBIC)
psc(p3c.nosts, 4L, FALSE, penaltyBIC)
