n.iterations.used = 2500L

#compute.simulated.variances <- function(modelp3, theta.locations) apply(modelp3$sf2[,1,-theta.locations], 1, var)
compute.simulated.variances <- function(modelp3, theta.locations) apply(modelp3$sf2[1:n.iterations.used,1,-(1:8)], 1, var)
simulated.statistics <- function(modelp3, theta.locations) cbind(modelp3$sf2[1:n.iterations.used,1,theta.locations], compute.simulated.variances(modelp3, theta.locations))
#observed.statistics  <- function(modelp3, theta.locations) c(modelp3$targets2[theta.locations],var(modelp3$targets2[-theta.locations]))
observed.statistics  <- function(modelp3, theta.locations) c(modelp3$targets2[theta.locations],var(modelp3$targets2[-(1:8)]))
simulated.deviations <- function(modelp3, theta.locations) t(apply(simulated.statistics(modelp3, theta.locations), 1, function(x) x - observed.statistics(modelp3, theta.locations)))
computation.simstat.statistics <-  function(modelp3, statistic = mean, theta.locations) apply(simulated.statistics(modelp3, theta.locations), 2, statistic)

simulated.statistics(ph3.stand, 1:7)[1:3,] # last is simulated variance
computation.simstat.statistics(ph3.stand, mean, 1:7)
computation.simstat.statistics(ph3.stand, var, 1:7)

c(estimate.stand$Sigma^2, estimate.notrt$Sigma^2, estimate.nosts$Sigma^2)

compute.scores.fixed.parameters <- function(modelp3, theta.locations, random.parameters = NULL, sigma2.hat) {
  if(!is.null(random.parameters)) scores.sigma2 <- rowSums(random.parameters[1:n.iterations.used,] * modelp3$ssc[1:n.iterations.used,1,-(1:8)])/(2*sigma2.hat)
  #if(!is.null(random.parameters)) scores.sigma2 <- rowSums(random.parameters * modelp3$ssc[,1,-(1:8)])
  else scores.sigma2 <- numeric(0)
  scores   <- cbind(modelp3$ssc[1:n.iterations.used,1,theta.locations], scores.sigma2)
  #simstats <- simulated.statistics(modelp3, theta.locations)
  simstats <- simulated.deviations(modelp3, theta.locations)
  derivative <- matrix(0, nrow = ncol(simstats), ncol = ncol(scores))
  for(i in 1:nrow(simstats)) derivative <- derivative + simstats[i,] %*% scores[i,,drop=FALSE]
  return(derivative/nrow(simstats))
}



J.stand.r <- compute.scores.fixed.parameters(ph3.stand, c(1:4, 6:8), theta.stand[,-(1:8)], estimate.stand$Sigma)
S.stand.r <- var(simulated.statistics(ph3.stand, c(1:4, 6:8)))
C.stand.r <- solve(J.stand.r) %*% S.stand.r %*% t(solve(J.stand.r))
sqrt(diag(C.stand.r))


J.stand.c <- compute.scores.fixed.parameters(p3c.stand, c(1:4, 6:8), NULL)
S.stand.c <- var(simulated.statistics(p3c.stand, c(1:4, 6:8)))
C.stand.c <- solve(J.stand.c[-8,]) %*% S.stand.c[1:7,1:7] %*% t(solve(J.stand.c[-8,]))
sqrt(diag(C.stand.c))



J.fullt.c <- compute.scores.fixed.parameters(p3c.fullt, 1:8, NULL)
S.fullt.c <- var(simulated.statistics(p3c.fullt, 1:8))
C.fullt.c <- solve(J.fullt.c[-9,]) %*% S.fullt.c[1:8,1:8] %*% t(solve(J.fullt.c[-9,]))
sqrt(diag(C.fullt.c))




J.notrt.r <- compute.scores.fixed.parameters(ph3.notrt, c(1:3,6:8), theta.notrt[,-(1:8)], estimate.notrt$Sigma)
S.notrt.r <- var(simulated.statistics(ph3.notrt, c(1:3,6:8)))
C.notrt.r <- solve(J.notrt.r) %*% S.notrt.r %*% t(solve(J.notrt.r))
sqrt(diag(C.notrt.r))

J.notrt.c <- compute.scores.fixed.parameters(p3c.notrt, c(1:3,6:8), NULL)
S.notrt.c <- var(simulated.statistics(p3c.notrt, c(1:3,6:8)))
C.notrt.c <- solve(J.notrt.c[-7,]) %*% S.notrt.c[1:6,1:6] %*% t(solve(J.notrt.c[-7,]))
sqrt(diag(C.notrt.c))



J.nosts.r <- compute.scores.fixed.parameters(ph3.nosts, 1:4, theta.nosts[,-(1:8)], estimate.nosts$Sigma)
S.nosts.r <- var(simulated.statistics(ph3.nosts, 1:4))
C.nosts.r <- solve(J.nosts.r) %*% S.nosts.r %*% t(solve(J.nosts.r))
sqrt(diag(C.nosts.r))

J.nosts.c <- compute.scores.fixed.parameters(p3c.nosts, 1:4, NULL)
S.nosts.c <- var(simulated.statistics(p3c.nosts, 1:4))
C.nosts.c <- solve(J.nosts.c[-5,]) %*% S.nosts.c[1:4,1:4] %*% t(solve(J.nosts.c[-5,]))
sqrt(diag(C.nosts.c))
model.nosts$se




