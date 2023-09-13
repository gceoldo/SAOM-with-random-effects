burn.in <- 600 # total 1500 iterations for tail average
# if short chains (30 iterations) are used:
# burn.in <- 20


computation.tail.statistic <- function(chain, burn.in, statistic = mean, type_chain = "parameters") {
  return(list(theta = apply(chain[[type_chain]]$theta[,-(1:burn.in)], 1, statistic), Sigma = statistic(chain[[type_chain]]$Sigma[-(1:burn.in)])))
}


n.iterations.phase3 <- 2500
# short phase 3:
# n.iterations.phase3 <- 100

algorithm.ph3 <- sienaAlgorithmCreate(projname = 'ph3', cond = FALSE,
                                      useStdInits = FALSE, nsub = 0, n3 = n.iterations.phase3, simOnly = TRUE)


estimate.stand <- computation.tail.statistic(ph2.stand, burn.in, mean)
estimate.stand$theta <- c(estimate.stand$theta[1:4],0,estimate.stand$theta[5:7])
theta.stand <- matrix(0, n.iterations.phase3, 8 + 39)
for(i in 1:n.iterations.phase3) {
  theta.stand[i,1:length(estimate.stand$theta)] <- estimate.stand$theta
  theta.stand[i,-(1:length(estimate.stand$theta))] <- sqrt(estimate.stand$Sigma)*rnorm(39)
}
ph3.stand <- siena07(algorithm.ph3, data=data, effects=effects.ranef.fullt, thetaValues = theta.stand)

# not estimated as no convergence
#estimate.fullt <- computation.tail.statistic(ph2.fullt, burn.in, mean)
#theta.fullt <- matrix(0, n.iterations.phase3, length(estimate.fullt$theta) + 39)
#for(i in 1:n.iterations.phase3) {
#  theta.fullt[i,1:length(estimate.fullt$theta)] <- estimate.fullt$theta
#  theta.fullt[i,-(1:length(estimate.fullt$theta))] <- sqrt(estimate.fullt$Sigma)*rnorm(39)
#}
#ph3.fullt <- siena07(algorithm.ph3, data=data, effects=effects.ranef.fullt, thetaValues = theta.fullt)

estimate.notrt <- computation.tail.statistic(ph2.notrt, burn.in, mean)
estimate.notrt$theta <- c(estimate.notrt$theta[1:3],0,0,estimate.notrt$theta[4:6])
theta.notrt <- matrix(0, n.iterations.phase3, 8 + 39)
for(i in 1:n.iterations.phase3) {
  theta.notrt[i,1:length(estimate.notrt$theta)] <- estimate.notrt$theta
  theta.notrt[i,-(1:length(estimate.notrt$theta))] <- sqrt(estimate.notrt$Sigma)*rnorm(39)
}
ph3.notrt <- siena07(algorithm.ph3, data=data, effects=effects.ranef.fullt, thetaValues = theta.notrt)

estimate.nosts <- computation.tail.statistic(ph2.nosts, burn.in, mean)
estimate.nosts$theta <- c(estimate.nosts$theta, 0,0,0,0)
theta.nosts <- matrix(0, n.iterations.phase3, 8 + 39)
for(i in 1:n.iterations.phase3) {
  theta.nosts[i,1:length(estimate.nosts$theta)] <- estimate.nosts$theta
  theta.nosts[i,-(1:length(estimate.nosts$theta))] <- sqrt(estimate.nosts$Sigma)*rnorm(39)
}
ph3.nosts <- siena07(algorithm.ph3, data=data, effects=effects.ranef.fullt, thetaValues = theta.nosts)



estcntrl.stand <- computation.tail.statistic(p2c.stand, burn.in, mean)
estcntrl.stand$theta <- c(estcntrl.stand$theta[1:4],0,estcntrl.stand$theta[5:7])
th.stand <- matrix(0, n.iterations.phase3, 8 + 39)
for(i in 1:n.iterations.phase3) {
  th.stand[i,1:length(estcntrl.stand$theta)] <- estcntrl.stand$theta
  th.stand[i,-(1:length(estcntrl.stand$theta))] <- sqrt(estcntrl.stand$Sigma)*rnorm(39)
}
p3c.stand <- siena07(algorithm.ph3, data=data, effects=effects.ranef.fullt, thetaValues = th.stand)
rm(th.stand)

estcntrl.fullt <- computation.tail.statistic(p2c.fullt, burn.in, mean)
th.fullt <- matrix(0, n.iterations.phase3, 8 + 39)
for(i in 1:n.iterations.phase3) {
  th.fullt[i,1:length(estcntrl.fullt$theta)] <- estcntrl.fullt$theta
  th.fullt[i,-(1:length(estcntrl.fullt$theta))] <- sqrt(estcntrl.fullt$Sigma)*rnorm(39)
}
p3c.fullt <- siena07(algorithm.ph3, data=data, effects=effects.ranef.fullt, thetaValues = th.fullt)
rm(th.fullt)

estcntrl.notrt <- computation.tail.statistic(p2c.notrt, burn.in, mean)
estcntrl.notrt$theta <- c(estcntrl.notrt$theta[1:3],0,0,estcntrl.notrt$theta[4:6])
th.notrt <- matrix(0, n.iterations.phase3, 8 + 39)
for(i in 1:n.iterations.phase3) {
  th.notrt[i,1:length(estcntrl.notrt$theta)] <- estcntrl.notrt$theta
  th.notrt[i,-(1:length(estcntrl.notrt$theta))] <- sqrt(estcntrl.notrt$Sigma)*rnorm(39)
}
p3c.notrt <- siena07(algorithm.ph3, data=data, effects=effects.ranef.fullt, thetaValues = th.notrt)
rm(th.notrt)

estcntrl.nosts <- computation.tail.statistic(p2c.nosts, burn.in, mean)
estcntrl.nosts$theta <- c(estcntrl.nosts$theta, 0,0,0,0)
th.nosts <- matrix(0, n.iterations.phase3, 8 + 39)
for(i in 1:n.iterations.phase3) {
  th.nosts[i,1:length(estcntrl.nosts$theta)] <- estcntrl.nosts$theta
  th.nosts[i,-(1:length(estcntrl.nosts$theta))] <- sqrt(estcntrl.nosts$Sigma)*rnorm(39)
}
p3c.nosts <- siena07(algorithm.ph3, data=data, effects=effects.ranef.fullt, thetaValues = th.nosts)
rm(th.nosts)



