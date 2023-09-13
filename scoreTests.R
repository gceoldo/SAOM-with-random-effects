n.iterations.used = 2500L

compute.simulated.variances <- function(modelp3, theta.locations) apply(modelp3$sf2[1:n.iterations.used,1,-(1:8)], 1, var)
simulated.statistics <- function(modelp3, theta.locations) cbind(modelp3$sf2[1:n.iterations.used,1,theta.locations], compute.simulated.variances(modelp3, theta.locations))
computation.simstat.statistics <-  function(modelp3, statistic = mean, theta.locations) apply(simulated.statistics(modelp3, theta.locations), 2, statistic)
observed.statistics <- function(modelp3) c(modelp3$targets2[1:8],var(modelp3$targets2[-(1:8)]))

compute.scores.fixed.parameters <- function(modelp3, theta.locations, random.parameters = NULL, sigma2.hat) {
  if(!is.null(random.parameters)) scores.sigma2 <- rowSums(random.parameters[1:n.iterations.used,] * modelp3$ssc[1:n.iterations.used,1,-(1:8)])/(2*sigma2.hat)
  #if(!is.null(random.parameters)) scores.sigma2 <- rowSums(random.parameters * modelp3$ssc[,1,-(min(theta.locations):max(theta.locations))])
  else scores.sigma2 <- numeric(0)
  scores   <- cbind(modelp3$ssc[1:n.iterations.used,1,theta.locations], scores.sigma2)
  simstats <- simulated.statistics(modelp3, 1:8)
  obsstats <- observed.statistics(modelp3)
  derivative <- matrix(0, nrow = ncol(simstats), ncol = ncol(scores))
  for(i in 1:nrow(simstats)) derivative <- derivative + (simstats[i,]-obsstats) %*% scores[i,,drop=FALSE]
  return(derivative/nrow(simstats))
}


# standard model no ranef, H0: sigma2 = 0
J.ranef.in.stnct <- compute.scores.fixed.parameters(p3c.stand, c(1:4,6:8), NULL)[-5,]
S.ranef.in.stnct <- var(simulated.statistics(p3c.stand, c(1:8)))[-5,-5]
g.ranef.in.stnct <- simulated.statistics(p3c.stand, c(1:4,6:8)) # summary(g.ranef.in.stnct)
G.ranef.in.stnct <- J.ranef.in.stnct[8,, drop = FALSE] %*% solve(J.ranef.in.stnct[1:7,])
X.ranef.in.stnct <- as.numeric(S.ranef.in.stnct[8,8] - 2 * G.ranef.in.stnct %*% S.ranef.in.stnct[1:7,8] + G.ranef.in.stnct %*% S.ranef.in.stnct[1:7,1:7] %*% t(G.ranef.in.stnct))
Y.ranef.in.stnct <- g.ranef.in.stnct[,8] - apply(g.ranef.in.stnct[,1:7], 1, function(x) G.ranef.in.stnct %*% x) # summary(Y.ranef.in.stnct)
y.ranef.in.stnct <- var(p3c.stand$targets2[9:47]) - as.numeric(G.ranef.in.stnct %*% p3c.stand$targets2[c(1:4,6:8)])
1-pnorm((y.ranef.in.stnct-mean(Y.ranef.in.stnct))/sqrt(X.ranef.in.stnct)) # 5.75e-05 < 4e-04
mean(Y.ranef.in.stnct >= y.ranef.in.stnct) # 0 < 4e-04


# standard model no ranef, H0: beta4 = 0
J.activ.in.stnct <- compute.scores.fixed.parameters(p3c.stand, c(1:4, 6:8), NULL)[-9,]
S.activ.in.stnct <- var(simulated.statistics(p3c.stand, c(1:8)))[-9,-9]
g.activ.in.stnct <- simulated.statistics(p3c.stand, c(1:4, 5, 6:8))[,-9] # summary(g.activ.in.stnct)
G.activ.in.stnct <- J.activ.in.stnct[5,, drop = FALSE] %*% solve(J.activ.in.stnct[c(1:4, 6:8),])
X.activ.in.stnct <- as.numeric(S.activ.in.stnct[5,5] - 2 * G.activ.in.stnct %*% S.activ.in.stnct[c(1:4, 6:8),5] + G.activ.in.stnct %*% S.activ.in.stnct[c(1:4, 6:8),c(1:4, 6:8)] %*% t(G.activ.in.stnct))
Y.activ.in.stnct <- g.activ.in.stnct[,5] - apply(g.activ.in.stnct[,c(1:4, 6:8)], 1, function(x) G.activ.in.stnct %*% x) # summary(Y.activ.in.stnct)
y.activ.in.stnct <- p3c.stand$targets2[5] - as.numeric(G.activ.in.stnct %*% p3c.stand$targets2[c(1:4,6:8)])
z2 <- ((y.activ.in.stnct - mean(Y.activ.in.stnct))^2)/X.activ.in.stnct
1-pchisq(z2, df = 1) 
#2*min(pnorm((y.activ.in.stnct-mean(Y.activ.in.stnct))/sqrt(X.activ.in.stnct)), 1-pnorm((y.activ.in.stnct-mean(Y.activ.in.stnct))/sqrt(X.activ.in.stnct)))
zsim <- ((Y.activ.in.stnct-mean(Y.activ.in.stnct))^2)/X.activ.in.stnct
mean(zsim > z2) 


# full model no ranef, H0: sigma2 = 0
J.ranef.in.fulct <- compute.scores.fixed.parameters(p3c.fullt, c(1:4, 5, 6:8), NULL)[,]
S.ranef.in.fulct <- var(simulated.statistics(p3c.fullt, c(1:8)))[,]
g.ranef.in.fulct <- simulated.statistics(p3c.fullt, c(1:4, 5, 6:8)) # summary(g.ranef.in.fulct)
G.ranef.in.fulct <- J.ranef.in.fulct[9,, drop = FALSE] %*% solve(J.ranef.in.fulct[1:8,])
X.ranef.in.fulct <- as.numeric(S.ranef.in.fulct[9,9] - 2 * G.ranef.in.fulct %*% S.ranef.in.fulct[1:8,9] + G.ranef.in.fulct %*% S.ranef.in.fulct[1:8,1:8] %*% t(G.ranef.in.fulct))
Y.ranef.in.fulct <- g.ranef.in.fulct[,9] - apply(g.ranef.in.fulct[,1:8], 1, function(x) G.ranef.in.fulct %*% x) # summary(Y.ranef.in.fulct)
y.ranef.in.fulct <- var(p3c.fullt$targets2[9:47]) - as.numeric(G.ranef.in.fulct %*% p3c.fullt$targets2[c(1:4, 5, 6:8)])
1-pnorm((y.ranef.in.fulct-mean(Y.ranef.in.fulct))/sqrt(X.ranef.in.fulct)) # 0.23
mean(Y.ranef.in.fulct >= y.ranef.in.fulct) # 0.1908


# standard model with ranef, H0: beta4 = 0
J.activ.in.stand <- compute.scores.fixed.parameters(ph3.stand, c(1:4, 6:8), theta.stand[,-(1:8)], estimate.stand$Sigma)[,] # summary(apply(theta.stand[,-(1:8)], 1, var))
S.activ.in.stand <- var(simulated.statistics(ph3.stand, c(1:8)))[,]
g.activ.in.stand <- simulated.statistics(ph3.stand, c(1:4, 5, 6:8))[,] # summary(g.activ.in.stand)
G.activ.in.stand <- J.activ.in.stand[5,, drop = FALSE] %*% solve(J.activ.in.stand[c(1:4, 6:9),])
X.activ.in.stand <- as.numeric(S.activ.in.stand[5,5] - 2 * G.activ.in.stand %*% S.activ.in.stand[c(1:4, 6:9),5] + G.activ.in.stand %*% S.activ.in.stand[c(1:4, 6:9),c(1:4, 6:9)] %*% t(G.activ.in.stand))
Y.activ.in.stand <- g.activ.in.stand[,5] - apply(g.activ.in.stand[,c(1:4, 6:9)], 1, function(x) G.activ.in.stand %*% x) # summary(Y.activ.in.stand)
y.activ.in.stand <- ph3.stand$targets2[5] - as.numeric(G.activ.in.stand %*% c(ph3.stand$targets2[c(1:4,6:8)], var(ph3.stand$targets2[9:47])))
z2 <- ((y.activ.in.stand - mean(Y.activ.in.stand))^2)/X.activ.in.stand
1-pchisq(z2, df = 1) # 0.73
2*min(pnorm((y.activ.in.stand-mean(Y.activ.in.stand))/sqrt(X.activ.in.stand)), 1-pnorm((y.activ.in.stand-mean(Y.activ.in.stand))/sqrt(X.activ.in.stand)))
zsim <- ((Y.activ.in.stand-mean(Y.activ.in.stand))^2)/X.activ.in.stand
mean(zsim > z2) # 0.7108


# no transitivity model with ranef, H0: beta3 = 0
J.trans.in.notrt <- compute.scores.fixed.parameters(ph3.notrt, c(1:3, 6:8), theta.notrt[,-(1:8)], estimate.notrt$Sigma)[-5,] # summary(apply(theta.stand[,-(1:8)], 1, var))
S.trans.in.notrt <- var(simulated.statistics(ph3.notrt, c(1:8)))[-5,-5]
g.trans.in.notrt <- simulated.statistics(ph3.notrt, c(1:4, 6:8))[,] # summary(g.trans.in.notrt)
G.trans.in.notrt <- J.trans.in.notrt[4,, drop = FALSE] %*% solve(J.trans.in.notrt[c(1:3, 5:8),])
X.trans.in.notrt <- as.numeric(S.trans.in.notrt[4,4] - 2 * G.trans.in.notrt %*% S.trans.in.notrt[c(1:3, 5:8),4] + G.trans.in.notrt %*% S.trans.in.notrt[c(1:3, 5:8),c(1:3, 5:8)] %*% t(G.trans.in.notrt))
Y.trans.in.notrt <- g.trans.in.notrt[,4] - apply(g.trans.in.notrt[,c(1:3, 5:8)], 1, function(x) G.trans.in.notrt %*% x) # summary(Y.trans.in.notrt)
y.trans.in.notrt <- ph3.notrt$targets2[4] - as.numeric(G.trans.in.notrt %*% c(ph3.notrt$targets2[c(1:3,6:8)], var(ph3.notrt$targets2[9:47])))
z2 <- ((y.trans.in.notrt - mean(Y.trans.in.notrt))^2)/X.trans.in.notrt
1-pchisq(z2, df = 1) # 0.0775
2*min(pnorm((y.trans.in.notrt-mean(Y.trans.in.notrt))/sqrt(X.trans.in.notrt)), 1-pnorm((y.trans.in.notrt-mean(Y.trans.in.notrt))/sqrt(X.trans.in.notrt)))
zsim <- ((Y.trans.in.notrt-mean(Y.trans.in.notrt))^2)/X.trans.in.notrt
mean(zsim > z2) # 0.0756


# no status model with ranef, H0: beta5 = beta6 = beta7 = 0
J.stats.in.nosts <- compute.scores.fixed.parameters(ph3.nosts, c(1:4), theta.nosts[,-(1:8)], estimate.nosts$Sigma)[-5,] # summary(apply(theta.stand[,-(1:8)], 1, var))
S.stats.in.nosts <- var(simulated.statistics(ph3.nosts, c(1:8)))[-5,-5]
g.stats.in.nosts <- simulated.statistics(ph3.nosts, c(1:4, 6:8))[,] # summary(g.stats.in.nosts)
G.stats.in.nosts <- J.stats.in.nosts[5:7,, drop = FALSE] %*% solve(J.stats.in.nosts[c(1:4, 8),])
X.stats.in.nosts <- S.stats.in.nosts[5:7,5:7] - (G.stats.in.nosts %*% S.stats.in.nosts[c(1:4, 8),5:7] + t(S.stats.in.nosts[c(1:4, 8),5:7]) %*% t(G.stats.in.nosts)) + G.stats.in.nosts %*% S.stats.in.nosts[c(1:4, 8),c(1:4, 8)] %*% t(G.stats.in.nosts)
Y.stats.in.nosts <- t(apply(g.stats.in.nosts, 1, function(x) x[5:7] - G.stats.in.nosts %*% x[c(1:4, 8)])) # summary(Y.stats.in.nosts)
y.stats.in.nosts <- ph3.nosts$targets2[6:8] - G.stats.in.nosts %*% c(ph3.nosts$targets2[1:4], var(ph3.nosts$targets2[9:47]))
z2 <- as.numeric(t(y.stats.in.nosts - colMeans(Y.stats.in.nosts)) %*% solve(X.stats.in.nosts) %*% (y.stats.in.nosts - colMeans(Y.stats.in.nosts)))
1-pchisq(z2, df = 1) # 1.87e-07
zsim <- apply(Y.stats.in.nosts, 1, function(x) as.numeric(t(x - colMeans(Y.stats.in.nosts)) %*% solve(X.stats.in.nosts) %*% (x - colMeans(Y.stats.in.nosts))))
mean(zsim > z2) # 0 # summary(zsim)



######## other tests

# no-transitivity model no ranef, H0: sigma2 = 0
J.ranef.in.ntrct <- compute.scores.fixed.parameters(p3c.notrt, c(1:3,6:8), NULL)[-c(4,5),]
S.ranef.in.ntrct <- var(simulated.statistics(p3c.notrt, c(1:8)))[-c(4,5),-c(4,5)]
g.ranef.in.ntrct <- simulated.statistics(p3c.notrt, c(1:3,6:8)); summary(g.ranef.in.ntrct)
G.ranef.in.ntrct <- J.ranef.in.ntrct[7,, drop = FALSE] %*% solve(J.ranef.in.ntrct[1:6,])
X.ranef.in.ntrct <- as.numeric(S.ranef.in.ntrct[7,7] - 2 * G.ranef.in.ntrct %*% S.ranef.in.ntrct[1:6,7] + G.ranef.in.ntrct %*% S.ranef.in.ntrct[1:6,1:6] %*% t(G.ranef.in.ntrct))
Y.ranef.in.ntrct <- g.ranef.in.ntrct[,7] - apply(g.ranef.in.ntrct[,1:6], 1, function(x) G.ranef.in.ntrct %*% x) # summary(Y.ranef.in.ntrct)
y.ranef.in.ntrct <- var(p3c.notrt$targets2[9:47]) - as.numeric(G.ranef.in.ntrct %*% p3c.notrt$targets2[c(1:3,6:8)])
1-pnorm((y.ranef.in.ntrct-mean(Y.ranef.in.ntrct))/sqrt(X.ranef.in.ntrct)) # 1.44329e-15
mean(Y.ranef.in.ntrct >= y.ranef.in.ntrct) # 0

# no-status model no ranef, H0: sigma2 = 0
J.ranef.in.nstct <- compute.scores.fixed.parameters(p3c.nosts, c(1:3), NULL)[-c(4:8),]
S.ranef.in.nstct <- var(simulated.statistics(p3c.nosts, c(1:8)))[-c(4:8),-c(4:8)]
g.ranef.in.nstct <- simulated.statistics(p3c.nosts, c(1:3)); summary(g.ranef.in.nstct)
G.ranef.in.nstct <- J.ranef.in.nstct[4,, drop = FALSE] %*% solve(J.ranef.in.nstct[1:3,])
X.ranef.in.nstct <- as.numeric(S.ranef.in.nstct[4,4] - 2 * G.ranef.in.nstct %*% S.ranef.in.nstct[1:3,4] + G.ranef.in.nstct %*% S.ranef.in.nstct[1:3,1:3] %*% t(G.ranef.in.nstct))
Y.ranef.in.nstct <- g.ranef.in.nstct[,4] - apply(g.ranef.in.nstct[,1:3], 1, function(x) G.ranef.in.nstct %*% x) # summary(Y.ranef.in.nstct)
y.ranef.in.nstct <- var(p3c.nosts$targets2[9:47]) - as.numeric(G.ranef.in.nstct %*% p3c.nosts$targets2[c(1:3)])
1-pnorm((y.ranef.in.nstct-mean(Y.ranef.in.nstct))/sqrt(X.ranef.in.nstct)) # 2.305141e-05
mean(Y.ranef.in.nstct >= y.ranef.in.nstct) # 8e-04

# no status model without ranef, H0: beta5 = beta6 = beta7 = 0
J.stats.in.nstct <- compute.scores.fixed.parameters(p3c.nosts, c(1:4), NULL)[-c(5,9),] 
S.stats.in.nstct <- var(simulated.statistics(p3c.nosts, c(1:8)))[-c(5,9),-c(5,9)]
g.stats.in.nstct <- simulated.statistics(p3c.nosts, c(1:4, 6:8))[,-8]; summary(g.stats.in.nstct)
G.stats.in.nstct <- J.stats.in.nstct[5:7,, drop = FALSE] %*% solve(J.stats.in.nstct[c(1:4),])
X.stats.in.nstct <- S.stats.in.nstct[5:7,5:7] - (G.stats.in.nstct %*% S.stats.in.nstct[c(1:4),5:7] + t(S.stats.in.nstct[c(1:4),5:7]) %*% t(G.stats.in.nstct)) + G.stats.in.nstct %*% S.stats.in.nstct[c(1:4),c(1:4)] %*% t(G.stats.in.nstct)
Y.stats.in.nstct <- t(apply(g.stats.in.nstct, 1, function(x) x[5:7] - G.stats.in.nstct %*% x[c(1:4)])); summary(Y.stats.in.nstct)
y.stats.in.nstct <- p3c.nosts$targets2[6:8] - G.stats.in.nstct %*% c(p3c.nosts$targets2[1:4])
z2 <- as.numeric(t(y.stats.in.nstct - colMeans(Y.stats.in.nstct)) %*% solve(X.stats.in.nstct) %*% (y.stats.in.nstct - colMeans(Y.stats.in.nstct)))
1-pchisq(z2, df = 1) # 1.807621e-11
zsim <- apply(Y.stats.in.nstct, 1, function(x) as.numeric(t(x - colMeans(Y.stats.in.nstct)) %*% solve(X.stats.in.nstct) %*% (x - colMeans(Y.stats.in.nstct))))
mean(zsim > z2) # 0


# no transitivity model with ranef, H0: beta3 = 0
J.trans.in.ntrct <- compute.scores.fixed.parameters(p3c.notrt, c(1:3, 6:8), NULL)[-c(5,8),] 
S.trans.in.ntrct <- var(simulated.statistics(p3c.notrt, c(1:8)))[-c(5,8),-c(5,8)]
g.trans.in.ntrct <- simulated.statistics(p3c.notrt, c(1:4, 6:8))[,-8]; summary(g.trans.in.ntrct)
G.trans.in.ntrct <- J.trans.in.ntrct[4,, drop = FALSE] %*% solve(J.trans.in.ntrct[c(1:3, 5:7),])
X.trans.in.ntrct <- as.numeric(S.trans.in.ntrct[4,4] - 2 * G.trans.in.ntrct %*% S.trans.in.ntrct[c(1:3, 5:7),4] + G.trans.in.ntrct %*% S.trans.in.ntrct[c(1:3, 5:7),c(1:3, 5:7)] %*% t(G.trans.in.ntrct))
Y.trans.in.ntrct <- g.trans.in.ntrct[,4] - apply(g.trans.in.ntrct[,c(1:3, 5:7)], 1, function(x) G.trans.in.ntrct %*% x); summary(Y.trans.in.ntrct)
y.trans.in.ntrct <- p3c.notrt$targets2[4] - as.numeric(G.trans.in.ntrct %*% c(p3c.notrt$targets2[c(1:3,6:8)])); y.trans.in.ntrct
#bias <- rowMeans(apply(g.trans.in.ntrct[,c(1:3, 5:7)], 1, function(x) x - p3c.notrt$targets2[c(1:3,6:8)]))
#Y.trans.in.ntrct <- g.trans.in.ntrct[,4] - apply(g.trans.in.ntrct[,c(1:3, 5:7)], 1, function(x) G.trans.in.ntrct %*% (x - bias)); summary(Y.trans.in.ntrct)
#y.trans.in.ntrct <- p3c.notrt$targets2[4] - as.numeric(G.trans.in.ntrct %*% c(p3c.notrt$targets2[c(1:3,6:8)] - bias)); y.trans.in.ntrct
#Y.trans.in.ntrct <- g.trans.in.ntrct[,4]; summary(Y.trans.in.ntrct)
#y.trans.in.ntrct <- p3c.notrt$targets2[4]; y.trans.in.ntrct
z2 <- ((y.trans.in.ntrct - mean(Y.trans.in.ntrct))^2)/X.trans.in.ntrct
1-pchisq(z2, df = 1) # 0.00258
zsim <- ((Y.trans.in.ntrct-mean(Y.trans.in.ntrct))^2)/X.trans.in.ntrct
mean(zsim > z2) # 0.496 # summary(zsim)
