scalepdf <- 2
pdf("algorithmnew.pdf",width=6 * scalepdf,height=3 * scalepdf)
par(mfrow=c(1,2))

col_ph2 <- 1
lwd_ph2 <- 1.5
plot(1:n_iterations, ph2.stand$parameters$Sigma[1,],xlim=c(1,n_iterations), col = col_ph2, type="l", main="Chain variance standard model, Phase 2", xlab="iteration", ylab="sigma2")
points(2140, estimate.stand$Sigma, col=col_ph2, pch=19)
abline(v=c(100, 200, 400), lty = 2, col = col_ph2, lwd = lwd_ph2); abline(v=600, col = col_ph2, lwd = lwd_ph2)
#abline(v=c(100, 200, 400), col = 1); segments(600, estimate.stand$Sigma, 2140, estimate.stand$Sigma, col=col_ph2)

rm(col_ph2, lwd_ph2)

#abline(v = 1000, col=8) = grey

lwd_line <- 1.5
plot(0,0,type="n",xlim=c(0,35),ylim=c(0,.22), main="Distribution variance, Phase 3", xlab="simulated variance", ylab="density")
points(density(apply(ph3.stand$sf2[,,8+(1:39)],1,var)), type="l", col=1, lty=1, lwd = lwd_line)
points(density(apply(ph3.notrt$sf2[,,8+(1:39)],1,var)), type="l", col=2, lty=1, lwd = lwd_line)
points(density(apply(ph3.nosts$sf2[,,8+(1:39)],1,var)), type="l", col=3, lty=1, lwd = lwd_line)
#points(density(apply(ph3.fullt$sf2[,,8+(1:39)],1,var)), type="l", col=4, lty=1, lwd = lwd_line) # no convergence
points(density(apply(p3c.stand$sf2[,,8+(1:39)],1,var)), type="l", col=1, lty=2, lwd = lwd_line)
points(density(apply(p3c.notrt$sf2[,,8+(1:39)],1,var)), type="l", col=2, lty=2, lwd = lwd_line)
points(density(apply(p3c.nosts$sf2[,,8+(1:39)],1,var)), type="l", col=3, lty=2, lwd = lwd_line)
points(density(apply(p3c.fullt$sf2[,,8+(1:39)],1,var)), type="l", col=4, lty=2, lwd = lwd_line)
abline(v=var(ph3.stand$targets2[9:47]),col=8, lwd=2)

text(23,.21,"model\nwith\nrandom\nout-deg.", cex=.7)
text(27,.21,"model\nwithout\nrandom\nout-deg.", cex=.7)
text(29,.165,pos=4,"standard\nno-transit.\nno-status\nfull",cex=1)
left_seg_1 <- 21.5
left_seg_2 <- 25.5
length_seg <- 3
y_first_seg <- 0.186
distance_segments <- .0105
segments(left_seg_1,y_first_seg,                     left_seg_1+length_seg, lwd=1.5)
segments(left_seg_1,y_first_seg-distance_segments,   left_seg_1+length_seg, col=2, lwd=1.5)
segments(left_seg_1,y_first_seg-2*distance_segments, left_seg_1+length_seg, col=3, lwd=1.5)
#segments(left_seg_1,y_first_seg-3*distance_segments, left_seg_1+length_seg, col=4, lwd=1.5)

segments(left_seg_2,y_first_seg,                     left_seg_2+length_seg, lwd=1.5, lty=2)
segments(left_seg_2,y_first_seg-distance_segments,   left_seg_2+length_seg, col=2, lwd=1.5, lty=2)
segments(left_seg_2,y_first_seg-2*distance_segments, left_seg_2+length_seg, col=3, lwd=1.5, lty=2)
segments(left_seg_2,y_first_seg-3*distance_segments, left_seg_2+length_seg, col=4, lwd=1.5, lty=2)

dev.off()
