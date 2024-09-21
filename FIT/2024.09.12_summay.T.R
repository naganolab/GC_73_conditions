

# summarize the results of Takanari

setwd("Tak.result")

load("pred.T.GC")   # 400 x 474 (5-474: genes)
load("pred.T.fld")  # 400 x 474 (5-474: genes)
load("pred.T.mix")  # 400 x 474 (5-474: genes)


#======================================#
# Remove genes with prediction failure #
#======================================#

FG.GC  <- apply(pred.T.GC[,5:474], 2, function(x) sum(is.na(x))>=1)   # 90, 134, 242, 256, 270, 287, 295, 303, 334,           396, 429 
FG.fld <- apply(pred.T.fld[,5:474], 2, function(x) sum(is.na(x))>=1)  # 90,                270,      295, 303, 334, 353, 395,          440 
FG.mix <- apply(pred.T.mix[,5:474], 2, function(x) sum(is.na(x))>=1)  # 90, 134, 242,      270,      295, 303, 334,           396, 429

FG.v <- apply(rbind(FG.GC, FG.fld, FG.mix), 2, function(x) sum(x)>=1) # 456 genes (14 failed genes)

pred.T.GC.nna  <- pred.T.GC[,5:474][, !FG.v]
pred.T.fld.nna <- pred.T.fld[,5:474][, !FG.v]
pred.T.mix.nna <- pred.T.mix[,5:474][, !FG.v]

pred.T.GC.nna <- cbind(pred.T.GC[, c(1,2,4)], pred.T.GC.nna) # 400 x (3 + 456)
pred.T.GC.nna$MAE <- apply(pred.T.GC.nna[, 4:459], 1, mean)

pred.T.fld.nna <- cbind(pred.T.fld[, c(1,2,4)], pred.T.fld.nna)
pred.T.fld.nna$MAE <- apply(pred.T.fld.nna[, 4:459], 1, mean)

pred.T.mix.nna <- cbind(pred.T.mix[, c(1,2,4)], pred.T.mix.nna)
pred.T.mix.nna$MAE <- apply(pred.T.mix.nna[, 4:459], 1, mean)



#===================#
# Process the files #
#===================#

### GC ###

pred.m.T.GC.64  <- as.matrix(pred.T.GC.nna[1:100, 4:459])
pred.m.T.GC.128 <- as.matrix(pred.T.GC.nna[101:200, 4:459])
pred.m.T.GC.256 <- as.matrix(pred.T.GC.nna[201:300, 4:459])
pred.m.T.GC.512 <- as.matrix(pred.T.GC.nna[301:400, 4:459])

min(pred.m.T.GC.64) ; max(pred.m.T.GC.64)     # 0.02158717, 111131638
min(pred.m.T.GC.128) ; max(pred.m.T.GC.128)   # 0.02975229, 10797198
min(pred.m.T.GC.256) ; max(pred.m.T.GC.256)   # 0.03336428, 78048764
min(pred.m.T.GC.512) ; max(pred.m.T.GC.512)   # 0.03576249, 431932940

pred.m.T.GC.64[pred.m.T.GC.64 >= 5] <- 5
pred.m.T.GC.128[pred.m.T.GC.128 >= 5] <- 5
pred.m.T.GC.256[pred.m.T.GC.256 >= 5] <- 5
pred.m.T.GC.512[pred.m.T.GC.512 >= 5] <- 5


### Field ###

pred.m.T.fld.64  <- as.matrix(pred.T.fld.nna[1:100, 4:459])
pred.m.T.fld.128 <- as.matrix(pred.T.fld.nna[101:200, 4:459])
pred.m.T.fld.256 <- as.matrix(pred.T.fld.nna[201:300, 4:459])
pred.m.T.fld.512 <- as.matrix(pred.T.fld.nna[301:400, 4:459])

min(pred.m.T.fld.64) ; max(pred.m.T.fld.64)     # 0.02158717, 311840413
min(pred.m.T.fld.128) ; max(pred.m.T.fld.128)   # 0.02589991, 5.688812e+17
min(pred.m.T.fld.256) ; max(pred.m.T.fld.256)   # 0.03342755, 1.403667e+13
min(pred.m.T.fld.512) ; max(pred.m.T.fld.512)   # 0.04206121, 938712851513

pred.m.T.fld.64[pred.m.T.fld.64 >= 5] <- 5
pred.m.T.fld.128[pred.m.T.fld.128 >= 5] <- 5
pred.m.T.fld.256[pred.m.T.fld.256 >= 5] <- 5
pred.m.T.fld.512[pred.m.T.fld.512 >= 5] <- 5


### Mix ###

pred.m.T.mix.64  <- as.matrix(pred.T.mix.nna[1:100, 4:459])
pred.m.T.mix.128 <- as.matrix(pred.T.mix.nna[101:200, 4:459])
pred.m.T.mix.256 <- as.matrix(pred.T.mix.nna[201:300, 4:459])
pred.m.T.mix.512 <- as.matrix(pred.T.mix.nna[301:400, 4:459])

min(pred.m.T.mix.64) ; max(pred.m.T.mix.64)     # 0.02158717, 809315.1
min(pred.m.T.mix.128) ; max(pred.m.T.mix.128)   # 0.0232929, 37738.86
min(pred.m.T.mix.256) ; max(pred.m.T.mix.256)   # 0.02581438, 1320.268
min(pred.m.T.mix.512) ; max(pred.m.T.mix.512)   # 0.02897933, 28.01769

pred.m.T.mix.64[pred.m.T.mix.64 >= 5] <- 5
pred.m.T.mix.128[pred.m.T.mix.128 >= 5] <- 5
pred.m.T.mix.256[pred.m.T.mix.256 >= 5] <- 5
pred.m.T.mix.512[pred.m.T.mix.512 >= 5] <- 5



#==================#
# Make the heatmap #
#==================#

library(viridis)


### GC ###

svg("accuracy.heatmap.T.GC.svg", width=5, height=6)

layout(matrix(c(0,0,1,2,3,0,4,0,5,0), nrow=5, byrow=T),
       widths = c(5, 1), 
       heights = c(1,5,5,5,5))

max.value <- 5 # max of the matrix
min.value <- 0

the.M <- pred.m.T.GC.64
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("A", side=3, line=0.5, at=-5)

par(mar = c(1.6, 3, 3.7, 2))
colorLevels <- seq(min.value, max.value, length=256)
image(1, colorLevels,
      matrix(colorLevels, ncol = length(colorLevels), nrow = 1),
      col = viridis(256), xlab=" ", ylab="", axes=F, bty="n")
axis(2, at=c(0:5), labels=c("0", "1", "2", "3", "4", ">=5"), las=1)

the.M <- pred.m.T.GC.128
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("B", side=3, line=0.5, at=-5)

the.M <- pred.m.T.GC.256
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("C", side=3, line=0.5, at=-5)

the.M <- pred.m.T.GC.512
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("D", side=3, line=0.5, at=-5)

dev.off()



### Field ###


svg("accuracy.heatmap.T.fld.svg", width=5, height=6)

layout(matrix(c(0,0,1,2,3,0,4,0,5,0), nrow=5, byrow=T),
       widths = c(5, 1), 
       heights = c(1,5,5,5,5))

max.value <- 5 # max of the matrix
min.value <- 0

the.M <- pred.m.T.fld.64
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("A", side=3, line=0.5, at=-5)

par(mar = c(1.6, 3, 3.7, 2))
colorLevels <- seq(min.value, max.value, length=256)
image(1, colorLevels,
      matrix(colorLevels, ncol = length(colorLevels), nrow = 1),
      col = viridis(256), xlab=" ", ylab="", axes=F, bty="n")
axis(2, at=c(0:5), labels=c("0", "1", "2", "3", "4", ">=5"), las=1)

the.M <- pred.m.T.fld.128
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("B", side=3, line=0.5, at=-5)

the.M <- pred.m.T.fld.256
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("C", side=3, line=0.5, at=-5)

the.M <- pred.m.T.fld.512
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("D", side=3, line=0.5, at=-5)

dev.off()



### Mix ###


svg("accuracy.heatmap.T.mix.svg", width=5, height=6)

layout(matrix(c(0,0,1,2,3,0,4,0,5,0), nrow=5, byrow=T),
       widths = c(5, 1), 
       heights = c(1,5,5,5,5))

max.value <- 5 # max of the matrix
min.value <- 0

the.M <- pred.m.T.mix.64
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("A", side=3, line=0.5, at=-5)

par(mar = c(1.6, 3, 3.7, 2))
colorLevels <- seq(min.value, max.value, length=256)
image(1, colorLevels,
      matrix(colorLevels, ncol = length(colorLevels), nrow = 1),
      col = viridis(256), xlab=" ", ylab="", axes=F, bty="n")
axis(2, at=c(0:5), labels=c("0", "1", "2", "3", "4", ">=5"), las=1)

the.M <- pred.m.T.mix.128
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("B", side=3, line=0.5, at=-5)

the.M <- pred.m.T.mix.256
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("C", side=3, line=0.5, at=-5)

the.M <- pred.m.T.mix.512
par(mar = c(1, 4, 3, 2))
image(1:ncol(the.M),1:nrow(the.M),t(the.M),
      col=viridis(256),
      zlim=c(min.value,max.value),
      asp=1,
      bty="n", 
      xlab="", 
      ylab="",
      axes=F
)
mtext("D", side=3, line=0.5, at=-5)

dev.off()




#=======================#
# Visualize the results #
#=======================#


svg("FIT.evaluate.T.svg", width=6, height=7)

layout(matrix(c(1,2,3,4,5,6), nrow=2, byrow=T),
       widths = c(1, 1, 1), 
       heights = c(1,1))

## no. of outlier (MAE >= 5) genes

plot.pred <- pred.T.GC.nna[, 4:ncol(pred.T.GC.nna)]
plot(rep(log2(c(64, 128, 256, 512)), each=100) + runif(400,-0.1, 0.1),
     apply(plot.pred, 1, function(x) sum(x >= 5)),
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(0,140))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(0, 140, by=10), cex.axis=0.8, lwd=0.5,
     lab=c("0", "", "", "", "", "50", "", "", "", "", "100", "", "", "", ""), las=1)
lines(c(5.9, 6.1), rep(mean(apply(plot.pred[1:100, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(6.9, 7.1), rep(mean(apply(plot.pred[101:200, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(7.9, 8.1), rep(mean(apply(plot.pred[201:300, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(8.9, 9.1), rep(mean(apply(plot.pred[301:400, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
mtext("A", side=3, line=1, at=5.2)
mtext("Number of genes with poor prediction", side=2, line=2.7, at=80, cex=0.8)


plot.pred <- pred.T.fld.nna[, 4:ncol(pred.T.fld.nna)]
plot(rep(log2(c(64, 128, 256, 512)), each=100) + runif(400,-0.1, 0.1),
     apply(plot.pred, 1, function(x) sum(x >= 5)),
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(0,20))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(0, 20, by=5), cex.axis=0.8, lwd=0.5,
     lab=c("0", "", "10", "", "20"), las=1)
lines(c(5.9, 6.1), rep(mean(apply(plot.pred[1:100, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(6.9, 7.1), rep(mean(apply(plot.pred[101:200, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(7.9, 8.1), rep(mean(apply(plot.pred[201:300, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(8.9, 9.1), rep(mean(apply(plot.pred[301:400, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
mtext("B", side=3, line=1, at=5.2)


plot.pred <- pred.T.mix.nna[, 4:ncol(pred.T.mix.nna)]
plot(rep(log2(c(64, 128, 256, 512)), each=100) + runif(400,-0.1, 0.1),
     apply(plot.pred, 1, function(x) sum(x >= 5)),
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(0,20))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(0, 20, by=5), cex.axis=0.8, lwd=0.5,
     lab=c("0", "", "10", "", "20"), las=1)
lines(c(5.9, 6.1), rep(mean(apply(plot.pred[1:100, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(6.9, 7.1), rep(mean(apply(plot.pred[101:200, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(7.9, 8.1), rep(mean(apply(plot.pred[201:300, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(8.9, 9.1), rep(mean(apply(plot.pred[301:400, ], 1, function(x) sum(x >= 5))), 2), lwd=2)
mtext("C", side=3, line=1, at=5.2)



# Median of MAE calculated across genes


plot.pred <- pred.T.GC.nna
plot(log2(plot.pred$smpl.size)+runif(400,-0.1, 0.1), 
     apply(plot.pred, 1, median), 
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(1.4, 2.1))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(1.4, 2.1, by=0.1), lab=c("1.4", "", "1.6", "", "1.8", "", "2.0", ""), las=1, cex.axis=0.8, lwd=0.5)
lines(c(5.9, 6.1), rep(mean(apply(plot.pred, 1, median)[1:100]), 2), lwd=2)
lines(c(6.9, 7.1), rep(mean(apply(plot.pred, 1, median)[101:200]), 2), lwd=2)
lines(c(7.9, 8.1), rep(mean(apply(plot.pred, 1, median)[201:300]), 2), lwd=2)
lines(c(8.9, 9.1), rep(mean(apply(plot.pred, 1, median)[301:400]), 2), lwd=2)
mtext("D", side=3, line=1, at=5.2)
mtext("Median of prediction deviations across 456 genes", side=2, line=2.7, at=1.7, cex=0.8)
mtext("Training data size", side=1, line=3, at=7.5, cex=0.8)

plot.pred <- pred.T.fld.nna
plot(log2(plot.pred$smpl.size)+runif(400,-0.1, 0.1), 
     apply(plot.pred, 1, median), 
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(1.4, 2.1))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(1.4, 2.1, by=0.1), lab=c("1.4", "", "1.6", "", "1.8", "", "2.0", ""), las=1, cex.axis=0.8, lwd=0.5)
lines(c(5.9, 6.1), rep(mean(apply(plot.pred, 1, median)[1:100]), 2), lwd=2)
lines(c(6.9, 7.1), rep(mean(apply(plot.pred, 1, median)[101:200]), 2), lwd=2)
lines(c(7.9, 8.1), rep(mean(apply(plot.pred, 1, median)[201:300]), 2), lwd=2)
lines(c(8.9, 9.1), rep(mean(apply(plot.pred, 1, median)[301:400]), 2), lwd=2)
mtext("E", side=3, line=1, at=5.2)
mtext("Training data size", side=1, line=3, at=7.5, cex=0.8)

plot.pred <- pred.T.mix.nna
plot(log2(plot.pred$smpl.size)+runif(400,-0.1, 0.1), 
     apply(plot.pred, 1, median), 
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(1.4, 2.1))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(1.4, 2.1, by=0.1), lab=c("1.4", "", "1.6", "", "1.8", "", "2.0", ""), las=1, cex.axis=0.8, lwd=0.5)
lines(c(5.9, 6.1), rep(mean(apply(plot.pred, 1, median)[1:100]), 2), lwd=2)
lines(c(6.9, 7.1), rep(mean(apply(plot.pred, 1, median)[101:200]), 2), lwd=2)
lines(c(7.9, 8.1), rep(mean(apply(plot.pred, 1, median)[201:300]), 2), lwd=2)
lines(c(8.9, 9.1), rep(mean(apply(plot.pred, 1, median)[301:400]), 2), lwd=2)
mtext("F", side=3, line=1, at=5.2)
mtext("Training data size", side=1, line=3, at=7.5, cex=0.8)

dev.off()



#=========================================#
# Gene-specific environmental data choice #
#=========================================#

# make sure the reason of prediction failure
# i.e., FIT can "choose" an environmental variable but
#       the coefficient for the variable can be zero.
#       This analysis take this into account.


## GC

t <- proc.time()
env.v.M.GC <- NULL
coef1.v.M.GC <- NULL
coef2.v.M.GC <- NULL
for (i in 301:400) {
  
  load(paste("model.T.GC.", i, sep=""))
  
  env.v <- NULL
  coef1.v <- NULL
  coef2.v <- NULL
  for (j in 1:length(the.model)) {
    
    if (!FG.v[j]) {
      if (the.model[[j]]$coefs[8]==0 && the.model[[j]]$coefs[9]==0) {
        env.v <- c(env.v, NA_character_)
        coef1.v <- c(coef1.v, 0)
        coef2.v <- c(coef2.v, 0)
      } else {
        env.v <- c(env.v, the.model[[j]]$env)
        coef1.v <- c(coef1.v, the.model[[j]]$coefs[8])
        coef2.v <- c(coef2.v, the.model[[j]]$coefs[9])
      }
    } else {}
  }
  env.v.M.GC <- rbind(env.v.M.GC, env.v)
  coef1.v.M.GC <- rbind(coef1.v.M.GC, coef1.v)
  coef2.v.M.GC <- rbind(coef2.v.M.GC, coef2.v)
}

proc.time() - t

dim(env.v.M.GC) # 100 x 456
rownames(env.v.M.GC) <- 1:nrow(env.v.M.GC)
colnames(env.v.M.GC) <- names(FG.v)[!FG.v]

prob.env.choice.GC <- rbind(apply(env.v.M.GC, 2, function(x) sum(x=="temperature", na.rm=T)/length(x)),
                            apply(env.v.M.GC, 2, function(x) sum(x=="radiation", na.rm=T)/length(x)))
rownames(prob.env.choice.GC) <- c("temperature", "radiation")                            



## Field

t <- proc.time()
env.v.M.fld <- NULL
coef1.v.M.fld <- NULL
coef2.v.M.fld <- NULL
for (i in 301:400) {
  
  load(paste("model.T.fld.", i, sep=""))
  
  env.v <- NULL
  coef1.v <- NULL
  coef2.v <- NULL
  for (j in 1:length(the.model)) {
    
    if (!FG.v[j]) {
      if (the.model[[j]]$coefs[8]==0 && the.model[[j]]$coefs[9]==0) {
        env.v <- c(env.v, NA_character_)
        coef1.v <- c(coef1.v, 0)
        coef2.v <- c(coef2.v, 0)
        } else {
          env.v <- c(env.v, the.model[[j]]$env)
          coef1.v <- c(coef1.v, the.model[[j]]$coefs[8])
          coef2.v <- c(coef2.v, the.model[[j]]$coefs[9])
        }
    } else {}
  }
  env.v.M.fld <- rbind(env.v.M.fld, env.v)
  coef1.v.M.fld <- rbind(coef1.v.M.fld, coef1.v)
  coef2.v.M.fld <- rbind(coef2.v.M.fld, coef2.v)
}

proc.time() - t

dim(env.v.M.fld)
rownames(env.v.M.fld) <- 1:nrow(env.v.M.fld)
colnames(env.v.M.fld) <- names(FG.v)[!FG.v]


prob.env.choice.fld <- rbind(apply(env.v.M.fld, 2, function(x) sum(x=="temperature", na.rm=T)/length(x)),
                            apply(env.v.M.fld, 2, function(x) sum(x=="radiation", na.rm=T)/length(x)))
rownames(prob.env.choice.fld) <- c("temperature", "radiation")    


## mix

t <- proc.time()
env.v.M.mix <- NULL
coef1.v.M.mix <- NULL
coef2.v.M.mix <- NULL
for (i in 301:400) {
  
  load(paste("model.T.mix.", i, sep=""))
  
  env.v <- NULL
  coef1.v <- NULL
  coef2.v <- NULL
  for (j in 1:length(the.model)) {
    
    if (!FG.v[j]) {
      if (the.model[[j]]$coefs[8]==0 && the.model[[j]]$coefs[9]==0) {
        env.v <- c(env.v, NA_character_)
        coef1.v <- c(coef1.v, 0)
        coef2.v <- c(coef2.v, 0)
      } else {
        env.v <- c(env.v, the.model[[j]]$env)
        coef1.v <- c(coef1.v, the.model[[j]]$coefs[8])
        coef2.v <- c(coef2.v, the.model[[j]]$coefs[9])
      }
    } else {}
  }
  env.v.M.mix <- rbind(env.v.M.mix, env.v)
  coef1.v.M.mix <- rbind(coef1.v.M.mix, coef1.v)
  coef2.v.M.mix <- rbind(coef2.v.M.mix, coef2.v)
}

proc.time() - t

dim(env.v.M.mix)
rownames(env.v.M.mix) <- 1:nrow(env.v.M.mix)
colnames(env.v.M.mix) <- names(FG.v)[!FG.v]

prob.env.choice.mix <- rbind(apply(env.v.M.mix, 2, function(x) sum(x=="temperature", na.rm=T)/length(x)),
                             apply(env.v.M.mix, 2, function(x) sum(x=="radiation", na.rm=T)/length(x)))
rownames(prob.env.choice.mix) <- c("temperature", "radiation")    


prob.env.choice.GC_Tak <- prob.env.choice.GC
prob.env.choice.fld_Tak <- prob.env.choice.fld
prob.env.choice.mix_Tak <- prob.env.choice.mix

save(prob.env.choice.GC_Tak, file="prob.env.choice.GC_Tak")
save(prob.env.choice.fld_Tak, file="prob.env.choice.fld_Tak")
save(prob.env.choice.mix_Tak, file="prob.env.choice.mix_Tak")


### visualize

## band graph

order.v <- order(prob.env.choice.mix_Tak["radiation", ], decreasing=T)
line.w <- 1
temp.col <- "dodgerblue"
rad.col <- "seagreen1"
none.col <- "gray70"

GC.T.sort <- prob.env.choice.GC_Tak[, order.v]
fld.T.sort <- prob.env.choice.fld_Tak[, order.v]
mix.T.sort <- prob.env.choice.mix_Tak[, order.v]

GC.T.sort <- rbind(GC.T.sort, apply(GC.T.sort, 2, function(x) 1 - sum(x)))
fld.T.sort <- rbind(fld.T.sort, apply(fld.T.sort, 2, function(x) 1 - sum(x)))
mix.T.sort <- rbind(mix.T.sort, apply(mix.T.sort, 2, function(x) 1 - sum(x)))

row.names(GC.T.sort) <- c("temperature", "radiation", "none")
row.names(fld.T.sort) <- c("temperature", "radiation", "none")
row.names(mix.T.sort) <- c("temperature", "radiation", "none")


svg("env.v.choice.band.T.svg", width=5, height=5)

par(mfrow=c(3,1))
par(mar=c(2, 3, 2, 2))
par(xpd=T)

plot(1,1,type="n", xlim=c(0,457), ylim=c(0,1), axes=F, xlab=NA)

the.DF <- fld.T.sort
for(i in 1:456) {
  if (the.DF["radiation", i]==0) {
  } else {
    lines(c(i,i), c(0,the.DF["radiation",i]), 
          lwd=line.w, col=rad.col)
  }
  if (the.DF["temperature", i]==0) {
  } else {
    lines(c(i,i), c(1, 1 - the.DF["temperature",i]), lwd=line.w, col=temp.col)
  }
  if (the.DF["none", i]==0) {
  } else {
    lines(c(i,i), c(the.DF["radiation",i], 1 - the.DF["temperature",i]), lwd=line.w, col=none.col)
  }
}
mtext("Field", side=3, line=0.5, at=5, adj=0, cex=0.8)
axis(2, at=seq(0, 1, by=0.2), las=1, cex.axis=0.8, lwd=0.5, labels=seq(0, 100, by=20))
led.x <- 477
led.wd <- 2
lines(c(led.x, led.x), c(0, 1/3), lwd=led.wd, col=rad.col)
lines(c(led.x, led.x), c(1/3, 2/3), lwd=led.wd, col=none.col)
lines(c(led.x, led.x), c(2/3, 1), lwd=led.wd, col=temp.col)


plot(1,1,type="n", xlim=c(0,457), ylim=c(0,1), axes=F, xlab=NA)

the.DF <- GC.T.sort
for(i in 1:456) {
  if (the.DF["radiation", i]==0) {
  } else {
    lines(c(i,i), c(0,the.DF["radiation",i]), 
          lwd=line.w, col=rad.col)
  }
  if (the.DF["temperature", i]==0) {
  } else {
    lines(c(i,i), c(1, 1 - the.DF["temperature",i]), lwd=line.w, col=temp.col)
  }
  if (the.DF["none", i]==0) {
  } else {
    lines(c(i,i), c(the.DF["radiation",i], 1 - the.DF["temperature",i]), lwd=line.w, col=none.col)
  }
}
mtext("GC", side=3, line=0.5, at=5, adj=0, cex=0.8)
axis(2, at=seq(0, 1, by=0.2), las=1, cex.axis=0.8, lwd=0.5, labels=seq(0, 100, by=20))



plot(1,1,type="n", xlim=c(0,457), ylim=c(0,1), axes=F, xlab=NA)

the.DF <- mix.T.sort
for(i in 1:456) {
  if (the.DF["radiation", i]==0) {
  } else {
    lines(c(i,i), c(0,the.DF["radiation",i]), 
          lwd=line.w, col=rad.col)
  }
  if (the.DF["temperature", i]==0) {
  } else {
    lines(c(i,i), c(1, 1 - the.DF["temperature",i]), lwd=line.w, col=temp.col)
  }
  if (the.DF["none", i]==0) {
  } else {
    lines(c(i,i), c(the.DF["radiation",i], 1 - the.DF["temperature",i]), lwd=line.w, col=none.col)
  }
}
mtext("Mix", side=3, line=0.5, at=5, adj=0, cex=0.8)
axis(2, at=seq(0, 1, by=0.2), las=1, cex.axis=0.8, lwd=0.5, labels=seq(0, 100, by=20))

dev.off()


## histogram

svg("env.v.choice.T.svg", width=7, height=3.2)

par(mfcol=c(2,3))
line <- par(lwd=0.5)

L.mar <- c(4.5, 5, 4, 0.5)   # left panel margin
M.mar <- c(4.5, 3, 4, 2.5)
R.mar <- c(4.5, 2.5, 4, 3)
PL.line <- 1.7               # panel label line
PL.x    <- -0.35              # panel label x
xlab.line <- 2.5

par(mar=L.mar)
hist(prob.env.choice.GC_Tak["temperature", ], col="dodgerblue",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("A", side=3, line=PL.line, at=PL.x)
mtext("Training data: smart GC", side=3, line=PL.line, at=PL.x + 0.15, adj=0, cex=0.8)
mtext("Temperature chosen", side=2, line=3, at=150, cex=0.7)
mtext("Number of subsamplings", side=1, line=xlab.line, at=0.5, cex=0.7)

par(mar=L.mar)
hist(prob.env.choice.GC_Tak["radiation", ],  col="seagreen1",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("Radiation chosen", side=2, line=3, at=150, cex=0.7)
mtext("Number of subsamplings", side=1, line=xlab.line, at=0.5, cex=0.7)


par(mar=M.mar)
hist(prob.env.choice.fld_Tak["temperature", ], col="dodgerblue",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("B", side=3, line=PL.line, at=PL.x)
mtext("Training data: field", side=3, line=PL.line, at=PL.x + 0.15, adj=0, cex=0.8)
mtext("Number of subsamplings", side=1, line=xlab.line, at=0.5, cex=0.7)

par(mar=M.mar)
hist(prob.env.choice.fld_Tak["radiation", ], col="seagreen1",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("Number of subsamplings", side=1, line=xlab.line, at=0.5, cex=0.7)


par(mar=R.mar)
hist(prob.env.choice.mix_Tak["temperature", ], col="dodgerblue",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("C", side=3, line=PL.line, at=PL.x)
mtext("Training data: mix", side=3, line=PL.line, at=PL.x + 0.15, adj=0, cex=0.8)
mtext("Number of subsamplings", side=1, line=xlab.line, at=0.5, cex=0.7)

par(mar=R.mar)
hist(prob.env.choice.mix_Tak["radiation", ], col="seagreen1",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("Frequency radiation chosen", side=1, line=xlab.line, at=0.5, cex=0.7)

dev.off()

