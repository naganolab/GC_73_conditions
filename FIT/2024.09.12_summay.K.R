

# summarize the results of Koshihikari

setwd("Kos.result")

load("pred.K.GC")   # 400 x 474 (5-474: genes)
load("pred.K.fld")  # 400 x 474 (5-474: genes)
load("pred.K.mix")  # 400 x 474 (5-474: genes)


#======================================#
# Remove genes with prediction failure #
#======================================#

FG.GC  <- apply(pred.K.GC[,5:474], 2, function(x) sum(is.na(x))>=1)   # 299, 353
FG.fld <- apply(pred.K.fld[,5:474], 2, function(x) sum(is.na(x))>=1)  #      353, 395, 440
FG.mix <- apply(pred.K.mix[,5:474], 2, function(x) sum(is.na(x))>=1)  #      353, 395

FG.v <- apply(rbind(FG.GC, FG.fld, FG.mix), 2, function(x) sum(x)>=1) # 470 genes

pred.K.GC.nna  <- pred.K.GC[,5:474][, !FG.v]
pred.K.fld.nna <- pred.K.fld[,5:474][, !FG.v]
pred.K.mix.nna <- pred.K.mix[,5:474][, !FG.v]

pred.K.GC.nna <- cbind(pred.K.GC[, c(1,2,4)], pred.K.GC.nna) # 400 x (3 + 466)
pred.K.GC.nna$MAE <- apply(pred.K.GC.nna[, 4:469], 1, mean)

pred.K.fld.nna <- cbind(pred.K.fld[, c(1,2,4)], pred.K.fld.nna)
pred.K.fld.nna$MAE <- apply(pred.K.fld.nna[, 4:469], 1, mean)

pred.K.mix.nna <- cbind(pred.K.mix[, c(1,2,4)], pred.K.mix.nna)
pred.K.mix.nna$MAE <- apply(pred.K.mix.nna[, 4:469], 1, mean)



#===================#
# Process the files #
#===================#

### GC ###

pred.m.K.GC.64  <- as.matrix(pred.K.GC.nna[1:100, 4:469])
pred.m.K.GC.128 <- as.matrix(pred.K.GC.nna[101:200, 4:469])
pred.m.K.GC.256 <- as.matrix(pred.K.GC.nna[201:300, 4:469])
pred.m.K.GC.512 <- as.matrix(pred.K.GC.nna[301:400, 4:469])

min(pred.m.K.GC.64) ; max(pred.m.K.GC.64)     # 0.1226628, 29323438
min(pred.m.K.GC.128) ; max(pred.m.K.GC.128)   # 0.1499764, 18631212
min(pred.m.K.GC.256) ; max(pred.m.K.GC.256)   # 0.1586678, 10708499
min(pred.m.K.GC.512) ; max(pred.m.K.GC.512)   # 0.202486,  7426871

pred.m.K.GC.64[pred.m.K.GC.64 >= 5] <- 5
pred.m.K.GC.128[pred.m.K.GC.128 >= 5] <- 5
pred.m.K.GC.256[pred.m.K.GC.256 >= 5] <- 5
pred.m.K.GC.512[pred.m.K.GC.512 >= 5] <- 5


### Field ###

pred.m.K.fld.64  <- as.matrix(pred.K.fld.nna[1:100, 4:469])
pred.m.K.fld.128 <- as.matrix(pred.K.fld.nna[101:200, 4:469])
pred.m.K.fld.256 <- as.matrix(pred.K.fld.nna[201:300, 4:469])
pred.m.K.fld.512 <- as.matrix(pred.K.fld.nna[301:400, 4:469])

min(pred.m.K.fld.64) ; max(pred.m.K.fld.64)     # 0.1226628, 1.221456e+12
min(pred.m.K.fld.128) ; max(pred.m.K.fld.128)   # 0.1262708, 130609.6
min(pred.m.K.fld.256) ; max(pred.m.K.fld.256)   # 0.1378236, 16513.17
min(pred.m.K.fld.512) ; max(pred.m.K.fld.512)   # 0.1567659, 6.005579e+14

pred.m.K.fld.64[pred.m.K.fld.64 >= 5] <- 5
pred.m.K.fld.128[pred.m.K.fld.128 >= 5] <- 5
pred.m.K.fld.256[pred.m.K.fld.256 >= 5] <- 5
pred.m.K.fld.512[pred.m.K.fld.512 >= 5] <- 5


### Mix ###

pred.m.K.mix.64  <- as.matrix(pred.K.mix.nna[1:100, 4:469])
pred.m.K.mix.128 <- as.matrix(pred.K.mix.nna[101:200, 4:469])
pred.m.K.mix.256 <- as.matrix(pred.K.mix.nna[201:300, 4:469])
pred.m.K.mix.512 <- as.matrix(pred.K.mix.nna[301:400, 4:469])

min(pred.m.K.mix.64) ; max(pred.m.K.mix.64)     # 0.1226628, 9.102094e+14
min(pred.m.K.mix.128) ; max(pred.m.K.mix.128)   # 0.1444141, 22687.78
min(pred.m.K.mix.256) ; max(pred.m.K.mix.256)   # 0.1429115, 571488659722
min(pred.m.K.mix.512) ; max(pred.m.K.mix.512)   # 0.165236,  625.3619

pred.m.K.mix.64[pred.m.K.mix.64 >= 5] <- 5
pred.m.K.mix.128[pred.m.K.mix.128 >= 5] <- 5
pred.m.K.mix.256[pred.m.K.mix.256 >= 5] <- 5
pred.m.K.mix.512[pred.m.K.mix.512 >= 5] <- 5



#==================#
# Make the heatmap #
#==================#

library(viridis)


### GC ###

svg("accuracy.heatmap.K.GC.svg", width=5, height=6)

layout(matrix(c(0,0,1,2,3,0,4,0,5,0), nrow=5, byrow=T),
       widths = c(5, 1), 
       heights = c(1,5,5,5,5))

max.value <- 5 # max of the matrix
min.value <- 0

the.M <- pred.m.K.GC.64
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

the.M <- pred.m.K.GC.128
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

the.M <- pred.m.K.GC.256
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

the.M <- pred.m.K.GC.512
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


svg("accuracy.heatmap.K.fld.svg", width=5, height=6)

layout(matrix(c(0,0,1,2,3,0,4,0,5,0), nrow=5, byrow=T),
       widths = c(5, 1), 
       heights = c(1,5,5,5,5))

max.value <- 5 # max of the matrix
min.value <- 0

the.M <- pred.m.K.fld.64
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

the.M <- pred.m.K.fld.128
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

the.M <- pred.m.K.fld.256
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

the.M <- pred.m.K.fld.512
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


svg("accuracy.heatmap.K.mix.svg", width=5, height=6)

layout(matrix(c(0,0,1,2,3,0,4,0,5,0), nrow=5, byrow=T),
       widths = c(5, 1), 
       heights = c(1,5,5,5,5))

max.value <- 5 # max of the matrix
min.value <- 0

the.M <- pred.m.K.mix.64
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

the.M <- pred.m.K.mix.128
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

the.M <- pred.m.K.mix.256
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

the.M <- pred.m.K.mix.512
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


svg("FIT.evaluate.K.svg", width=6, height=7)

layout(matrix(c(1,2,3,4,5,6), nrow=2, byrow=T),
       widths = c(1, 1, 1), 
       heights = c(1,1))

## no. of outlier (MAE >= 5) genes

plot.pred <- pred.K.GC.nna[, 4:ncol(pred.K.GC.nna)]
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


plot.pred <- pred.K.fld.nna[, 4:ncol(pred.K.fld.nna)]
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


plot.pred <- pred.K.mix.nna[, 4:ncol(pred.K.mix.nna)]
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


plot.pred <- pred.K.GC.nna
plot(log2(plot.pred$smpl.size)+runif(400,-0.1, 0.1), 
     apply(plot.pred, 1, median), 
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(1.4, 2))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(1.4, 2, by=0.1), lab=c("1.4", "", "1.6", "", "1.8", "", "2.0"), las=1, cex.axis=0.8, lwd=0.5)
lines(c(5.9, 6.1), rep(mean(apply(plot.pred, 1, median)[1:100]), 2), lwd=2)
lines(c(6.9, 7.1), rep(mean(apply(plot.pred, 1, median)[101:200]), 2), lwd=2)
lines(c(7.9, 8.1), rep(mean(apply(plot.pred, 1, median)[201:300]), 2), lwd=2)
lines(c(8.9, 9.1), rep(mean(apply(plot.pred, 1, median)[301:400]), 2), lwd=2)
mtext("D", side=3, line=1, at=5.2)
mtext("Median of prediction deviations across 466 genes", side=2, line=2.7, at=1.7, cex=0.8)
mtext("Training data size", side=1, line=3, at=7.5, cex=0.8)

plot.pred <- pred.K.fld.nna
plot(log2(plot.pred$smpl.size)+runif(400,-0.1, 0.1), 
     apply(plot.pred, 1, median), 
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(1.4, 2))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(1.4, 2, by=0.1), lab=c("1.4", "", "1.6", "", "1.8", "", "2.0"), las=1, cex.axis=0.8, lwd=0.5)
lines(c(5.9, 6.1), rep(mean(apply(plot.pred, 1, median)[1:100]), 2), lwd=2)
lines(c(6.9, 7.1), rep(mean(apply(plot.pred, 1, median)[101:200]), 2), lwd=2)
lines(c(7.9, 8.1), rep(mean(apply(plot.pred, 1, median)[201:300]), 2), lwd=2)
lines(c(8.9, 9.1), rep(mean(apply(plot.pred, 1, median)[301:400]), 2), lwd=2)
mtext("E", side=3, line=1, at=5.2)
mtext("Training data size", side=1, line=3, at=7.5, cex=0.8)

plot.pred <- pred.K.mix.nna
plot(log2(plot.pred$smpl.size)+runif(400,-0.1, 0.1), 
     apply(plot.pred, 1, median), 
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(1.4, 2))
axis(1, at=c(6:9), lab=c(64, 128, 256, 512), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(1.4, 2, by=0.1), lab=c("1.4", "", "1.6", "", "1.8", "", "2.0"), las=1, cex.axis=0.8, lwd=0.5)
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

# Specify the use of environmental data
# i.e., FIT can "choose" an environmental variable but
#       the coefficient for the variable can be zero.
#       This analysis take this into account.


## GC

t <- proc.time()
env.v.M.GC <- NULL
coef1.v.M.GC <- NULL
coef2.v.M.GC <- NULL
for (i in 301:400) {
  
  load(paste("model.K.GC.", i, sep=""))
  
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

dim(env.v.M.GC) # 100 x 466
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
  
  load(paste("model.K.fld.", i, sep=""))
  
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
  
  load(paste("model.K.mix.", i, sep=""))
  
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

prob.env.choice.GC_Kos <- prob.env.choice.GC
prob.env.choice.fld_Kos <- prob.env.choice.fld
prob.env.choice.mix_Kos <- prob.env.choice.mix

save(prob.env.choice.GC_Kos, file="prob.env.choice.GC_Kos")
save(prob.env.choice.fld_Kos, file="prob.env.choice.fld_Kos")
save(prob.env.choice.mix_Kos, file="prob.env.choice.mix_Kos")


### visualize

## band graph

order.v <- order(prob.env.choice.mix_Kos["radiation", ], decreasing=T)
line.w <- 1
temp.col <- "dodgerblue"
rad.col <- "seagreen1"
none.col <- "gray70"

GC.K.sort <- prob.env.choice.GC_Kos[, order.v]
fld.K.sort <- prob.env.choice.fld_Kos[, order.v]
mix.K.sort <- prob.env.choice.mix_Kos[, order.v]

GC.K.sort <- rbind(GC.K.sort, apply(GC.K.sort, 2, function(x) 1 - sum(x)))
fld.K.sort <- rbind(fld.K.sort, apply(fld.K.sort, 2, function(x) 1 - sum(x)))
mix.K.sort <- rbind(mix.K.sort, apply(mix.K.sort, 2, function(x) 1 - sum(x)))

row.names(GC.K.sort) <- c("temperature", "radiation", "none")
row.names(fld.K.sort) <- c("temperature", "radiation", "none")
row.names(mix.K.sort) <- c("temperature", "radiation", "none")


svg("env.v.choice.band.K.svg", width=5, height=5)

par(mfrow=c(3,1))
par(mar=c(2, 3, 2, 2))
par(xpd=T)

plot(1,1,type="n", xlim=c(0,467), ylim=c(0,1), axes=F, xlab=NA)

the.DF <- fld.K.sort
for(i in 1:466) {
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
axis(2, at=seq(0, 1, by=0.2), las=1, cex.axis=0.8, lwd=0.5,
     labels=seq(0, 100, by=20))
led.x <- 487
led.wd <- 2
lines(c(led.x, led.x), c(0, 1/3), lwd=led.wd, col=rad.col)
lines(c(led.x, led.x), c(1/3, 2/3), lwd=led.wd, col=none.col)
lines(c(led.x, led.x), c(2/3, 1), lwd=led.wd, col=temp.col)


plot(1,1,type="n", xlim=c(0,467), ylim=c(0,1), axes=F, xlab=NA)

the.DF <- GC.K.sort
for(i in 1:466) {
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



plot(1,1,type="n", xlim=c(0,467), ylim=c(0,1), axes=F, xlab=NA)

the.DF <- mix.K.sort
for(i in 1:466) {
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

svg("env.v.choice.K.svg", width=7, height=3.2)

par(mfcol=c(2,3))
line <- par(lwd=0.5)

L.mar <- c(4.5, 5, 4, 0.5)   # left panel margin
R.mar <- c(4.5, 2.5, 4, 3)
PL.line <- 1.7               # panel label line
PL.x    <- -0.35              # panel label x
xlab.line <- 2.5

par(mar=L.mar)
hist(prob.env.choice.GC_Kos["temperature", ], col="dodgerblue",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("A", side=3, line=PL.line, at=PL.x)
mtext("Training data: smart GC", side=3, line=PL.line, at=PL.x + 0.15, adj=0, cex=0.8)
mtext("Temperature chosen", side=2, line=3, at=150, cex=0.7)
mtext("Number of subsamplings", side=1, line=xlab.line, at=0.5, cex=0.7)

par(mar=L.mar)
hist(prob.env.choice.GC_Kos["radiation", ], col="seagreen1",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("Radiation chosen", side=2, line=3, at=150, cex=0.7)
mtext("Number of subsampling", side=1, line=xlab.line, at=0.5, cex=0.7)


par(mar=R.mar)
hist(prob.env.choice.fld_Kos["temperature", ], col="dodgerblue",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("B", side=3, line=PL.line, at=PL.x)
mtext("Training data: field", side=3, line=PL.line, at=PL.x + 0.15, adj=0, cex=0.8)
mtext("Number of subsampling", side=1, line=xlab.line, at=0.5, cex=0.7)

par(mar=R.mar)
hist(prob.env.choice.fld_Kos["radiation", ], col="seagreen1",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("Number of subsampling", side=1, line=xlab.line, at=0.5, cex=0.7)


par(mar=R.mar)
hist(prob.env.choice.mix_Kos["temperature", ],  col="dodgerblue",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("C", side=3, line=PL.line, at=PL.x)
mtext("Training data: mix", side=3, line=PL.line, at=PL.x + 0.15, adj=0, cex=0.8)
mtext("Number of subsampling", side=1, line=xlab.line, at=0.5, cex=0.7)

par(mar=R.mar)
hist(prob.env.choice.mix_Kos["radiation", ], col="seagreen1",
     breaks=seq(0, 1, by=0.1), ylim=c(0, 350), axes=F, xlab=NA, ylab=NA, main=NA)
axis(1, at=seq(0, 1, by=0.1), lab=c("0", "", "20", "", "40", "", "60", "", "80", "", "100"), cex.axis =0.8, lwd=0.5)
axis(2, at=seq(0, 350, by=50), lab=c("0", "", "100", "", "200", "", "300", ""), las=1, cex.axis=0.8, lwd=0.5)
mtext("Number of subsampling", side=1, line=xlab.line, at=0.5, cex=0.7)

dev.off()

