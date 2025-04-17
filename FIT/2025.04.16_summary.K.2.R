

# summarize the results of Koshihikari

setwd("Kos.result.mix.2")

load("pred.neg.mix.K")   # 100 x 474 (5-474: 470 genes)
load("pred.pos.mix.K")   # 100 x 474 (5-474: 470 genes)
load("pred.mix.mix.K")   # 100 x 474 (5-474: 470 genes)

#======================================#
# Remove genes with prediction failure #
#======================================#

FG.neg <- apply(pred.neg.mix.K[,5:474], 2, function(x) sum(is.na(x))>=1)   # 353
FG.pos <- apply(pred.pos.mix.K[,5:474], 2, function(x) sum(is.na(x))>=1)   # 353
FG.mix <- apply(pred.mix.mix.K[,5:474], 2, function(x) sum(is.na(x))>=1)   # 353

FG.v <- apply(rbind(FG.neg, FG.pos, FG.mix), 2, function(x) sum(x)>=1) # 469 genes

pred.K.neg.nna <- pred.neg.mix.K[,5:474][, !FG.v]
pred.K.pos.nna <- pred.pos.mix.K[,5:474][, !FG.v]
pred.K.mix.nna <- pred.mix.mix.K[,5:474][, !FG.v]

pred.K.neg.nna <- cbind(pred.neg.mix.K[, c(1,2,4)], pred.K.neg.nna) # 100 x 472 (3 + 469)
pred.K.neg.nna$MAE <- apply(pred.K.neg.nna[, 4:472], 1, mean)

pred.K.pos.nna <- cbind(pred.pos.mix.K[, c(1,2,4)], pred.K.pos.nna)
pred.K.pos.nna$MAE <- apply(pred.K.pos.nna[, 4:472], 1, mean)

pred.K.mix.nna <- cbind(pred.mix.mix.K[, c(1,2,4)], pred.K.mix.nna)
pred.K.mix.nna$MAE <- apply(pred.K.mix.nna[, 4:472], 1, mean)


#===================#
# Process the files #
#===================#

### negative ###
pred.m.K.neg  <- as.matrix(pred.K.neg.nna[ , 4:472])
min(pred.m.K.neg) ; max(pred.m.K.neg)     # 0.08860306, 151572009
pred.m.K.neg[pred.m.K.neg >= 5] <- 5

### positive ###
pred.m.K.pos  <- as.matrix(pred.K.pos.nna[ , 4:472])
min(pred.m.K.pos) ; max(pred.m.K.pos)     # 0.08172196, 3634.527
pred.m.K.pos[pred.m.K.pos >= 5] <- 5

### negative ###
pred.m.K.mix <- as.matrix(pred.K.mix.nna[ , 4:472])
min(pred.m.K.mix) ; max(pred.m.K.mix)     # 0.08509356, 12893.97
pred.m.K.mix[pred.m.K.mix >= 5] <- 5



#=======================#
# Visualize the results #
#=======================#

library(viridis)

svg("FIT.evaluate.K.2.svg", width=5, height=3)

par(mfrow=c(1,2))
par(mar=c(4, 4, 4, 1))

## no. of outlier (MAE >= 5) genes

plot(rep(c(1, 2, 3), each=100) + runif(300,-0.1, 0.1),
     c(apply(pred.K.neg.nna[, 4:ncol(pred.K.neg.nna)], 1, function(x) sum(x >= 5)),
       apply(pred.K.pos.nna[, 4:ncol(pred.K.pos.nna)], 1, function(x) sum(x >= 5)),
       apply(pred.K.mix.nna[, 4:ncol(pred.K.mix.nna)], 1, function(x) sum(x >= 5))),
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(0,10))
axis(1, at=c(1:3), lab=c("Negative", "Positive", "Mix"), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(0, 10, by=2), cex.axis=0.8, lwd=0.5,
     lab=seq(0, 10, by=2), las=1)
lines(c(0.9, 1.1), rep(mean(apply(pred.K.neg.nna[, 4:ncol(pred.K.neg.nna)], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(1.9, 2.1), rep(mean(apply(pred.K.pos.nna[, 4:ncol(pred.K.pos.nna)], 1, function(x) sum(x >= 5))), 2), lwd=2)
lines(c(2.9, 3.1), rep(mean(apply(pred.K.mix.nna[, 4:ncol(pred.K.mix.nna)], 1, function(x) sum(x >= 5))), 2), lwd=2)

mtext("A", side=3, line=1, at=0.8)
mtext("Number of genes with\n poor prediction", side=2, line=2.3, at=5, cex=0.8)


# Median of MAE calculated across genes


plot(rep(c(1, 2, 3), each=100) + runif(300,-0.1, 0.1), 
     c(apply(pred.K.neg.nna[, 4:ncol(pred.K.neg.nna)], 1, median),
       apply(pred.K.pos.nna[, 4:ncol(pred.K.pos.nna)], 1, median),
       apply(pred.K.mix.nna[, 4:ncol(pred.K.mix.nna)], 1, median)),
     pch=19, col=rgb(0, 1, 1, alpha=0.3), axes=F, xlab=NA, ylab=NA, ylim=c(1.4, 2))
axis(1, at=c(1,2,3), lab=c("Negative", "Positive", "Mix"), cex.axis=0.8, lwd=0.5)
axis(2, at=seq(1.4, 2, by=0.1), lab=c("1.4", "", "1.6", "", "1.8", "", "2.0"), las=1, cex.axis=0.8, lwd=0.5)
lines(c(0.9, 1.1), rep(mean(apply(pred.K.neg.nna[, 4:ncol(pred.K.neg.nna)], 1, median)), 2), lwd=2)
lines(c(1.9, 2.1), rep(mean(apply(pred.K.pos.nna[, 4:ncol(pred.K.pos.nna)], 1, median)), 2), lwd=2)
lines(c(2.9, 3.1), rep(mean(apply(pred.K.mix.nna[, 4:ncol(pred.K.mix.nna)], 1, median)), 2), lwd=2)
mtext("B", side=3, line=1, at=0.8)
mtext("Median of prediction deviations\n across 469 genes", side=2, line=2.3, at=1.7, cex=0.8)

dev.off()

