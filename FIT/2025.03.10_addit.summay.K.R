

# summarize the results of Koshihikari

setwd("Kos.result.GC")
load("pred.K.GC")   # 400 x 474 (5-474: genes)

setwd("../Kos.result.fld")
load("pred.K.fld")  # 400 x 474 (5-474: genes)

setwd("../Kos.result.mix")
load("pred.K.mix")  # 400 x 474 (5-474: genes)


#======================================#
# Remove genes with prediction failure #
#======================================#

FG.GC  <- apply(pred.K.GC[,5:474], 2, function(x) sum(is.na(x))>=1)   # 299, 353
FG.fld <- apply(pred.K.fld[,5:474], 2, function(x) sum(is.na(x))>=1)  #      353, 395, 440
FG.mix <- apply(pred.K.mix[,5:474], 2, function(x) sum(is.na(x))>=1)  #      353, 395

FG.v <- apply(rbind(FG.GC, FG.fld, FG.mix), 2, function(x) sum(x)>=1) # 470 genes

pred.K.GC.nna  <- pred.K.GC[301:400, 5:474][, !FG.v]
pred.K.fld.nna <- pred.K.fld[301:400, 5:474][, !FG.v]
pred.K.mix.nna <- pred.K.mix[301:400, 5:474][, !FG.v]

#=============================#
# Gene-wise median difference #
#=============================#

pred.K.GC.med  <- apply(pred.K.GC.nna,  2, median) # 466 genes
min(pred.K.GC.med) ; max(pred.K.GC.med) # 0.2241169, 9920.206
pred.K.GC.med[pred.K.GC.med >= 5] <- 5

pred.K.fld.med <- apply(pred.K.fld.nna, 2, median) # 466 genes
min(pred.K.fld.med) ; max(pred.K.fld.med) # 0.1815851, 3.155332
pred.K.fld.med[pred.K.fld.med >= 5] <- 5

pred.K.mix.med <- apply(pred.K.mix.nna, 2, median) # 466 genes
min(pred.K.mix.med) ; max(pred.K.mix.med) # 0.203237, 3.228288
pred.K.mix.med[pred.K.mix.med >= 5] <- 5


#===========#
# Visualize #
#===========#

setwd("..")
svg("genewise.medians.K.svg", width=6.9, height=2.7)

par(mfrow=c(1,3))

plot(pred.K.mix.med, pred.K.GC.med, pch=19, asp=1, xlim=c(0, 5), ylim=c(0, 5),
     col=rgb(0, 0.7, 0.7, alpha=0.2), axes=F,
     xlab="Mix-model gene-wise medians", ylab="GC-model gene-wise medians")
lines(c(0, 5), c(0, 5), lty=2)
axis(1, at=c(0:5), lab=c(0:5), cex.axis=0.8, lwd=0.5)
axis(2, at=c(0:5), cex.axis=0.8, lwd=0.5,
     lab=c(0:5), las=1)

plot(pred.K.mix.med, pred.K.fld.med, pch=19, asp=1, xlim=c(0, 5), ylim=c(0, 5),
     col=rgb(0, 0.7, 0.7, alpha=0.2), axes=F,
     xlab="Mix-model gene-wise medians", ylab="Field-model gene-wise medians")
lines(c(0, 5), c(0, 5), lty=2)
axis(1, at=c(0:5), lab=c(0:5), cex.axis=0.8, lwd=0.5)
axis(2, at=c(0:5), cex.axis=0.8, lwd=0.5,
     lab=c(0:5), las=1)

plot(pred.K.fld.med, pred.K.GC.med, pch=19, asp=1, xlim=c(0, 5), ylim=c(0, 5),
     col=rgb(0, 0.7, 0.7, alpha=0.2), axes=F,
     xlab="Field-model gene-wise medians", ylab="GC-model gene-wise medians")
lines(c(0, 5), c(0, 5), lty=2)
axis(1, at=c(0:5), lab=c(0:5), cex.axis=0.8, lwd=0.5)
axis(2, at=c(0:5), cex.axis=0.8, lwd=0.5,
     lab=c(0:5), las=1)


dev.off()

save(pred.K.mix.med, file="pred.K.mix.med")
save(pred.K.fld.med, file="pred.K.fld.med")
save(pred.K.GC.med,  file="pred.K.GC.med")


