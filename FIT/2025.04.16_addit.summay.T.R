

# summarize the results of Koshihikari

setwd("Tak.result.GC")
load("pred.T.GC")   # 400 x 474 (5-474: genes)

setwd("../Tak.result.fld")
load("pred.T.fld")  # 400 x 474 (5-474: genes)

setwd("../Tak.result.mix")
load("pred.T.mix")  # 400 x 474 (5-474: genes)


#======================================#
# Remove genes with prediction failure #
#======================================#

FG.GC  <- apply(pred.T.GC[,5:474], 2, function(x) sum(is.na(x))>=1)   # 90, 134, 242, 256, 270, 287, 295, 303, 334,           396, 429
FG.fld <- apply(pred.T.fld[,5:474], 2, function(x) sum(is.na(x))>=1)  # 90,                270,      295, 303, 334, 353, 395,           440
FG.mix <- apply(pred.T.mix[,5:474], 2, function(x) sum(is.na(x))>=1)  # 90, 134, 242,      270,      295, 303, 334,           396, 429

FG.v <- apply(rbind(FG.GC, FG.fld, FG.mix), 2, function(x) sum(x)>=1) # 470 genes (14 genes failed)

pred.T.GC.nna  <- pred.T.GC[301:400,  5:474][, !FG.v]
pred.T.fld.nna <- pred.T.fld[301:400, 5:474][, !FG.v]
pred.T.mix.nna <- pred.T.mix[301:400, 5:474][, !FG.v]

#=============================#
# Gene-wise median difference #
#=============================#

pred.T.GC.med  <- apply(pred.T.GC.nna,  2, median) # 456 genes
min(pred.T.GC.med) ; max(pred.T.GC.med) # 0.04775849, 17520.48
pred.T.GC.med[pred.T.GC.med >= 5] <- 5

pred.T.fld.med <- apply(pred.T.fld.nna, 2, median) # 456 genes
min(pred.T.fld.med) ; max(pred.T.fld.med) # 0.05622028, 3.072547
pred.T.fld.med[pred.T.fld.med >= 5] <- 5

pred.T.mix.med <- apply(pred.T.mix.nna, 2, median) # 456 genes
min(pred.T.mix.med) ; max(pred.T.mix.med) # 0.0505758, 3.054804
pred.T.mix.med[pred.T.mix.med >= 5] <- 5


#===========#
# Visualize #
#===========#

setwd("..")
svg("genewise.medians.T.svg", width=6.9, height=2.7)

par(mfrow=c(1,3))

plot(pred.T.mix.med, pred.T.GC.med, pch=19, asp=1, xlim=c(0, 5), ylim=c(0, 5),
     col=rgb(0, 0.7, 0.7, alpha=0.2), axes=F,
     xlab="Mix-model gene-wise medians", ylab="GC-model gene-wise medians")
lines(c(0, 5), c(0, 5), lty=2)
axis(1, at=c(0:5), lab=c(0:5), cex.axis=0.8, lwd=0.5)
axis(2, at=c(0:5), cex.axis=0.8, lwd=0.5,
     lab=c(0:5), las=1)

plot(pred.T.mix.med, pred.T.fld.med, pch=19, asp=1, xlim=c(0, 5), ylim=c(0, 5),
     col=rgb(0, 0.7, 0.7, alpha=0.2), axes=F,
     xlab="Mix-model gene-wise medians", ylab="Field-model gene-wise medians")
lines(c(0, 5), c(0, 5), lty=2)
axis(1, at=c(0:5), lab=c(0:5), cex.axis=0.8, lwd=0.5)
axis(2, at=c(0:5), cex.axis=0.8, lwd=0.5,
     lab=c(0:5), las=1)

plot(pred.T.fld.med, pred.T.GC.med, pch=19, asp=1, xlim=c(0, 5), ylim=c(0, 5),
     col=rgb(0, 0.7, 0.7, alpha=0.2), axes=F,
     xlab="Field-model gene-wise medians", ylab="GC-model gene-wise medians")
lines(c(0, 5), c(0, 5), lty=2)
axis(1, at=c(0:5), lab=c(0:5), cex.axis=0.8, lwd=0.5)
axis(2, at=c(0:5), cex.axis=0.8, lwd=0.5,
     lab=c(0:5), las=1)


dev.off()

save(pred.T.mix.med, file="pred.T.mix.med")
save(pred.T.fld.med, file="pred.T.fld.med")
save(pred.T.GC.med,  file="pred.T.GC.med")


