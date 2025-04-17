
library("FIT")
library("foreach")
library("doParallel")
library("dplyr")

load("genes") # character vector of length 474
load("log.rpm.K")

load("weather.01")
load("attribute.all.K")

genes <- genes[genes %in% colnames(log.rpm.K)] #470 genes

setwd("Kos.result.mix.2")
smpl.size <- rep(144*2, each=100)


#==================#
# define functions #
#==================#

# Function to find prediction failures
failed.gene.extractor <- function(x, n, genes.v) { # x: list of FIT prediction, n: expected number of prediction points
  failed.gene.id <- NULL
  for (i in 1:length(genes.v)) {
    pred.subset <- x[[i]]
    if(length(pred.subset)!=n) {
      failed.gene.id <- c(failed.gene.id, i) } else {}
  }
  failed.gene.id
}


# Function to remove predictions for some genes and to unlist the FIT prediction
unlist.subset <- function(x, v, genes.v) { # x: list, v: vector giving components to be removed in integer
  res.v <- NULL
  for (i in 1:length(genes.v)) {
    pred.subset <- x[[i]]
    if (i %in% v) { 
    } else {
      res.v <- c(res.v, pred.subset)
    }
  }
  res.v
}


#==================#
# make predictions #
#==================#

oldest.K <- min(as.POSIXct(paste(paste(attribute.all.K$year, attribute.all.K$month, attribute.all.K$day, sep="-"),
                paste(attribute.all.K$hour, attribute.all.K$min, "00", sep=":"), sep=" ")))
weather.K <- weather.01[weather.01$time >= oldest.K + as.difftime(-72, units="hours"),]
oldest.wt.K   <- min(weather.K$time)
time.K   <- seq(0, (nrow(weather.K)-1)*10, by=10)
weather.K <- cbind(time.K, weather.K[,c("temperature","radiation")])
colnames(weather.K)   <- c("time", "temperature", "radiation")
train.weather.K   <- FIT::convert.weather(weather.K,   c("temperature", "radiation"))

K.tst <- 1390:2004
attribute.tst.K <- attribute.all.K[K.tst, ]
attribute.tst.K$time <- as.numeric(as.POSIXct(paste(paste(attribute.tst.K$year, attribute.tst.K$month, attribute.tst.K$day, sep="-"),
                                                    paste(attribute.tst.K$hour, attribute.tst.K$min, "00", sep=":"), 
                                                    sep=" "))-oldest.wt.K, units="mins")
tst.attribute.K <- FIT::convert.attribute(attribute.tst.K)


#==========#
# negative #
#==========#

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=6)  # Choose the number of threads you'd like to use

t <- proc.time()
result <- foreach(i = 1:100, .combine="rbind", .inorder=FALSE, .packages = "dplyr") %dopar% {
  
  load(paste("model.neg.mix.K.", i, sep=""))

  pred.K <- FIT::predict(the.model, tst.attribute.K, train.weather.K)
  # genes with prediction failure
  FGID.K   <- failed.gene.extractor(pred.K, length(K.tst), genes)
  FGID <- FGID.K

  pred.v.K <- unlist.subset(pred.K, FGID, genes)
  pred.mat.K <- matrix(pred.v.K,   ncol=length(genes)-length(FGID))
  colnames(pred.mat.K) <- if(length(FGID)!=0) {genes[-FGID]} else {genes}
  rownames(pred.mat.K) <- rownames(log.rpm.K[K.tst, ])
  pred.df.K <- as.data.frame(pred.mat.K)
  
  obs.K <- log.rpm.K[K.tst, if(length(FGID)!=0) {genes[-FGID]} else {genes}]
  mae.K <- abs(pred.df.K - obs.K)
  mae.K <- apply(mae.K, 2, mean)
  
  if (length(FGID)==0) {} else {
    for (j in 1:length(FGID)) {
      mae.K   <- append(mae.K, NA_real_, after=FGID[j]-1)
    }
  }
  names(mae.K) <- genes

  c(i, smpl.size[i], length(FGID), mean(mae.K, na.rm=T), mae.K)
}
proc.time() - t
stopCluster(cl)

# 
colnames(result) <- c("model.ID", "smpl.size", "no.FG", "MAE", genes)
pred.neg.mix.K <- result
pred.neg.mix.K <- as.data.frame(pred.neg.mix.K)
save(pred.neg.mix.K, file="pred.neg.mix.K")



#==========#
# positive #
#==========#

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=6)  # Choose the number of threads you'd like to use

t <- proc.time()
result <- foreach(i = 1:100, .combine="rbind", .inorder=FALSE, .packages = "dplyr") %dopar% {
  
  load(paste("model.pos.mix.K.", i, sep=""))
  
  pred.K <- FIT::predict(the.model, tst.attribute.K, train.weather.K)
  # genes with prediction failure
  FGID.K   <- failed.gene.extractor(pred.K, length(K.tst), genes)
  FGID <- FGID.K
  
  pred.v.K <- unlist.subset(pred.K, FGID, genes)
  pred.mat.K <- matrix(pred.v.K,   ncol=length(genes)-length(FGID))
  colnames(pred.mat.K) <- if(length(FGID)!=0) {genes[-FGID]} else {genes}
  rownames(pred.mat.K) <- rownames(log.rpm.K[K.tst, ])
  pred.df.K <- as.data.frame(pred.mat.K)
  
  obs.K <- log.rpm.K[K.tst, if(length(FGID)!=0) {genes[-FGID]} else {genes}]
  mae.K <- abs(pred.df.K - obs.K)
  mae.K <- apply(mae.K, 2, mean)
  
  if (length(FGID)==0) {} else {
    for (j in 1:length(FGID)) {
      mae.K   <- append(mae.K, NA_real_, after=FGID[j]-1)
    }
  }
  names(mae.K) <- genes
  
  c(i, smpl.size[i], length(FGID), mean(mae.K, na.rm=T), mae.K)
}
proc.time() - t
stopCluster(cl)

# 
colnames(result) <- c("model.ID", "smpl.size", "no.FG", "MAE", genes)
pred.pos.mix.K <- result
pred.pos.mix.K <- as.data.frame(pred.pos.mix.K)
save(pred.pos.mix.K, file="pred.pos.mix.K")


#=====================#
# negative + positive #
#=====================#

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=6)  # Choose the number of threads you'd like to use

t <- proc.time()
result <- foreach(i = 1:100, .combine="rbind", .inorder=FALSE, .packages = "dplyr") %dopar% {
  
  load(paste("model.mix.mix.K.", i, sep=""))
  
  pred.K <- FIT::predict(the.model, tst.attribute.K, train.weather.K)
  # genes with prediction failure
  FGID.K   <- failed.gene.extractor(pred.K, length(K.tst), genes)
  FGID <- FGID.K
  
  pred.v.K <- unlist.subset(pred.K, FGID, genes)
  pred.mat.K <- matrix(pred.v.K,   ncol=length(genes)-length(FGID))
  colnames(pred.mat.K) <- if(length(FGID)!=0) {genes[-FGID]} else {genes}
  rownames(pred.mat.K) <- rownames(log.rpm.K[K.tst, ])
  pred.df.K <- as.data.frame(pred.mat.K)
  
  obs.K <- log.rpm.K[K.tst, if(length(FGID)!=0) {genes[-FGID]} else {genes}]
  mae.K <- abs(pred.df.K - obs.K)
  mae.K <- apply(mae.K, 2, mean)
  
  if (length(FGID)==0) {} else {
    for (j in 1:length(FGID)) {
      mae.K   <- append(mae.K, NA_real_, after=FGID[j]-1)
    }
  }
  names(mae.K) <- genes
  
  c(i, smpl.size[i], length(FGID), mean(mae.K, na.rm=T), mae.K)
}
proc.time() - t
stopCluster(cl)

# 
colnames(result) <- c("model.ID", "smpl.size", "no.FG", "MAE", genes)
pred.mix.mix.K <- result
pred.mix.mix.K <- as.data.frame(pred.mix.mix.K)
save(pred.mix.mix.K, file="pred.mix.mix.K")

