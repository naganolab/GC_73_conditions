
library("FIT")
library("foreach")
library("doParallel")
library("dplyr")

load("genes") # character vector of length 474
load("log.rpm.T")

load("weather.01")
load("attribute.all.T")

genes <- genes[genes %in% colnames(log.rpm.T)] #470 genes

setwd("Tak.result")
smpl.size <- rep(c(64, 128, 256, 512), each=100)


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

oldest.T <- min(as.POSIXct(paste(paste(attribute.all.T$year, attribute.all.T$month, attribute.all.T$day, sep="-"),
                paste(attribute.all.T$hour, attribute.all.T$min, "00", sep=":"), sep=" ")))
weather.T <- weather.01[weather.01$time >= oldest.T + as.difftime(-72, units="hours"),]
oldest.wt.T   <- min(weather.T$time)
time.T   <- seq(0, (nrow(weather.T)-1)*10, by=10)
weather.T <- cbind(time.T, weather.T[,c("temperature","radiation")])
colnames(weather.T)   <- c("time", "temperature", "radiation")
train.weather.T   <- FIT::convert.weather(weather.T,   c("temperature", "radiation"))

T.tst <- 1298:1836
attribute.tst.T <- attribute.all.T[T.tst, ]
attribute.tst.T$time <- as.numeric(as.POSIXct(paste(paste(attribute.tst.T$year, attribute.tst.T$month, attribute.tst.T$day, sep="-"),
                                                    paste(attribute.tst.T$hour, attribute.tst.T$min, "00", sep=":"), 
                                                    sep=" "))-oldest.wt.T, units="mins")
tst.attribute.T <- FIT::convert.attribute(attribute.tst.T)


#=======#
# Field #
#=======#

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=6)  # Choose the number of threads you'd like to use

t <- proc.time()
result <- foreach(i = 1:400, .combine="rbind", .inorder=FALSE, .packages = "dplyr") %dopar% {
  
  load(paste("model.T.fld.", i, sep=""))

  pred.T <- FIT::predict(the.model, tst.attribute.T, train.weather.T)
  # genes with prediction failure
  FGID.T   <- failed.gene.extractor(pred.T, length(T.tst), genes)
  FGID <- FGID.T

  
  pred.v.T <- unlist.subset(pred.T, FGID, genes)
  pred.mat.T <- matrix(pred.v.T,   ncol=length(genes)-length(FGID))
  colnames(pred.mat.T) <- if(length(FGID)!=0) {genes[-FGID]} else {genes}
  rownames(pred.mat.T) <- rownames(log.rpm.T[T.tst, ])
  pred.df.T <- as.data.frame(pred.mat.T)
  
  obs.T <- log.rpm.T[T.tst, if(length(FGID)!=0) {genes[-FGID]} else {genes}]
  mae.T <- abs(pred.df.T - obs.T)
  mae.T <- apply(mae.T, 2, mean)
  
  if (length(FGID)==0) {} else {
    for (j in 1:length(FGID)) {
      mae.T   <- append(mae.T, NA_real_, after=FGID[j]-1)
    }
  }
  names(mae.T) <- genes

  c(i, smpl.size[i], length(FGID), mean(mae.T, na.rm=T), mae.T)
}
proc.time() - t
stopCluster(cl)

# 
colnames(result) <- c("model.ID", "smpl.size", "no.FG", "MAE", genes)
pred.T.fld <- result
pred.T.fld <- as.data.frame(pred.T.fld)
save(pred.T.fld, file="pred.T.fld")




#======#
#  GC  #
#======#

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=5)  # Choose the number of threads you'd like to use (20 for sumii?)

t <- proc.time()
result <- foreach(i = 1:400, .combine="rbind", .inorder=FALSE, .packages = "dplyr") %dopar% {
  
  load(paste("model.T.GC.", i, sep=""))
  
  pred.T <- FIT::predict(the.model, tst.attribute.T, train.weather.T)
  # genes with prediction failure
  FGID.T   <- failed.gene.extractor(pred.T, length(T.tst), genes)
  FGID <- FGID.T
  
  
  pred.v.T <- unlist.subset(pred.T, FGID, genes)
  pred.mat.T <- matrix(pred.v.T,   ncol=length(genes)-length(FGID))
  colnames(pred.mat.T) <- if(length(FGID)!=0) {genes[-FGID]} else {genes}
  rownames(pred.mat.T) <- rownames(log.rpm.T[T.tst, ])
  pred.df.T <- as.data.frame(pred.mat.T)
  
  obs.T <- log.rpm.T[T.tst, if(length(FGID)!=0) {genes[-FGID]} else {genes}]
  mae.T <- abs(pred.df.T - obs.T)
  mae.T <- apply(mae.T, 2, mean)
  
  if (length(FGID)==0) {} else {
    for (j in 1:length(FGID)) {
      mae.T   <- append(mae.T, NA_real_, after=FGID[j]-1)
    }
  }
  names(mae.T) <- genes
  
  c(i, smpl.size[i], length(FGID), mean(mae.T, na.rm=T), mae.T)
}
proc.time() - t
stopCluster(cl)

# 
colnames(result) <- c("model.ID", "smpl.size", "no.FG", "MAE", genes)
pred.T.GC <- result
pred.T.GC <- as.data.frame(pred.T.GC)
save(pred.T.GC, file="pred.T.GC")





#=======#
#  mix  #
#=======#

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=5)  # Choose the number of threads you'd like to use (20 for sumii?)

t <- proc.time()
result <- foreach(i = 1:400, .combine="rbind", .inorder=FALSE, .packages = "dplyr") %dopar% {
  
  load(paste("model.T.mix.", i, sep=""))
  
  pred.T <- FIT::predict(the.model, tst.attribute.T, train.weather.T)
  # genes with prediction failure
  FGID.T   <- failed.gene.extractor(pred.T, length(T.tst), genes)
  FGID <- FGID.T
  
  
  pred.v.T <- unlist.subset(pred.T, FGID, genes)
  pred.mat.T <- matrix(pred.v.T,   ncol=length(genes)-length(FGID))
  colnames(pred.mat.T) <- if(length(FGID)!=0) {genes[-FGID]} else {genes}
  rownames(pred.mat.T) <- rownames(log.rpm.T[T.tst, ])
  pred.df.T <- as.data.frame(pred.mat.T)
  
  obs.T <- log.rpm.T[T.tst, if(length(FGID)!=0) {genes[-FGID]} else {genes}]
  mae.T <- abs(pred.df.T - obs.T)
  mae.T <- apply(mae.T, 2, mean)
  
  if (length(FGID)==0) {} else {
    for (j in 1:length(FGID)) {
      mae.T   <- append(mae.T, NA_real_, after=FGID[j]-1)
    }
  }
  names(mae.T) <- genes
  
  c(i, smpl.size[i], length(FGID), mean(mae.T, na.rm=T), mae.T)
}
proc.time() - t
stopCluster(cl)

# 
colnames(result) <- c("model.ID", "smpl.size", "no.FG", "MAE", genes)
pred.T.mix <- result
pred.T.mix <- as.data.frame(pred.T.mix)
save(pred.T.mix, file="pred.T.mix")

