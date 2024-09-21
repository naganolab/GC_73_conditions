# Train FIT models with Koshihikari data

num_cores <- 6 # Choose the number of threads you'd like to use

library("FIT")
library("foreach")
library("doParallel")

load("genes") # character vector of length 474

load("weather.01")
load("attribute.all.T")
load("weight.all.T")
load("log.rpm.T")

setwd("Tak.result")

#================#
# shared objects #
#================#

grid.coords <- list(  # grid search parameter
  env.temperature.threshold = c(10, 15, 20, 25, 30),
  env.temperature.amplitude = c(-100/30, -1/30, 1/30, 100/30),
  env.temperature.period = c(10, 30, 90, 270, 720, 1440, 1440*3),
  gate.temperature.phase = seq(0, 23*60, 1*60),
  gate.temperature.threshold = -1,
  gate.temperature.amplitude = c(-5, 5),
  
  env.radiation.threshold = c(1, 10, 20, 30, 40),
  env.radiation.amplitude = c(-100/80, -1/80, 1/80, 100/80),
  env.radiation.period = c(10, 30, 90, 270, 720, 1440, 1440*3),
  gate.radiation.phase = seq(0, 23*60, 1*60),
  gate.radiation.threshold = -1,
  gate.radiation.amplitude = c(-5, 5)
)

recipe <- FIT::make.recipe(c('temperature','radiation'), 
                           init = 'gridsearch',
                           optim = c('lm'),
                           fit = 'fit.lasso',
                           init.data = grid.coords,
                           time.step = 10, 
                           gate.open.min = 479)

genes <- genes[genes %in% colnames(log.rpm.T)] #470 genes

#=======================================#
# Run the bootstrap for Koshihikari out # i.e., field
#=======================================#

iteration <- 100
sample.sizes <- c(64,128,256,512)
sample.size.v <- rep(sample.sizes, each=iteration)
set.seed(1729)  # Choose any number. 1729 is known as the taxi number, btw.

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=num_cores)

t <- proc.time()
result <- foreach(i = 1:length(sample.size.v),.combine="rbind",.inorder=FALSE,.packages = "dplyr") %dopar% {
  
  #====================#
  # make training data #
  #====================#

  ################
  ### Takanari ###
  
  # sample type  /  n/rows
  # 73-conditions 583 rows 1    - 583
  # 2015          418 rows 584  - 1001
  # 2016          296 rows 1002 - 1297
  # 2017          539 rows 1298 - 1836
  
  trn.data.size <- sample.size.v[i]  # !!! Choose even number!
  GC.sample.seq  <- sample(1:583, trn.data.size, replace=F)    # !!! cultivar-specific; do not use 2017 data
  fld.sample.seq <- sample(584:1297, trn.data.size, replace=F) # !!! cultivar-specific; do not use 2017 data
  mix.sample.seq <- c(GC.sample.seq[1:(trn.data.size/2)], fld.sample.seq[1:(trn.data.size/2)])
 
  # use of the the 3 potions:
  #trn.row <- sort(GC.sample.seq)
  trn.row <- sort(fld.sample.seq)
  #trn.row <- sort(mix.sample.seq)

  attrib <- attribute.all.T[trn.row, ] # !!! cultivar-specific
  oldest.smpl <- min(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                            paste(attrib$hour, attrib$min, "00", sep=":"), sep=" ")))
  wt.sub <- weather.01[weather.01$time >= oldest.smpl + as.difftime(-72, units="hours"),]
  oldest.wt <- min(wt.sub$time)
  time.label <- seq(0, (nrow(wt.sub)-1)*10, by=10)
  wt.sub <- cbind(time.label, wt.sub[, c("temperature","radiation")])
  colnames(wt.sub) <- c("time", "temperature", "radiation")
  train.wt <- FIT::convert.weather(wt.sub, c("temperature", "radiation"))
  
  #======================#
  # Parameter estimation #
  #======================#
  
  attrib$time <- as.numeric(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                             paste(attrib$hour, attrib$min, "00", sep=":"), 
                                             sep=" "))-oldest.wt, units="mins")
  train.attrib  <- FIT::convert.attribute(attrib)

  
  # choose genes you'd like to predict
  train.expression <- FIT::convert.expression(as.matrix(log.rpm.T[trn.row, ]), genes)
  train.weight     <- FIT::convert.weight(weight.all.T[trn.row, ], genes)
  the.model <- FIT::train(train.expression,
                          train.attrib,
                          train.wt,
                          recipe,
                          train.weight)
  the.model <- unlist(the.model)

  # save the models and the training data
  save(the.model, file=paste("model.T.fld", i, sep="."))
  
  trn.smpl <- rep(NA_real_, 512)
  trn.smpl[1:length(trn.row)] <- trn.row
  c(i, trn.data.size, trn.smpl)
}
proc.time() - t
stopCluster(cl)

colnames(result) <- c("run.ID", "smpl.size", paste("trn.smpl.", 1:512, sep=""))
trn.smpl.T.fld <- as.data.frame(result)
save(trn.smpl.T.fld, file="trn.smpl.T.fld")




#======================================#
# Run the bootstrap for Koshihikari in # i.e., GC
#======================================#

iteration <- 100
sample.sizes <- c(64,128,256,512)
sample.size.v <- rep(sample.sizes, each=iteration)
set.seed(1729)  # Choose any number. 1729 is known as the taxi number, btw.

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=num_cores)

t <- proc.time()
result <- foreach(i = 1:length(sample.size.v),.combine="rbind",.inorder=FALSE,.packages = "dplyr") %dopar% {
  
  #====================#
  # make training data #
  #====================#
  
  ################
  ### Takanari ###
  
  # data type    /  n/row
  # 73-conditions 583 rows 1    - 583
  # 2015          418 rows 584  - 1001
  # 2016          296 rows 1002 - 1297
  # 2017          539 rows 1298 - 1836
  
  trn.data.size <- sample.size.v[i]  # !!! Choose even number!
  GC.sample.seq  <- sample(1:583, trn.data.size, replace=F)    # !!! cultivar-specific; do not use 2017 data
  fld.sample.seq <- sample(584:1297, trn.data.size, replace=F) # !!! cultivar-specific; do not use 2017 data
  mix.sample.seq <- c(GC.sample.seq[1:(trn.data.size/2)], fld.sample.seq[1:(trn.data.size/2)])
  
  # use of the the 3 potions:
  trn.row <- sort(GC.sample.seq)
  #trn.row <- sort(fld.sample.seq)
  #trn.row <- sort(mix.sample.seq)
  
  attrib <- attribute.all.T[trn.row, ] # !!! cultivar-specific
  oldest.smpl <- min(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                      paste(attrib$hour, attrib$min, "00", sep=":"), sep=" ")))
  wt.sub <- weather.01[weather.01$time >= oldest.smpl + as.difftime(-72, units="hours"),]
  oldest.wt <- min(wt.sub$time)
  time.label <- seq(0, (nrow(wt.sub)-1)*10, by=10)
  wt.sub <- cbind(time.label, wt.sub[, c("temperature","radiation")])
  colnames(wt.sub) <- c("time", "temperature", "radiation")
  train.wt <- FIT::convert.weather(wt.sub, c("temperature", "radiation"))
  
  #======================#
  # Parameter estimation #
  #======================#
  
  attrib$time <- as.numeric(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                             paste(attrib$hour, attrib$min, "00", sep=":"), 
                                             sep=" "))-oldest.wt, units="mins")
  train.attrib  <- FIT::convert.attribute(attrib)
  
  
  # choose genes you'd like to predict
  train.expression <- FIT::convert.expression(as.matrix(log.rpm.T[trn.row, ]), genes)
  train.weight     <- FIT::convert.weight(weight.all.T[trn.row, ], genes)
  the.model <- FIT::train(train.expression,
                          train.attrib,
                          train.wt,
                          recipe,
                          train.weight)
  the.model <- unlist(the.model)
  
  # save the models and the training data
  save(the.model, file=paste("model.T.GC", i, sep="."))
  
  trn.smpl <- rep(NA_real_, 512)
  trn.smpl[1:length(trn.row)] <- trn.row
  c(i, trn.data.size, trn.smpl)
}
proc.time() - t
stopCluster(cl)

colnames(result) <- c("run.ID", "smpl.size", paste("trn.smpl.", 1:512, sep=""))
trn.smpl.T.GC <- as.data.frame(result)
save(trn.smpl.T.GC, file="trn.smpl.T.GC")




#=======================================#
# Run the bootstrap for Koshihikari mix #
#=======================================#

iteration <- 100
sample.sizes <- c(64, 128, 256, 512)
sample.size.v <- rep(sample.sizes, each=iteration)
set.seed(1729)  # Choose any number. 1729 is known as the taxi number, btw.

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=num_cores)  # Choose the number of threads you'd like to use (20 for sumii?)

t <- proc.time()
result <- foreach(i = 1:length(sample.size.v),.combine="rbind",.inorder=FALSE,.packages = "dplyr") %dopar% {
  
  #====================#
  # make training data #
  #====================#
  
  ################
  ### Takanari ###
  
  # sample type  /  n/row
  # 73-conditions 583 rows 1    - 583
  # 2015          418 rows 584  - 1001
  # 2016          296 rows 1002 - 1297
  # 2017          539 rows 1298 - 1836
  
  trn.data.size <- sample.size.v[i]  # !!! Choose even number!
  GC.sample.seq  <- sample(1:583, trn.data.size, replace=F)    # !!! cultivar-specific; do not use 2017 data
  fld.sample.seq <- sample(584:1297, trn.data.size, replace=F) # !!! cultivar-specific; do not use 2017 data
  mix.sample.seq <- c(GC.sample.seq[1:(trn.data.size/2)], fld.sample.seq[1:(trn.data.size/2)])
  
  # use of the the 3 potions:
  #trn.row <- sort(GC.sample.seq)
  #trn.row <- sort(fld.sample.seq)
  trn.row <- sort(mix.sample.seq)
  
  attrib <- attribute.all.T[trn.row, ] # !!! cultivar-specific
  oldest.smpl <- min(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                      paste(attrib$hour, attrib$min, "00", sep=":"), sep=" ")))
  wt.sub <- weather.01[weather.01$time >= oldest.smpl + as.difftime(-72, units="hours"),]
  oldest.wt <- min(wt.sub$time)
  time.label <- seq(0, (nrow(wt.sub)-1)*10, by=10)
  wt.sub <- cbind(time.label, wt.sub[, c("temperature","radiation")])
  colnames(wt.sub) <- c("time", "temperature", "radiation")
  train.wt <- FIT::convert.weather(wt.sub, c("temperature", "radiation"))
  
  #======================#
  # Parameter estimation #
  #======================#
  
  attrib$time <- as.numeric(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                             paste(attrib$hour, attrib$min, "00", sep=":"), 
                                             sep=" "))-oldest.wt, units="mins")
  train.attrib  <- FIT::convert.attribute(attrib)
  
  # choose genes you'd like to predict
  train.expression <- FIT::convert.expression(as.matrix(log.rpm.T[trn.row, ]), genes)
  train.weight     <- FIT::convert.weight(weight.all.T[trn.row, ], genes)
  the.model <- FIT::train(train.expression,
                          train.attrib,
                          train.wt,
                          recipe,
                          train.weight)
  the.model <- unlist(the.model)
  
  # save the models and the training data
  save(the.model, file=paste("model.T.mix", i, sep="."))
  
  trn.smpl <- rep(NA_real_, 512)
  trn.smpl[1:length(trn.row)] <- trn.row
  c(i, trn.data.size, trn.smpl)
}
proc.time() - t
stopCluster(cl)

colnames(result) <- c("run.ID", "smpl.size", paste("trn.smpl.", 1:512, sep=""))
trn.smpl.T.mix <- as.data.frame(result)
save(trn.smpl.T.mix, file="trn.smpl.T.mix")
