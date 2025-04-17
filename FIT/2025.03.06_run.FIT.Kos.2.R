# Train FIT models with Koshihikari data

Num_cores <- 10

library("FIT")
library("foreach")
library("doParallel")

load("genes") # character vector of length 474

load("weather.01")
load("attribute.all.K")
load("weight.all.K")
load("log.rpm.K")

setwd("Kos.result.mix.2")

#================#
# shared objects #
#================#

grid.coords <- list(  # parameters for the grid search
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

genes <- genes[genes %in% colnames(log.rpm.K)] #470 genes


#=======================================#
# Run the bootstrap for Koshihikari mix #
#=======================================#

div.mat <- NULL # The matrix that defines the diving points of Env.Setting
row.counter <- 1
for(i in 2:584) {   # 584: rows of attribute.all.K
  time.prev <- attribute.all.K$time[i - 1]
  time.now  <- attribute.all.K$time[i]
  if(time.now - time.prev == 180) {
    row.counter <- row.counter + 1
  } else {
    div.mat <- rbind(div.mat, c((i - 1) - (row.counter - 1), row.counter))
    row.counter <- 1
  }
  if(i==584) {div.mat <- rbind(div.mat, c((i - 1) - (row.counter - 1), row.counter))}
}


pos.DIF <- NULL # logical vector
neg.DIF <- NULL # logical vector

the.ES <- NULL
for(i in 1:nrow(div.mat)) {
  the.ES <- attribute.all.K[(div.mat[i, 1]):(div.mat[i, 1] + div.mat[i, 2] - 1), ]
  the.dates <- as.POSIXct(paste(paste(the.ES$year, the.ES$month, the.ES$day, sep="-"),
                                paste(the.ES$hour, the.ES$min, "00", sep=":"), sep=" "))
  the.weather <- weather.01[weather.01$time %in% the.dates, ]
  
  # judge the positive/negative DIF
  if(var(the.weather$radiation)==0) {
    pos.DIF <- c(pos.DIF, rep(FALSE, nrow(the.weather)))
    neg.DIF <- c(neg.DIF, rep(FALSE, nrow(the.weather)))
  } else {
    dark.row  <- which.min(the.weather$radiation)
    light.row <- which.max(the.weather$radiation)
    if(the.weather$temperature[dark.row] < the.weather$temperature[light.row]) {
      pos.DIF <- c(pos.DIF, rep(TRUE,  nrow(the.weather)))
      neg.DIF <- c(neg.DIF, rep(FALSE, nrow(the.weather)))
    } else {
      pos.DIF <- c(pos.DIF, rep(FALSE,  nrow(the.weather)))
      neg.DIF <- c(neg.DIF, rep(TRUE, nrow(the.weather)))
    }
  }
}



iteration <- 100
sample.sizes <- 144*2 # sum(neg.DIF)
sample.size.v <- rep(sample.sizes, each=iteration)

cl <- makeCluster(detectCores())
registerDoParallel(cl, cores=Num_cores)  # Choose the number of threads you'd like to use (20 for sumii?)

t <- proc.time()
result <- foreach(i = 1:length(sample.size.v),.combine="rbind",.inorder=FALSE,.packages = "dplyr") %dopar% {
  
  set.seed(i %% 100)  # Choose any number. 1729 is known as the taxi number, btw.
    
  #====================#
  # make training data #
  #====================#

  ###################
  ### Koshihikari ###

  # Sample type/    n/row/
  # 73-conditions 584 rows 1    - 584
  # 2015          462 rows 585  - 1046
  # 2016          343 rows 1047 - 1389
  # 2017          615 rows 1390 - 2004
  
  trn.data.size <- sample.size.v[i]  # !!! Choose even number!
  negative.GC.seq  <- sample(which(neg.DIF), trn.data.size/2, replace=F)
  positive.GC.seq  <- sample(which(pos.DIF), trn.data.size/2, replace=F)
  fld.seq <- sample(585:1389, trn.data.size/2, replace=F) # !!! cultivar-specific; do not use 2017 data
  negative.mix.seq <- c(negative.GC.seq, fld.seq)
  positive.mix.seq <- c(positive.GC.seq, fld.seq)
  mix.mix.seq <- c(negative.GC.seq[1:(trn.data.size/4)],
                   positive.GC.seq[1:(trn.data.size/4)], fld.seq)
  
  #==========#
  # negative #
  #==========#
  
  # use of the training data:
  trn.row <- sort(negative.mix.seq)
  
  attrib <- attribute.all.K[trn.row, ]
  oldest.smpl <- min(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                            paste(attrib$hour, attrib$min, "00", sep=":"), sep=" ")))
  wt.sub <- weather.01[weather.01$time >= oldest.smpl + as.difftime(-72, units="hours"),]
  oldest.wt <- min(wt.sub$time)
  time.label <- seq(0, (nrow(wt.sub)-1)*10, by=10)
  wt.sub <- cbind(time.label, wt.sub[, c("temperature","radiation")])
  colnames(wt.sub) <- c("time", "temperature", "radiation")
  train.wt <- FIT::convert.weather(wt.sub, c("temperature", "radiation"))
  
  # Parameter estimation
  attrib$time <- as.numeric(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                             paste(attrib$hour, attrib$min, "00", sep=":"), 
                                             sep=" "))-oldest.wt, units="mins")
  train.attrib  <- FIT::convert.attribute(attrib)

  # choose genes you'd like to predict
  train.expression <- FIT::convert.expression(as.matrix(log.rpm.K[trn.row, ]), genes)
  train.weight     <- FIT::convert.weight(weight.all.K[trn.row, ], genes)
  the.model <- FIT::train(train.expression,
                          train.attrib,
                          train.wt,
                          recipe,
                          train.weight)
  the.model <- unlist(the.model)

  # save the models and the training data
  save(the.model, file=paste("model.neg.mix.K", i, sep="."))
  
  
  #==========#
  # positive #
  #==========#
  
  # use of the training data:
  trn.row <- sort(positive.mix.seq)
  
  attrib <- attribute.all.K[trn.row, ]
  oldest.smpl <- min(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                      paste(attrib$hour, attrib$min, "00", sep=":"), sep=" ")))
  wt.sub <- weather.01[weather.01$time >= oldest.smpl + as.difftime(-72, units="hours"),]
  oldest.wt <- min(wt.sub$time)
  time.label <- seq(0, (nrow(wt.sub)-1)*10, by=10)
  wt.sub <- cbind(time.label, wt.sub[, c("temperature","radiation")])
  colnames(wt.sub) <- c("time", "temperature", "radiation")
  train.wt <- FIT::convert.weather(wt.sub, c("temperature", "radiation"))
  
  # Parameter estimation
  attrib$time <- as.numeric(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                             paste(attrib$hour, attrib$min, "00", sep=":"), 
                                             sep=" "))-oldest.wt, units="mins")
  train.attrib  <- FIT::convert.attribute(attrib)
  
  # choose genes you'd like to predict
  train.expression <- FIT::convert.expression(as.matrix(log.rpm.K[trn.row, ]), genes)
  train.weight     <- FIT::convert.weight(weight.all.K[trn.row, ], genes)
  the.model <- FIT::train(train.expression,
                          train.attrib,
                          train.wt,
                          recipe,
                          train.weight)
  the.model <- unlist(the.model)
  
  # save the models and the training data
  save(the.model, file=paste("model.pos.mix.K", i, sep="."))
  
  #=====================#
  # positive + negative #
  #=====================#
  
  # use of the training data:
  trn.row <- sort(mix.mix.seq)
  
  attrib <- attribute.all.K[trn.row, ]
  oldest.smpl <- min(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                      paste(attrib$hour, attrib$min, "00", sep=":"), sep=" ")))
  wt.sub <- weather.01[weather.01$time >= oldest.smpl + as.difftime(-72, units="hours"),]
  oldest.wt <- min(wt.sub$time)
  time.label <- seq(0, (nrow(wt.sub)-1)*10, by=10)
  wt.sub <- cbind(time.label, wt.sub[, c("temperature","radiation")])
  colnames(wt.sub) <- c("time", "temperature", "radiation")
  train.wt <- FIT::convert.weather(wt.sub, c("temperature", "radiation"))
  
  # Parameter estimation
  attrib$time <- as.numeric(as.POSIXct(paste(paste(attrib$year, attrib$month, attrib$day, sep="-"),
                                             paste(attrib$hour, attrib$min, "00", sep=":"), 
                                             sep=" "))-oldest.wt, units="mins")
  train.attrib  <- FIT::convert.attribute(attrib)
  
  # choose genes you'd like to predict
  train.expression <- FIT::convert.expression(as.matrix(log.rpm.K[trn.row, ]), genes)
  train.weight     <- FIT::convert.weight(weight.all.K[trn.row, ], genes)
  the.model <- FIT::train(train.expression,
                          train.attrib,
                          train.wt,
                          recipe,
                          train.weight)
  the.model <- unlist(the.model)
  
  # save the models and the training data
  save(the.model, file=paste("model.mix.mix.K", i, sep="."))


  ## record training data
  trn.smpl <- c(negative.mix.seq, positive.mix.seq)
  c(i, trn.data.size, trn.smpl)
}
proc.time() - t
stopCluster(cl)

colnames(result) <- c("run.ID", "smpl.size", 
                      paste("neg.trn.smpl.", 1:(144*2), sep=""),
                      paste("pos.trn.smpl.", 1:(144*2), sep=""))
trn.smpl.K.mix <- as.data.frame(result)
save(trn.smpl.K.mix, file="trn.smpl.K.mix")

# mix.mix model uses the first halves of positive and negative data sets


