

### load files ------

load("OsTM0_1_rpm") ; rpm.1 <- rpm
load("OsTM0_2_rpm") ; rpm.2 <- rpm
load("OsTM0_3_rpm") ; rpm.3 <- rpm
load("OsTM0_4_rpm") ; rpm.4 <- rpm
load("OsTM0_5_rpm") ; rpm.5 <- rpm
load("OsTM0_6_rpm") ; rpm.6 <- rpm
load("OsTM0_7_rpm") ; rpm.7 <- rpm
load("OsTM0_8_rpm") ; rpm.8 <- rpm

load("180604_des") # the object is "des"


### combine data ------

rpm.1 <- rpm.1[,-1] #  1st column is "blank" ; T577 - T736
rpm.2 <- rpm.2[,-1] #  1st column is "blank" ; T737 - T896
rpm.3 <- rpm.3[,-1] #  1st column is "blank" ; T897 - T1056
rpm.4 <- rpm.4[,-1] #  1st column is "blank" ; T1057 - T1168/ T1 - T48
rpm.5 <- rpm.5[,-1] #  1st column is "blank" ; T49 - T208
rpm.6 <- rpm.6[,-1] #  1st column is "blank" ; T209 - T368
rpm.7 <- rpm.7[,-1] #  1st column is "blank" ; T369 - T528
rpm.8 <- rpm.8[,-1] #  1st column is "blank" ; T529 - T576 + "TK" + "TT" + "T861"
### remove TK and TT
rpm.8 <- rpm.8[,-c(49,50)]

rpm <- cbind(rpm.1, rpm.2, rpm.3, rpm.4, rpm.5, rpm.6, rpm.7, rpm.8)


### Process data ------

rpm.d <- rpm[des[,"NormalizationGroup"]=="data",]

# remove non-detected genes
rpm.d.ex <- rpm.d[apply(rpm.d, 1, sum)!=0,]

## log-transform
l.rpm <- log2(rpm.d.ex + 1)

# remove genes with low expression
log2(8 + 1)*ncol(l.rpm) # 3699.302
sum.log.rpm <- apply(l.rpm, 1, sum)
l.rpm.an <- l.rpm[sum.log.rpm >= 3699.302,]

dim(l.rpm.an) # 9804 genes, 1168 samples


### associate RNAseq with leaf.temp etc.

comb.d <- t(l.rpm.an)

leaf <- read.csv("combined.mod.csv",header=T,sep=",")
annotation <- read.table("sample.annotation.txt", header=T)

rownames(comb.d) <- substr(rownames(comb.d),2,nchar(rownames(comb.d)))
rownames(comb.d) <- as.integer(rownames(comb.d))

comb.d <- cbind(comb.d, annotation[row.names(comb.d), c("EnvID","PotID","Time","Cultivar")])
comb.d$Cultivar <- as.character(comb.d$Cultivar)
leaf$Cultivar   <- as.character(leaf$Cultivar)

leaf.sub <- NULL
for (i in 1:nrow(comb.d)) {
  temp <- subset(leaf, Env_ID==comb.d$EnvID[i]&Time==comb.d$Time[i]&Cultivar==comb.d$Cultivar[i])[,c("average","Tmp","LD")]
  leaf.sub <- rbind(leaf.sub, temp)
}

comb.d <- cbind(comb.d,leaf.sub) # this is the data file for analysis


### ------------------------------------------------------------
### With data excluding Env_ID 4-8,14-18 (seemingly exchanged plates), and 64,30,34,19, and 26 (randomly chosen),
### perform LASSO to predict the temperature of the smart GC,
### then predict the temperature of the smart GC for the samples not used to make the model

library("glmnet")
library("ggplot2")
library("gridExtra")

comb.d$Cultivar <- ifelse(comb.d$Cultivar=="Kos",1,0) # Kos = 1, Tak = 0
comb.d$LD       <- ifelse(comb.d$LD=="L",1,0)         # L = 1, D = 0
comb.d <- as.matrix(comb.d)

# make levels vector for EnvID
EnvID.mat <- matrix(c(37:73, 1:36, rep(NA_integer_, 2)), byrow = T, ncol=5)
for (i in 1:15) {
  assign(paste("levels.pl.", i, sep=""), as.factor(as.character(EnvID.mat[i,])))
}


for(j in c(1,6,11)) {
  plate.v <- j:(j + 4)
  for (i in plate.v[1]:plate.v[length(plate.v)]) {

    # split the data into target (test) data and training data
    tar.data.row <- annotation[annotation[,"PlateID"]==i, "SampleID"]

    tar.data  <- comb.d[rownames(comb.d) %in% as.character(tar.data.row), ]
    trn.data  <- comb.d[!(rownames(comb.d) %in% as.character(tar.data.row)), ]

    ### perform LASSO for variable selection
    lasso.cv <- cv.glmnet(x=trn.data[,1:(ncol(trn.data)-7)], y=trn.data[,"Tmp"], family="gaussian",alpha=1)

    ####################################
    ## fit to the training data
    lam.min <- lasso.cv$lambda.min
    lasso.min <- glmnet(x=trn.data[,1:(ncol(trn.data)-7)], y=trn.data[,"Tmp"], family="gaussian", lambda = lam.min, alpha=1)
    est.min <- predict(lasso.min, newx= trn.data[,1:(ncol(trn.data)-7)], s=lam.min, type="response")
    trn.min.AE <- abs(trn.data[,"Tmp"]-est.min)
    trn.min.MAE <- mean(trn.min.AE)

    ###################################
    ## fit to the target data & draw figures
    est.min.tar <- predict(lasso.min, newx= tar.data[,1:(ncol(tar.data)-7)], s=lam.min, type="response")
    tar.min.AE  <- abs(tar.data[,"Tmp"]-est.min.tar)
    tar.min.MAE <- mean(tar.min.AE)
    tar.min.pl  <- cbind(tar.data[,(ncol(tar.data)-6):ncol(tar.data)], est.min.tar)

    hist.min.df <-  data.frame(AE <- c(trn.min.AE, tar.min.AE),
                               data <- c(rep("trn",length(trn.min.AE)), rep("tar",length(tar.min.AE))))
    colnames(hist.min.df) <- c("AE","data")
    hist.min.df$data <- factor(hist.min.df$data, levels = c("trn","tar"))

    hist.min <- ggplot(hist.min.df, aes(x=AE,fill=data)) +
      geom_histogram(aes(y=..density..),position="identity", alpha=0.3, binwidth = 2) +
      scale_fill_manual(values=c("grey","red")) +
      geom_vline(xintercept=trn.min.MAE,colour="darkgrey") +
      geom_vline(xintercept=tar.min.MAE,colour="red") +
      theme_classic() +
      xlab("Absolute error") +
      xlim(c(-1,21)) +
      theme(legend.position = "none")

    tar.min.pl <- as.data.frame(tar.min.pl)
    colnames(tar.min.pl)[8] <- "est"
    tar.min.pl$Cultivar <- ifelse(tar.min.pl$Cultivar==1,"Kos","Tak") # Kos = 1, Tak = 0
    tar.min.pl$LD <- ifelse(tar.min.pl$LD==1,"L","D") # Kos = 1, Tak = 0
    tar.min.pl$EnvID <- factor(tar.min.pl$EnvID, levels = get(paste("levels.pl.", i, sep="")))

    plot.min <- ggplot(tar.min.pl) +
      geom_point(aes(x=Time - 24*floor(Time/20), y=est, color=Cultivar), size=3) +
      geom_point(aes(x=Time - 24*floor(Time/20), y=max(est)+2, color=LD), shape=8) +
      scale_color_manual(values = c("darkgray","blue","orange","red")) +
      geom_line(aes(x=Time - 24*floor(Time/20), y=Tmp)) +
      theme(legend.position = "none") +
      xlab("Hours") +
      ylab("Estimated temp.") +
      facet_wrap(~EnvID, ncol=5)
  
    assign(paste("hist.min.", i, sep=""), hist.min)
    assign(paste("plot.min.", i, sep=""), plot.min)
  }
}

pdf(file="1st.rnd.pl.01-05.pdf", width=8.27, height=11.69)
grid.arrange(hist.min.1, plot.min.1, 
             hist.min.2, plot.min.2,
             hist.min.3, plot.min.3, 
             hist.min.4, plot.min.4,
             hist.min.5, plot.min.5, 
             widths=c(1,3), ncol=2)
dev.off()

pdf(file="1st.rnd.pl.06-10.pdf", width=8.27, height=11.69)
grid.arrange(hist.min.6, plot.min.6, 
             hist.min.7, plot.min.7,
             hist.min.8, plot.min.8, 
             hist.min.9, plot.min.9,
             hist.min.10, plot.min.10, 
             widths=c(1,3), ncol=2)
dev.off()

pdf(file="1st.rnd.pl.11-15.pdf", width=8.27, height=11.69)
grid.arrange(hist.min.11, plot.min.11, 
             hist.min.12, plot.min.12,
             hist.min.13, plot.min.13, 
             hist.min.14, plot.min.14,
             hist.min.15, plot.min.15, 
             widths=c(1,3), ncol=2)
dev.off()


###------------------------------------------------
### exchange pl. 9 <-> pl.11, and pl.10 <-> pl.12

RNA.pl.9  <- comb.d[comb.d[,"EnvID"]%in%c(4:8),   1:(ncol(comb.d)-7)]
att.pl.9  <- comb.d[comb.d[,"EnvID"]%in%c(4:8),   (ncol(comb.d)-6):ncol(comb.d)]
RNA.pl.11 <- comb.d[comb.d[,"EnvID"]%in%c(14:18), 1:(ncol(comb.d)-7)]
att.pl.11 <- comb.d[comb.d[,"EnvID"]%in%c(14:18), (ncol(comb.d)-6):ncol(comb.d)]
RNA.pl.10 <- comb.d[comb.d[,"EnvID"]%in%c(9:13),  1:(ncol(comb.d)-7)]
att.pl.10 <- comb.d[comb.d[,"EnvID"]%in%c(9:13),  (ncol(comb.d)-6):ncol(comb.d)]
RNA.pl.12 <- comb.d[comb.d[,"EnvID"]%in%c(19:23), 1:(ncol(comb.d)-7)]
att.pl.12 <- comb.d[comb.d[,"EnvID"]%in%c(19:23), (ncol(comb.d)-6):ncol(comb.d)]

others <- comb.d[!(comb.d[,"EnvID"]%in%c(4:8, 14:18, 9:13, 19:23)), ]

cr.pl.11 <- cbind(RNA.pl.9,  att.pl.11)
cr.pl.9  <- cbind(RNA.pl.11, att.pl.9)
cr.pl.12 <- cbind(RNA.pl.10, att.pl.12)
cr.pl.10 <- cbind(RNA.pl.12, att.pl.10)

rownames(cr.pl.11) <- rownames(att.pl.11)
rownames(cr.pl.9)  <- rownames(att.pl.9)
rownames(cr.pl.12)  <- rownames(att.pl.12)
rownames(cr.pl.10)  <- rownames(att.pl.10)

comb.d.cr <- rbind(cr.pl.9, cr.pl.10, cr.pl.11, cr.pl.12, others)


for(j in c(1,6,11)) {
  plate.v <- j:(j + 4)
  for (i in plate.v[1]:plate.v[length(plate.v)]) {
  
    # split the data into target (test) data and training data
    tar.data.row <- annotation[annotation[,"PlateID"]==i, "SampleID"]
  
    tar.data  <- comb.d.cr[rownames(comb.d.cr) %in% as.character(tar.data.row), ]
    trn.data  <- comb.d.cr[!(rownames(comb.d.cr) %in% as.character(tar.data.row)), ]
  
    ### perform LASSO for variable selection
    lasso.cv <- cv.glmnet(x=trn.data[,1:(ncol(trn.data)-7)], y=trn.data[,"Tmp"], family="gaussian",alpha=1)
  
    ####################################
    ## fit to the training data
    lam.min <- lasso.cv$lambda.min
    lasso.min <- glmnet(x=trn.data[,1:(ncol(trn.data)-7)], y=trn.data[,"Tmp"], family="gaussian", lambda = lam.min, alpha=1)
    est.min <- predict(lasso.min, newx= trn.data[,1:(ncol(trn.data)-7)], s=lam.min, type="response")
    trn.min.AE <- abs(trn.data[,"Tmp"]-est.min)
    trn.min.MAE <- mean(trn.min.AE)
  
    ###################################
    ## fit to the target data & draw figures
    est.min.tar <- predict(lasso.min, newx= tar.data[,1:(ncol(tar.data)-7)], s=lam.min, type="response")
    tar.min.AE  <- abs(tar.data[,"Tmp"]-est.min.tar)
    tar.min.MAE <- mean(tar.min.AE)
    tar.min.pl  <- cbind(tar.data[,(ncol(tar.data)-6):ncol(tar.data)], est.min.tar)
   
    hist.min.df <-  data.frame(AE <- c(trn.min.AE, tar.min.AE),
                               data <- c(rep("trn",length(trn.min.AE)), rep("tar",length(tar.min.AE))))
    colnames(hist.min.df) <- c("AE","data")
    hist.min.df$data <- factor(hist.min.df$data, levels = c("trn","tar"))
  
    hist.min <- ggplot(hist.min.df, aes(x=AE,fill=data)) +
      geom_histogram(aes(y=..density..),position="identity", alpha=0.3, binwidth = 2) +
      scale_fill_manual(values=c("grey","red")) +
      geom_vline(xintercept=trn.min.MAE,colour="darkgrey") +
      geom_vline(xintercept=tar.min.MAE,colour="red") +
      theme_classic() +
      xlab("Absolute error") +
      xlim(c(-1,21)) +
      theme(legend.position = "none")
  
    tar.min.pl <- as.data.frame(tar.min.pl)
    colnames(tar.min.pl)[8] <- "est"
    tar.min.pl$Cultivar <- ifelse(tar.min.pl$Cultivar==1,"Kos","Tak") # Kos = 1, Tak = 0
    tar.min.pl$LD <- ifelse(tar.min.pl$LD==1,"L","D") # Kos = 1, Tak = 0
    tar.min.pl$EnvID <- factor(tar.min.pl$EnvID, levels = get(paste("levels.pl.", i, sep="")))
  
    plot.min <- ggplot(tar.min.pl) +
      geom_point(aes(x=Time - 24*floor(Time/20), y=est, color=Cultivar), size=3) +
      geom_point(aes(x=Time - 24*floor(Time/20), y=max(est)+2, color=LD), shape=8) +
      scale_color_manual(values = c("darkgray","blue","orange","red")) +
      geom_line(aes(x=Time - 24*floor(Time/20), y=Tmp)) +
      theme(legend.position = "none") +
      xlab("Hours") +
      ylab("Estimated temp.") +
      facet_wrap(~EnvID, ncol=5)
  
    assign(paste("hist.min.", i, sep=""), hist.min)
    assign(paste("plot.min.", i, sep=""), plot.min)
  }
}

pdf(file="2nd.rnd.pl.01-05.pdf", width=8.27, height=11.69)
grid.arrange(hist.min.1, plot.min.1, 
             hist.min.2, plot.min.2,
             hist.min.3, plot.min.3, 
             hist.min.4, plot.min.4,
             hist.min.5, plot.min.5, 
             widths=c(1,3), ncol=2)
dev.off()

pdf(file="2nd.rnd.pl.06-10.pdf", width=8.27, height=11.69)
grid.arrange(hist.min.6, plot.min.6, 
             hist.min.7, plot.min.7,
             hist.min.8, plot.min.8, 
             hist.min.9, plot.min.9,
             hist.min.10, plot.min.10, 
             widths=c(1,3), ncol=2)
dev.off()

pdf(file="2nd.rnd.pl.11-15.pdf", width=8.27, height=11.69)
grid.arrange(hist.min.11, plot.min.11, 
             hist.min.12, plot.min.12,
             hist.min.13, plot.min.13, 
             hist.min.14, plot.min.14,
             hist.min.15, plot.min.15, 
             widths=c(1,3), ncol=2)
dev.off()

