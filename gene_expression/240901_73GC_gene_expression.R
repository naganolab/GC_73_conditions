### Author: Yoichi Hashida
### library ----
library("viridis") #options: viridis, magma, plasma and inferno
library("tidyverse")
library("coin")
library("ggbeeswarm")
library("cowplot")
library("fields")
library("Rtsne")
library("RColorBrewer")

color5 <- c("#ff4b00","#4dc4ff","#03af7a","#fff100","#005aff") #yellow
source("../script/ng.Colors.R")

### GO analysis
# load the GO analysis functions
source("../script/GOanalysis/GOanalysis_functions.R")

# load the table of Gene ID and GO (as ulg)
load("../script/GOanalysis/ulg.RAP_190219")

### KEGG analysis 
# load the GO analysis functions
source("../script/KEGGanalysis/200107_KEGG_function.R")

# load kegg_rice
load("../script/KEGGanalysis/kegg_rice")
load("../script/KEGGanalysis/kegg_description")

#memory.size(50000)

dir.create("../output")

### load files ----
load("../input/20180511_Araport11_genes.201606.transcript.rep_ERCC_Virus7457.desc")   # object name: des
des.A.thal <- des
load("../input/180604_des")
load("../input/at.lab.additional")
load("../input/log.rpm")
rownames(log.rpm) <- as.integer(rownames(log.rpm))
load("../input/rawcnt.lab")

# leaf <- read.csv("../input/20190722_combined.mod.csv",header=T,sep=",")
# leaf[565,"Tmp"] <- 30 ; leaf[565,"LD"] <- "D" ; leaf[1149,"Tmp"] <- 30 ; leaf[1149,"LD"] <- "D"

annotation <- read.table("../input/sample.annotation.txt", header=T)
env.settings <- read.table("../input/env.settings.txt")

# log.rpm       1167 x 39772 (39765 genes + 7 attributes)
# env.settings  73 x 7

mean.exp <- apply(log.rpm[, 1:39765], 2, mean)
var.log.rpm <- apply(log.rpm[, 1:39765], 2, var)
comb.d   <- log.rpm                                 # the name is to use the old code

exp.ind.ratio.calculator <- function(x) {
  no.ind <- length(x)
  no.exp.ind <- sum(x!=0)
  no.exp.ind/no.ind
}
exp.ind.ratio <- apply(log.rpm[, 1:39765], 2, exp.ind.ratio.calculator)


plot(mean.exp, var.log.rpm, pch=16, col=rgb(exp.ind.ratio,0,1 - exp.ind.ratio,alpha=0.2))

sum((mean.exp<10)&(var.log.rpm<11)) # 39600 genes are within this range

local.dens.mat <- matrix(NA_real_, nrow=10, ncol=11)
for (i in 1:10) {
  for (j in 1:11) {
    local.dens.mat[i,j] <- sum((mean.exp>(i-1))&(mean.exp<i)&(var.log.rpm>(j-1))&(var.log.rpm<j))/39600
  }
}

local.sample.dens.mat <- exp(-local.dens.mat)/sum(exp(-local.dens.mat))

IDs.for.mv.plot <- c(1:39765)[(mean.exp>10)|(var.log.rpm>11)]
for (i in 1:10) {
  for (j in 1:11) {
    if(local.dens.mat[i,j]==0) {} else {
      id.selected <- sample(c(1:39765)[(mean.exp>(i-1))&(mean.exp<i)&(var.log.rpm>(j-1))&(var.log.rpm<j)],
                            min(local.sample.dens.mat[i,j]*3000, sum((mean.exp>(i-1))&(mean.exp<i)&(var.log.rpm>(j-1))&(var.log.rpm<j))), replace=F)
      IDs.for.mv.plot <- c(IDs.for.mv.plot, id.selected)
    }
  }
}


layout.mat <- matrix(c(1,  1,1,1,    1,1,1,    1,1,1,    1,2,2,    2,2,2,
                       3,  8,9,10,   11,12,13, 14,15,16, 17,17,17, 18,19,20,
                       4,  21,22,23, 24,25,26, 27,27,27, 28,29,30, 31,32,33,
                       5,  34,35,36, 37,37,37, 38,39,40, 41,42,43, 44,45,46,
                       6,  47,47,47, 48,49,50, 51,52,53, 54,55,56, 57,58,59,
                       7,  60,61,62, 63,64,65, 66,67,68, 69,70,71, 72,73,74,
                       0,  0,0,75,   0,0,76,   0,0,77,   0,0,78,   0,0,79),byrow=T,nrow=7,ncol=16)

### geneplot ----
geneplot <- function(glist){
  tmp <- 1:39765
  names(tmp) <- colnames(log.rpm)[1:39765]
  rnd.seed <- tmp[as.character(glist)]
  rnd.seed <- rnd.seed[!is.na(rnd.seed)]
  
  for (i in 1:length(rnd.seed)) {
    
    the.gene <- comb.d[ ,c(rnd.seed[i], 39766:39769)] # 1168 x 5 (RNA, EnvID, PotID, Time, Cultivar)
    #the.gene <- comb.d[, c(which(colnames(comb.d)=="Os01g0182600"), 39766:39769)]
    
    #pdf(paste(sprintf("%04d",i),".",colnames(the.gene)[1], ".pdf", sep=""))
    
    par(oma=c(3,3,3,3))
    layout(layout.mat,
           widths = c(2.7, rep(c(2,2,2.7), 5)), heights = c(1.2,rep(1,6)))
    
    #===========#
    # text part #
    #===========#
    par(mar=c(0,0,0,0))
    plot.new() ; plot.window(xlim=c(0,1), ylim=c(0,1))
    text(0,0.95,colnames(the.gene[1]),cex=1, adj=0)
    #text(0.2,0.95,sugargene_list$Gene.name[sugargene_list$Gene.ID..RAP.DB.==colnames(the.gene[1])], cex=1, adj=0)
    #text(0.2,0.95,clockgene_list$name[clockgene_list$ID==colnames(the.gene[1])], cex=1, adj=0)
    #text(0.2,0.95, sprintf("bsLASSO_No.%s",i), cex=1, adj=0)
    #text(0.2,0.95, sprintf("%s (%s)",photo_list$Protein[i], photo_list$Pathway[i]), cex=1, adj=0)
    #text(0.2,0.95, dreb_list$description[i], cex=0.8, adj=0)
    text(0,0.8, labels=paste("Description: ", des[rownames(des)==colnames(the.gene)[1], "Description"], sep=""), adj=c(0,1), cex=0.7)
    A.thal.gene <-des[rownames(des)==colnames(the.gene)[1], "A.thaliana_gene_id"]
    A.thal.short.des <- if(substr(as.character(A.thal.gene), 1, 2)=="AT") {
      if(nchar(des.A.thal[rownames(des.A.thal)==A.thal.gene, "computational_description"])>=2) { #priority: comp.des > Note > curator_sum
        des.A.thal[rownames(des.A.thal)==A.thal.gene, "computational_description"] } else {
          if(nchar(des.A.thal[rownames(des.A.thal)==A.thal.gene, "Note"])>=2) {
            des.A.thal[rownames(des.A.thal)==A.thal.gene, "Note"] } else {
              if(nchar(des.A.thal[rownames(des.A.thal)==A.thal.gene, "curator_summary"])>=2) {
                des.A.thal[rownames(des.A.thal)==A.thal.gene, "curator_summary"] } else {
                  des.A.thal[rownames(des.A.thal)==A.thal.gene, "curator_summary"] }}}} else {NA}
    text(0,0.65, labels=paste("The most similar A. thaliana gene: ", 
                              des[rownames(des)==colnames(the.gene)[1], "A.thaliana_gene_id"], " | ",
                              des[rownames(des)==colnames(the.gene)[1], "symbol"], "\n",
                              "Short description: ", A.thal.short.des, sep=""), adj=c(0,1), cex=0.6)
    text(c(0, 0.17, 0.34), 0.3, labels=c("Plot colour: ", "Koshihikari ","Takanari") ,col=c(1,4,2), adj=0, cex=0.8)
    
    #====================#
    # mean-variance plot #
    #====================#
    par(mar=c(5,5,0,1), mgp=c(2, 0.7, 0))
    x.v <- c(mean.exp[IDs.for.mv.plot], mean.exp[rnd.seed[i]])
    y.v <- c(var.log.rpm[IDs.for.mv.plot], var.log.rpm[rnd.seed[i]])
    pch.v <- c(rep(16, length(IDs.for.mv.plot)), 4)
    cex.v <- c(rep(1, length(IDs.for.mv.plot)), 2)
    red.v <- c(exp.ind.ratio[IDs.for.mv.plot], 0)
    blue.v <- c(1 - exp.ind.ratio[IDs.for.mv.plot], 0)
    alpha.v <- c(rep(0.2, length(IDs.for.mv.plot)), 1)
    plot(x.v, y.v, pch=pch.v, col=rgb(red.v,0,blue.v,alpha=alpha.v),cex=cex.v,
         xlab="mean(log2(rpm+1))", ylab="var(log2(rpm+1))", cex.lab=0.6, cex.axis=0.6)
    
    #==================================#
    # expression time series: all dark #
    #==================================#
    par(mar=c(1,0,0,1), mgp=c(1.5, 0.7, 0))
    xx <- seq(-2,19,length=50)
    y.lim.v <- c(min(comb.d[, rnd.seed]), max(comb.d[, rnd.seed]))
    
    ##########################
    ### L = na, D = 35-15C ###
    for (j in c(35,30,25,20,15)) {
      env.param <- env.settings[env.settings$L==0&env.settings$D.Tmp==j,]
      fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID,]
      col.v <- c("blue", "red")
      plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
           ylim=y.lim.v, xaxt=if(j==15) {"s"} else{"n"} , cex=1)
      
      sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                          fig.df[,1][fig.df$Cultivar=="Kos"])
      sp.tak <- splinefun(fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"]),
                          fig.df[,1][fig.df$Cultivar=="Tak"])
      kos.est <- sp.kos(xx)
      tak.est <- sp.tak(xx)
      polygon(c(20,20,-3,-3), c(y.lim.v, y.lim.v[2], y.lim.v[1]), col = rgb(0,0,0,alpha=0.1),border = F)
      lines(xx,kos.est,col="blue")
      lines(xx,tak.est,col="red")
    }
    
    ###############
    ### L and D ###
    for (j in c(35,30,25,20,15)) {  # dark temp
      for (k in c(20,25,30,35,40)) { # light temp
        
        if(j==k) {
          par(mar=c(0,0,0,0))
          plot.new() ; plot.window(xlim=c(0,1), ylim=c(0,1))
          arrows(0, 0.7, 0.2, 0.7, length=0.05, code=3)
          arrows(0.1, 0.55, 0.1, 0.35, length=0.05, code=3)
          text(0.25, 0.7, labels=paste("dark.tmp= ", j, sep=""), cex=1, adj=0)
          text(0.25, 0.45, labels=paste("light.tmp= ", k, sep=""), cex=1, adj=c(0))
          
        } else {
          env.param <- env.settings[env.settings$L.Tmp==k&env.settings$D.Tmp==j&env.settings$L!=0&env.settings$D!=0,]
          
          for (l in c(8,12,16)) {
            par(mar=if(l==16) {c(1,0,0,1)} else {c(1,0,0,0)})
            fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID[env.param$L==l],]
            plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
                 ylim=y.lim.v, cex=1, xaxt=if((j==15&&l==8)||(j==15&&l==12)) {"s"} else {"n"}, yaxt="n")
            sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                                fig.df[,1][fig.df$Cultivar=="Kos"])
            sp.tak <- splinefun(fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"]),
                                fig.df[,1][fig.df$Cultivar=="Tak"])
            kos.est <- sp.kos(xx)
            tak.est <- sp.tak(xx)
            polygon(c(21 - l, 21 - l,-3,-3), c(y.lim.v, y.lim.v[2], y.lim.v[1]), col = rgb(0,0,0,alpha=0.1),border = F)
            lines(xx,kos.est,col="blue")
            lines(xx,tak.est,col="red")
          }
        }
      }
    }
    
    #################
    ### all light ###
    for (j in c(20,25,30,35,40)) {
      env.param <- env.settings[env.settings$L.Tmp==j&env.settings$D==0,]
      fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID,]
      plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
           ylim=y.lim.v, cex=1)
      sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                          fig.df[,1][fig.df$Cultivar=="Kos"])
      sp.tak <- splinefun(fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"]),
                          fig.df[,1][fig.df$Cultivar=="Tak"])
      kos.est <- sp.kos(xx)
      tak.est <- sp.tak(xx)
      lines(xx,kos.est,col="blue")
      lines(xx,tak.est,col="red")
    }
  }
}

### Fig.1b ----
geneplot2 <- function(glist){
  tmp <- 1:39765
  names(tmp) <- colnames(log.rpm)[1:39765]
  rnd.seed <- tmp[as.character(glist)]
  rnd.seed <- rnd.seed[!is.na(rnd.seed)]
  
  for (i in 1:length(rnd.seed)) {
    
    the.gene <- comb.d[ ,c(rnd.seed[i], 39766:39769)] # 1168 x 5 (RNA, EnvID, PotID, Time, Cultivar)
    #the.gene <- comb.d[, c(which(colnames(comb.d)=="Os01g0182600"), 39766:39769)]
    
    #pdf(paste(sprintf("%04d",i),".",colnames(the.gene)[1], ".pdf", sep=""))
    
    par(oma=c(3,3,3,3))
    layout(layout.mat,
           widths = c(2.7, rep(c(2,2,2.7), 5)), heights = c(1.2,rep(1,6)))
    
    #===========#
    # text part #
    #===========#
    par(mar=c(0,0,0,0))
    plot.new() ; plot.window(xlim=c(0,1), ylim=c(0,1))
    text(0,0.95,colnames(the.gene[1]),cex=1, adj=0)
    #text(0.2,0.95,sugargene_list$Gene.name[sugargene_list$Gene.ID..RAP.DB.==colnames(the.gene[1])], cex=1, adj=0)
    #text(0.2,0.95,clockgene_list$name[clockgene_list$ID==colnames(the.gene[1])], cex=1, adj=0)
    #text(0.2,0.95, sprintf("bsLASSO_No.%s",i), cex=1, adj=0)
    #text(0.2,0.95, sprintf("%s (%s)",photo_list$Protein[i], photo_list$Pathway[i]), cex=1, adj=0)
    #text(0.2,0.95, dreb_list$description[i], cex=0.8, adj=0)
    text(0,0.8, labels=paste("Description: ", des[rownames(des)==colnames(the.gene)[1], "Description"], sep=""), adj=c(0,1), cex=0.7)
    A.thal.gene <-des[rownames(des)==colnames(the.gene)[1], "A.thaliana_gene_id"]
    A.thal.short.des <- if(substr(as.character(A.thal.gene), 1, 2)=="AT") {
      if(nchar(des.A.thal[rownames(des.A.thal)==A.thal.gene, "computational_description"])>=2) { #priority: comp.des > Note > curator_sum
        des.A.thal[rownames(des.A.thal)==A.thal.gene, "computational_description"] } else {
          if(nchar(des.A.thal[rownames(des.A.thal)==A.thal.gene, "Note"])>=2) {
            des.A.thal[rownames(des.A.thal)==A.thal.gene, "Note"] } else {
              if(nchar(des.A.thal[rownames(des.A.thal)==A.thal.gene, "curator_summary"])>=2) {
                des.A.thal[rownames(des.A.thal)==A.thal.gene, "curator_summary"] } else {
                  des.A.thal[rownames(des.A.thal)==A.thal.gene, "curator_summary"] }}}} else {NA}
    text(0,0.65, labels=paste("The most similar A. thaliana gene: ", 
                              des[rownames(des)==colnames(the.gene)[1], "A.thaliana_gene_id"], " | ",
                              des[rownames(des)==colnames(the.gene)[1], "symbol"], "\n",
                              "Short description: ", A.thal.short.des, sep=""), adj=c(0,1), cex=0.6)
    text(c(0, 0.17, 0.34), 0.3, labels=c("Plot colour: ", "Koshihikari ","Takanari") ,col=c("black",color5[1:2]), adj=0, cex=0.8)
    
    #====================#
    # mean-variance plot #
    #====================#
    par(mar=c(5,5,0,1), mgp=c(2, 0.7, 0))
    x.v <- c(mean.exp[IDs.for.mv.plot], mean.exp[rnd.seed[i]])
    y.v <- c(var.log.rpm[IDs.for.mv.plot], var.log.rpm[rnd.seed[i]])
    pch.v <- c(rep(16, length(IDs.for.mv.plot)), 4)
    cex.v <- c(rep(1, length(IDs.for.mv.plot)), 2)
    red.v <- c(exp.ind.ratio[IDs.for.mv.plot], 0)
    blue.v <- c(1 - exp.ind.ratio[IDs.for.mv.plot], 0)
    alpha.v <- c(rep(0.2, length(IDs.for.mv.plot)), 1)
    plot(x.v, y.v, pch=pch.v, col=rgb(red.v,0,blue.v,alpha=alpha.v),cex=cex.v,
         xlab="mean(log2(rpm+1))", ylab="var(log2(rpm+1))", cex.lab=0.6, cex.axis=0.6)
    
    #==================================#
    # expression time series: all dark #
    #==================================#
    par(mar=c(1,0,0,1), mgp=c(1.5, 0.7, 0))
    xx <- seq(-2,19,length=50)
    y.lim.v <- c(min(comb.d[, rnd.seed]), max(comb.d[, rnd.seed]))
    
    ##########################
    ### L = na, D = 35-15C ###
    for (j in c(35,30,25,20,15)) {
      env.param <- env.settings[env.settings$L==0&env.settings$D.Tmp==j,]
      fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID,]
      col.v <- color5[1:2]
      plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
           ylim=y.lim.v, xaxt=if(j==15) {"s"} else{"n"} , cex=1)
      
      sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                          fig.df[,1][fig.df$Cultivar=="Kos"])
      sp.tak <- splinefun(fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"]),
                          fig.df[,1][fig.df$Cultivar=="Tak"])
      kos.est <- sp.kos(xx)
      tak.est <- sp.tak(xx)
      polygon(c(20,20,-3,-3), c(y.lim.v, y.lim.v[2], y.lim.v[1]), col = rgb(0,0,0,alpha=0.1),border = F)
      lines(xx,kos.est,col=color5[1])
      lines(xx,tak.est,col=color5[2])
    }
    
    ###############
    ### L and D ###
    for (j in c(35,30,25,20,15)) {  # dark temp
      for (k in c(20,25,30,35,40)) { # light temp
        
        if(j==k) {
          par(mar=c(0,0,0,0))
          plot.new() ; plot.window(xlim=c(0,1), ylim=c(0,1))
          arrows(0, 0.7, 0.2, 0.7, length=0.05, code=3)
          arrows(0.1, 0.55, 0.1, 0.35, length=0.05, code=3)
          text(0.25, 0.7, labels=paste("dark.tmp= ", j, sep=""), cex=1, adj=0)
          text(0.25, 0.45, labels=paste("light.tmp= ", k, sep=""), cex=1, adj=c(0))
          
        } else {
          env.param <- env.settings[env.settings$L.Tmp==k&env.settings$D.Tmp==j&env.settings$L!=0&env.settings$D!=0,]
          
          for (l in c(8,12,16)) {
            par(mar=if(l==16) {c(1,0,0,1)} else {c(1,0,0,0)})
            fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID[env.param$L==l],]
            plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
                 ylim=y.lim.v, cex=1, xaxt=if((j==15&&l==8)||(j==15&&l==12)) {"s"} else {"n"}, yaxt="n")
            sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                                fig.df[,1][fig.df$Cultivar=="Kos"])
            sp.tak <- splinefun(fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"]),
                                fig.df[,1][fig.df$Cultivar=="Tak"])
            kos.est <- sp.kos(xx)
            tak.est <- sp.tak(xx)
            polygon(c(21 - l, 21 - l,-3,-3), c(y.lim.v, y.lim.v[2], y.lim.v[1]), col = rgb(0,0,0,alpha=0.1),border = F)
            lines(xx,kos.est,col=color5[1])
            lines(xx,tak.est,col=color5[2])
          }
        }
      }
    }
    
    #################
    ### all light ###
    for (j in c(20,25,30,35,40)) {
      env.param <- env.settings[env.settings$L.Tmp==j&env.settings$D==0,]
      fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID,]
      plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
           ylim=y.lim.v, cex=1)
      sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                          fig.df[,1][fig.df$Cultivar=="Kos"])
      sp.tak <- splinefun(fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"]),
                          fig.df[,1][fig.df$Cultivar=="Tak"])
      kos.est <- sp.kos(xx)
      tak.est <- sp.tak(xx)
      lines(xx,kos.est,col=color5[1])
      lines(xx,tak.est,col=color5[2])
    }
  }
}

pdf("../output/Fig.1b.pdf",height=5,width=5)
geneplot2("Os01g0182600") #OsGI
dev.off()

### gene filtering ----
data <- t(log.rpm[,158:39765])
dat <- data[apply(data, 1, mean) > 1,] # 17742 genes
gene_17742 <- rownames(dat)

### total_read_count_sample ----
a <- rawcnt.lab[names(mean.exp)[158:39765],]
b <- log10(apply(a, 2, sum))

exp <- factor(rep("Exp", 1167))
d <- data.frame(exp, reads=b)

g <- ggplot(d, aes(x = reads, fill=exp))+
  geom_histogram(color="black", linewidth=0.25)+
  scale_fill_manual(values = c("white","gray"))+
  xlim(4.8,7)+
  labs(x="Mapped read count (log10)", y="Frequency")+
  theme_cowplot(font_size = 7, line_size=0.25)+
  theme(plot.title = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")
p <- plot_grid(g)

ggsave(p, file="../output/Fig.S2a.pdf", # Fig.S2a
       width = 60, height = 40, units="mm")

### mean_read_count_gene ----
a <- rawcnt.lab[names(mean.exp)[158:39765],]
e <- apply(a, 1, mean)
e <- e[e>0]

d <- data.frame(exp=log10(e))

g <- ggplot(d, aes(x = exp))+
  geom_histogram(bins=50, color="black", linewidth=0.25)+
  geom_vline(xintercept = 1, color=color5[1], size =0.5)+
  scale_fill_manual(values = "gray")+
  xlim(-6,5)+
  ylim(0,3000)+
  labs(x="Mean read count (log10)", y="Frequency")+
  theme_cowplot(font_size = 7, line_size=0.25)+
  theme(plot.title = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "top",
        legend.justification = "center")
p <- plot_grid(g)

ggsave(p, file="../output/Fig.S2c.pdf",# Fig.S2c
       width = 60, height = 40, units="mm")

### Fig. 1c-e t-SNE ----
tmp <- env.settings[,c(1,2)]
colnames(tmp) <- c("EnvID", "day_length")
ttmp <- left_join(at.lab.additional,tmp, by="EnvID")

set.seed(1)
tsne <- Rtsne(t(dat), check_duplicates = FALSE, verbose=F)
#save(tsne, file="../output/240906_tsne")

#load("../output/240906_tsne")
d <- data.frame(tsne$Y, ttmp)
d$Time <- factor(d$Time, levels = c("22","1","4","7","10","13","16","19"))
d$day_length <- factor(d$day_length,
                       levels = c("0L24D","8L16D","12L12D","16L8D", "24L0D"))
d$day_length_LD <- factor(paste(d$day_length,d$LD, sep = "_"),
                          levels = c("0L24D_D","8L16D_L","8L16D_D","12L12D_L",
                                     "12L12D_D","16L8D_L", "16L8D_D", "24L0D_L"))

gg.pca.daylength <- ggplot(data=d,aes(x=d[,1], y=d[,2],
                                      col=day_length, fill=day_length_LD,
                                      shape=Cultivar))+
  geom_point(size=0.8, stroke =0.6) +
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=color5)+
  scale_fill_manual(values=c(color5[1],"transparent",color5[2],"transparent",color5[3],
                             "transparent",color5[4],color5[5]))+
  labs(x="tSNE1", y="tSNE2")+
  theme_cowplot(font_size = 7, line_size = 0.5)+
  panel_border(colour=1, size=1)+
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7))#,
#legend.position = "none")

color8 <- brewer.pal(8,"Accent")
d$Time_LD <- factor(paste(d$Time,d$LD, sep = "_"),
                    levels = c("22_D","22_L","1_D","1_L",
                               "4_D","4_L","7_D","7_L",
                               "10_D","10_L","13_D","13_L",
                               "16_D","16_L","19_D","19_L"))

gg.pca.samplingtime <- ggplot(data=d,aes(x=d[,1], y=d[,2],
                                         col=Time, fill=Time_LD, shape=Cultivar))+
  geom_point(size=0.8, stroke =0.6) +
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=color8)+
  scale_fill_manual(values=c(color8[1],"transparent",color8[2],"transparent",
                             color8[3],"transparent",color8[4],"transparent",
                             color8[5],"transparent",color8[6],"transparent",
                             color8[7],"transparent",color8[8],"transparent"))+
  labs(x="tSNE1", y="tSNE2")+
  theme_cowplot(font_size = 8, line_size = 0.5)+
  panel_border(colour=1, size=1)+
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7))#,
#legend.position = "none")

d$Tmp <- factor(d$Tmp, levels = c("15","20","25","30","35","40"))
d$Tmp_LD <- factor(paste(d$Tmp,d$LD, sep = "_"),
                   levels = c("15_D","15_L","20_D","20_L",
                              "25_D","25_L","30_D","30_L",
                              "35_D","35_L","40_D","40_L"))
color6 <- brewer.pal(6,"Oranges")

gg.pca.temp <- ggplot(data=d,aes(x=d[,1], y=d[,2],
                                 col=Tmp, fill=Tmp_LD, shape=Cultivar))+
  geom_point(size=0.8, stroke =0.6) +
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=color6)+
  scale_fill_manual(values=c(color6[1],color6[2],"transparent",
                             color6[3],"transparent",color6[4],"transparent",
                             color6[5],"transparent","transparent"))+
  labs(x="tSNE1", y="tSNE2")+
  theme_cowplot(font_size = 8, line_size = 0.5)+
  panel_border(colour=1, size=1)+
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7))#,
#legend.position = "none")

gg.dl <- gg.pca.daylength + theme(legend.position = "none")
gg.sm <- gg.pca.samplingtime + theme(legend.position = "none")
gg.tm <- gg.pca.temp + theme(legend.position = "none")

p <- plot_grid(gg.dl, gg.sm, gg.tm, ncol=1)
ggsave(p, file=sprintf("../output/Fig.1c-e_nolegend.pdf"), 
       height = 130, width = 46, units = "mm")

pp <- plot_grid(gg.pca.daylength)
ggsave(pp, file=sprintf("../output/Fig.1c_legend.pdf"), 
       height = 60, width = 75, units = "mm")

pp <- plot_grid(gg.pca.samplingtime)
ggsave(pp, file=sprintf("../output/Fig.1d_legend.pdf"), 
       height = 105, width = 90, units = "mm")

pp <- plot_grid(gg.pca.temp)
ggsave(pp, file=sprintf("../output/Fig.1e_legend.pdf"), 
       height = 75, width = 90, units = "mm")

### amplitude ----
kos.amp <- matrix(NA, nrow=73, ncol=length(gene_17742))
colnames(kos.amp) <- gene_17742
tak.amp <- matrix(NA, nrow=73, ncol=length(gene_17742))
colnames(tak.amp) <- gene_17742

xx <- seq(-3,21,1.5) # every 90 min

for (j in 1:73){
  env.param <- env.settings[j,]
  
  for (i in 1:length(gene_17742)){
    the.gene <- cbind(log.rpm[,gene_17742[i]],log.rpm[,39766:39769])
    fig.df <- the.gene[the.gene$EnvID==env.param$Env_ID,]
    
    kos.time <- fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"])+0.5
    sp.kos <- smooth.spline(kos.time,fig.df[,1][fig.df$Cultivar=="Kos"],
                            spar = 0.3)
    kos.est <- predict(sp.kos, xx)$y
    
    tak.time <- fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"])+0.5
    sp.tak <- smooth.spline(tak.time,fig.df[,1][fig.df$Cultivar=="Tak"],
                            spar = 0.3)
    tak.est <- predict(sp.tak, xx)$y
    
    kos.amp[j,i] <- max(kos.est)-min(kos.est)
    tak.amp[j,i] <- max(tak.est)-min(tak.est)
    
  }
  cat(sprintf("%s\n",j))  
}

save(kos.amp, file="../output/kos.amp_240906") #every 90 min
save(tak.amp, file="../output/tak.amp_240906")

### amplitude fig ----
#load("../output/kos.amp_240906")
#load("../output/tak.amp_240906")

kos.amp.mean <- apply(kos.amp, 2, mean)
tak.amp.mean <- apply(tak.amp, 2, mean)

gene_amplitude <- gene_17742[(kos.amp.mean > 2) & (tak.amp.mean > 2)] # 12770 genes

kos.amp2 <- kos.amp[,gene_amplitude]
tak.amp2 <- tak.amp[,gene_amplitude]

kos.amp.mean2 <- apply(kos.amp2, 1, mean)
tak.amp.mean2 <- apply(tak.amp2, 1, mean)

label_Ltemp <- c("20K","20T","25K","25T","30K","30T","35K","35T","40K","40D")
label_Dtemp <- c("15K","15T","20K","20T","25K","25T","30K","30T","35K","35T")

# amp L temp
data <- data.frame(env=rep(env.settings$L.Tmp,2),
                   value=c(kos.amp.mean2, tak.amp.mean2),
                   geno=rep(c("K","T"),each=73))[c(6:68,79:141),]

g_amp_L <- ggplot(data, aes(x = paste(env,geno), y = value, color=geno)) +
  geom_boxplot(outlier.shape = NA, color="black") +
  geom_beeswarm(size = 0.5) +
  scale_color_manual(values = color5[c(1,2)])+
  scale_x_discrete(labels=label_Ltemp)+
  scale_y_continuous(breaks=seq(0.5,4.5,length=5),limits=c(1.5,4.5))+
  xlab("Light_period temperature")+
  ylab("Average amplitude of\ndiurnal oscillation")+
  theme_cowplot() +
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

# amp D temp
data <- data.frame(env=rep(env.settings$D.Tmp,2), 
                   value=c(kos.amp.mean2, tak.amp.mean2),
                   geno=rep(c("K","T"),each=73))[c(6:68,79:141),]

g_amp_D <- ggplot(data, aes(x = paste(env,geno), y = value, color=geno)) +
  geom_boxplot(outlier.shape = NA, color="black") +
  geom_beeswarm(size = 0.5) +
  scale_color_manual(values = color5[c(1,2)])+
  scale_x_discrete(labels=label_Dtemp)+
  scale_y_continuous(breaks=seq(0.5,4.5,length=5),limits=c(1.5,4.5))+
  xlab("Dark period temperature")+
  ylab("Average amplitude of\ndiurnal oscillation")+
  theme_cowplot() +
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

# amp day length
data <- data.frame(env=factor(paste(rep(env.settings$L,2),rep(c("K","T"),each=73),sep = "_"),
                              levels = c("0_K","0_T","8_K","8_T","12_K","12_T",
                                         "16_K","16_T","24_K","24_T")), 
                   value=c(kos.amp.mean2, tak.amp.mean2),
                   geno=rep(c("K","T"),each=73))

g_amp_DL <- ggplot(data, aes(x=env, y = value, color=geno)) +
  geom_boxplot(outlier.shape = NA, color="black") +
  geom_beeswarm(size = 0.5) +
  scale_color_manual(values = color5[c(1,2)])+
  scale_y_continuous(breaks=seq(0.5,4.5,length=5),limits=c(1.5,4.5))+
  xlab("Light period length")+
  ylab("Average amplitude of\ndiurnal oscillation")+
  theme_cowplot() +
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

### cultivar specific genes ----
kos.mean <- matrix(NA, nrow=73, ncol=length(gene_17742))
colnames(kos.mean) <- gene_17742
tak.mean <- matrix(NA, nrow=73, ncol=length(gene_17742))
colnames(tak.mean) <- gene_17742

for (j in 1:73){
  env.param <- env.settings[j,]
  
  for (i in 1:length(gene_17742)){
    the.gene <- cbind(log.rpm[,gene_17742[i]],log.rpm[,39766:39769])
    fig.df <- the.gene[the.gene$EnvID==env.param$Env_ID,]
    kos.mean[j,i] <- mean(fig.df$`log.rpm[, gene_17742[i]]`[1:8],na.rm = T)
    tak.mean[j,i] <- mean(fig.df$`log.rpm[, gene_17742[i]]`[9:16],na.rm = T)
  }
  cat(sprintf("%s\n",j))  
}

kos.mean.73 <- apply(kos.mean,2,mean)
tak.mean.73 <- apply(tak.mean,2,mean)

kos_specific <- gene_17742[tak.mean.73 < 0.5 & ((kos.mean.73-tak.mean.73) > 2)] # 483 genes
tak_specific <- gene_17742[kos.mean.73 < 0.5 & ((tak.mean.73-kos.mean.73) > 2)] # 110 genes

### t-test without specific genes ----
ttest.pvalue <- matrix(NA, nrow=73, ncol=17742)
colnames(ttest.pvalue) <- gene_17742

for (j in c(1:31,33:73)){
  env.param <- env.settings[j,]
  
  for (i in 1:17742){
    the.gene <- cbind(log.rpm[,gene_17742[i]],log.rpm[,39766:39769])
    fig.df <- the.gene[the.gene$EnvID==env.param$Env_ID,]
    ko <- fig.df$`log.rpm[, gene_17742[i]]`[1:8]
    ta <- fig.df$`log.rpm[, gene_17742[i]]`[9:16]
    ttest.pvalue[j,i] <- t.test(ko,ta, paired = T)$p.value
  }
  cat(sprintf("%s: %s\n",j, Sys.time()))  
}

# EnvID=32 
j =32
env.param <- env.settings[j,]

for (i in 1:17742){
  the.gene <- cbind(log.rpm[,gene_17742[i]],log.rpm[,39766:39769])
  fig.df <- the.gene[the.gene$EnvID==env.param$Env_ID,]
  ko <- fig.df$`log.rpm[, gene_17742[i]]`[c(1,2,4,5,6,7,8)]
  ta <- fig.df$`log.rpm[, gene_17742[i]]`[9:15]
  ttest.pvalue[j,i] <- t.test(ko,ta, paired = T)$p.value
}

gene_not_specific <- gene_17742[-which(gene_17742 %in% c(kos_specific,tak_specific))] #17149 genes

ttest.qvalue <- apply(ttest.pvalue,1, function(x) p.adjust(x,method = "BH"))

aa <- function(x){
  length(which(x < 0.05))
}

ttest.signif <- apply(ttest.qvalue, 2, aa)
ttest.signif.gene <- ttest.qvalue < 0.05

tmp <- t((kos.mean[,gene_not_specific]-tak.mean[,gene_not_specific]) > 0) #kos > tak
tmp2 <- t((kos.mean[,gene_not_specific]-tak.mean[,gene_not_specific]) < 0) #tak > kos

ttest.signif.gene.kos <- ttest.signif.gene[gene_not_specific,] & tmp
ttest.signif.gene.tak <- ttest.signif.gene[gene_not_specific,] & tmp2

ttest.signif.kos <- apply(ttest.signif.gene.kos, 2, sum)
ttest.signif.tak <- apply(ttest.signif.gene.tak, 2, sum)

### ttest fig ----
label_Ltemp <- c("20K","20T","25K","25T","30K","30T","35K","35T","40K","40D")
label_Dtemp <- c("15K","15T","20K","20T","25K","25T","30K","30T","35K","35T")

# ttest L temp
data <- data.frame(env=rep(env.settings$L.Tmp,2), 
                   value=c(ttest.signif.kos, ttest.signif.tak),
                   geno=rep(c("K","T"),each=73))
data <- data[!is.na(data$env),]

g_amp_L_ttest <- ggplot(data, aes(x = paste(env,geno), y = value, color=geno)) +
  geom_boxplot(outlier.shape = NA, color="black") +
  geom_beeswarm(size = 0.5) +
  scale_color_manual(values = color5[c(1,2)])+
  scale_x_discrete(labels=label_Ltemp)+
  scale_y_continuous(breaks=seq(0,6000,length=4),limits=c(0,6000))+
  #ggtitle("Amplitude Day air temperature") +
  xlab("Light period temperature")+
  ylab("The number of DEGs")+
  theme_cowplot() +
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

# amp D temp
data <- data.frame(env=rep(env.settings$D.Tmp,2), 
                   value=c(ttest.signif.kos, ttest.signif.tak),
                   geno=rep(c("K","T"),each=73))
data <- data[!is.na(data$env),]

g_amp_D_ttest <- ggplot(data, aes(x = paste(env,geno), y = value, color=geno)) +
  geom_boxplot(outlier.shape = NA, color="black") +
  geom_beeswarm(size = 0.5) +
  scale_color_manual(values = color5[c(1,2)])+
  scale_x_discrete(labels=label_Dtemp)+
  scale_y_continuous(breaks=seq(0,6000,length=4),limits=c(0,6000))+
  #ggtitle("Amplitude Night air temperature") +
  xlab("Dark period temperature")+
  ylab("The number of DEGs")+
  theme_cowplot() +
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

### day length
data <- data.frame(env=factor(paste(rep(env.settings$L,2),rep(c("K","T"),each=73),sep = "_"),
                              levels = c("0_K","0_T","8_K","8_T","12_K","12_T",
                                         "16_K","16_T","24_K","24_T")), 
                   value=c(ttest.signif.kos, ttest.signif.tak),
                   geno=rep(c("K","T"),each=73))

g_amp_DL_ttest <- ggplot(data, aes(x=env, y = value, color=geno)) +
  geom_boxplot(outlier.shape = NA, color="black") +
  geom_beeswarm(size = 0.5) +
  scale_color_manual(values = color5[c(1,2)])+
  scale_y_continuous(breaks=seq(0,6000,length=4),limits=c(0,6000))+
  #ggtitle("Amplitude: day length") +
  xlab("Light period length")+
  ylab("The number of DEGs")+
  theme_cowplot() +
  theme(plot.title = element_text(size=7, hjust=0.5),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position = "none")

### Fig.1f-k amplitude and t-test ----
p <- plot_grid(g_amp_DL, g_amp_L,g_amp_D, 
               g_amp_DL_ttest,g_amp_L_ttest,g_amp_D_ttest,
               ncol=3, align = "hv")

ggsave(p, file=sprintf("../output/Fig.1f-k.pdf"), # Fig.1f-k 
       height = 70, width = 180, units = "mm")

### Table.S6-7, Fig.2h correlation with air temperature ----
correlation <- rep(NA, nrow(dat))
for(i in 1:nrow(dat)){
  correlation[i] <- cor(at.lab.additional$Tmp,dat[i,])
}

names(correlation) <- rownames(dat)
correlation <- sort(correlation)

correlation_high <- round(correlation[abs(correlation) >= 0.5],digits = 2)

a <- des[names(correlation_high),]
b <- cbind(correlation_high,a)

write.table(b, "../output/Table.S6-7.tsv", # 
            sep="\t", row.names = F)

# figure
data <- data.frame(xx=correlation)
gg.temp <- ggplot(data, aes(x=xx)) +
  geom_histogram(binwidth = 0.1,colour="black",fill="white") +
  xlim(-1,1)+
  xlab("r")+
  ylab("Gene frequency")+
  theme_cowplot() +
  theme(plot.title = element_text(size=6, hjust=0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")  
print(gg.temp)

ggsave(gg.temp, file="../output/Fig.2h.pdf",
       width = 45, height = 40, units = "mm")

### Table S8 representative 500 genes ----
load("../input/T.to.R.genes.80_Kos") # 92
load("../input/T.to.R.genes.80_Tak") # 101

des_genelist_T.to.R.genes.80_kos <- des[setdiff(T.to.R.genes.80_Kos,T.to.R.genes.80_Tak),c(2:3,8)]
des_genelist_T.to.R.genes.80_tak <- des[setdiff(T.to.R.genes.80_Tak,T.to.R.genes.80_Kos),c(2:3,8)]
des_genelist_T.to.R.genes.80_kos_tak <- des[intersect(T.to.R.genes.80_Kos,T.to.R.genes.80_Tak),c(2:3,8)]

### Table 8 
write.table(des_genelist_T.to.R.genes.80_kos,"../output/Table.S8_kos.tsv",
            sep="\t", row.names = F)
write.table(des_genelist_T.to.R.genes.80_tak,"../output/Table.S8_tak.tsv",
            sep="\t", row.names = F)
write.table(des_genelist_T.to.R.genes.80_kos_tak,"../output/Table.S8_kos_tak.tsv",
            sep="\t", row.names = F)

# The number of genes belong to the cluster of T.to.R.genes.80_kos and T.to.R.genes.80_tak genes
load("../input/gene.cluster")
id_converter <- read.delim("../input/230728_all_genes_id_conversion.tsv")
rownames(id_converter) <- id_converter$microarray

exemplars_new <- exemplars
names(exemplars_new) <- id_converter[id_converter$microarray %in% names(exemplars),"rnaseq"]
exemplars_new <- exemplars_new[!is.na(names(exemplars_new))]

names(exemplars[exemplars %in% exemplars_new[T.to.R.genes.80_Kos]]) #3636 genes
names(exemplars[exemplars %in% exemplars_new[T.to.R.genes.80_Tak]]) #3988 genes
names(exemplars[exemplars %in% exemplars_new[intersect(T.to.R.genes.80_Kos,T.to.R.genes.80_Tak)]]) #1630 genes
