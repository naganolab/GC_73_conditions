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
library("NSM3")
library("reshape2")
library("multcompView")
library("vioplot")

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


### Fig.2a ----
layout.mat2 <- matrix(c(
  1,  0, 6,7,8,    0, 9,10,11,  0, 12,13,14, 0, 15,15,15, 0, 16,17,18,0,
  2,  0, 19,20,21, 0, 22,23,24, 0, 25,25,25, 0, 26,27,28, 0, 29,30,31,0,
  3,  0, 32,33,34, 0, 35,35,35, 0, 36,37,38, 0, 39,40,41, 0, 42,43,44,0,
  4,  0, 45,45,45, 0, 46,47,48, 0, 49,50,51, 0, 52,53,54, 0, 55,56,57,0,
  5,  0, 58,59,60, 0, 61,62,63, 0, 64,65,66, 0, 67,68,69, 0, 70,71,72,0,
  0,  0, 0,0,73,   0, 0,0,74,   0, 0,0,75,   0, 0,0,76,   0, 0,0,77,0), byrow=T,nrow=6,ncol=22)

geneplot2 <- function(glist){
  tmp <- 1:39765
  names(tmp) <- colnames(log.rpm)[1:39765]
  rnd.seed <- tmp[as.character(glist)]
  rnd.seed <- rnd.seed[!is.na(rnd.seed)]
  
  for (i in 1:length(rnd.seed)) {
    
    the.gene <- comb.d[ ,c(rnd.seed[i], 39766:39769)] # 1168 x 5 (RNA, EnvID, PotID, Time, Cultivar)
    #the.gene <- comb.d[, c(which(colnames(comb.d)=="Os01g0182600"), 39766:39769)]
    
    #pdf(paste(sprintf("%04d",i),".",colnames(the.gene)[1], ".pdf", sep=""))
    
    par(oma=c(1,3,3,0))
    par(tcl = -0.3)
    layout(layout.mat2,
           widths = c(2, rep(c(0.5,2,2,2), 5),0.5), heights = c(rep(1,6)))
    
    #==================================#
    # expression time series: all dark #
    #==================================#
    par(mar=c(0.4,0,0,0), mgp=c(1.5, 0.2, 0))
    xx <- seq(-2,19,length=50)
    y.lim.v <- c(min(comb.d[, rnd.seed]), max(comb.d[, rnd.seed]))
    
    ##########################
    ### L = na, D = 35-15C ###
    for (j in c(35,30,25,20,15)) {
      env.param <- env.settings[env.settings$L==0&env.settings$D.Tmp==j,]
      fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID,]
      col.v <- color5[1:2]
      plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
           ylim=y.lim.v, xaxt="n", yaxt = "n", cex=0.5)
      # # Y axis
      # axis(2, at = c(0, 4, 8), cex.axis = 0.5)
      # # X axis
      # if (j == 15){
      #   axis(1, at = c(0, 15), labels = FALSE, cex.axis = 0.5)
      # }
      axis(1,at = c(0, 15),labels = FALSE, tck = -0.1)
      axis(2,at = c(0, 4, 8),labels = FALSE, tck = -0.1)
      
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
          # arrows(0, 0.7, 0.2, 0.7, length=0.05, code=3)
          # arrows(0.1, 0.55, 0.1, 0.35, length=0.05, code=3)
          # text(0.25, 0.7, labels=paste("dark.tmp= ", j, sep=""), cex=1, adj=0)
          # text(0.25, 0.45, labels=paste("light.tmp= ", k, sep=""), cex=1, adj=c(0))
          
        } else {
          env.param <- env.settings[env.settings$L.Tmp==k&env.settings$D.Tmp==j&env.settings$L!=0&env.settings$D!=0,]
          
          for (l in c(8,12,16)) {
            par(mar = c(0.4,0,0,0))
            fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID[env.param$L==l],]
            plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
                 ylim=y.lim.v, cex=0.5, xaxt="n", yaxt="n")
            # # X axis
            # if((j==15&&l==8)||(j==15&&l==12)){
            #   axis(1, at = c(0, 15), cex.axis = 0.5)
            axis(1, at = c(0, 15), labels = FALSE, tck = -0.1)
            if(l==8){axis(2,at = c(0, 4, 8),labels = FALSE, tck = -0.1)}
            
            sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                                fig.df[,1][fig.df$Cultivar=="Kos"])
            sp.tak <- splinefun(fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"]),
                                fig.df[,1][fig.df$Cultivar=="Tak"])
            kos.est <- sp.kos(xx)
            tak.est <- sp.tak(xx)
            polygon(c(21 - l, 21 - l,-3,-3), c(y.lim.v, y.lim.v[2], y.lim.v[1]), col = rgb(0,0,0,alpha=0.15),border = F)
            lines(xx,kos.est,col=color5[1])
            lines(xx,tak.est,col=color5[2])
          }
        }
      }
    }
    
    #################
    ### all light ###
    for (j in c(20,25,30,35,40)) {
        par(mar=c(0.4,0,0,0))
      env.param <- env.settings[env.settings$L.Tmp==j&env.settings$D==0,]
      fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID,]
      plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
           ylim=y.lim.v, cex=0.5, xaxt = "n", yaxt = "n")
      # axis(1, at = c(0, 15), cex.axis = 0.5, tcl = -0.3)
      # axis(2, at = c(0,4,8), cex.axis = 0.5)
      axis(1,at = c(0, 15),labels = FALSE, tck = -0.1)
      axis(2,at = c(0, 4, 8),labels = FALSE, tck = -0.1)
      sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                          fig.df[,1][fig.df$Cultivar=="Kos"])
      sp.tak <- splinefun(fig.df$Time[fig.df$Cultivar=="Tak"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Tak"]),
                          fig.df[,1][fig.df$Cultivar=="Tak"])
      kos.est <- sp.kos(xx)
      tak.est <- sp.tak(xx)
      lines(xx,kos.est,col=color5[1])
      lines(xx,tak.est,col=color5[2])
    }
      # label
  # left
  #mtext(expression("Gene expression (" * log[2](rpm) + 1 * ")"), side=2, outer=TRUE, line=3, cex=0.5)
  mtext("Dark temperature (??C)", side=2, outer=TRUE, line=1.8, cex=0.5)
  mtext(c("Continuous\nlight","15", "20", "25", "30", "35"), side=2, outer=TRUE, line=0.8, at=seq(1/12, 11/12, length.out=6), cex=0.5, las=3)
  
  # upper
  mtext("Light temperature (??C)", side=3, outer=TRUE, line=1.8, cex=0.5)
  mtext(c("Continuous\ndark","20", "25", "30", "35", "40"), side=3, outer=TRUE, line=0.8, at=c(1/35,seq(5.5/35, 31.5/35, length.out=5)), cex=0.5)
  mtext(c("8L", "12L", "16L", "8L", "12L", "16L","8L", "12L", "16L","8L", "12L", "16L"),
        side=3, outer=TRUE, line=0, at=c(3.5/35,5.5/35,7.5/35,
                                         10/35,12/35,14/35,
                                         16.5/35, 18.5/35,20.5/35, 
                                         23/35,25/35,27/35,
                                         29.5/35,31.5/35,33.5/35
        ),
        cex=0.5)
  
  # bottom
  #mtext("time of day (h)", side=1, outer=TRUE, line=0, cex=0.5)
    
  }

}

pdf("../output/Fig.2a.pdf", width=5.3, height=2.8)
geneplot2("Os01g0182600") #OsGI
dev.off()

pdf("../output/Fig.S6.pdf", width=5.3, height=2.8)
geneplot2("Os01g0840100") #Hsp70
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

ggsave(p, file="../output/Fig.S1a.pdf", 
       width = 60, height = 40, units="mm")

### mean_read_count_gene ----
a <- log.rpm[,names(mean.exp)[158:39765]]
e <- apply(a, 2, mean)
e <- e[e>0]

d <- data.frame(exp=e)

g <- ggplot(d, aes(x = exp)) +
  geom_histogram(bins = 14, boundary = 1, color = "black", fill = "gray", linewidth = 0.25) +
  geom_vline(xintercept = 1, color = color5[1], size = 0.5) +
  coord_cartesian(xlim = c(0, 14)) +
  labs(x = "Mean expression (log2(rpm+1))", y = "Frequency") +
  theme_cowplot(font_size = 7, line_size = 0.25) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position = "top",
    legend.justification = "center"
  )
p <- plot_grid(g)
p

ggsave(p, file="../output/Fig.S1c.pdf",
       width = 50, height = 35, units="mm")

### PCA 73 conditions ----
data <- t(log.rpm[,158:39765])
dat <- data[apply(data, 1, mean) > 1,] # 17742 genes
gene_17742 <- rownames(dat)

pca <- prcomp(t(dat), scale = TRUE)
pc <- pca$x

pca_var <- pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_percent <- (pca_var / sum(pca_var)) * 100  # Convert to percentage

# Print the variance explained by PC1 and PC2
cat("PC1 Contribution:", round(pca_var_percent[1], 2), "%\n")
cat("PC2 Contribution:", round(pca_var_percent[2], 2), "%\n")

### PCA function
plot_pca_base_env <- function(env_param, hide_axis_labels = FALSE, hide_yaxis = FALSE) {
  
  # Fixed dataset variables
  d <- data.frame(
    Cultivar = at.lab.additional$Cultivar[at.lab.additional$EnvID == env_param],
    Time = at.lab.additional$Time[at.lab.additional$EnvID == env_param],
    LD = at.lab.additional$LD[at.lab.additional$EnvID == env_param],
    PC1 = pc[at.lab.additional$EnvID == env_param, "PC1"],
    PC2 = pc[at.lab.additional$EnvID == env_param, "PC2"]
  )
  
  # Color and shape settings
  color_palette <- setNames(c("#ff4b00", "#4dc4ff"), unique(d$Cultivar))
  shapes <- setNames(c(21, 22), unique(d$Cultivar))
  
  # Axis label and tick options
  xlab_text <- if (hide_axis_labels) "" else "PC1"
  ylab_text <- if (hide_axis_labels) "" else "PC2"
  xaxt_opt <- if (hide_axis_labels) "n" else "s"
  yaxt_opt <- if (hide_axis_labels) "n" else "s"
  
  # Plot
  plot(d$PC1, d$PC2, type = "n",
       xlab = xlab_text, ylab = ylab_text,
       xlim = c(-130, 150), ylim = c(-150, 60),
       xaxt = xaxt_opt, yaxt = yaxt_opt)
  
  # Add tick marks only (without labels) if axis labels are hidden
  if (hide_axis_labels) {
    axis(1, labels = FALSE, tck = -0.05)  # longer x-axis ticks
    if (hide_yaxis==F){
    axis(2, labels = FALSE, tck = -0.05)  # longer y-axis ticks
    }
  }
  # Draw lines by Cultivar
  for (cultivar in unique(d$Cultivar)) {
    subset_data <- d[d$Cultivar == cultivar, ]
    if (nrow(subset_data) > 1) {
      lines(subset_data$PC1, subset_data$PC2,
            col = color_palette[cultivar],
            lwd = 0.5)
    }
  }
  
  # Add points with alpha by LD
  for (i in 1:nrow(d)) {
    alpha_value <- ifelse(d$LD[i] == "D", 0.3, 1.0)
    point_color <- adjustcolor(color_palette[d$Cultivar[i]], alpha.f = alpha_value)
    
    points(d$PC1[i], d$PC2[i],
           pch = shapes[d$Cultivar[i]],
           col = point_color,
           bg = point_color,
           cex = 0.5)
  }
}

### Fig. 2b ----
pdf("../output/Fig.2b.pdf", width=5.3, height=2.8)

par(oma=c(1,3,3,0))
#layout(layout.mat2,
#       widths = c(2, rep(c(2,2,2), 5)), heights = c(rep(1,6)))

layout(layout.mat2,
       widths = c(2, rep(c(0.5,2,2,2), 5),0.5), heights = c(rep(1,6)))

### L = na, D = 35-15C ###
par(mar=c(0.3,0,0,0), mgp=c(1.5, 0.2, 0))

for (j in c(35,30,25,20,15)) {
  env.param <- env.settings[env.settings$L==0 & env.settings$D.Tmp==j, "Env_ID"]
  plot_pca_base_env(env_param = env.param, hide_axis_labels = T,hide_yaxis = F)
}

### L and D ###
for (j in c(35,30,25,20,15)) {  # dark temp
  for (k in c(20,25,30,35,40)) { # light temp
    
    if(j == k) {
      par(mar=c(0,0,0,0))
      plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1))
      
    } else {
      for (l in c(8,12,16)) {
        par(mar=c(0.3,0,0,0))
        env.param <- env.settings[env.settings$L.Tmp==k & env.settings$D.Tmp==j & env.settings$L==l,"Env_ID"]
        if(l==8){
        plot_pca_base_env(env_param = env.param, hide_axis_labels = T ,hide_yaxis = F)
        }else{
        plot_pca_base_env(env_param = env.param, hide_axis_labels = T ,hide_yaxis = T)
        }
      }
    }
  }
}

### all light ###
for (j in c(20,25,30,35,40)) {
  par(mar=c(0.3,0,0,0))
  env.param <- env.settings[env.settings$L.Tmp==j & env.settings$D==0,"Env_ID"]
  plot_pca_base_env(env_param = env.param, hide_axis_labels = T, hide_yaxis = F)
}

# label
mtext("Dark temperature (??C)", side=2, outer=TRUE, line=1.8, cex=0.5)
mtext("Light temperature (??C)", side=3, outer=TRUE, line=1.8, cex=0.5)

mtext(c("Continuous\nlight","15", "20", "25", "30", "35"), side=2, outer=TRUE, line=0.5, at=seq(1/12, 11/12, length.out=6), cex=0.5, las=3)
mtext(c("Continuous\ndark","20", "25", "30", "35", "40"), side=3, outer=TRUE, line=0.8, at=c(1/35,seq(5.5/35, 31.5/35, length.out=5)), cex=0.5)
mtext(c("8L", "12L", "16L", "8L", "12L", "16L","8L", "12L", "16L","8L", "12L", "16L"),
      side=3, outer=TRUE, line=0, at=c(3.5/35,5.5/35,7.5/35,
                                       10/35,12/35,14/35,
                                       16.5/35, 18.5/35,20.5/35, 
                                       23/35,25/35,27/35,
                                       29.5/35,31.5/35,33.5/35
      ),
      cex=0.5)

dev.off()

### Fig. S3 t-SNE ----
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

tol5 <- c(
  "#332288",  # Navy
  "#88CCEE",  # Sky Blue
  "#117733",  # Dark Green
  "#DDCC77",  # Sand Yellow
  "#882255"   # Burgundy
)

gg.pca.daylength <- ggplot(data=d,aes(x=d[,1], y=d[,2],
                                      col=day_length, fill=day_length_LD,
                                      shape=Cultivar))+
  geom_point(size=0.8, stroke =0.6) +
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=tol5)+
  scale_fill_manual(values=c(tol5[1],"transparent",tol5[2],"transparent",tol5[3],
                             "transparent",tol5[4],tol5[5]))+
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

color8 <- viridis::viridis(8)

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
ggsave(p, file=sprintf("../output/Fig.S3_nolegend.pdf"), 
       height = 130, width = 46, units = "mm")

pp <- plot_grid(gg.pca.daylength)
ggsave(pp, file=sprintf("../output/Fig.S3a_legend.pdf"), 
       height = 60, width = 75, units = "mm")

pp <- plot_grid(gg.pca.samplingtime)
ggsave(pp, file=sprintf("../output/Fig.S3b_legend.pdf"), 
       height = 105, width = 90, units = "mm")

pp <- plot_grid(gg.pca.temp)
ggsave(pp, file=sprintf("../output/Fig.S3c_legend.pdf"), 
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

### Kruskal-Wallis test and Steel-Dwass test function ----
create_p_value_matrix <- function(data, ordered_groups) {
  # Split 'comparison' column into two groups
  split_pairs <- strsplit(as.character(data$comparison), " - ")
  data$group1 <- sapply(split_pairs, `[`, 1)
  data$group2 <- sapply(split_pairs, `[`, 2)
  data$comparison <- NULL
  data$p_value <- as.numeric(data$p_value)
  
  # Duplicate rows in reverse order to ensure matrix symmetry
  reversed_data <- data
  reversed_data$group1 <- data$group2
  reversed_data$group2 <- data$group1
  full_data <- rbind(data, reversed_data)
  
  # Create all possible group combinations
  expanded <- expand.grid(group1 = ordered_groups,
                          group2 = ordered_groups,
                          stringsAsFactors = FALSE)
  
  # Merge pairwise comparisons with all combinations
  merged_data <- merge(expanded, full_data, by = c("group1", "group2"), all.x = TRUE)
  
  # Convert to matrix format
  pmat <- acast(merged_data, group1 ~ group2, value.var = "p_value", drop = FALSE)
  
  # Fill lower triangle to ensure symmetry
  pmat[lower.tri(pmat)] <- t(pmat)[lower.tri(pmat)]
  diag(pmat) <- NA  # Set diagonal to NA
  
  return(pmat)
}

assign_significance_letters <- function(data, pairwise_p, group_cols = c("env", "geno"), value_col = "value") {
  # Create group label by combining specified grouping columns
  data <- data %>%
    mutate(group_label = do.call(paste, c(across(all_of(group_cols)), sep = "_")))
  
  # Calculate median value per group
  label <- data %>%
    group_by(group_label) %>%
    summarize(median_value = median(.data[[value_col]]), .groups = "drop") %>%
    arrange(desc(median_value))
  
  # Extract ordered group labels
  ordered_labels <- label$group_label
  
  # Create p-value matrix
  p_matrix <- create_p_value_matrix(pairwise_p, ordered_labels)
  
  # Reorder matrix based on descending mean values
  p_matrix_ordered <- p_matrix[ordered_labels, ordered_labels]
  
  # Assign significance letters
  letters_result <- multcompLetters(p_matrix_ordered)
  
  # Merge letter assignments
  label$letter <- letters_result$Letters[label$group_label]
  
  # Return result table
  return(label)
}

### amp DL ----
data <- data.frame(env=rep(env.settings$L,2), 
                   value=c(kos.amp.mean2, tak.amp.mean2),
                   geno=rep(c("K","T"),each=73))

data$env_geno <- factor(paste(data$env,data$geno ,sep = "_"))

kruskal.test(value ~ env_geno, data = data)
set.seed(1)
c <- pSDCFlig(x = data$value, g = data$env_geno, method = "Monte Carlo")
#save(c, file = "../output/Steel-Dwass_DL")

# Load Steel-Dwass result data
#load("Steel-Dwass_DL")

# Create pairwise comparison table
pairwise_p <- data.frame(
  comparison = c$labels,
  p_value = c$p.val
)

result_table <- assign_significance_letters(
  data = data,
  pairwise_p = pairwise_p,
  group_cols = c("env", "geno"), 
  value_col = "value"
)

print(result_table)
write.table(result_table, "../output/Amp_DL_Steel_Dwass.txt", sep="\t",row.names = F)

### amp Ltemp ----
data <- data.frame(env=rep(env.settings$L.Tmp,2), 
                   value=c(kos.amp.mean2, tak.amp.mean2),
                   geno=rep(c("K","T"),each=73))[c(6:73,79:146),]

data$env_geno <- factor(paste(data$env,data$geno ,sep = "_"))

kruskal.test(value ~ env_geno, data = data)
set.seed(1)
#a <- pSDCFlig(x = data$value, g = data$env_geno, method = "Monte Carlo")
#save(a, file="../output/Steel-Dwass_Ltemp")
#load("../output/Steel-Dwass_Ltemp")

# Create pairwise comparison table
pairwise_p <- data.frame(
  comparison = a$labels,
  p_value = a$p.val
)

result_table <- assign_significance_letters(
  data = data,
  pairwise_p = pairwise_p,
  group_cols = c("env", "geno"), 
  value_col = "value"
)

print(result_table)
write.table(result_table, "../output/Amp_Ltemp_Steel_Dwass.txt", sep="\t",row.names = F)

### amp Dtemp ----
data <- data.frame(env=rep(env.settings$D.Tmp,2), 
                   value=c(kos.amp.mean2, tak.amp.mean2),
                   geno=rep(c("K","T"),each=73))[c(1:68,74:141),]

data$env_geno <- factor(paste(data$env,data$geno ,sep = "_"))

kruskal.test(value ~ env_geno, data = data)
set.seed(1)
#b <- pSDCFlig(x = data$value, g = data$env_geno, method = "Monte Carlo")
#save(b, file = "../output/Steel-Dwass_Dtemp")
#load("../output/Steel-Dwass_Dtemp")
# Create pairwise comparison table
pairwise_p <- data.frame(
  comparison = b$labels,
  p_value = b$p.val
)

result_table <- assign_significance_letters(
  data = data,
  pairwise_p = pairwise_p,
  group_cols = c("env", "geno"), 
  value_col = "value"
)

print(result_table)
write.table(result_table, "../output/Amp_Dtemp_Steel_Dwass.txt", sep="\t",row.names = F)

### amplitude violin plot 73 conditions -----
pdf("../output/Fig.S4.pdf", width=6.5, height=4)

par(oma=c(1,3,3,0))
#layout(layout.mat2,
#       widths = c(2, rep(c(2,2,2), 5)), heights = c(rep(1,6)))

layout(layout.mat2,
       widths = c(2, rep(c(0.5, 2,2,2), 5),0.5), heights = c(rep(1,6)))

### L = na, D = 35-15C ###
par(mar=c(0.4,0,0,0), mgp=c(1.5, 0.2, 0))

for (j in c(35,30,25,20,15)) {
  env.param <- env.settings[env.settings$L==0 & env.settings$D.Tmp==j, "Env_ID"]
  vioplot(kos.amp2[env.param,], tak.amp2[env.param,],col= color5[1:2], border="black",
          xlab="", ylab="",  xaxt = "n", yaxt = "n",ylim=c(0,15))
  axis(2, at = c(0, 5, 10, 15), labels = FALSE, tck = -0.1)
}

### L and D ###
for (j in c(35,30,25,20,15)) {  # dark temp
  for (k in c(20,25,30,35,40)) { # light temp
    
    if(j == k) {
      par(mar=c(0,0,0,0))
      plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1))
      
    } else {
      for (l in c(8,12,16)) {
        par(mar=c(0.4,0,0,0))
        env.param <- env.settings[env.settings$L.Tmp==k & env.settings$D.Tmp==j & env.settings$L==l,"Env_ID"]
        vioplot(kos.amp2[env.param,], tak.amp2[env.param,], col= color5[1:2], border="black",
                xlab="", ylab="", ylim=c(0,15),
                xaxt = "n", yaxt = "n")
        if(l==8){axis(2, at = c(0, 5, 10, 15), labels = FALSE, tck = -0.1)}

      }
    }
  }
}

### all light ###
for (j in c(20,25,30,35,40)) {
  par(mar=c(0.4,0,0,0))
  env.param <- env.settings[env.settings$L.Tmp==j & env.settings$D==0,"Env_ID"]
  vioplot(kos.amp2[env.param,], tak.amp2[env.param,],col= color5[1:2], border="black",
          xlab="", ylab="",
          xaxt = "n", yaxt = "n",ylim=c(0,15))
  axis(2, at = c(0, 5, 10, 15),labels = FALSE, tck = -0.1)
}

# label
mtext("Dark temperature (??C)", side=2, outer=TRUE, line=1.8, cex=0.5)
mtext("Light temperature (??C)", side=3, outer=TRUE, line=1.8, cex=0.5)

mtext(c("Continuous\nlight","15", "20", "25", "30", "35"), side=2, outer=TRUE, line=0.5, at=seq(1/12, 11/12, length.out=6), cex=0.5, las=3)
mtext(c("Continuous\ndark","20", "25", "30", "35", "40"), side=3, outer=TRUE, line=0.8, at=c(1/35,seq(5.5/35, 31.5/35, length.out=5)), cex=0.5)
mtext(c("8L", "12L", "16L", "8L", "12L", "16L","8L", "12L", "16L","8L", "12L", "16L"),
      side=3, outer=TRUE, line=0, at=c(3.5/35,5.5/35,7.5/35,
                                       10/35,12/35,14/35,
                                       16.5/35, 18.5/35,20.5/35, 
                                       23/35,25/35,27/35,
                                       29.5/35,31.5/35,33.5/35
      ),
      cex=0.5)

dev.off()

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

### t-test Kruskal-Wallis test and Steel-Dwass test ----
### t-test DL ----
data <- data.frame(env=rep(env.settings$L,2), 
                   value=c(ttest.signif.kos, ttest.signif.tak),
                   geno=rep(c("K","T"),each=73))

data$env_geno <- factor(paste(data$env,data$geno ,sep = "_"))

kruskal.test(value ~ env_geno, data = data)
set.seed(1)
#c <- pSDCFlig(x = data$value, g = data$env_geno, method = "Monte Carlo")
#save(c, file = "../output/Steel-Dwass_DL_t-test")

# Load Steel-Dwass result data
#load("../output/Steel-Dwass_DL_t-test")

# Create pairwise comparison table
pairwise_p <- data.frame(
  comparison = c$labels,
  p_value = c$p.val
)

result_table <- assign_significance_letters(
  data = data,
  pairwise_p = pairwise_p,
  group_cols = c("env", "geno"), 
  value_col = "value"
)

print(result_table)
write.table(result_table, "../output/t-test_DL_Steel_Dwass.txt", sep="\t",row.names = F)

### t-test Ltemp ----
data <- data.frame(env=rep(env.settings$L.Tmp,2), 
                   value=c(ttest.signif.kos, ttest.signif.tak),
                   geno=rep(c("K","T"),each=73))[c(6:73,79:146),]

data$env_geno <- factor(paste(data$env,data$geno ,sep = "_"))

kruskal.test(value ~ env_geno, data = data)
set.seed(1)
#a <- pSDCFlig(x = data$value, g = data$env_geno, method = "Monte Carlo")
#save(a, file="Steel-Dwass_Ltemp_t-test")

# Load Steel-Dwass result data
load("../output/Steel-Dwass_Ltemp_t-test")

# Create pairwise comparison table
pairwise_p <- data.frame(
  comparison = a$labels,
  p_value = a$p.val
)

result_table <- assign_significance_letters(
  data = data,
  pairwise_p = pairwise_p,
  group_cols = c("env", "geno"), 
  value_col = "value"
)

print(result_table)
write.table(result_table, "../output/t-test_Ltemp_Steel_Dwass.txt", sep="\t",row.names = F)

### t-test Dtemp ----
data <- data.frame(env=rep(env.settings$D.Tmp,2), 
                   value=c(ttest.signif.kos, ttest.signif.tak),
                   geno=rep(c("K","T"),each=73))[c(1:68,74:141),]

data$env_geno <- factor(paste(data$env,data$geno ,sep = "_"))

kruskal.test(value ~ env_geno, data = data)
set.seed(1)
#b <- pSDCFlig(x = data$value, g = data$env_geno, method = "Monte Carlo")
#save(b, file = "Steel-Dwass_Dtemp_t-test")

# Load Steel-Dwass result data
load("../output/Steel-Dwass_Dtemp_t-test")

# Create pairwise comparison table
pairwise_p <- data.frame(
  comparison = b$labels,
  p_value = b$p.val
)

result_table <- assign_significance_letters(
  data = data,
  pairwise_p = pairwise_p,
  group_cols = c("env", "geno"), 
  value_col = "value"
)

print(result_table)
write.table(result_table, "../output/t-test_Dtemp_Steel_Dwass.txt", sep="\t",row.names = F)

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

### Fig.2c-h amplitude and t-test ----
p <- plot_grid(g_amp_DL, g_amp_L,g_amp_D, 
               g_amp_DL_ttest,g_amp_L_ttest,g_amp_D_ttest,
               ncol=3, align = "hv")

ggsave(p, file=sprintf("../output/Fig.2c-h.pdf"),  
       height = 70, width = 170, units = "mm")

### Table.S6-7, Fig.3h correlation with air temperature ----
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

### representative 500 genes ---- 
# The number of genes belong to the cluster of T.to.R.genes.80_kos and T.to.R.genes.80_tak genes
load("../input/gene.cluster")

# exemplars["Os01g0182600"] #OsGI
# exemplars[4227] #Os04g0595100 
# names(exemplars[exemplars==4227])
# names(exemplars)[unique(exemplars)]

id_converter <- read.delim("../input/230728_all_genes_id_conversion.tsv")
rownames(id_converter) <- id_converter$microarray
id_converter_500 <- id_converter[names(exemplars)[unique(exemplars)],]
aa <- exemplars[names(exemplars) %in% rownames(id_converter_500)] 
id_converter_500$number <- aa[rownames(id_converter_500)]
representative_474_rnaseq <- id_converter_500[id_converter_500$flag!=3,"number"]
names(representative_474_rnaseq) <- id_converter_500[id_converter_500$flag!=3,"rnaseq"]

exemplars_new <- exemplars
names(exemplars_new) <- id_converter[names(exemplars),"rnaseq"]

exemplars_new <- exemplars_new[exemplars_new %in% as.numeric(representative_474_rnaseq)]
exemplars_new <- exemplars_new[!is.na(names(exemplars_new))] # 15907 genes

exemplars_new %in% as.numeric(representative_474_rnaseq[T.to.R.genes.80_Kos_rnaseq]) %>% sum # 2887
exemplars_new %in% as.numeric(representative_474_rnaseq[T.to.R.genes.80_Tak_rnaseq]) %>% sum # 3142
exemplars_new %in% as.numeric(representative_474_rnaseq[intersect(T.to.R.genes.80_Kos_rnaseq,T.to.R.genes.80_Tak_rnaseq)]) %>% sum # 1167

group_474 <- list()
for(i in 1:474){
  id <- (exemplars_new %>% table %>% names)[i]
  group_474[[i]] <- names(exemplars_new)[exemplars_new %>% as.character() %in% id]
}

### choice of predictor ----
load("../input/prob.env.choice.fld_Kos")
load("../input/prob.env.choice.GC_Kos")
load("../input/prob.env.choice.mix_Kos")
load("../input/prob.env.choice.fld_Tak")
load("../input/prob.env.choice.GC_Tak")
load("../input/prob.env.choice.mix_Tak")

### radiation > temperature
## kos
# Field
(prob.env.choice.fld_Kos[2,]-prob.env.choice.fld_Kos[1,] > 0) %>% sum #128
# GC
(prob.env.choice.GC_Kos[2,]-prob.env.choice.GC_Kos[1,] > 0) %>% sum #308
# mix
(prob.env.choice.mix_Kos[2,]-prob.env.choice.mix_Kos[1,] > 0) %>% sum #378

## Tak
# Field
(prob.env.choice.fld_Tak[2,]-prob.env.choice.fld_Tak[1,] > 0) %>% sum #86
# GC
(prob.env.choice.GC_Tak[2,]-prob.env.choice.GC_Tak[1,] > 0) %>% sum #310
# mix
(prob.env.choice.mix_Tak[2,]-prob.env.choice.mix_Tak[1,] > 0) %>% sum #359

### Table S8 radiation or temperature kos ----
fld <- rep(NA,466)
names(fld) <- colnames(prob.env.choice.fld_Kos)

for(i in 1:466){
  tmp <- prob.env.choice.fld_Kos[2,i]-prob.env.choice.fld_Kos[1,i]
  
  if(tmp > 0){
    fld[i] <- "radiation"
  }else{
    if(tmp < 0){
    fld[i] <- "temperature"
  }else{
    fld[i] <- "radiation = temperature"
  }
  }  
}

GC <- rep(NA,466)
names(GC) <- colnames(prob.env.choice.GC_Kos)

for(i in 1:466){
  tmp <- prob.env.choice.GC_Kos[2,i]-prob.env.choice.GC_Kos[1,i]
  
  if(tmp > 0){
    GC[i] <- "radiation"
  }else{
    if(tmp < 0){
      GC[i] <- "temperature"
    }else{
      GC[i] <- "radiation = temperature"
    }
  }  
}

mix <- rep(NA,466)
names(mix) <- colnames(prob.env.choice.mix_Kos)

for(i in 1:466){
  tmp <- prob.env.choice.mix_Kos[2,i]-prob.env.choice.mix_Kos[1,i]
  
  if(tmp > 0){
    mix[i] <- "radiation"
  }else{
    if(tmp < 0){
      mix[i] <- "temperature"
    }else{
      mix[i] <- "radiation = temperature"
    }
  }  
}

predictor <- data.frame(fld,GC,mix,des[names(fld),c(2,3,8)])
predictor <- predictor %>% arrange(Locus_ID)

write.table(predictor,file="../output/TableS8.txt",sep="\t",row.names=F)

### Table S10radiation or temperature  Tak ----
fld <- rep(NA,456)
names(fld) <- colnames(prob.env.choice.fld_Tak)

for(i in 1:456){
  tmp <- prob.env.choice.fld_Tak[2,i]-prob.env.choice.fld_Tak[1,i]
  
  if(tmp > 0){
    fld[i] <- "radiation"
  }else{
    if(tmp < 0){
      fld[i] <- "temperature"
    }else{
      fld[i] <- "radiation = temperature"
    }
  }  
}

GC <- rep(NA,456)
names(GC) <- colnames(prob.env.choice.GC_Tak)

for(i in 1:456){
  tmp <- prob.env.choice.GC_Tak[2,i]-prob.env.choice.GC_Tak[1,i]
  
  if(tmp > 0){
    GC[i] <- "radiation"
  }else{
    if(tmp < 0){
      GC[i] <- "temperature"
    }else{
      GC[i] <- "radiation = temperature"
    }
  }  
}

mix <- rep(NA,456)
names(mix) <- colnames(prob.env.choice.mix_Tak)

for(i in 1:456){
  tmp <- prob.env.choice.mix_Tak[2,i]-prob.env.choice.mix_Tak[1,i]
  
  if(tmp > 0){
    mix[i] <- "radiation"
  }else{
    if(tmp < 0){
      mix[i] <- "temperature"
    }else{
      mix[i] <- "radiation = temperature"
    }
  }  
}

predictor <- data.frame(fld,GC,mix,des[names(fld),c(2,3,8)])
predictor <- predictor %>% arrange(Locus_ID)

write.table(predictor,file="../output/TableS10.txt",sep="\t",row.names=F)

### Table S9 and S11 GO and KEGG radiation > temperature ----
pp <- function(genelist){names(exemplars_new)[exemplars_new %in% as.numeric(representative_474_rnaseq[genelist])]}

fld_kos_rep_gene <- colnames(prob.env.choice.fld_Kos)[(prob.env.choice.fld_Kos[2,]-prob.env.choice.fld_Kos[1,] > 0)]
fld_kos_rep_gene_all <- pp(fld_kos_rep_gene)

GC_kos_rep_gene <- colnames(prob.env.choice.GC_Kos)[(prob.env.choice.GC_Kos[2,]-prob.env.choice.GC_Kos[1,] > 0)]
GC_kos_rep_gene_all <- pp(GC_kos_rep_gene)

mix_kos_rep_gene <- colnames(prob.env.choice.mix_Kos)[(prob.env.choice.mix_Kos[2,]-prob.env.choice.mix_Kos[1,] > 0)]
mix_kos_rep_gene_all <- pp(mix_kos_rep_gene)

fld_tak_rep_gene <- colnames(prob.env.choice.fld_Tak)[(prob.env.choice.fld_Tak[2,]-prob.env.choice.fld_Tak[1,] > 0)]
fld_tak_rep_gene_all <- pp(fld_tak_rep_gene)

GC_tak_rep_gene <- colnames(prob.env.choice.GC_Tak)[(prob.env.choice.GC_Tak[2,]-prob.env.choice.GC_Tak[1,] > 0)]
GC_tak_rep_gene_all <- pp(GC_tak_rep_gene)

mix_tak_rep_gene <- colnames(prob.env.choice.mix_Tak)[(prob.env.choice.mix_Tak[2,]-prob.env.choice.mix_Tak[1,] > 0)]
mix_tak_rep_gene_all <- pp(mix_tak_rep_gene)

dir.create("../output/GO")

# GO
use.genelist <- rownames(dat)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))

glist <- list(fld_kos_rep_gene_all=fld_kos_rep_gene_all,
              GC_kos_rep_gene_all=GC_kos_rep_gene_all,
              mix_kos_rep_gene_all=mix_kos_rep_gene_all,
              fld_tak_rep_gene_all=fld_tak_rep_gene_all,
              GC_tak_rep_gene_all=GC_tak_rep_gene_all,
              mix_tak_rep_gene_all=mix_tak_rep_gene_all)

for (i in 1:6){

  fn <- sprintf("../output/GO/TableS9_11_GO_%s_radiation.csv",names(glist[i]))
  gl <-  glist[[i]] 
  
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  adp <- p.adjust(result[,"p.value"], method = "BH")
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)

}

dir.create("../output/KEGG")

# KEGG
use.genelist <- rownames(dat)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]

for (i in 1:6){
  
  fn <- sprintf("../output/KEGG/TableS9_11_KEGG_%s_radiation.csv",names(glist[i]))
  gl <-  glist[[i]] 
  
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)

}

### Table S9 and S11 GO and KEGG radiation < temperature ----
pp <- function(genelist){names(exemplars_new)[exemplars_new %in% as.numeric(representative_474_rnaseq[genelist])]}

fld_kos_rep_gene <- colnames(prob.env.choice.fld_Kos)[(prob.env.choice.fld_Kos[2,]-prob.env.choice.fld_Kos[1,] < 0)]
fld_kos_rep_gene_all <- pp(fld_kos_rep_gene)

GC_kos_rep_gene <- colnames(prob.env.choice.GC_Kos)[(prob.env.choice.GC_Kos[2,]-prob.env.choice.GC_Kos[1,] < 0)]
GC_kos_rep_gene_all <- pp(GC_kos_rep_gene)

mix_kos_rep_gene <- colnames(prob.env.choice.mix_Kos)[(prob.env.choice.mix_Kos[2,]-prob.env.choice.mix_Kos[1,] < 0)]
mix_kos_rep_gene_all <- pp(mix_kos_rep_gene)

fld_tak_rep_gene <- colnames(prob.env.choice.fld_Tak)[(prob.env.choice.fld_Tak[2,]-prob.env.choice.fld_Tak[1,] < 0)]
fld_tak_rep_gene_all <- pp(fld_tak_rep_gene)

GC_tak_rep_gene <- colnames(prob.env.choice.GC_Tak)[(prob.env.choice.GC_Tak[2,]-prob.env.choice.GC_Tak[1,] < 0)]
GC_tak_rep_gene_all <- pp(GC_tak_rep_gene)

mix_tak_rep_gene <- colnames(prob.env.choice.mix_Tak)[(prob.env.choice.mix_Tak[2,]-prob.env.choice.mix_Tak[1,] < 0)]
mix_tak_rep_gene_all <- pp(mix_tak_rep_gene)

# GO
use.genelist <- rownames(dat)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))

glist <- list(fld_kos_rep_gene_all=fld_kos_rep_gene_all,
              GC_kos_rep_gene_all=GC_kos_rep_gene_all,
              mix_kos_rep_gene_all=mix_kos_rep_gene_all,
              fld_tak_rep_gene_all=fld_tak_rep_gene_all,
              GC_tak_rep_gene_all=GC_tak_rep_gene_all,
              mix_tak_rep_gene_all=mix_tak_rep_gene_all)

for (i in 1:6){
  
  fn <- sprintf("../output/GO/TableS9_11_GO_%s_temperature.csv",names(glist[i]))
  gl <-  glist[[i]] 
  
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  adp <- p.adjust(result[,"p.value"], method = "BH")
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

# KEGG
use.genelist <- rownames(dat)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]

for (i in 1:6){
  
  fn <- sprintf("../output/KEGG/TableS9_11_KEGG_%s_temperature.csv",names(glist[i]))
  gl <-  glist[[i]] 
  
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

### Table S12 from temperature to radiation ----
T.to.R.genes.80_Kos <- intersect(colnames(prob.env.choice.fld_Kos)[prob.env.choice.fld_Kos[1,] >= 0.8],
                                 colnames(prob.env.choice.mix_Kos)[prob.env.choice.mix_Kos[2,] >= 0.8])
T.to.R.genes.80_Tak <- intersect(colnames(prob.env.choice.fld_Tak)[prob.env.choice.fld_Tak[1,] >= 0.8],
                                 colnames(prob.env.choice.mix_Tak)[prob.env.choice.mix_Tak[2,] >= 0.8])

des_genelist_T.to.R.genes.80_kos <- des[setdiff(T.to.R.genes.80_Kos,T.to.R.genes.80_Tak),c(2:3,8)]
des_genelist_T.to.R.genes.80_tak <- des[setdiff(T.to.R.genes.80_Tak,T.to.R.genes.80_Kos),c(2:3,8)]
des_genelist_T.to.R.genes.80_kos_tak <- des[intersect(T.to.R.genes.80_Kos,T.to.R.genes.80_Tak),c(2:3,8)]

write.table(des_genelist_T.to.R.genes.80_kos,"../output/Table.S12_kos.tsv",
            sep="\t", row.names = F)
write.table(des_genelist_T.to.R.genes.80_tak,"../output/Table.S12_tak.tsv",
            sep="\t", row.names = F)
write.table(des_genelist_T.to.R.genes.80_kos_tak,"../output/Table.S12_kos_tak.tsv",
            sep="\t", row.names = F)

### Table S13 from radiation to radiation ----
R.to.R.genes.80_Kos <- intersect(colnames(prob.env.choice.fld_Kos)[prob.env.choice.fld_Kos[2,] >= 0.8],
                                 colnames(prob.env.choice.mix_Kos)[prob.env.choice.mix_Kos[2,] >= 0.8]) #18
R.to.R.genes.80_Tak <- intersect(colnames(prob.env.choice.fld_Tak)[prob.env.choice.fld_Tak[2,] >= 0.8],
                                 colnames(prob.env.choice.mix_Tak)[prob.env.choice.mix_Tak[2,] >= 0.8]) #11

des_genelist_R.to.R.genes.80_kos <- des[setdiff(R.to.R.genes.80_Kos,R.to.R.genes.80_Tak),c(2:3,8)] #14
des_genelist_R.to.R.genes.80_tak <- des[setdiff(R.to.R.genes.80_Tak,R.to.R.genes.80_Kos),c(2:3,8)] #7
des_genelist_R.to.R.genes.80_kos_tak <- des[intersect(R.to.R.genes.80_Kos,R.to.R.genes.80_Tak),c(2:3,8)] #4

write.table(des_genelist_R.to.R.genes.80_kos,"../output/Table.S13_kos.tsv",
            sep="\t", row.names = F)
write.table(des_genelist_R.to.R.genes.80_tak,"../output/Table.S13_tak.tsv",
            sep="\t", row.names = F)
write.table(des_genelist_R.to.R.genes.80_kos_tak,"../output/Table.S13_kos_tak.tsv",
            sep="\t", row.names = F)

### Table S14 from temperature to temperature ----
T.to.T.genes.80_Kos <- intersect(colnames(prob.env.choice.fld_Kos)[prob.env.choice.fld_Kos[1,] >= 0.8],
                                 colnames(prob.env.choice.mix_Kos)[prob.env.choice.mix_Kos[1,] >= 0.8]) #4
T.to.T.genes.80_Tak <- intersect(colnames(prob.env.choice.fld_Tak)[prob.env.choice.fld_Tak[1,] >= 0.8],
                                 colnames(prob.env.choice.mix_Tak)[prob.env.choice.mix_Tak[1,] >= 0.8]) #19

des_genelist_T.to.T.genes.80_kos <- des[setdiff(T.to.T.genes.80_Kos,T.to.T.genes.80_Tak),c(2:3,8)] #3
des_genelist_T.to.T.genes.80_tak <- des[setdiff(T.to.T.genes.80_Tak,T.to.T.genes.80_Kos),c(2:3,8)] #18
des_genelist_T.to.T.genes.80_kos_tak <- des[intersect(T.to.T.genes.80_Kos,T.to.T.genes.80_Tak),c(2:3,8)] #1

write.table(des_genelist_T.to.T.genes.80_kos,"../output/Table.S14_kos.tsv",
            sep="\t", row.names = F)
write.table(des_genelist_T.to.T.genes.80_tak,"../output/Table.S14_tak.tsv",
            sep="\t", row.names = F)
write.table(des_genelist_T.to.T.genes.80_kos_tak,"../output/Table.S14_kos_tak.tsv",
            sep="\t", row.names = F)

### Table S15 T to R genes GO and KEGG pathway analysis ----
use.genelist <- rownames(dat)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))

glist <- list(
  Kos = pp(T.to.R.genes.80_Kos),
  Tak = pp(T.to.R.genes.80_Tak),
  Kos_Tak = pp(intersect(T.to.R.genes.80_Kos,T.to.R.genes.80_Tak))
)

for(i in 1:3){
  
  gl <-  glist[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("../output/GO/TableS15_T.to.R.genes.80_%s_GO.csv",names(glist)[i]) 
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

### KEGG
use.genelist <- rownames(dat)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]

for (i in 1:3){
  
  gl <- glist[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("../output/KEGG/TableS15_T.to.R.genes.80_%s_KEGG.csv", names(glist)[i]) 
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

### Table S16 R to R genes GO and KEGG pathway analysis ----
use.genelist <- rownames(dat)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))

glist <- list(
  Kos = pp(R.to.R.genes.80_Kos),
  Tak = pp(R.to.R.genes.80_Tak),
  Kos_Tak = pp(intersect(pp(R.to.R.genes.80_Kos),pp(R.to.R.genes.80_Tak)))
)

for(i in 1:3){
  
  gl <-  glist[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("../output/GO/TableS16_R.to.R.genes.80_%s_GO.csv",names(glist)[i]) 
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

### KEGG
use.genelist <- rownames(dat)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]

for (i in 1:3){
  
  gl <- glist[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("../output/KEGG/TableS16_R.to.R.genes.80_%s_KEGG.csv", names(glist)[i]) 
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

### Table S17 T to T genes GO and KEGG pathway analysis ----
use.genelist <- rownames(dat)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))

glist <- list(
  Kos = pp(T.to.T.genes.80_Kos),
  Tak = pp(T.to.T.genes.80_Tak),
  Kos_Tak = pp(intersect(pp(T.to.T.genes.80_Kos),pp(T.to.T.genes.80_Tak)))
)

for(i in 1:3){
  
  gl <-  glist[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("../output/GO/TableS17_T.to.T.genes.80_%s_GO.csv",names(glist)[i]) 
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

### KEGG
use.genelist <- rownames(dat)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]

for (i in 1:3){
  
  gl <- glist[[i]]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("../output/KEGG/TableS17_T.to.T.genes.80_%s_KEGG.csv", names(glist)[i]) 
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

### Table S18 and 19 high performance in GC than field or in mix than field GO analysis Koshihikari ----
load("../input/pred.K.fld.med")
load("../input/pred.K.GC.med")
load("../input/pred.K.mix.med") 

GC_than_fld_0 <- exemplars_new[names(pred.K.fld.med)[(pred.K.fld.med - pred.K.GC.med > 0 )]] #202
GC_than_fld_0_all <- names(exemplars_new)[exemplars_new %in% as.numeric(GC_than_fld_0)] #7059

mix_than_fld_0 <- exemplars_new[names(pred.K.fld.med)[(pred.K.fld.med - pred.K.mix.med > 0 )]] #360
mix_than_fld_0_all <- names(exemplars_new)[exemplars_new %in% as.numeric(mix_than_fld_0)] #12109



GC_mix_than_fld_0 <- intersect(names(GC_than_fld_0),names(mix_than_fld_0)) #192 genes

des_GC_mix_than_fld_0 <- des[GC_mix_than_fld_0, c(2,3,8)] #192
des_GC_than_fld_0 <- des[setdiff(names(GC_than_fld_0),GC_mix_than_fld_0), c(2,3,8)] # 10 genes
des_mix_than_fld_0 <- des[setdiff(names(mix_than_fld_0),GC_mix_than_fld_0), c(2,3,8)] # 168 genes

write.table(des_GC_mix_than_fld_0,file="../output/GO/TableS18_mix_GC_than_fld_genes.txt",sep="\t",row.names = F) 
write.table(des_GC_than_fld_0,file="../output/GO/TableS18_GC_than_fld_genes.txt",sep="\t",row.names = F) 
write.table(des_mix_than_fld_0,file="../output/GO/TableS18_mix_than_GC_genes.txt",sep="\t", row.names = F) 

## GO
glist <- list(GC_than_fld = GC_minus_fld_0_all,
              mix_than_fld = mix_than_fld_0_all)

for(i in 1:2){
  
  gl <-  glist[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("../output/GO/TableS19_high_performance_in_%s.csv",names(glist)[i])

  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

## KEGG
use.genelist <- rownames(dat)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]

for (i in 1:2){
  
  fn <- sprintf("../output/KEGG/TableS19_high_performance_in_%s.csv",names(glist[i]))
  gl <-  glist[[i]] 
  
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

### Table S20 and 21 high performance in GC than field or in mix than field GO analysis Takanari ----
load("../input/pred.T.fld.med")
load("../input/pred.T.GC.med")
load("../input/pred.T.mix.med") 

GC_than_fld_0 <- exemplars_new[names(pred.T.fld.med)[(pred.T.fld.med - pred.T.GC.med > 0 )]] #179
GC_than_fld_0_all <- names(exemplars_new)[exemplars_new %in% as.numeric(GC_than_fld_0)] #6250

mix_than_fld_0 <- exemplars_new[names(pred.T.fld.med)[(pred.T.fld.med - pred.T.mix.med > 0 )]] #335
mix_than_fld_0_all <- names(exemplars_new)[exemplars_new %in% as.numeric(mix_than_fld_0)] #11215


GC_mix_than_fld_0 <- intersect(names(GC_than_fld_0),names(mix_than_fld_0)) #168 genes

des_GC_mix_than_fld_0 <- des[GC_mix_than_fld_0, c(2,3,8)] #168
des_GC_than_fld_0 <- des[setdiff(names(GC_than_fld_0),GC_mix_than_fld_0), c(2,3,8)] # 11 genes
des_mix_than_fld_0 <- des[setdiff(names(mix_than_fld_0),GC_mix_than_fld_0), c(2,3,8)] # 167 genes

write.table(des_GC_mix_than_fld_0,file="../output/GO/TableS20_Takanari_mix_GC_than_fld_genes.txt",sep="\t",row.names = F) 
write.table(des_GC_than_fld_0,file="../output/GO/TableS20_Takanari_GC_than_fld_genes.txt",sep="\t",row.names = F) 
write.table(des_mix_than_fld_0,file="../output/GO/TableS20_Takanari_mix_than_GC_genes.txt",sep="\t", row.names = F) 

## GO
glist <- list(GC_than_fld = GC_minus_fld_0_all,
              mix_than_fld = mix_than_fld_0_all)

for(i in 1:2){
  
  gl <-  glist[[i]]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("../output/GO/TableS21_Takanari_high_performance_in_%s.csv",names(glist)[i])
  
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}

## KEGG
use.genelist <- rownames(dat)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]

for (i in 1:2){
  
  fn <- sprintf("../output/KEGG/TableS21_Takanari_high_performance_in_%s.csv",names(glist[i]))
  gl <-  glist[[i]] 
  
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  adp <- p.adjust(result[,"p.value"], method = "BH")
  
  tmp.id <- rownames(result)[adp < 0.05]
  tmp.adp <- adp[adp < 0.05]
  tmp.description <- kegg_description[tmp.id]
  tmp.xnn <- result[adp < 0.05, c("xtt", "xtn", "xnt", "xnn")]
  if (length(tmp.id)==1){
    out <- t(as.matrix(c(tmp.adp, tmp.id, tmp.description, tmp.xnn)))
  }else{
    out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  }
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  write.csv(out, fn, row.names = F)
  
}
