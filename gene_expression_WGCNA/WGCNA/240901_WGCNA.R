library("tidyverse")
library("coin")
library("cowplot")
library("fields")
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/WGCNA/WGCNA_1.70-3.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
library("WGCNA")

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

### load des, at, env, setting ----
load("../input/180604_des")
load("../input/at.lab.additional")
env.settings <- read.table("../input/env.settings.txt")

### gene filtering ----
load("../input/log.rpm")

rownames(log.rpm) <- as.integer(rownames(log.rpm))
data <- t(log.rpm[,158:39765])
dat <- data[apply(data, 1, mean) > 1,] # 17742 genes
gene_17742 <- rownames(dat)

### WGCNA (lecture document by Yasunori Ichihashi 190604) ----
# https://shiokoji11235.com/wgcna-analysis-part2
### TOM ----
#adjacency = adjacency(t(dat),power=14) ## 12 conservative value for experiments over 40 samples
#TOM = TOMsimilarity(adjacency, TOMType = "signed")
#save(adjacency, file="240906_adjacency")
#save(TOM, file="240906_TOM")

load("240906_adjacency")
load("240906_TOM")

dissTOM = 1-TOM
id=rownames(dat)
colnames(TOM)=id
rownames(TOM)=id

### soft thresholding ----
# powers <- c(1:30)
# sft <- pickSoftThreshold(t(dat), powerVector=powers,
#                          networkType = "signed",
#                          RsquaredCut=0.9, verbose=5,
#                          corFnc = "bicor",
#                          corOptions = list(maxPOutliers = 0.05))
# 
# tiff("SoftThresholding.tif", width=16, height=8, unit="in",compression="lzw",res=100)
# par(mfrow = c(1,2))
# cex1 = 0.8
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.9,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()
# 
# sft$powerEstimate
# 
# # check scale free topology
# k=as.vector(apply(adjacency, 2, sum, na.rm=T))
# hist(k)
# scaleFreePlot(k, main="Check scale free topology\n")
# 
### hierarchical clustering ----
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04)

### module detection ----
deepSplit <- 4
minModuleSize <- 20

dynamicMods <- cutreeDynamic(dendro = geneTree, 
                             distM = dissTOM,
                             deepSplit = deepSplit, 
                             pamStage = TRUE,
                             pamRespectsDendro = TRUE,
                             minClusterSize = minModuleSize)
print(table(dynamicMods))

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(dendro = geneTree,
                    colors = dynamicColors,
                    groupLabels = "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

moduleColors <- dynamicColors #58+grey
colorOrder <- c("grey", standardColors(58))
moduleLabels <- match(moduleColors, colorOrder)

### merge module ----
MEList <- moduleEigengenes(t(dat), colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres <- 0.2
abline(h=MEDissThres, col = "red")

merge <- mergeCloseModules(t(dat), dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
MEs <- merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(58))
moduleLabels <- match(moduleColors, colorOrder)
n_module <- length(unique(moduleLabels))-1

### eigengenes ----
datExpr <- t(dat)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")) #kME

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

modNames = substring(names(MEs), 3)
names(geneModuleMembership) = modNames
names(MMPvalue) = paste("p.MM", modNames, sep="")

module_id <- match(modNames, colorOrder)

MES <- MEs
colnames(MES) <- module_id
MES <- MES[,order(module_id)]

tmp <- 1:39765
names(tmp) <- colnames(log.rpm)[1:39765]
d_eigengenes <- cbind(MES,log.rpm[,39766:39769])

write.table(d_eigengenes,file="d_eigengenes.txt")

### Fig.2b differences in expression levels between Koshihikari and Takanari ----
dif_MES <- rep(NA, ncol(MES))
names(dif_MES) <- colnames(MES)
for(i in 1:ncol(MES)){
  dif_MES[i] <- MES[at.lab.additional$Cultivar=="Tak",i] %>% mean - MES[at.lab.additional$Cultivar=="Kos",i] %>% mean
}

data <- data.frame(xx=dif_MES[1:22])
gg.dif.module <- ggplot(data, aes(x=xx)) +
  geom_histogram(binwidth = 0.01,colour="black",fill="white") +
  scale_y_continuous(breaks = seq(0, 12, by = 4))+
  xlim(-0.1,0.1)+
  xlab("Tak - Kos")+
  ylab("Module frequency")+
  theme_cowplot() +
  theme(plot.title = element_text(size=6, hjust=0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")  
print(gg.dif.module)

### Fig.2c correlation between expression of eigengenes and trait ----
d_traits <- data.frame(at.lab.additional$average,at.lab.additional$Tmp)
cor_trait_eigengenes <- cor(d_traits, MES[1:22]) # 22 modules

cor_trait_eigengenes_table <- t(cor_trait_eigengenes)[order(as.numeric(colnames(cor_trait_eigengenes))),]
colnames(cor_trait_eigengenes_table) <- c("measured temperature","setting temperature")

# air temperature
data <- data.frame(xx=cor_trait_eigengenes[2,])
gg.temp.module <- ggplot(data, aes(x=xx)) +
  geom_histogram(binwidth = 0.1, colour="black",fill="white") +
  xlim(-1,1)+
  xlab("r")+
  ylab("Module frequency")+
  theme_cowplot() +
  theme(plot.title = element_text(size=6, hjust=0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")  
print(gg.temp.module) #Fig.2c

### Fig.2b,c ----
gggg <- plot_grid(gg.dif.module, gg.temp.module, ncol=1, align = "hv") # Fig.2b,c
ggsave(gggg, file="../output/Fig.2bc.pdf",
       width = 45, height = 80, units = "mm")

### Fig.2d module genes ----
layout.mat2 <- matrix(c(
  1,  6,7,8,   9,10,11, 12,13,14, 15,15,15, 16,17,18,
  2,  19,20,21, 22,23,24, 25,25,25, 26,27,28, 29,30,31,
  3,  32,33,34, 35,35,35, 36,37,38, 39,40,41, 42,43,44,
  4,  45,45,45, 46,47,48, 49,50,51, 52,53,54, 55,56,57,
  5,  58,59,60, 61,62,63, 64,65,66, 67,68,69, 70,71,72,
  0,  0,0,73,   0,0,74,   0,0,75,   0,0,76,   0,0,77), byrow=T,nrow=6,ncol=16)

  the.gene <- d_eigengenes[ ,c(i+1,(n_module+2):(n_module+5))] # 1168 x 5 (RNA, EnvID, PotID, Time, Cultivar)
  
  par(oma=c(3,3,3,3))
  layout(layout.mat2,
         widths = c(2.7, rep(c(2,2,2.7), 5)), heights = rep(1,6))
  
  #==================================#
  # expression time series: all dark #
  #==================================#
  par(mar=c(1,0,0,1), mgp=c(1.5, 0.7, 0))

  y.lim.v <- c(min(d_eigengenes[,2:23]), max(d_eigengenes[,2:23]))
  #y.lim.v <- c(-0.2,0.2)

glist <- "Os01g0840100"
tmp <- 1:39765
names(tmp) <- colnames(log.rpm)[1:39765]
rnd.seed <- tmp[as.character(glist)]
rnd.seed <- rnd.seed[!is.na(rnd.seed)]
the.gene <- log.rpm[ ,c(rnd.seed, 39766:39769)]
y.lim.v <- range(the.gene$Os01g0840100)
xx <- seq(-2,19,length=50)

the.gene2 <- log.rpm[c(rnd.seed, 39771,39772)]
the.gene2$LD <- as.numeric(the.gene2$LD) #dark=1, light=2

### fig
#j=20,k=35
pdf("../output/Fig.2d_Os01g0840100_35_20.pdf",width=2.3, height=1.3)
layout.mat2 <- matrix(c(1,2,3),byrow=F,nrow=1,ncol=3)
par(oma=c(3,3,3,3))
layout(layout.mat2,
widths = c(2,2,2), heights = 1)
j=20 # dark temp
k=35 # light temp
col.v <- color5[1:2]
env.param <- env.settings[env.settings$L.Tmp==k&env.settings$D.Tmp==j&env.settings$L!=0&env.settings$D!=0,]

for (l in c(8,12,16)) {
        par(mar=if(l==16) {c(1,0,0,0)} else {c(1,0,0,0)})
        fig.df    <- the.gene[(the.gene$EnvID==env.param$Env_ID[env.param$L==l]) & the.gene$Cultivar=="Kos",]
        plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
             ylim=y.lim.v, cex=1, yaxt=if(l==8) {"s"} else {"n"})
        
        sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                            fig.df[,1][fig.df$Cultivar=="Kos"])
        kos.est <- sp.kos(xx)
        polygon(c(21 - l, 21 - l,-3,-3), c(y.lim.v, y.lim.v[2], y.lim.v[1]), col = rgb(0,0,0,alpha=0.1),border = F)
        lines(xx,kos.est,col=color5[1])
}
dev.off()

### fig
#j=35,k=20
pdf("../output/Fig.2d_Os01g0840100_20_35.pdf",width=2.3, height=1.3)
layout.mat2 <- matrix(c(1,2,3),byrow=F,nrow=1,ncol=3)
par(oma=c(3,3,3,3))
layout(layout.mat2,
       widths = c(2,2,2), heights = 1)
j=35 # dark temp
k=20 # light temp
col.v <- color5[1:2]
env.param <- env.settings[env.settings$L.Tmp==k&env.settings$D.Tmp==j&env.settings$L!=0&env.settings$D!=0,]

for (l in c(8,12,16)) {
  par(mar=if(l==16) {c(1,0,0,0)} else {c(1,0,0,0)})
  fig.df    <- the.gene[(the.gene$EnvID==env.param$Env_ID[env.param$L==l]) & the.gene$Cultivar=="Kos",]
  plot(fig.df$Time - 24*floor(0.05*fig.df$Time) - 0.5, fig.df[,1], col=col.v[as.factor(fig.df$Cultivar)],
       ylim=y.lim.v, cex=1, yaxt=if(l==8) {"s"} else {"n"})
  
  sp.kos <- splinefun(fig.df$Time[fig.df$Cultivar=="Kos"]-24*floor(0.05*fig.df$Time[fig.df$Cultivar=="Kos"]),
                      fig.df[,1][fig.df$Cultivar=="Kos"])
  kos.est <- sp.kos(xx)
  polygon(c(21 - l, 21 - l,-3,-3), c(y.lim.v, y.lim.v[2], y.lim.v[1]), col = rgb(0,0,0,alpha=0.1),border = F)
  lines(xx,kos.est,col=color5[1])
}
dev.off()

### Fig.2e scatter plot ----
the.gene2$LD <- factor(the.gene2$LD,levels = c(2,1))

#alpha
a <- ggplot(the.gene2, aes(x = Tmp, y = Os01g0840100)) +
     geom_point(size=0.5,alpha=0.2,shape=21,fill="black",color="black") +
     labs(x = "air temperature (C)", y = "Os01g0840100\nlog2(rpm+1)") +
     theme_cowplot() +
     theme(plot.title = element_text(size=6, hjust=0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")  

b <- ggplot(the.gene2, aes(x = LD, y = Os01g0840100)) +
  geom_point(size=0.5,alpha=0.1,shape=21,fill="black",color="black") +
  scale_x_discrete(labels=c("Light","Dark"))+
  labs(x = "Light/Dark", y = "Os01g0840100\nlog2(rpm+1)") +
  theme_cowplot() +
  theme(plot.title = element_text(size=6, hjust=0.5),
        axis.title.x = element_text(size=6),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.position = "none")

g <- plot_grid(a, b, ncol=2, align = "hv",rel_widths = c(2,1.2))
ggsave(g, file="../output/Fig.2e.pdf", # Fig.2e
       width = 45, height = 45, units = "mm")

### Fig.S4 eigengenes ----
layout.mat2 <- matrix(c(
                       1,  6,7,8,   9,10,11, 12,13,14, 15,15,15, 16,17,18,
                       2,  19,20,21, 22,23,24, 25,25,25, 26,27,28, 29,30,31,
                       3,  32,33,34, 35,35,35, 36,37,38, 39,40,41, 42,43,44,
                       4,  45,45,45, 46,47,48, 49,50,51, 52,53,54, 55,56,57,
                       5,  58,59,60, 61,62,63, 64,65,66, 67,68,69, 70,71,72,
                       0,  0,0,73,   0,0,74,   0,0,75,   0,0,76,   0,0,77), byrow=T,nrow=6,ncol=16)

for (i in 1:(n_module)){
    
  the.gene <- d_eigengenes[ ,c(i+1,(n_module+2):(n_module+5))] # 1168 x 5 (RNA, EnvID, PotID, Time, Cultivar)
  
  pdf(sprintf("../output/Fig.S4_eigengene_module%s.pdf",i))
  
  par(oma=c(3,3,3,3))
  layout(layout.mat2,
         widths = c(2.7, rep(c(2,2,2.7), 5)), heights = rep(1,6))
  
  #==================================#
  # expression time series: all dark #
  #==================================#
  par(mar=c(1,0,0,1), mgp=c(1.5, 0.7, 0))
  xx <- seq(-2,19,length=50)
  y.lim.v <- c(min(d_eigengenes[,2:23]), max(d_eigengenes[,2:23]))
  #y.lim.v <- c(-0.2,0.2)
  ##########################
  ### L = na, D = 35-15C ###
  for (j in c(35,30,25,20,15)) {
    env.param <- env.settings[env.settings$L==0&env.settings$D.Tmp==j,]
    fig.df    <- the.gene[the.gene$EnvID==env.param$Env_ID,]
    #col.v <- c("blue", "red")
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
    #lines(xx,kos.est,col="blue")
    #lines(xx,tak.est,col="red")
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
          #lines(xx,kos.est,col="blue")
          #lines(xx,tak.est,col="red")
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
    #lines(xx,kos.est,col="blue")
    #lines(xx,tak.est,col="red")
    lines(xx,kos.est,col=color5[1])
    lines(xx,tak.est,col=color5[2])
  }
  dev.off()
}

### Table.S3 module GO analysis ----
use.genelist <- rownames(dat)
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]
term <- ng.GetGOTerms(unique(ulg2[,"GOid"]))

dir.create("GO")

for (i in sort(module_id)){
  gl <-  rownames(dat)[moduleLabels==i]
  result <- yh.mft(cgt=ulg2, gn.test=gl)[(term != "NA"),]
  
  fn <- sprintf("GO/GO_module%s.csv", i)
  
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
  
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in sort(module_id)){
  a <- read.csv(sprintf("GO/GO_module%s.csv",  i))
  if (is.na(a[1,1])){
    next
  }else{
    module <- rep(i, nrow(a))
    a <- cbind(module, a)
    tmp <- rbind(tmp, a)
  }
}
write.csv(tmp, file= "module_GO.csv", row.names = F)

### import module_id_convert.csv
module_id_convert <- read.csv("../input/module_id_convert.csv")

tmp <- read.csv("module_GO.csv")
tmp[,1] <- as.vector(module_id_convert$new[match(tmp$module, module_id_convert$old)])
tmp[,2] <- sprintf('%.2e', tmp[,2])
tmp <- tmp[!is.na(tmp$module),]

write.csv(tmp, file="../output/Table.S3.csv", row.names = F)

### Table.S4 module KEGG passway analysis ----
use.genelist <- rownames(dat)
kegg <- kegg_rice[is.element(kegg_rice[,"locus"], use.genelist),]
dir.create("KEGG")

for (i in sort(module_id)){
  gl <- rownames(dat)[moduleLabels==i]
  result <- kegg.mft(cgt=kegg, gn.test=gl)
  
  fn <- sprintf("KEGG/KEGG_%s.csv", i)
  
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
  cat(sprintf("%s / %s\n", i, Sys.time()))
}

# merge csv file
tmp <- NULL
for (i in sort(module_id)){
  a <- read.csv(sprintf("KEGG/KEGG_%s.csv",  i))
  if (is.na(a[1,1])){
    next
  }else{
    module <- rep(i, nrow(a))
    a <- cbind(module, a)
    tmp <- rbind(tmp, a)
  }
  #cat(sprintf("%s / %s\n", i, Sys.time()))
}

write.csv(tmp, file= "module_KEGG.csv", row.names = F)

tmp <- read.csv("module_KEGG.csv")
tmp[,1] <- as.vector(module_id_convert$new[match(tmp$module, module_id_convert$old)])
tmp[,2] <- sprintf('%.2e', tmp[,2])
tmp <- tmp[!is.na(tmp$module),]

write.csv(tmp, file="../output/Table.S4.csv", row.names = F)

### Fig.2f module GO figure ----
# adapted from http://clonecad.net/protocol_db/article/3697
# Import the table containing the enriched GO terms
GO_all <- read.csv("GO/GO_module25.csv",header=T,stringsAsFactors = T)
KEGG_all <-  read.csv("KEGG/KEGG_25.csv",header=T,stringsAsFactors = T)
GO_KEGG <- rbind(GO_all, KEGG_all)
GO_KEGG <- GO_KEGG[nrow(GO_KEGG):1,]

GO_data <- data.frame(GO_KEGG = paste(GO_KEGG$ID, GO_KEGG$Description, sep = " "),
           Gene_number = GO_KEGG$B,
           qvalue = GO_KEGG$Adjusted.P.value,
           log_qvalue = -(log10(GO_KEGG$Adjusted.P.value)),
           type = rep("GO",nrow(GO_KEGG))
           )

a <- ggplot(GO_data, aes(x = GO_KEGG, y = type, alpha=1)) +
  #geom_hline(yintercept = 1, linetype="dashed", 
  #           color = "azure4", size=.5)+
  geom_point(data=GO_data, aes(x=GO_KEGG, y=type, size = Gene_number, colour = log_qvalue, alpha=.7))+
  #scale_x_discrete(limits= GO_data$GO_KEGG)+
  scale_color_gradient(low="cyan",high="magenta",limits=c(0, NA))+
  coord_flip()+
  theme_bw(base_size=7)+
  theme(legend.direction = "horizontal",
        legend.position = "right")+
#        legend.justification = c(0,0))+
  #theme(axis.ticks.length=unit(-0.1, "cm"),
  #      axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
  #      axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
  #      axis.text = element_text(color = "black"),
  #      panel.grid.minor = element_blank(),
  #      legend.title.align=0.5)+
  xlab("")+
  ylab("")+
  labs(color="-log10(qvalue)", size="Number\nof genes")+
  guides(y = guide_axis(order=2),
         colour = guide_colourbar(order=1),
         alpha = "none")
         
ggsave(a, file="../output/Fig.2f.pdf", # Fig.2f
       units = "mm", width=150, height=55)

### export to cytoscape module2-59 (new module1-22)----
# Select modules
modules = colorOrder[2:59]
# Select module probes
probes = gene_17742
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule] # 5889 genes

modTOM = (adjacency > 0.05) * TOM
modTOM = modTOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "CytoscapeInput-edges_selected.txt",
                               nodeFile = "CytoscapeInput-nodes_selected.txt",
                               weighted = TRUE,
                               threshold = 0.05,
                               nodeNames = modProbes,
                               altNodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])

# color
a <- t(col2rgb(standardColors(22)))
rownames(a) <- standardColors(22)
#write.csv(a, file="240906_rgb_new.csv")
color_module <- rownames(a)

a <- read.delim("CytoscapeInput-nodes_selected.txt")
aa <- moduleLabels
names(aa) <- probes
a$module <- aa[a$nodeName]

converter <- data.frame(module = module_id_convert$old, 
                        module_new = module_id_convert$new,
                        color = c("grey",color_module))

aaa <- left_join(a,converter,by="module")

write.table(aaa, file = "CytoscapeInput-nodes2_selected.txt",
            sep= "\t", quote = F, row.names = F)

### Fig.2g export to cytoscape module25 ----
modules = colorOrder[25]
# Select module probes
probes = gene_17742
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]

modTOM = (adjacency > 0.05) * TOM
modTOM = modTOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "CytoscapeInput-edges_module25_all.txt",
                               nodeFile = "CytoscapeInput-nodes_module25_all.txt",
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               altNodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]) # input of Fig.2g

a <- des[modProbes,]
write.csv(a,"../output/Table.S5_module25.csv", row.names = F) # Table.S5

### detecting hubgene by cytoscape ----
# loading the output from cytoscape
a <- read.csv("network_analyzer_module25.csv")
rownames(a) <- a$name

hubgenes_module25 <- a$name[a$BetweennessCentrality > 0.3 & a$Degree >= 4] 
#"Os08g0525600","Os02g0115900","Os03g0113700" # Fig.2g

### overlap with module 5 and r > 0.5 ----
a <- read.delim("../output/Table.S6-7.tsv")
ddd <- gene_17742[moduleLabels == 25]

intersect(a$Locus_ID[61:151],ddd) # 38 genes
