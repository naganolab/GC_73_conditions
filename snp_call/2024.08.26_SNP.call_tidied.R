
library(vcfR)
library(ggplot2)

load("kostak_specific_snp")
vcf.raw = read.vcfR("output_maxmissing01_minDP1.recode.vcf", verbose = FALSE, limit = 100000000000000 )

# convert vcf to gt
gt.raw <- extract.gt(vcf.raw, return.alleles = TRUE) # transform vcf -> genotype matrix
#save(gt.raw, file = "gt_raw_TM")
load("gt_raw_TM") # 247701 x 1536

# remove raws of gt.raw without SNP info
gt.flt <- NULL
gt.kos.sp.mini <- NULL
gt.tak.sp.mini <- NULL

for (i in 1:length(gt.kos.specific)) {
  locus <- names(gt.kos.specific)[i]
  if (sum(rownames(gt.raw)==locus)==1) {
    gt.flt <- rbind(gt.flt, gt.raw[rownames(gt.raw)==locus,])
    gt.kos.sp.mini <- c(gt.kos.sp.mini, gt.kos.specific[i])
    gt.tak.sp.mini <- c(gt.tak.sp.mini, gt.tak.specific[i])
  } else{}
}

#gt.flt 26 746 x 1536 matrix
#gt.***.sp.mini 26 746 character vector

cultivar.identifier <- function(x) {  # x: a column of gt.flt
  ifelse(x==".", ".", ifelse(x==gt.kos.sp.mini, "kos", ifelse(x==gt.tak.sp.mini, "tak", "unknown")))
}

kt.mat <- apply(gt.flt, 2, cultivar.identifier)
rownames(kt.mat) <- names(gt.kos.sp.mini)

#save(kt.mat, file = "kt.mat.rda")
load("kt.mat.rda")


#------------------
# load data for annotation

annotation <- read.table("sample.annotation.txt", header=T)
index      <- read.csv("IndexPrimer_ID.csv", header=T, sep=",")

index$name <- as.character(index$name)
index$name <- as.integer(substr(index$name,2,nchar(index$name)))
index$lib.grp <- as.character(index$lib.grp)

sample.ID <- NULL         # sample IDs for annotation
true.sample.col <- NULL   # kt.mat columns of samples (not blank), which correspond to "sample.ID"
for (i in 1:ncol(kt.mat)) {
  lib <- substr(colnames(kt.mat)[i], 1, 7)
  ind <- as.integer(substr(colnames(kt.mat)[i], 14, 16))
  temp.ind <- index[(index$lib.grp==lib)&(index$index==ind), ]
  
  if (nrow(temp.ind)==0) {} else {
    if (is.na(temp.ind$name)) {} else {
    sample.ID       <- c(sample.ID, temp.ind$name)
    true.sample.col <- c(true.sample.col, i)
    }
  }
}  
  
kt.mat.id <- kt.mat[, true.sample.col]
colnames(kt.mat.id) <- sample.ID  

# make a data.frame of kosishness, takishness, and cultivar
cultivar <- ifelse(colnames(kt.mat.id)%in%annotation[annotation$Cultivar=="Kos",1], "kos", "tak")

kosishness.calculator <- function(x) { # x is a vector of loci with cultivar identity
  no.loci        <- sum(x!=".")
  no.kosish.loci <- sum(x=="kos")
  no.kosish.loci/no.loci 
}

takishness.calculator <- function(x) {
  no.loci        <- sum(x!=".")
  no.takish.loci <- sum(x=="tak")
  no.takish.loci/no.loci 
}

kosishness <- apply(kt.mat.id, 2, kosishness.calculator)
takishness <- apply(kt.mat.id, 2, takishness.calculator)


SNP.call.df <- data.frame(cultivar, kosishness, takishness)

scat.plot <- ggplot(SNP.call.df) +
  geom_point(aes(x=kosishness, y=takishness, color=cultivar), shape=19, alpha=0.3, size=2) +
  xlim(c(0,1)) + ylim(c(0,1)) +
  theme_classic() + theme(aspect.ratio=1) +
  xlab("Proportion of loci with Koshihikari genotype") + 
  ylab("Proportion of loci with Takanari genotype") +
  scale_color_discrete(name = "Nominal\ncultivar") +
  geom_vline(xintercept=0.7, linetype="dashed", color="grey") +
  geom_hline(yintercept=0.7, linetype="dashed", color="grey")
scat.plot

# fig output
pdf("Fig.S2b.pdf",width = 3.8, height = 3.8)
scat.plot
dev.off()

# which sample is weird?
SNP.call.df[(SNP.call.df$kosishness < 0.7)&(SNP.call.df$takishness < 0.7), ]   # 507
annotation[annotation$SampleID==507,]

# this sample has above-average total read (~2 million, while the average is ~1 million)


#=======================================================#
# how many loci show SNP after excluding blank samples? #
#=======================================================#

gt.sample <- gt.raw[,true.sample.col]

gt.counter <- function(x) {
  length(levels(as.factor(x)))
}

no.genotype <- apply(gt.sample, 1, gt.counter)
sum(no.genotype==1)
hist(no.genotype,ylim=c(0,100))
nrow(gt.sample)-56 # 247 645 loci show SNP (i.e., polymorphic)


# no. of sequenced loci/sample
seq.loci.counter <- function(x) {
  sum(x!=".")
}

no.seq.loci <- apply(kt.mat.id, 2, seq.loci.counter)

min(no.seq.loci) # 4005
max(no.seq.loci) # 20953
mean(no.seq.loci) # 12 786.78
median(no.seq.loci) # 12717.5
hist(no.seq.loci,breaks=50)
