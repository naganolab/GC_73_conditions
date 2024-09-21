library("KEGGREST")
library("tidyverse")

# make KEGG passway vs gene ID list ----
kegg_dosa <- keggLink("pathway", "dosa")
locus <- gsub("dosa:", "", names(kegg_dosa))
passway <- gsub("path:", "", kegg_dosa)

for(i in 1:length(locus)){
  locus[i] <- str_replace(locus[i], pattern="t", replacement="g") %>% 
              str_sub(.,1,12) 
}

description <- rep(NA, length(unique(passway)))

for (i in 1:length(unique(passway))){
  description[i] <- keggGet(unique(passway)[i])
}

kegg_description <- rep(NA, length(description))

for (i in 1:length(kegg_description)){
  kegg_description[i] <- description[[i]]$PATHWAY_MAP
}
names(kegg_description) <- unique(passway)


kegg_rice <- cbind(locus, passway, description)
rownames(kegg_rice) <- locus
save(kegg_rice, file = "kegg_rice")
save(kegg_description, file = "kegg_description")

# gene id not duplicated
#for (i in unique(passway)){
#  a <- length(unique(locus[(passway == i)]))
#  b <- length(locus[(passway == i)])
#  if (a == b){
#    cat("")
#  }else{
#    cat(sprintf("%s",i))
#  }
#}

# function ----
kegg.mft <- function(
  cgt, #kegg_rice, [,1]:"locus", [,2]:"passway"
  gn.test, #contig names for test
  alternative="greater"
){
  
  #cat(sprintf("%s\n", Sys.time()))
  
  gid.u <- unique(cgt[,"passway"])
  
  ft.in <- matrix(0, nrow=length(gid.u), ncol=9)
  colnames(ft.in) <- c("xtt", "xft", "xtf", "xff", "xnt", "xnf", "xtn", "xfn", "xnn")
  rownames(ft.in) <- gid.u
  
  #               gn.test
  #             TRUE FALSE
  #Group  TRUE   xtt   xft   xnt
  #      FALSE   xtf   xff   xnf
  #              xtn   xfn   xnn
  
  ft.in[,"xnn"] <- length(unique(cgt[, "locus"]))
  
  gn.pp.gid <- table(cgt[, "passway"])
  ft.in[names(gn.pp.gid), "xnt"] <- gn.pp.gid
  ft.in[,"xnf"] <- ft.in[,"xnn"] - ft.in[,"xnt"]
  
  ft.in[,"xtn"] <- length(intersect(gn.test, unique(cgt[, "locus"])))
  ft.in[,"xfn"] <- ft.in[,"xnn"] - ft.in[,"xtn"]
  
  gsea.test <- cgt[is.element(cgt[,"locus"], gn.test), ]
  gn.test.gid <- table(gsea.test[, "passway"])
  ft.in[names(gn.test.gid), "xtt"] <- gn.test.gid
  
  ft.in[,"xtf"] <- ft.in[,"xtn"] - ft.in[,"xtt"]
  ft.in[,"xft"] <- ft.in[,"xnt"] - ft.in[,"xtt"]
  ft.in[,"xff"] <- ft.in[,"xnf"] - ft.in[,"xtf"]
  
  #cat(sprintf("%s\n", Sys.time()))
  
  #Fisher's exact test.  8? sec
  fr <- rep(1, nrow(ft.in))
  dt <- rep(1, nrow(ft.in))
  for(i in 1:nrow(ft.in)){
    start <- Sys.time()
    if(ft.in[i,"xtn"] > 1 && ft.in[i,"xnt"] > 1){ 
      contable <- matrix(ft.in[i, 1:4], ncol=2)
      tmp <- fisher.test(contable, alternative = alternative)
      fr[i] <- tmp$p.value
    } else {
    }
    end <- Sys.time()
    dt[i] <- end - start
  }
  
  out <- cbind(fr, ft.in, dt)
  colnames(out) <- c("p.value", colnames(ft.in), "time")
  rownames(out) <- rownames(ft.in)
  
  #cat(sprintf("%s\n", Sys.time()))
  
  return(out)
  
}

# how to use ----
rpm2.ave <- rpm.ave[raw_10,]
use.genelist <- rownames(rpm2.ave)
kegg <- kegg_rice[is.element(kegg_rice$locus, use.genelist),]

gl <- UNREP_high
result <- kegg.mft(cgt=kegg, gn.test=gl)

fn <- sprintf("%s/%s_KEGG_Unrep_high.csv", dir.output,exec.date)

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

