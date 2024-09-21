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
  if(is.null(nrow(gsea.test))==T){
    gn.test.gid <- table(gsea.test["passway"])
  }else{
  gn.test.gid <- table(gsea.test[, "passway"])
  }
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
