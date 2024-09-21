
# 以下のbioconductorのライブラリがインストールされていることが必要
# GO.db

# set the directory of GO analysis files
dir.go <- "GOanalysis_170821"

# load the GO analysis functions
fn <- sprintf("%s/GOanalysis_functions.R", dir.go)
source(fn)

# load the table of Gene ID and GO (as ulg)
fn <- sprintf("%s/ulg.Ahg_TAIR_140828", dir.go) # for A. halleri
load(fn) 





####
# 解析例
####

# load an example gene list (as gl)
# Flowering genes by Coupland lab.
fn <- sprintf("%s/example.genelist", dir.go)
load(fn) 

# fisher's exact test for all GO
result <- ng.mft(cgt=ulg, gn.test=gl)

# output results as a csv file
fn <- sprintf("%s/output.csv", dir.go)
write.csv(ng.prepGOtestOutTable(result), file=fn)





####
# 分母になる遺伝子セットを絞る場合
####

# load the gene list passed criteria for seasonal analysis 
fn <- sprintf("%s/use.genelist", dir.go)
load(fn) 

# make subset of ulg
ulg2 <- ulg[is.element(ulg[,"locus"], use.genelist),]

# fisher's exact test for all GO
result2 <- ng.mft(cgt=ulg2, gn.test=gl)

# output results as a csv file
fn <- sprintf("%s/output2.csv", dir.go)
write.csv(ng.prepGOtestOutTable(result), file=fn)



####
# glの中からあるGOを持つものを抽出する
####

withgo <- ulg[ulg[,"GOid"]=="GO:0019748", "locus"]
intersect(gl, withgo)
