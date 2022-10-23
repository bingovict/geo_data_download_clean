############################################################################===
#### 15. 读取Agilent芯片数据
############################################################################===


#### 1 使用limma包读取Agilent芯片数据 GSE41657 ----
# GPL6480, Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41657

#control type:0是对照，1就是基因
#glsaboveBG：信号强度是否高于背景强度，如果是的话就是1，说明这个值是可靠的


## Targets
library(limma)
Targets <- limma::readTargets(path = "agilent/GSE41657", row.names = "sample_id")
Targets[1:10, ]
table(Targets$group)

## read.maimages
GSE41657_raw <- read.maimages(files = Targets$FileName,
                              source = "agilent",#source代表是经过哪种程序得到的，有些Agilent芯片是通过genepix处理的
                              path = "agilent/GSE41657/GSE41657_RAW/",
                              names = Targets$sample_id,
                              other.columns = "gIsWellAboveBG",#读取进去是为了判断是否高于背景值
                              green.only = T)#代表是个单色芯片，默认是false

#GSE41657_raw的E代表expression,EB代表background

## add targets info
GSE41657_raw$targets <- Targets

## probe type summary
table(GSE41657_raw$genes$ControlType)
# -1     0     1   0说明不是control，需要检测这43376个基因，1是阳性对照，-1是阴性对照
# 153 43376  1486

## view data
head(GSE41657_raw$E)[, 1:5]
dim(GSE41657_raw)

## save data
save(GSE41657_raw, file = "agilent/GSE41657/GSE41657_raw.Rdata")





#### 2 使用AgiMicroRna包读取Agilent miRNA芯片数据 GSE108724 ----
# GPL20712 Agilent-070156 Human miRNA [miRNA version]
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108724
# SurePrint Human miRNA Microarrays G4872A-070156 
# Human miRNA Microarray, Release 21.0, 8 x 60K
#Agilent的miRNA芯片很受欢迎

## Targets file 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AgiMicroRna")

library(AgiMicroRna)
targets.micro <- readTargets(infile = "agilent/GSE108724/Targets1.txt")
rownames(targets.micro) <- paste0(targets.micro$Treatment, 1:7) 

## readMicroRnaAFE
GSE108724_raw <- readMicroRnaAFE(targets.micro)
dd.micro <- GSE108724_raw 

dim(dd.micro)
print(names(dd.micro))

## save data
save(targets.micro, dd.micro, GSE108724_raw, 
     file = "agilent/GSE108724/GSE108724_import.Rdata")

#### END ----

