#################################################===
##### 18. 提取GEO注释信息 
#################################################===


#### 1 annotation from GEOquery ----
library(GEOquery)
GPL6480 <- getGEO(GEO = "GPL6480", destdir = "geo")
GPL6480_2 <- getGEO(filename = "geo/GPL6480.soft")
identical(GPL6480, GPL6480_2)#判断两个是不是一样的
# [1] TRUE

Meta(GPL6480)$title#meta函数提取平台的注释信息
# [1] "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)"

#两种方法提取需要的数据，一是table，二是提取子集的方式
anno_GPL6480_geoquery <- Table(GPL6480)
anno_GPL6480_geoquery <- GPL6480@dataTable@table

#### 2 annotation from GPL SOFT file ----
library(readr)
anno_GPL6480_geosoft <- read_tsv("geo/GPL6480.soft", comment = "!", skip = 43866)#不推荐这种方式，因为要看skip了多少行

#### 3 annotation from microarray manufacture website ----
# https://earray.chem.agilent.com/earray/catalogGeneLists.do?action=displaylist

#安捷伦芯片需要知道设计ID号才方便在官网上找到信息，Agilent-014850 ，014850 就是ID号
anno_GPL6480_agilent <- read_tsv("geo/014850_hs1_1586063945255.zip")#tsv的优点是可以直接读取zip文件
all(anno_GPL6480_agilent$ProbeID %in% anno_GPL6480_geoquery$ID)


#### 4 annotation from Bioconductor ----
 # https://bioconductor.org/packages/release/BiocViews.html#___AgilentChip
library(hgug4112a.db)
hgug4112a()

anno_GPL6480_bioc <- as.data.frame(hgug4112aSYMBOL)
anno_GPL6480_bioc <- toTable(hgug4112aSYMBOL)
anno_GPL6480_bioc2 <- toTable(hgug4112aENTREZID)

all(anno_GPL6480_bioc$probe_id %in% anno_GPL6480_geoquery$ID)


#### END ----