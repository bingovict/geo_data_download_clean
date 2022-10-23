############################################################################===
##### 14. 读取Illumina芯片数据 
############################################################################===


#### 1 使用beadarray包读取illumina芯片数据 ----
  ###  1.1 beadarray::readIllumina, txt GSE140882 ----
# GPL10558
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140882
# bead-level data without TIFF

## 文件所在目录
dir <- "geo/GSE140882/GSE140882_RAW"
# dir <- "geo/GSE140882/GSE140882_RAW/"  # 不能加/

## 读取目录中的文件, bead-level
library(beadarray)
GSE140882_BLdata <- readIllumina(dir = dir)

## add annotation
annotation(GSE140882_BLdata) <-  suggestAnnotation(GSE140882_BLdata)

## view data
head(GSE140882_BLdata[[10]])
slotNames(GSE140882_BLdata)
str(GSE140882_BLdata)
dim(GSE140882_BLdata)

## save data 
save(GSE140882_BLdata, 
     file = "geo/GSE140882/GSE140882_beadlevel1.Rda")

load(file = "geo/GSE140882/GSE140882_beadlevel1.Rda")

  ###  1.2 beadarray::readIllumina, txt tiff GSE59183 ----
# GPL10558, GSE59183
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59183
# bead-level data with TIFF

## 文件所在目录
dir = "illumina/GSE59183/GSE59183_RAW"

## 读取目录中的文件, bead-level
library(beadarray)

# read raw data by readIllumina
# without TIFF
data_ntiff <- readIllumina(dir = dir, useImages = F, dec = ".", 
                           illuminaAnnotation = "Humanv4")
head(data_ntiff[[1]])

# with TIFF
data_tiff <- readIllumina(dir, useImages = T, dec = ".", 
                          illuminaAnnotation = "Humanv4")
head(data_tiff[[1]])

### read TIFF 
tif_files <- list.files("illumina/GSE59183/GSE59183_RAW/", 
                        pattern = "tif", full.names = T)
TIFF <- readTIFF(tif_files[1])

dim(TIFF)
class(TIFF)

data <- data_ntiff
xcoords <- getBeadData(data, array = 1, what = "GrnX")
ycoords <- getBeadData(data, array = 1, what = "GrnY")

# 
bg <- medianBackground(TIFF, cbind(xcoords, ycoords))
bg2 <- illuminaBackground(TIFF, cbind(xcoords, ycoords))
data <- insertBeadData(data, array = 1, what = "GrnRB", bg)
data <- insertBeadData(data, array = 1, what = "GrnRB2", bg2)
head(data[[1]])

# apply Illumina's image filtering
TIFF2 <- illuminaSharpen(TIFF)

# calculate foreground values 前景值-背景值=信号值
# 对于illumina数据而言，就不需要矫正背景数据了，就没有了背景矫正这个环节
fg <- illuminaForeground(TIFF2, cbind(xcoords, ycoords))
data <- insertBeadData(data, array = 1, what = "GrnF", fg)
head(data[[1]])

#bg
data <- backgroundCorrectSingleSection(data, array = 1, fg = "GrnF", bg = "GrnRB", newName = "GrnR")
head(data[[1]])
head(data_tiff[[1]])
head(data_ntiff[[1]])

# add annotation
suggestAnnotation(data_ntiff, verbose = T)
annotation(data) <- suggestAnnotation(data, verbose = TRUE)

# structure of data 
slotNames(data)
data@sectionData$Targets$greenImage
data@sectionData$Targets$textFile
head(data[[1]])
colnames(data[[1]])
sectionNames(data)
numBeads(data)

  ###  1.3 beadarray::readIdatFiles,idat GSE113440 ----
# GPL10558
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113440
# summary-level data from idat files

## 文件所在目录
dir <- "illumina/GSE113440/GSE113440_RAW"
idat_files <- list.files(dir, pattern = "idat", full.names = T)

## 读取idat files
library(beadarray)
GSE113440_idat <- readIdatFiles(idatFiles = idat_files)

## raw expression matrix 
GSE113440_expr_idat <- GSE113440_idat@assayData$exprs

## probe status
GSE113440_probe_status <- GSE113440_idat@featureData@data
table(GSE113440_probe_status$Status)

## save data
save(GSE113440_idat, GSE113440_expr_idat, GSE113440_probe_status,
     file = "illumina/GSE113440/GSE113440_idat.Rdata")






#### 2 使用lumi包读取illumina芯片数据 GSE108369 ----
# GPL10558, GSE108369
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108369
# summary-level data from txt files

## 文件所在目录
dir <- "illumina/GSE108369/GSE108369_RAW"
files <- list.files(dir, pattern = "GSM", full.names = T)
files

## 整理targets
library(tidyverse)
files
targets <- data.frame(FileName = files,
                      sample_id = str_match(files, "GSM\\d+")[, 1],
                      sample_name = str_match(files, "eum_(.*).txt")[, 2],
                      group = str_match(files, "eum_(.*)_")[, 2],
                      row.names = str_match(files, "GSM\\d+")[, 1],
                      stringsAsFactors = F)
targets

## preview files
preview <- read_tsv(files[1])

## revise file format
files
for (i in 1:length(files)) {
  # read_tsv读入，dplyr::select选择第1列和第3到第6列
  tmp <- readr::read_tsv(files[i]) %>% dplyr::select(1,3:6)
  
  # 改列名
  colnames(tmp) <- paste0(c("ProbeID", "AVG_SIGNAL","BEAD_STD",'Avg_NBEADS','DETECTION'),
                          "_", targets$sample_id[i])
  # 写出
  readr::write_tsv(tmp, path = paste0("illumina/GSE108369/", basename(files[i])))
  # write.table(tmp, file = basename(files[i]), sep = "\t", col.names = T,row.names = F, quote = F)
}

## 应用 lumi读取文件 
library(lumi)
GSE108369_lumi <- lumiR.batch(fileList = paste0("illumina/GSE108369/", basename(files)), 
                              detectionTh = 0.01,
                              sampleInfoFile = targets)

## exprssion matrix
GSE108369_expr_lumi <- exprs(GSE108369_lumi)
dim(GSE108369_expr_lumi)
sum(is.na(GSE108369_expr_lumi))#判断是否有缺失值，如果有的话就返回1，再计算总和就可以知道有几个缺失值

## save data 
save(GSE108369_lumi, GSE108369_expr_lumi, 
     file = "illumina/GSE108369/GSE108369_lumi.Rdata")





#### 3 使用limma包读取illumina芯片数据 ----
  ###  3.1 GSE16997 ----
## GPL6884	Illumina HumanWG-6 v3.0 expression beadchip
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16997

## 文件所在目录
dir <- "illumina/GSE16997/GSE16697_RAW"
files <- list.files(dir, full.names = T)
files

## read.ilmn (with control)
library(limma)
# read.ilmn(files=NULL, ctrlfiles=NULL, path=NULL, ctrlpath=NULL, probeid="Probe",
#           annotation=c("TargetID", "SYMBOL"), expr="AVG_Signal",
#           other.columns="Detection", sep="\t", quote="\"", verbose=TRUE, ...)

GSE16997_with_ctrl <- read.ilmn(files = files[3],#files参数是指probe profile
                                ctrlfiles = files[1],#ctrlfiles参数是指control probe profile
                                other.columns = "Detection")#后续的处理的时候需要Detection

x <- GSE16997_with_ctrl
str(x)#看数据结构
table(x$genes$Status)#看探针的总结，有哪几种探针类型，包括对照的、管家基因。这个数据可以用于质控
head(x$E)
dim(x$E)#看数据维度
boxplot(log2(x$E), range = 0, ylab = "log2 intensity")

## read.ilmn (without control)，没有control的时候就可以只读原始信息就够了，GSE16997_raw
preview <- read_tsv(files[2])
colnames(preview)

GSE16997_no_ctrl <- read.ilmn(files = files[2],
                              expr = "Sample",#告诉函数表达量函数是以sample开头的
                              other.columns = "Detection",
                              probeid = "ID_REF")#不告知函数哪一列是需要的话就和函数默认的命名不能匹配上

head(GSE16997_no_ctrl$E)
dim(GSE16997_no_ctrl$E)

## save data
save(GSE16997_with_ctrl, GSE16997_no_ctrl, 
     file = "illumina/GSE16997/GSE16697_limma.Rdata")

  ###  3.2 GSE42242 ----

## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42242
preview <- read_tsv("GSE42242_non-normalized.txt")
colnames(preview)

GSE42242_non_norm <- read.ilmn(files = "GSE42242_non-normalized.txt",
                               expr = "SAMPLE",#这里要根据数据实际的样子来改，有时候又变成了sample全
                               other.columns = "Detection",
                               probeid = "ID_REF")

head(GSE42242_non_norm$E)
dim(GSE42242_non_norm$E)

## save data 
save(GSE42242_non_norm, file = "geo/GSE42242/GSE42242_non-normalized.Rdata")



#### END ----




####总结
#如果数据量很大，每个gse整百万行的那种，一般都是bead-level的数据，用beadarray::readIllumina
#对summary-level来说，如果是idat格式的，可以用beadarray::readIdatFiles来读取
#看到数据有四列，"ProbeID", "AVG_SIGNAL","BEAD_STD",'Avg_NBEADS','DETECTION'，需要把熟做转换，就是改个名字+样本名字
#limma包就是一个sample一个detection这种




