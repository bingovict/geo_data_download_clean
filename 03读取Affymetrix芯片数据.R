############################################################################===
##### 04. 读取Affymetrix芯片数据 
############################################################################===


#### 1 使用affy包读取CEL文件 ----
# GPL570, GSE29450
# [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29450

### 读入Targets
setwd("D:/R/lianxi")
library(limma)
Targets <- readTargets("geo/GSE29450/Targets.txt")  # 不指定行名所在的列
Targets <- readTargets("geo/GSE29450/Targets.txt",
                       row.names = "sample_id")            # 指定行名所在的列

### 读入CEL
library(affy)
cel_files <- affy::list.celfiles("geo/GSE29450/GSE29450_RAW/", 
                                 full.name = T)
cel_files

## affy包中有2个读取CEL文件的函数，推荐ReadAffy(),更灵活易用
# cel <- read.affybatch(filenames = cel_files)
GSE29450_cel <- ReadAffy(filenames = Targets$FileName,
                         celfile.path = "geo/GSE29450/GSE29450_RAW",
                         phenoData = Targets)

## 使用widget选择读入
cel_widget <- ReadAffy(widget = T)

cel <- GSE29450_cel
# 表达矩阵 exprs
head(exprs(cel))[, 1:5]
# 样本信息 pData
GSE29450_targets <- pData(cel)
# 样本名称 sampeNames
sampleNames(cel)
# ScanDate
cel@protocolData@data$ScanDate
# expression matrix
GSE29450_expr_cel <- exprs(cel)

### 保存数据，备下一步分析用
save(GSE29450_cel, file = "geo/GSE29450/GSE29450_cel.Rdata")

#### 2 使用oligo包读取CEL文件 ----
# GPL6244, GSE24129 
# [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24129

library(oligo)
cel_files <- list.celfiles("geo/GSE24129/GSE24129_RAW/", 
                           full.name = T, 
                           listGzipped = T)
cel_files 

GSE24129_cel <- read.celfiles(cel_files)

cel <- GSE24129_cel

# exprs
head(exprs(cel))[, 1:5]
# pData
pData(cel)
# sampeNames
sampleNames(cel)

## new cel data 
## targets for GSE24129
cel_files

library(stringr)
FileName <- cel_files
sample_id <- str_extract(cel_files, "GSM\\d+")

group <- str_match(cel_files, "GSM\\d+_(.*)_\\d_Hu")[, 2]

sample_name <- str_match(cel_files, "GSM\\d+_(.*_\\d)_Hu")[, 2]

Targets <- data.frame(FileName = FileName,
                      sample_id = sample_id,
                      sample_name = sample_name,
                      group = group,
                      row.names = sample_id,
                      stringsAsFactors = F)

## new cel
GSE24129_cel <- read.celfiles(filenames = Targets$FileName,
                              phenoData = AnnotatedDataFrame(data = Targets), 
                              sampleNames = Targets$sample_id)
cel <- GSE24129_cel
# exprs
expr <- exprs(cel)
head(expr)[, 1:5]
# pData
pData(cel)
# sampeNames
sampleNames(cel)
# featureNames
tail(featureNames(cel))

### 保存数据，备下一步分析用
save(GSE24129_cel, file = "geo/GSE24129/GSE24129_cel.Rdata")

#### END ----
