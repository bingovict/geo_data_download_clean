############################################################################===
##### 23. Illumina芯片数据标准化GSE140882 -  Bead-level data
############################################################################===

#### GPL10558平台数据集 GSE140882 
# GPL10558, GSE140882 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140882

#### 1 GSE140882 BLdata import ----
###  1.1 load data ----
library(beadarray)
library(tidyverse)

load_input <- load("geo/GSE140882/GSE140882_beadlevel.Rda")
load_input

###  1.2 BLdata ----
BLdata <- GSE140882_BLdata
targets <- BLdata@sectionData$Targets#提取targets，但是这个targets仅仅是一部分信息，还需要自己整理
annotation(BLdata)

### boxplot for raw data
# samples <- str_sub(names(GSE140882_BL@beadData), 1, 10) 
# df <- list()
# for (i in 1:27) {
#   df[[i]] <- data.frame(sample_id = samples[i],
#                         probe_id = GSE140882_BL@beadData[[i]]$ProbeID$ProbeID,
#                         value = GSE140882_BL@beadData[[i]]$Grn$Grn)
# }
# df <- bind_rows(df)
# ggplot(data = df, aes(x = value, col = sample_id)) +
#     geom_density()

###  1.2 GSE140882_targets ----
library(GEOquery)
library(tidyverse)

GSE140882_sm <- getGEO(filename = "illumina/GSE140882/GSE140882_series_matrix.txt.gz",
                            getGPL = F) 

GSE140882_targets <- pData(GSE140882_sm) %>%
  dplyr::select(2, 1) %>%#选出两列，title 这一列其实包含了四列内容
  separate(title, into = paste0("x", 1:4), sep = " ") %>%#把title 这一列拆成四列，根据空格来分割
  mutate(bead_num = names(GSE140882_BLdata@beadData),#mutate就是加列
         sample_id = geo_accession,
         group = str_c(x1, x2, x3, sep = "_"),#str_c将多个东西连在一起
         sample_name = str_c(x1, x2, x3, x4, sep = "_"),
         cell_type = x1,
         treatment = x2,
         time = x3,
         rep = x4) %>%
  dplyr::select(bead_num:rep)#从bead_num这一列到rep这一列选择出来

#### 2 QC using beadarray ----

# library(hwriter)
# expressionQCPipeline(BLData = BLdata,
#                      qcDir = "illumina/GSE140882_QC",
#                      overWrite = T)
# 耗内存，耗时长
# qcReport <- makeQCTable (BLData =BLdata)


#### 3 data summary ----
#把重复的数据总结成bead type，从bead level变成summary level
datasumm <- beadarray::summarize(BLData = BLdata )
# error:there is no package called ‘illuminaHumanv4.db’
BiocManager::install("illuminaHumanv4.db")

dim(datasumm)
# Features  Samples Channels 
# 48107O       27        1 

GSE140882_expr_BLsumm <- exprs(datasumm)#提取表达矩阵
#修改列名，取1到10
colnames(GSE140882_expr_BLsumm) <- str_sub(colnames(GSE140882_expr_BLsumm), 1, 10)

exprs(datasumm)[1:10 , 1:2]
se.exprs(datasumm)[1:10 , 1:2]#表达量的标准差

boxplot(exprs(datasumm),#表达箱线图
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)

boxplot (nObservations(datasumm),#beads的分布图，每个探针重复的beads数量
         ylab = " number of beads ",
         las = 2,
         outline = FALSE)
#Error in nObservations(datasumm) : 没有"nObservations"这个函数


det <- calculateDetection (datasumm)#lumi包有个detection列，也可以加上这一列
head (det)#detection的值越小，信号值越可靠
Detection (datasumm) <- det


#### 4 QC using arrayQualityMrtics ----
#对summary level的数据进行质控
library(arrayQualityMetrics)
pData(GSE140882_Summdata) <- data.frame(pData(GSE140882_Summdata), GSE140882_targets)
arrayQualityMetrics(GSE140882_Summdata, 
                    outdir = "illumina/GSE140882/GSE140882_QC_aqm", 
                    do.logtransform = FALSE,
                    force = TRUE,
                    intgroup = "group")

#### 5 PCA plot ----
PCA_new(GSE140882_expr_BLsumm, group = GSE140882_targets$group)#可以看到分成了两陀细胞
PCA_new(GSE140882_expr_BLsumm[, 1:15], group = GSE140882_targets$group[1:15])#先看HBL1细胞
PCA_new(GSE140882_expr_BLsumm[, 16:27], group = GSE140882_targets$group[16:27])#再看HT细胞


#### 6 data normalization ----
GSE140882_expr_BLsummnorm <- limma::normalizeBetweenArrays(GSE140882_expr_BLsumm)


#### 7 Plots after normalization ----
## boxplot after normalization
boxplot(GSE140882_expr_BLsummnorm,#箱线图看一下标准化后的数据
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)

## PCA after normalization
PCA_new(GSE140882_expr_BLsummnorm, group = GSE140882_targets$group)
PCA_new(GSE140882_expr_BLsummnorm[, 1:15], group = GSE140882_targets$group[1:15])
PCA_new(GSE140882_expr_BLsummnorm[, 16:27], group = GSE140882_targets$group[16:27])
#可以看到不同颜色的差别越来越远
#PC1 78%，PC1可以解释78%的数据差异变化

####  save data
GSE140882_Summdata <- datasumm
ls(pattern = "^GSE140882")
save(list = ls(pattern = "^GSE140882")[2:6], 
     file = "illumina/GSE140882/GSE140882_beadarray_processed.Rda")

#### End ----

