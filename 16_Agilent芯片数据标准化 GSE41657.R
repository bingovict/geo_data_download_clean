############################################################################===
##### 27. Agilent芯片数据标准化 GSE41657
############################################################################===


#### GPL6480平台数据集 GSE41657
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41657

library(limma)
library(GEOquery)
library(tidyverse) 

#### 1 GSE41657 data import ----
###  1.1 load data ----
load_input <- load("agilent/GSE41657/GSE41657_raw.Rdata")
load_input

dim(GSE41657_raw)

###  1.2 GSE41657_targets ----
GSE41657_targets <- GSE41657_raw$targets

GSE41657_sm <- getGEO(filename = "agilent/GSE41657/GSE41657_series_matrix.txt.gz",
                      getGPL = F) 
GSE41657_pd <- pData(GSE41657_sm) 

#### 2 Plots before QC ----

boxplot(log2(GSE41657_raw$E),
        ylab = expression(log[2](intensity)),
        las = 2,
        col = factor(GSE41657_targets$group),
        outline = FALSE)

PCA_new(log2(GSE41657_raw$E), 
        ntop = nrow(GSE41657_raw$E),
        group = GSE41657_targets$tissue_patho,
        show_name = F)

#### 3 QC using arrayQualityMrtics ----
library(arrayQualityMetrics)

#arrayQualityMetrics不接受ElistRaw，所以要构建ExpressionSet
rownames(GSE41657_targets) <- colnames(GSE41657_raw$E) <- GSE41657_targets$sample_id
eset <- ExpressionSet(assayData = GSE41657_raw$E, 
                      phenoData = AnnotatedDataFrame(data = GSE41657_targets))

arrayQualityMetrics(eset, 
                    outdir = "agilent/GSE41657/GSE41657_QC_raw", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()

#### 4 data normalization ----
#背景矫正是个在各自芯片水平上进行处理，标准化是在芯片之间的水平进行处理
GSE41657_bgc <-  backgroundCorrect(RG = GSE41657_raw, 
                                   method = "normexp",#推荐normexp进行背景矫正
                                   offset = 50,#补偿值50
                                   normexp.method = "mle")#limma包推荐的mle

GSE41657_norm <- normalizeBetweenArrays(GSE41657_bgc, #normalizeBetweenArrays适用于单色芯片
                                        method = "quantile")#normalizeinArrays适用于双色芯片


GSE41657_expr_norm <- GSE41657_norm$E

sum(is.na(GSE41657_expr_norm))

#### 5 Plots after normalization ----
## boxplot after normalization
boxplot(GSE41657_expr_norm,
        ylab = expression(log[2](intensity)),
        las = 2,
        col = factor(GSE41657_targets$group),
        outline = FALSE)

## PCA after normalization
PCA_new(GSE41657_expr_norm, 
        ntop = nrow(GSE41657_expr_norm),
        group = GSE41657_targets$tissue_patho)

##  GSE41657_sm
boxplot(exprs(GSE41657_sm),
        ylab = expression(log[2](intensity)),
        las = 2,
        col = factor(GSE41657_targets$group),
        outline = FALSE)
#可以看到series matrix的箱线图并不好，还不如从rawdata开始分析，灵活性更好，采用的方法更多
#从原始数据开始可以每一步都把关
PCA_new(exprs(GSE41657_sm), 
        ntop = nrow(exprs(GSE41657_sm)),
        group = GSE41657_targets$tissue_patho)

#### 6 QC after normalization ----
library(arrayQualityMetrics)

eset2 <- ExpressionSet(assayData = GSE41657_norm$E, 
                      phenoData = AnnotatedDataFrame(data = GSE41657_targets))
arrayQualityMetrics(eset2, 
                    outdir = "agilent/GSE41657/GSE41657_QC_norm", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()

#### 7 save data ----
ls(pattern = "^GSE41657")
save(list = ls(pattern = "^GSE41657")[c(-1, -5)], 
     file = "agilent/GSE41657/GSE41657_limma_processed.Rda")

#### End ----


