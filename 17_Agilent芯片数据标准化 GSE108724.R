############################################################################===
#####  28. 读取Agilent芯片数据 GSE108724 miRNA data
############################################################################===


#### GPL6480平台数据集 GSE108724
# GPL20712 Agilent-070156 Human miRNA [miRNA version]
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108724

# AgiMicroRna
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037903/

library(AgiMicroRna)#专门为agilent的miRNA芯片打造
library(GEOquery)
library(tidyverse) 

#### 1 GSE108724 data import ----
###  1.1 load data ----
load_input <- load("agilent/GSE108724/GSE108724_import.Rdata")
load_input

identical(GSE108724_raw, dd.micro)

dim(GSE108724_raw)
GSE108724_expr_raw <- GSE108724_raw$meanS#原始的miRNA数据里有很多重复的，最开始数据是保存再meanS里的

head(dd.micro$TGS)[, 1:5]# total gene signal
head(dd.micro$meanS)[, 1:5]


###  1.2 GSE108724_targets ----
GSE108724_sm <- getGEO(filename = "agilent/GSE108724/GSE108724_series_matrix.txt.gz",
                      getGPL = F) 
GSE108724_pd <- pData(GSE108724_sm) 
GSE108724_targets <- targets.micro %>% 
  mutate(sample_id = str_extract(FileName, pattern = "GSM\\d+"),
         group = Treatment)
  
#### 2 Plots before normallization  ----
#3 箱线图
mi_raw <- dd.micro
boxplotMicroRna(log2(mi_raw$meanS),
                maintitle='log2 Mean Signal',
                colorfill= 'orange')
## 密度图
plotDensityMicroRna(log2(mi_raw$meanS),
                    maintitle='log2 Mean Signal')
## hclust
hierclusMicroRna(log2(mi_raw$meanS),
                 GErep,#为什么这里是GErep
                 methdis = "euclidean",
                 methclu = "complete",
                 sel = TRUE,
                 size = 300)

#
RleMicroRna(log2(mi_raw$meanS),
            maintitle = 'log2 Mean Signal - RLE')

# PCA
PCA_new(expr = log2(GSE108724_expr_raw),
        ntop = nrow(GSE108724_expr_raw),
        group = GSE108724_targets$group,
        show_name = T)



#### 3 QC using arrayQualityMertics ----

colnames(GSE108724_expr_raw) <- rownames(GSE108724_targets) <- GSE108724_targets$sample_id
eset <- ExpressionSet(assayData = GSE108724_expr_raw, 
                      phenoData = AnnotatedDataFrame(data = GSE108724_targets))

library(arrayQualityMetrics)
arrayQualityMetrics(eset, 
                    outdir = "agilent/GSE108724/GSE108724_QC_raw", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()

#### 4 data normalization ----
ddTGS <- tgsMicroRna(mi_raw,
                     half = TRUE,
                     makePLOT = FALSE,
                     verbose = FALSE)
#从744016变成了35980，每个miRNA有20个重复，并且把重复归纳了
#相当于affy芯片先做summary 后做normalization
ddNORM <- tgsNormalization(ddTGS,
                           "quantile",
                           makePLOTpre = FALSE,
                           makePLOTpost = FALSE,
                           targets.micro,
                           verbose = TRUE)

#要么是上面两步，要么是一步rma，借鉴的是affy包的rma方法
ddTGS.rma=rmaMicroRna(mi_raw,
                      normalize=TRUE,
                      background=TRUE)


#### 5 Plots after normalization ----
## boxplot after normalization
GSE108724_expr_norm <- ddNORM$TGS#因为已经归纳过了，就从TGS来提取

boxplot(GSE108724_expr_norm,
        ylab = expression(log[2](intensity)),
        las = 2,
        col = factor(GSE108724_targets$group),
        outline = FALSE)

## PCA after normalization
PCA_new(GSE108724_expr_norm, 
        ntop = nrow(GSE108724_expr_norm),
        group = GSE108724_targets$group)

##  GSE108724_sm
boxplot(exprs(GSE108724_sm),
        ylab = expression(log[2](intensity)),
        las = 2,
        col = factor(GSE108724_targets$group),
        outline = FALSE)
PCA_new(exprs(GSE108724_sm), 
        ntop = nrow(exprs(GSE108724_sm)),
        group = GSE108724_targets$group)

#### 6 QC after normalization ----
library(arrayQualityMetrics)


colnames(GSE108724_expr_norm) <- rownames(GSE108724_targets) <- GSE108724_targets$sample_id

eset2 <- ExpressionSet(assayData = GSE108724_expr_norm, 
                       phenoData = AnnotatedDataFrame(data = GSE108724_targets))
arrayQualityMetrics(eset2, 
                    outdir = "agilent/GSE108724/GSE108724_QC_norm", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()

#### 7 save data ----
GSE108724_norm <- ddTGS
GSE108724_norm_rma <- ddTGS.rma
  
ls(pattern = "^GSE108724")
save(list = ls(pattern = "^GSE108724"), 
     file = "agilent/GSE108724/GSE108724_agiMicroRna_processed.Rda")

#### End ----
