############################################################################===
##### 25. Illumina芯片数据标准化 GSE108369-  Summary-level data
############################################################################===


#### GPL10558平台数据集 GSE108369
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108369
# summary-level data from txt files

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("lumi")

library(lumi)
library(tidyverse)

#### 1 GSE108369 data import ----
###  1.1 load data ----
load_input <- load("illumina/GSE108369/GSE108369_lumi.Rdata")
load_input

###  1.2 GSE108369_targets ----
library(GEOquery)
library(tidyverse)

GSE108369_sm <- getGEO(filename = "illumina/GSE108369/GSE108369_series_matrix.txt.gz",
                       getGPL = F) 

GSE108369_pd <- pData(GSE108369_sm) 

GSE108369_targets <- GSE108369_pd %>%
  dplyr::select(2, 1) %>%
  separate(title, into = paste0("x", 1:4), sep = " ") %>%
  mutate(sample_id = geo_accession,
         sample_name = x4,
         group = str_sub(x4 ,1, -3)) %>%
  dplyr::select(sample_id:group)

#### 3 Plots before QC ----
identical(colnames(GSE108369_expr_lumi),GSE108369_targets$sample_id)

boxplot(log2(GSE108369_expr_lumi),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
PCA_new(log2(GSE108369_expr_lumi), 
        ntop = nrow(GSE108369_expr_lumi),
        group = GSE108369_targets$group,
        show_name = F)

#### 4 QC using arrayQualityMrtics ----
library(arrayQualityMetrics)

citation(package = "lumi")
citation(package = "arrayQualityMetrics")

rownames(GSE108369_targets) <- colnames(GSE108369_expr_lumi)

#因为arrayQualityMetrics不支持lumi格式的数据，此时需要先构建ExpressionSet
eset <- ExpressionSet(assayData = GSE108369_expr_lumi, 
                      phenoData = AnnotatedDataFrame(data = GSE108369_targets))

arrayQualityMetrics(eset, 
                    outdir = "illumina/GSE108369/GSE108369_QC_aqm", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()

#### 6 data normalization ----
GSE108369_lumiN <-  lumiB(GSE108369_lumi, method = 'none') %>% #backgroud
                    lumiT(method = "vst") %>% #transfomation 
                    lumiN(method = "quantile") %>% #normaliztion
                    lumiQ()#quality control

GSE108369_lumiN2 <- lumiExpresso(GSE108369_lumi)#lumiExpresso将上述几步一步运行完

GSE108369_expr_lumiN <- exprs(GSE108369_lumiN2)


#### 7 Plots after normalization ----
## boxplot after normalization
boxplot(GSE108369_expr_lumiN,
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)

## PCA after normalization
PCA_new(GSE108369_expr_lumiN, 
        ntop = nrow(GSE108369_expr_lumiN),
        group = GSE108369_targets$group)

##  GSE108369_sm
boxplot(exprs(GSE108369_sm),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
PCA_new(exprs(GSE108369_sm), 
        ntop = nrow(GSE108369_expr_lumiN),
        group = GSE108369_targets$group)

####  save data
ls(pattern = "^GSE108369")
save(list = ls(pattern = "^GSE108369"), 
     file = "illumina/GSE108369/GSE108369_lumi_processed.Rda")

#### End ----
## 推荐VST转换，先用beadarry包把beadlevel的数据变成summary level的数据
#归纳之后用limi包的流程来走比较好


