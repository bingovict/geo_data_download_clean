############################################################################===
#####  16. GEOquery下载读取数据
############################################################################===


#### 1 加载包 ----
library(readxl)
library(tidyverse)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)
# 需引用的文献
citation(package = "GEOquery") 

#### 2 getGEO 含单个数据集 GSE5327 ----
GSE5327_1 <- getGEO(GEO = "GSE5327", destdir = "geo", getGPL = F)#getGPL = F是因为gpl一般很大
GSE5327_11 <- GSE5327_1[[1]]

#当已经下载好了series_matrix文件的时候可以用如下代码读取，避免了一步提取的步骤
GSE5327_2 <- getGEO(filename = "geo/GSE5327_series_matrix.txt.gz", getGPL = F)

pd_GSE5327 <- pData(GSE5327_2)
#观察临床信息中的data processing,Microarray suite,MAS 5.0，即标准化方法，如果已经经过了MAS，一般情况下不需要再标准化了，避免矫枉过正

#### 3 getGEO 含两个数据集 GSE3494 ----

GSE3494_geo <- getGEO(GEO = "GSE3494", destdir = "geo", getGPL = F)
names(GSE3494_geo)

expr_GSE3494_1 <- exprs(GSE3494_geo[[1]])
expr_GSE3494_2 <- exprs(GSE3494_geo[[2]])

# GSE3494_geo_2 <- getGEO(filename = "geo/GSE3494_family.soft.gz", getGPL = F)
#读取soft格式的文件也不常用
#当确实没有series matrix格式的时候，只有soft格式的时候就不得不用soft格式

GSMs <- names(GSE3494_geo_2@gsms)

# tmp <- lapply(1:251, function(i) {
#   GSE3494_geo_2@gsms %>%
#     .[[i]] %>%
#     .@dataTable %>%
#     .@table %>%
#     select(1, 2) %>%
#     mutate(GSM = GSMs[i])
#   }) %>%
#   bind_rows() %>%
#   pivot_wider(names_from = GSM,
#               values_from = VALUE) %>%
#   arrange(ID_REF) %>%
#   column_to_rownames("ID_REF") %>%
#   as.matrix()

head(tmp)[, 1:6]

#### 4 getGSEDataTables ----
#用于下载geo页面的表格
GSE3494_table <- getGSEDataTables("GSE3494")

summary(GSE3494_table)

#### 5 直接读取 GSE169267 ----

GSE169267 <- read_tsv("geo/GSE119267_processed_data.txt")

GSE169267_2 <- read_xlsx("geo/GSE119267_processed_data.xlsx")#尽量不要用excel来处理数据，会出问题
# 基因名变化

#### 6 getGEO GPL6480 ----
  ##　AnnotGPL = T，用的是.gz的文件
GPL6480_1 <- getGEO("GPL6480", destdir = "geo", AnnotGPL = T)
GPL6480_11 <- Table(GPL6480_1)

  ## AnnotGPL = F，用的是.soft的文件
GPL6480_2 <- getGEO("GPL6480", destdir = "geo", AnnotGPL = F)
GPL6480_21 <- Table(GPL6480_2)

#### 7 getGEOSuppFiles ----
## 直接下载文件到当前工作目录下
getGEOSuppFiles("GSE29450", makeDirectory = F, baseDir = "geo")

## 仅获取下载链接
url_GSE29450 <-getGEOSuppFiles("GSE29450", fetch_files = F, makeDirectory = F)
url_GSE29450$url

#### END ----

