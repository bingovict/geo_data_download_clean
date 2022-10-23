###############################################################################===
#####  17. 提取GEO表达矩阵 
###############################################################################===


#### 1 load packages ----
library(GEOquery)
library(readxl)
library(tidyverse)

#### 2 exprs from GEOquery ----
gse <- getGEO(filename = "geo/GSE5327_series_matrix.txt.gz", destdir = ".", getGPL = F )
expr <- exprs(gse)
head(expr)[, 1:6]

#### 3 expr from Series matrix file ----
expr1 <- read_tsv(file = "geo/GSE5327_series_matrix.txt.gz", comment = "!") #把前面有感叹号的行给去除掉
head(expr1)[, 1:6]
## A tibble: 6 x 6 tibble是tidyverse定义的一个新的数据框格式，可以显示每一列的数据类型，chr是字符，dbl是双精度的数值


#### 4 expr 从excel中读取的注意事项，基因名可能改变 ----
  #避免从excel中读取
expr2 <- read_xlsx("geo/GSE119267_processed_data.xlsx")%>% 
  tibble::column_to_rownames(var = "ID_REF")#把ID_REF这一列变成行名，也就是column_to_rownames
head(expr2)[, 1:6]

expr3 <- read_tsv("geo/GSE119267_processed_data.txt") %>% 
  tibble::column_to_rownames(var = "ID_REF")
head(expr3)[, 1:6]


#### 5 提取原始数据中的表达矩阵 ----
# from Affybatch
# from oligo
# from illumina
# from Agilent
# from limma
# 前面已经演示，不再重复

#### END ----

