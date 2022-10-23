#################################################===
##### 19. 提取GEO芯片样本及临床信息 
#################################################===


#### 1 clinical data from pData ----
library(GEOquery)
# gse <- getGEO("GSE5327", destdir = "geo", getGPL = F)
# pd <- phenoData(gse[[1]]) phenoData用法不常见

gse1 <- getGEO(filename = "geo/GSE5327_series_matrix.txt.gz",  getGPL = F )
pd1 <- pData(gse1)

#### 2 targets data from pData ----
GSE29450_sm <- getGEO(filename = "affymetrix/GSE29450/GSE29450_series_matrix.txt.gz",
                           getGPL = F)

library(tidyverse)
GSE29450_targets <- pData(GSE29450_sm) %>%
  dplyr::select(sample_id = geo_accession,# dplyr::select选中提取，并且重命名列名
                sample_name = title,
                tissue_type = source_name_ch1) %>%
  mutate(group = str_sub(tissue_type,-6,-1))#从tissue_type里找到cancer和normal,从倒数第一个到倒数第六个，mutate就是增加一列，列名叫group

#### 3 clinical data from csv file  ----
library(readr)
clin1 <- read_csv("geo/GSE124535_clinical_info.csv")

#### 4 clinical data from Excel file ----
library(readxl)
clin2 <- read_xls("geo/GSE10141_clinical_data_addition_081611.xls", sheet = 1)#sheet 表一还是表二

#### END ----