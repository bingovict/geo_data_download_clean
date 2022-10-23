############################################################################===
#####  12. 整理Targets数据
############################################################################===


### 1 加载R包 ----
library(tidyverse)
library(limma)

### 2 查看Targets文件说明 ----
readTargets(file="Targets.txt", 
            path=NULL, 
            sep="\t", 
            row.names=NULL, 
            quote="\"",...)

### 3 对象及列命名的几种习惯 ----
# 可读性、简洁、统一
# snake_case 
# camelCase
# CamelCase
# 一般推荐snake_case, 因可读性强，经典范例，stringr::str_系列👍
stringr::str_

### 4 创建Targets数据，以GSE29450为例 ----
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29450

FileName <- list.files(path = "geo/GSE29450/GSE29450_RAW/")
FileName
sample_id <- str_sub(FileName, 1, 9)
tissue_type <- rep(c("CCOC", "HOSE"), each = 10)
group <- rep(c("cancer", "normal"), each = 10)

# 从GEO网页复制samples信息
tmp <- "GSM729034	CCOC-1411G
        GSM729035	CCOC-1609
        GSM729036	CCOC-1638
        GSM729037	CCOC-2316
        GSM729038	CCOC-2348G
        GSM729039	CCOC-2501
        GSM729040	CCOC-451G
        GSM729041	CCOC-671G
        GSM729042	CCOC-805G
        GSM729043	CCOC-810G
        GSM729044	HOSE2008
        GSM729045	HOSE2061
        GSM729046	HOSE2064
        GSM729047	HOSE2085
        GSM729048	HOSE2225
        GSM729049	HOSE2226
        GSM729050	HOSE2228
        GSM729051	HOSE2230
        GSM729052	HOSE2234
        GSM729053	HOSE2237"

tmp <- str_split(tmp, pattern = "\\s+", simplify = T)

sample_name <- tmp[seq(2, 40, 2)]

Targets <- data.frame(FileName = FileName,
                      sample_id = sample_id,
                      sample_name = sample_name,
                      tissue_type = tissue_type,
                      group = group,
                      row.names = sample_id,
                      stringsAsFactors = F)

# options()$stringsAsFactors
head(Targets)
str(Targets)

## 写出Targets文件
write.table(x = Targets,
            file = "geo/GSE29450/Targets1.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write_tsv(x = Targets, 
          path = "geo/GSE29450/Targets.txt")

## 读取Targets文件
Targets_in <- readTargets(file = "geo/GSE29450/Targets.txt",
                          row.names = "sample_id")

Targets_in <- readTargets(path = "geo/GSE29450/",
                          row.names = "sample_id")

### END ----

save(Targets, file = "geo/GSE29450/GSE29450_targets.Rdata")
