############################################################################===
#####  11. 下载原始芯片数据
############################################################################===

#### 使用GEOquery包下载芯片原始文件

### 1 加载包 ----
library(GEOquery)

### 2 GSE数据，以GSE29450为例 ----
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29450

##  2.1 直接下载文件到当前工作目录下 ----

getGEOSuppFiles("GSE29450", baseDir ="geo", makeDirectory = T)  # 新建一个文件夹
getGEOSuppFiles("GSE29450", baseDir = "geo", makeDirectory = F) # 不会新建对应的文件夹
# 可能下载速度较慢

##  2.2 仅获取下载链接 ----
supp_url <- getGEOSuppFiles("GSE29450", fetch_files = F, makeDirectory = F)
url <- as.character(supp_url$url)
url

##  1.3.3 使用R函数download.file下载（getGEOSuppFiles中应用的函数）
download.file(url, destfile = "gse29450.tar")

### 3 GPL数据，以GPL6480为例 ----
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6480

## 获取下载链接
supp_url2 <- getGEOSuppFiles("GPL6480", fetch_files = F)
url2 <- as.character(supp_url2$url)
url2

library(tidyverse)

GSEs <- paste0("GSE", 29450:29459)
urls <- lapply(GSEs, function(x) {getGEOSuppFiles(x, makeDirectory = F, fetch_files = F)})
urls2 <- urls %>% 
  bind_rows() %>% 
  .[[2]] %>% 
  write_lines("urls.txt")

## GEO下载GPL注释表格

# GPL list
GPL6480 <- getGEO("GPL6480", destdir = "geo", AnnotGPL = T)

# GPL annotation table
GPL6480_table <- Table(GPL6480)

### END ----
