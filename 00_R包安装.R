############################################################################===
#####  01. R包安装 
############################################################################===


#### 0 开始之前，查看文件：01 Readme_First.txt ----
readLines("01 Readme_First.txt")  
readr::read_lines("01 Readme_First.txt")

#### 1 设置镜像 ----
###  1.1 查看当前的CRAN及Bioconductor镜像 ----
options()$repos
options()$BioC_mirror
BiocManager::repositories()

###  1.2 手动选择镜像 ----
chooseCRANmirror()
chooseBioCmirror()

###  1.3 直接设定镜像 ----
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options("BioC_mirror" = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

#### 2 安装R包的一个自定义函数pkgs_in() ----
pkgs_in <- function(pkgs) {
  # 设置为国内清华镜像
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options("BioC_mirror" = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
  
  # 首先安装BiocManager，若已安装则跳过
  if (!requireNamespace("BiocManager", quietly = TRUE)) { 
    install.packages("BiocManager",ask = F, update = F)
  }
  
  # 安装stringr，若已安装则跳过
  if (!requireNamespace("stringr", quietly = TRUE)) { 
    install.packages("stringr",ask = F, update = F)
  }
  
  # 去重，识别有无github格式的安装包
  pkgs <- unique(pkgs)
  pkgs2 <- pkgs
  logi <- stringr::str_detect(pkgs2, "/")
  pkgs2[logi] <- stringr::str_match(pkgs2[logi], ".*/(.*)$")[,2]
  
  # 安装pkgs中尚未安装的包
  new <- !(sapply(pkgs2, requireNamespace, quietly = T))
  
  # 显示需安装的包
  if (sum(new) > 0) {
    cat("pkgs to install: ", pkgs[new], "\n")
  } else {
    cat("All pkgs already installed \n")
  }
  
  # install pkgs
  if(any(new)) BiocManager::install(pkgs[new], ask = F, update = F)
}

#### 3 需要安装的pkgs ----
pkgs <- c("tidyverse", "limma", "affy", "oligo", "lumi",
          "beadarray", "GEOquery", "simpleaffy", "gcrma", "readxl",
          "impute", "genefilter", "pd.hugene.1.0.st.v1", "pd.hg.u133.plus.2",
          "tkWidgets", "illuminaHumanv4.db", "AnnotationDbi", "org.Hs.eg.db",
          "hgug4112a.db", "AgiMicroRna", "sva", "DESeq2", "edgeR",
          "lumiHumanIDMapping", "remotes", "pheatmap", "shiny", "aggregation",
          "tidyverse/dplyr", "limma", "hwriter", "devtools")

#### 4 安装pkgs中的R包 ----
###  pkgs_in()函数的特点 
# 可随时添加一个或多个安装包，若已安装的则不会再重复安装
# 支持安装github包，输入格式如“tidyverse/dplyr”，中间包含 / 
# 确认输入的R包名字是正确的，注意区分字母大小写
# 可反复运行，安装不成功时重启R session,再次安装
# 重启R session快捷键：Crtl + Shift + F10
# 反复运行，依然安装不成功，再次检查输入是否正确；也可能由于网络不佳，可换用其他安装方法

pkgs_in(pkgs)  
pkgs_in(pkgs)
pkgs_in(pkgs)

# 最后运行pkgs_in()时应无提示，表明安装成功；
# library()，不报error，提示安装成功。
# 若仍有安装不成功的包，可参考以下的方法。

#### 5 常规安装方法，用于上述方法未能成功安装的R包  ----
###  5.1 安装CRAN来源的R包 ----
# 替换为需要安装的R包名字, 如 --> "BiocManager" <--
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", ask = F, update = F)
}

# or
# install.packages("BiocManager")

###  5.2 安装Bioconductor来源的R包 ----
# 首先确认 BiocManager 包是否已安装成功
library(BiocManager)

## 安装Bioconductor R包 
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", ask = F, update = F)
}

# or
# BiocManager::install("limma")

## 若为GitHub来源的包, pkgs_in函数也支持；但依赖于网络状况, 有时下载速度慢
library(sleuth)
# BiocManager::install("sleuth")
# install.packages("sleuth")
# pkgs_in("pachterlab/sleuth")
# 
# devtools::install_github("pachterlab/sleuth")
# BiocManager::install("pachterlab/sleuth")


###  5.3 从下载的R包压缩文件安装 ----
# 以limma包为例
# https://bioconductor.org/packages/release/bioc/html/limma.html

# 重启一次 R session   Crtl + Shift + F10

# Windows Binary
remove.packages("limma")
library(limma)
install.packages("packages/limma_3.42.2.zip", repos = NULL, type = "win.binary")
library(limma)
# 
# Source Package
remove.packages("limma")
library(limma)
install.packages("packages/limma_3.42.2.tar.gz", repos = NULL, type = "source")
library(limma)

###  5.4 从GitHub网站下载的已解压的源文件安装 ----
# 以sleuth为例
# https://github.com/pachterlab/sleuth

# remove.packages("sleuth")
# install.packages("packages/sleuth-master/", repos = NULL, type = "source")
library(sleuth)
sleuth::
# 注意查看error报错的信息，缺少什么包就安装什么包

# 以上的包安装方法能够解决绝大多数关于R包安装的问题；但仍不足以 🔑one-fit-all🔒
# 遇到问题，可网上搜索解决办法
#### END ----  
