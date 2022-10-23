############################################################################===
##### 21. Affymetrix芯片数据标准化 - affy读取
############################################################################===


#### GPL570平台数据集 GSE29450

#### 1 加载Affymetrix芯片GSE29450数据 ----
## load data
load_input <- load("affymetrix/GSE29450/GSE29450_cel.Rdata")
load_input

library(affy)
GSE29450_expr_cel <- exprs(GSE29450_cel)#提取表达矩阵
head(GSE29450_expr_cel)[, 1:5]

GSE29450_targets <- pData(GSE29450_cel) #提取targets（样本分组）数据，在phenoData中的data

#### 2 芯片数据质量评估 基于affy/simpleaffy ----
  ###   2.1 查看数据----
affy::pm(GSE29450_cel)[1:5, 1:5]# perfect matches 
affy::mm(GSE29450_cel)[1:5, 1:5]# mismatches 这是一个背景数据，没有完全匹配的数据
GSE29450_expr_cel[1:5, 1:5]

affy::pm(GSE29450_cel)[1, 1:5]#369707这是id号
GSE29450_expr_cel["369707", 1:5]#定位到369707这一行，发现两者是一致的

## 查看某个探针集的表达情况 
library(tidyverse)#%>% 管道符就是tidyverse里的

GSE29450_pm_1007 <- affy::pm(GSE29450_cel, "1007_s_at") %>% 
  as.data.frame() %>% 
  rownames_to_column("probe_id") %>% 
  pivot_longer(cols = starts_with("GSM"),
               names_to = "sample_id",
               values_to = "pm_value") %>% 
  left_join(GSE29450_targets %>% dplyr::select(sample_id, group), 
            by = "sample_id") %>% 
  mutate(probe_id = factor(probe_id, levels = c(paste0("1007_s_at", 1:16))))

ggplot(GSE29450_pm_1007, aes(x = probe_id, y = pm_value, col = group)) +
  geom_point() +
  labs(x = "Probe set 1007_s_at", y = "Intensity of PM values") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  ###   2.2 CELs to images ----
image(GSE29450_cel[, 6], main = GSE29450_targets$sample_id[6])
image(GSE29450_cel[, 12], main = GSE29450_targets$sample_id[12])

dir.create("affymetrix/GSE29450/GSE29450_plots/")
for (i in 1:20) {
  tiff(filename = paste0("affymetrix/GSE29450/GSE29450_plots/",
                         GSE29450_targets$sample_id[i], "_", i,".tif"),
       height = 10,
       width = 10,
       units = "in",
       res = 300)
  image(GSE29450_cel[, i], main = GSE29450_targets$sample_id[i])
  dev.off()
}

# 芯片扫描日期可作为一种重要的批次效应因素
GSE29450_targets$batch <- GSE29450_cel@protocolData@data$ScanDate %>% 
  str_sub(1, 8)#提取1到8个字符

table(GSE29450_targets$batch)
# #03/12/04 04/09/04 05/25/06 09/30/04 11/02/06 
# 1        3        6        6        4 
#20个样本分了5批，还是有明显的批次效应的，这种数据集质量堪忧

  ###   2.3 RNA degradation  ---- 
  #评估RNA的讲解程度
RNAdeg <- AffyRNAdeg(GSE29450_cel)
summaryAffyRNAdeg(RNAdeg)
#RNA降解从5‘端开始降解，所5'信号值要比3'信号值低一些，所以假如低得多了，斜率就大了，说明降解得很厉害



#  RNA degradation plot
cols <- rainbow(nrow(GSE29450_targets))
plotAffyRNAdeg(RNAdeg, cols = cols)
legend("topleft",
       ncol = 2,
       legend = sampleNames(GSE29450_cel), 
       lty = 1, 
       lwd = 2,
       cex = 0.8,
       box.lty=0,
       bg = "transparent",
       col = cols)
box()#这里画图出了问题，legend的放置位置错了

  ###   2.4 simpleaffy包的质量评估函数qc ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("simpleaffy")

library(simpleaffy)
GSE29450_qc <- simpleaffy::qc(GSE29450_cel)
plot(GSE29450_qc)
# P/M/A:present/ marginal present/ absent
# actin 3/5 - 3, gapdh 3/5 - 1


  ###   2.5 boxplot and histogram ----
#箱线图
boxplot(GSE29450_cel, las = 2, col = rep(c("blue", "red"), each = 10))

#密度图
hist(GSE29450_cel, lty = 1:3, col = cols)
legend("topright",
       legend = sampleNames(GSE29450_cel),
       lty = 1:3,
       cex = 0.8,
       col = cols,
       box.col = "transparent",
       xpd = TRUE)
box()

  ###   2.6 PCA_new ----
## new function for PCA 
#参考了Deseq2包来写的PCA函数
PCA_new <- function(expr, ntop = 500, group, show_name = F){
  library(ggplot2)
  library(ggrepel)
  rv <- genefilter::rowVars(expr)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(expr[select, ]))#最核心的代码
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[, 1], 
                  PC2 = pca$x[, 2], 
                  group = group, 
                  name = colnames(expr))
  attr(d, "percentVar") <- percentVar[1:2]
  if (show_name) {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      geom_text_repel(aes(label = name),
                      size = 3,
                      segment.color = "black",
                      show.legend = FALSE )
  } else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2",color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
  }
}

PCA_new(log2(GSE29450_expr_cel), 
        nrow(GSE29450_expr_cel), 
        group = GSE29450_targets$group,
        show_name = T)
#很明显的看到分成三批，正常来讲normal和cancer不会离得这么远，所以更多的是批次效应而不是生物学差异

#### 3 芯片数据质量评估 基于arrayQulitMetrics ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("arrayQualityMetrics")

library(arrayQualityMetrics)
arrayQualityMetrics(GSE29450_cel, 
                    outdir = "affymetrix/GSE29450/GSE29450_QC_raw", #通过该文件夹下的index来看结果
                    force = TRUE,#如果本身有个文件夹的话就把文件夹刷新了
                    intgroup = "group",#GSE29450_cel@phenoData@data[["group"]]
                    do.logtransform = TRUE)
dev.off()#dev.off报错了也无妨，就是为了让它报错

#### 4 背景校正 ----
GSE29450_bgc <- affy::bg.correct(GSE29450_cel, method = "rma")#最常用rma方法
bgcorrect.methods()

#背景矫正后作图
boxplot(GSE29450_bgc, las = 2, col = rep(c("blue", "red"), each = 10))


PCA_new(log2(exprs(GSE29450_bgc)), 
        nrow(GSE29450_bgc), 
        group = GSE29450_targets$group,
        show_name = T)
#可以看到PCA图还是没什么改善

  #### 5 芯片数据标准化 ----
  ###对于affy芯片来说，rma或者gcrma已经把背景矫正整合进去了，所以没有必要单独做背景矫正

  ###   5.1 RMA ----
GSE29450_rma <- affy::rma(GSE29450_cel)#用的是原始cel数据，不是背景矫正后的数据

boxplot(GSE29450_rma, las = 2, col = rep(c("blue", "red"), each = 10))
#看上去中位线就在一条线上了，不过依然可以看出异常值，比如第6个，之前就看出它降解最大

PCA_new(exprs(GSE29450_rma), 
        nrow(GSE29450_rma), 
        group = GSE29450_targets$group)

#芯片数据标准化之后看数据质量是否改善
library(arrayQualityMetrics)
arrayQualityMetrics(GSE29450_rma, 
                    outdir = "affymetrix/GSE29450/GSE29450_QC_rma", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = TRUE)

  ###   5.2 GCRMA ----
library(gcrma)
GSE29450_gcrma <- gcrma(GSE29450_cel)
GSE29450_expr_gcrma <- exprs(GSE29450_gcrma)

boxplot(GSE29450_gcrma, las = 2, col = rep(c("blue", "red"), each = 10))


# dev.off()
PCA_new(exprs(GSE29450_gcrma), 
        nrow(GSE29450_gcrma), 
        group = GSE29450_targets$group,
        show_name = T)
#可以看到标准化后PCA图没有明显的改善，说明还是批次效应造成的；还可以看到gsm729039这个明显的离群值
#PCA图和箱线图两个图之间得到相互印证

library(arrayQualityMetrics)
arrayQualityMetrics(GSE29450_gcrma, 
                    outdir = "affymetrix/GSE29450/GSE29450_QC_gcrma", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = TRUE)
#arrayQualityMetrics也显示gsm729039这个样本有问题




#### 6 缺失值补充 ----
#首先看一下有没有缺失值
sum(is.na(exprs(GSE29450_gcrma)))

## 模拟缺失值
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

library(impute)

test <- GSE29450_expr_gcrma

test[1:10, 5] <- NA#人工造出NA哈哈
sum(is.na(test))

test2 <- impute.knn(test, maxp = 30000)#maxp默认的是1500，好像会经常报错
sum(is.na(test2))

#比较一下原始的和补充之后的数值
GSE29450_expr_gcrma[1:10, 5]
test2$data[1:10, 5]

#### 7 标准化结果可视化 ----
##  已介绍，用箱线图和PCA图

#### 8 该数据集的GEO series matrix 数据 ----
library(GEOquery)
GSE29450_sm <- getGEO(filename = "affymetrix/GSE29450/GSE29450_series_matrix.txt.gz", 
                      getGPL = F)
GSE29450_expr_sm <- exprs(GSE29450_sm)#提取表达矩阵
GSE29450_pd <- pData(GSE29450_sm)#提取临床信息

#看数据的处理方式
GSE24950_pd$data_processing[1]

boxplot(GSE29450_sm, las = 2, col = rep(c("blue", "red"), each = 10))

PCA_new(exprs(GSE29450_sm), 
        nrow(GSE29450_sm),
        group = GSE29450_targets$group)

#### 9 save data ----

save(GSE29450_expr_gcrma, GSE29450_targets, 
     file = "affymetrix/GSE29450/GSE29450_after_gcrma.Rdata")

#### End ----
