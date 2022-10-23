############################################################################===
##### 22. Affymetrix芯片数据标准化 - oligo读取
############################################################################===


#### GPL6244平台数据集 GSE24129 

#### 1 加载Affymetrix芯片数据2 ----
## load data
load_input <- load("affymetrix/GSE24129/GSE24129_cel.Rdata")
load_input

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligo")

library(oligo)
GSE24129_expr_cel <- exprs(GSE24129_cel)
head(GSE24129_expr_cel)[, 1:5]
GSE24129_targets <- pData(GSE24129_cel) 

#### 2 芯片数据质量评估 基于oligo ----
###  2.1 查看数据----
oligo::pm(GSE24129_cel)[1:5, 1:5]
GSE24129_expr_cel[1:5, 1:5]

oligo::pm(GSE24129_cel)[1, 1:5]
GSE24129_expr_cel["1056", 1:5]

###  2.2 CEL to images ----
oligo::image(GSE24129_cel[, 6], main = GSE24129_targets$sample_id[6])
oligo::image(GSE24129_cel[, 12], main = GSE24129_targets$sample_id[12])

# dir.create("affymetrix/GSE24129/GSE24129_plots/")
# for (i in 1:24) {
#   tiff(filename = paste0("affymetrix/GSE24129/GSE24129_plots/",
#                          GSE24129_targets$sample_id[i], "_", i,".tif"),
#        height = 10,
#        width = 10,
#        units = "in",
#        res = 300)
#   image(GSE24129_cel[, i], main = GSE24129_targets$sample_id[i])
#   dev.off()
# }

# 芯片扫描日期可作为一种重要的批次效应因素
library(tidyverse)
library(dplyr)
GSE24129_targets$batch <- GSE24129_cel@protocolData@data$dates %>% 
  as.character() %>% 
  str_sub(1, 10)

table(GSE24129_targets$batch)

###  2.3 boxplot  ggplot ----

box_data <-  log2(GSE24129_expr_cel) %>% #对表达矩阵log2转换
  as.data.frame() %>% #tideverse包不接受矩阵，只接受数据框
  rownames_to_column("probe_id") %>% #行名变列，这一列叫probe_id
  slice_sample(n = 50000) %>%  #因为运行速度很慢，所以只挑了50000个数据来运行
  pivot_longer(cols = starts_with("GSM"),#对数据进行长宽转换，从宽数据到长数据，转换以gsm开头的列
               names_to = "sample_id",#从gsm变成sample_id
               values_to = "value") %>% #新建一列叫value
  left_join(GSE24129_targets %>% dplyr::select(sample_id, group), #left_join相当于取交集，以左边的为主
            by = "sample_id") #以sample_id为锚点取交集，多了一列就是group
 
#使用ggplot来作图，aes确定映射，col分组作为颜色
ggplot(data = box_data, aes(x = sample_id, y = value, col = group)) +
  geom_boxplot() + 
  ylab("log2 value") +
  xlab("samples") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))#X轴左边做90度旋转
#可以看到ggplot2的图颜值要高于boxplot哈哈哈


#密度图
ggplot(data = box_data, aes(x = value, col = sample_id)) +
  geom_density()
#可以看到有个也是突出来的

dev.off()

###  2.4 PCA_new ----
## new function for PCA 
library(tidyverse)

PCA_new(log2(GSE24129_expr_cel), 
        ntop = 5000,
        group = GSE24129_targets$group,
        show_name = T)


#### 3 芯片数据质量评估 基于arrayQulitMetrics ----
library(arrayQualityMetrics)
arrayQualityMetrics(GSE24129_cel, 
                    outdir = "affymetrix/GSE24129/GSE24129_QC_raw", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = TRUE)
dev.off()

#### 4 背景校正 ----
backgroundCorrectionMethods()#针对oligo包读取有三种背景校正方法
# LESN - Low End Signal is Noise Background corrections 很少用LESN这种方法
# affyPLM::bg.correct.LESN()

GSE24129_bgc <- oligo::backgroundCorrect(GSE24129_cel, method = "rma")
#there is no package called ‘pd.hugene.1.0.st.v1’
BiocManager::install("pd.hugene.1.0.st.v1")


PCA_new(log2(exprs(GSE24129_bgc)), 
        nrow(GSE24129_bgc), 
        group = GSE24129_targets$group,
        show_name = T)
#背景矫正后就完全把gse59052完全孤立出去了

boxplot(log2(exprs(GSE24129_bgc)), 
              las = 2, 
              outline = FALSE,
              col = as.factor(GSE24129_targets$group))

#### RMA标准化
GSE24129_rma <- rma(GSE24129_cel)
#gcRMA不支持oligo读取的数据，gcRMA只支持affybatch


boxplot(exprs(GSE24129_rma), 
              las = 2, 
              outline = FALSE,
              col = as.factor(GSE24129_targets$group))
#可以看到箱线图变好看了

PCA_new(exprs(GSE24129_rma), 
        nrow(GSE24129_rma), 
        group = GSE24129_targets$group,
        show_name = T)

library(arrayQualityMetrics)
arrayQualityMetrics(GSE24129_rma, 
                    outdir = "affymetrix/GSE24129/GSE24129_QC_rma", 
                    force = TRUE,
                    intgroup = "group",
                    do.logtransform = TRUE)
dev.off()

#### 5 缺失值补充 ----
sum(is.na(exprs(GSE24129_rma)))

#### 6 该数据集的GEO series matrrix 数据情况 ----
library(GEOquery)

GSE24129_sm <- getGEO(filename = "affymetrix/GSE24129/GSE24129_series_matrix.txt.gz", 
                      getGPL = F)
GSE24129_expr_sm <- exprs(GSE24129_sm)
GSE24129_pd <- pData(GSE24129_sm)

boxplot(exprs(GSE24129_sm), 
              las = 2, 
              outline = FALSE,
              col = as.factor(GSE24129_targets$group))

PCA_new(exprs(GSE24129_sm), 
        nrow(GSE24129_sm),
        group = GSE24129_targets$group)

#### 7 save data ----
save(GSE24129_rma, GSE24129_targets, 
     file = "affymetrix/GSE24129/GSE24129_after_rma.Rda")

#### End ----
