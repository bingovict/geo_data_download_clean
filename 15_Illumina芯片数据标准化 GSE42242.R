############################################################################===
##### 26. Illumina芯片数据标准化 GSE42242 - Summary-level data
############################################################################===


#### GPL10558平台数据集 GSE42242
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42242
# summary-level data from txt files

library(limma)
library(tidyverse) 

#### 1 GSE42242 data import ----
###  1.1 load data ----
load_input <- load("geo/GSE42242/GSE42242_non-normalized.Rdata")
load_input

###  1.2 GSE42242_targets ----
library(GEOquery)

GSE42242_sm <- getGEO(filename = "geo/GSE42242/GSE42242_series_matrix.txt.gz",
                       getGPL = F) 

GSE42242_pd <- pData(GSE42242_sm) 

GSE42242_targets <- GSE42242_pd %>%
  dplyr::select(2, 1) %>%
  separate(title, into = paste0("x", 1:4), sep = " ") %>%
  mutate(sample_id = geo_accession,
         sample_name = paste0(x1, "_", 1:4),
         group = x1) %>%
  dplyr::select(sample_id:group)

#### 2 Plots before QC ----

boxplot(log2(GSE42242_non_norm$E),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
PCA_new(log2(GSE42242_non_norm$E), 
        ntop = nrow(GSE42242_non_norm$E),
        group = GSE42242_targets$group,
        show_name = F)

#### 3 QC using arrayQualityMrtics ----
library(arrayQualityMetrics)

rownames(GSE42242_targets) <- colnames(GSE42242_non_norm$E) <- GSE42242_targets$sample_id

eset <- ExpressionSet(assayData = GSE42242_non_norm$E, 
                      phenoData = AnnotatedDataFrame(data = GSE42242_targets))

arrayQualityMetrics(eset, 
                    outdir = "geo/GSE42242/GSE42242_QC_aqm", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()

#### 4 data normalization ----
GSE42242_neqc <- neqc(GSE42242_non_norm, #neqc专门为illumina设计，又做背景矫正，又做标准化
                      detection.p="Detection")

GSE42242_expr_neqc <- GSE42242_neqc$E

sum(is.na(GSE42242_expr_neqc))

#### 5 Plots after normalization ----
## boxplot after normalization
boxplot(GSE42242_expr_neqc,
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)

## PCA after normalization
PCA_new(GSE42242_expr_neqc, 
        ntop = nrow(GSE42242_expr_neqc),
        group = GSE42242_targets$group)

##  GSE42242_sm
boxplot(exprs(GSE42242_sm),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
PCA_new(exprs(GSE42242_sm), 
        ntop = nrow(exprs(GSE42242_sm)),
        group = GSE42242_targets$group)

####  save data
ls(pattern = "^GSE42242")
save(list = ls(pattern = "^GSE42242")[-3], 
     file = "illumina/GSE42242/GSE42242_neqc_processed.Rda")

#### End ----

