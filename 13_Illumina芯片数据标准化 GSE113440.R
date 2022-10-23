############################################################################===
##### 24. Illumina芯片数据标准化 GSE113440 -  Summary-level data
############################################################################===


#### GPL10558平台数据集 GSE113440 
# GPL10558, GSE113440 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113440
library(beadarray)
library(tidyverse)

#### 1 GSE113440 BSdata import ----
###  1.1 load data ----
load_input <- load("illumina/GSE113440/GSE113440_idat.Rdata")
load_input

###  1.2 BSdata ----
BSdata <- GSE113440_idat
annotation(BSdata)
det <- calculateDetection(BSdata)#加一Detection列
head (det)
Detection (BSdata) <- det

###  1.3 GSE113440_targets ----
library(GEOquery)
library(tidyverse)

GSE113440_sm <- getGEO(filename = "illumina/GSE113440/GSE113440_series_matrix.txt.gz",
                       getGPL = F) 

GSE113440_pd <- pData(GSE113440_sm) 

GSE113440_targets <- GSE113440_pd %>%
  dplyr::select(2, 1) %>%
  separate(title, into = paste0("x", 1:2), sep = " ") %>%
  mutate(sample_id = geo_accession,
         sample_name = str_c(x1, x2, sep = "_"),
         group = x1,
         rep = x2) %>%
  dplyr::select(sample_id:rep)

#### 3 Plots before QC ----
boxplot(log2(exprs(BSdata)),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
boxplot (nObservations (BSdata),
         ylab = " number of beads ",
         las = 2,
         outline = FALSE)
PCA_new(log2(exprs(BSdata)), 
        group = GSE113440_targets$group,
        show_name = T)

#### 4 QC using arrayQualityMrtics ----
library(arrayQualityMetrics)
pData(BSdata) <- data.frame(pData(BSdata), GSE113440_targets)
arrayQualityMetrics(BSdata, 
                    outdir = "illumina/GSE113440/GSE113440_QC_aqm", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "group")
dev.off()

#### 6 data normalization ----
library(vsn)
GSE113440_norm <- beadarray::normaliseIllumina(BSData = BSdata,
                                               method = "quantile",
                                               transform = "log2")

GSE113440_expr_norm <- exprs(GSE113440_norm)

#### 7 Plots after normalization ----
## boxplot after normalization
boxplot(GSE113440_expr_norm,
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)

## PCA after normalization
PCA_new(GSE113440_expr_norm, group = GSE113440_targets$group)

##  GSE113440_sm
boxplot(exprs(GSE113440_sm),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
PCA_new(exprs(GSE113440_sm), group = GSE113440_targets$group)

####  save data
GSE113440_BSdata <- BSdata
ls(pattern = "^GSE113440")
save(list = ls(pattern = "^GSE113440")[c(-2, -4, -7)], 
     file = "illumina/GSE113440/GSE113440_beadarray_processed.Rda")

#### End ----

