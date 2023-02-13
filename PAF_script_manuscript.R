#PAF analysis - general clustering and DEG analysis 

#load in the libraries 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(SeuratDisk)
library(SeuratData)
library(RColorBrewer)
library(ArchR)
library(pheatmap)
library(Dune)
library(MAST)

#load in filtered matrices - the seurat object 

#samples included: P5, P6, P7, P8, P11, P13, P15, P16, P17, P18, P19, P20, P21



#------------begin to assess QC metrics----------------#

#add mito percent 
pf[["percent.mt"]] <- PercentageFeatureSet(pf, pattern = "^MT-")

#add "complexity" score
pf$log10GenesPerUMI <- log10(pf$nFeature_RNA) / log10(pf$nCount_RNA)

#visualize QC metrics by creating violin plots per sample 
VlnPlot(pf, features = c("nFeature_RNA"), pt.size = 0)
VlnPlot(pf, features = c("nCount_RNA"), pt.size = 0)
VlnPlot(pf, features = c("percent.mt"), pt.size = 0)

#create scatter plot comparing mito percent vs. gene 
FeatureScatter(pf, feature1 = 'nFeature_RNA', feature2 = 'percent.mt')

#cluster prior to QC - assess clustering for batch or sample specific clusters 
#follow Seurat's Pipeline


pf <- NormalizeData(pf)
pf <- FindVariableFeatures(pf, selection.method = "vst", nfeatures = 2000)
pf <- ScaleData(pf)
pf <- RunPCA(pf, features = VariableFeatures(object = pf))
#PCA Plot of data
DimPlot(pf, reduction = "pca")
#elbow plot to determine PC cutoff 
ElbowPlot(pf)
pf <- FindNeighbors(pf, dims = 1:15)
pf <- FindClusters(pf, resolution = 0.5)
pf <- RunUMAP(pf, dims = 1:15)
DimPlot(pf, reduction = "umap")

#----------begin to perform QC--------------#

#view "complexity" violin plot
VlnPlot(pf, features = "log10GenesPerUMI", pt.size = 0) + geom_hline(yintercept = 0.8)

#cell level filtering 
pf_filt <- subset(pf, subset = nCount_RNA >= 300 & nFeature_RNA >= 250 & nFeature_RNA <= 6500 & log10GenesPerUMI > 0.80 & percent.mt < 30)

#gene level filtering 

#output a logical vector for every gene on whether the more than 0 counts per cell and extract counts
counts <- GetAssayData(object = pf_filt, slot = "counts")

#output a logical vector for every gene on whether the more than zero counts per cell 
nonzero <- counts > 0

#sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

#only keeping those genes expressed in more than 10 cells 
filtered_counts <- counts[keep_genes, ]

#reassign to filtered seurat object 
pf_gfilt <- CreateSeuratObject(filtered_counts, meta.data = pf_filt@meta.data)

#assess distribution of data post QC 
VlnPlot(pf_gfilt, features = c("nFeature_RNA"), pt.size = 0)
VlnPlot(pf_gfilt, features = c("nCount_RNA"), pt.size = 0)
VlnPlot(pf_gift, features = c("percent.mt"), pt.size = 0)
FeatureScatter(pf_gfilt, feature1 = 'nFeature_RNA', feature2 = 'percent.mt')

#---------------perform clustering post QC----------------#
#follow Seurat's pipeline - determine if there are batch specific or sample specific clusters 
pf_gfilt <- NormalizeData(pf_gfilt)
pf_gfilt <- FindVariableFeatures(pf_gfilt, selection.method = "vst", nfeatures = 2000)
pf_gfilt <- ScaleData(pf_gfilt)
pf_gfilt <- RunPCA(pf_gfilt, features = VariableFeatures(object = pf_gfilt))
DimPlot(pf_gfilt, reduction = "pca")
ElbowPlot(pf_gfilt)
pf_gfilt <- FindNeighbors(pf_gfilt, dims = 1:15)
pf_gfilt <- FindClusters(pf_gfilt, resolution = 0.5)
pf_gfilt <- RunUMAP(pf_gfilt, dims = 1:15)
DimPlot(pf_gfilt, reduction = "umap")

#--------------Perform Batch effect correction - 3 methods: Harmony, CCA, and rPCA------------------#

#perform Harmony Batch Efefct Correction from Harmony Pipeline
pf_h <- pf_gfilt %>%
  RunHarmony("batch")

harmony_embeddings <- Embeddings(pf_h, 'harmony')

pf_h <- pf_h %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
#visualize clusters (can also visualize based on different metadata criteria)
DimPlot(pf_h, reduction = "umap")

#save "Harmony" Seurat Object - pf_h

#perform CCA Batch Effect Correction from Seurat Pipeline

#split the dataset into a list of two seurat objects; integrate across batch
pf.list <- SplitObject(pf_gfilt, split.by = "batch")

#normalize and identify variable features for each dataset independently 
pf.list <- lapply(X = pf.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
})

#select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = pf.list)

pf.anchors <- FindIntegrationAnchors(object.list = pf.list, anchor.features = features)

#create an integrated data assay 
pf.combined <- IntegrateData(anchorset = pf.anchors)

#perform integrated analysis

#create an integrated assay 
DefaultAssay(pf.combined) <- "integrated"

#standard workflow for visualization and clustering 
pf.combined <- ScaleData(pf.combined, verbose = FALSE)
pf.combined <- RunPCA(pf.combined, npcs = 15)
pf.combined <- RunUMAP(pf.combined, reduction = "pca", dims = 1:15)
pf.combined <- FindNeighbors(pf.combined, reduction = "pca", dims = 1:15)
pf.combined <- FindClusters(pf.combined, resolution = 0.5)
#visualize clusters (can also visualize based on different metadata criteria)
DimPlot(pf.combined, reduction = "umap")

#save "CCA" Seurat Object - pf_c

#perform rPCA Batch Effect Correction from Seurat Pipeline 

#split the dataset into a list of two seurat objects and integrate across batch 
pf.list <- SplitObject(pf_gfilt, split.by = "batch")

#normalize and identify variable features for each dataset independently 
pf.list <- lapply(X = pf.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
})

#select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = pf.list)
pf.list <- lapply(X = pf.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})  


pf.anchors <- FindIntegrationAnchors(object.list = pf.list, anchor.features = features, reduction = "rpca")

#create an integrated data assay 
pf.combined <- IntegrateData(anchorset = pf.anchors)

#perform integrated analysis

#create an integrated assay 
DefaultAssay(pf.combined) <- "integrated"

#standard workflow for visualization and clustering 
pf.combined <- ScaleData(pf.combined, verbose = FALSE)
pf.combined <- RunPCA(pf.combined, npcs = 15)
pf.combined <- RunUMAP(pf.combined, reduction = "pca", dims = 1:15)
pf.combined <- FindNeighbors(pf.combined, reduction = "pca", dims = 1:15)
pf.combined <- FindClusters(pf.combined, resolution = 0.5)
#visualize clusters (can also visualize based on different metadata criteria)
DimPlot(pf.combined, reduction = "umap")

#save "rPCA" Seurat object - pf_r

#----------------------Compute Rand Index to determine clustering stability---------------------------#

#load in three seurat objects that are batch-effect corrected: pf_h, pf_c, pf_r

#create dataframe 
df_pfc <- data.frame(paste0(pf_c$seurat_clusters))
df_pfr <- data.frame(paste0(pf_r$seurat_clusters))
df_pfh <- data.frame(paste0(pf_h$seurat_clusters))

#merge the dataframes 
ri <- cbind(df_pfr, df_pfh, df_pfc)

#compute the rand index using dune 
plotARIs(ri %>% select(paste0.pf_c.seurat_clusters., paste0.pf_h.seurat_clusters., paste0.pf_r.seurat_clusters.))

#compute the rand index using the merged clusters
merger <- Dune(clusMat = ri %>% select(paste0.pf_c.seurat_clusters., paste0.pf_h.seurat_clusters., paste0.pf_r.seurat_clusters.))
names(merger)

#ARI for the merged clusters
plotARIs(clusMat = merger$currentMat)


#--------------------Create Confusion Matrix Heatmaps to Assess Clustering Similarity-----------------#

cM_rpca_cca <- confusionMatrix(paste0(pf_r$seurat_clusters), paste0(pf_c$seurat_clusters))
cM_rpca_cca <- cM_rpca_cca / Matrix::rowSums(cM_rpca_cca)
cM_rpca_cca_plot <- pheatmap::pheatmap(
  mat = as.matrix(cM_rpca_cca), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)

cM_rpca_harmony <- confusionMatrix(paste0(pf_r$seurat_clusters), paste0(pf_h$seurat_clusters))
cM_rpca_harmony <- cM_rpca_harmony / Matrix::rowSums(cM_rpca_harmony)
cM_rpca_harmony_plot <- pheatmap::pheatmap(
  mat = as.matrix(cM_rpca_harmony), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)

cM_cca_harmony <- confusionMatrix(paste0(pf_c$seurat_clusters), paste0(pf_h$seurat_clusters))
cM_cca_harmony <- cM_cca_harmony / Matrix::rowSums(cM_cca_harmony)
cM_cca_harmony_plot <- pheatmap::pheatmap(
  mat = as.matrix(cM_cca_harmony), 
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)

#-------------------Annotate Clusters using Marker genes from Literature-------------------#
#use rPCA Batch Effect Corrected Data for annotating clusters

#separate epithelial and immune cells based on annotations

#repeat rPCA batch effect correction (same as above) for epithelial and immune cells 

#annotate epithelial and immune cells using marker genes from literature 

#-----------------Perform DEG Analysis Using MAST------------------#

#prior to DEG analysis, remove samples with low number of cells (also removed the non-IBD control)

#example for goblet cells; comparing inflamed vs. non-inflamed goblet cells 

#create the MAST function

MAST_runfun <- function(object,assay = "SCT",slot = "data", features = NULL, ident.1 = NULL, ident.2 = NULL, cov = NULL,min.cells.group = 3, reduction = NULL,pseudocount.use = 1){
  ## cov: selected covariance variables for fitting in the hurdle model, which need to have the same name as 
  ##      in Seurat object meta data slot.
  
  ## Getting the cell ids in ident.1 and ident.2,respectively
  cells.1 = row.names(object@meta.data)[which(Idents(object) %in% ident.1)]
  cells.2 = row.names(object@meta.data)[which(Idents(object) %in% ident.2)]
  if(is.null(features)){
    features = row.names(object)
  }
  ## error checking
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
  
  
  ## constructing group information datasframe
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  
  ## Using objectled nCount_RNA as default covariance 
  group.info$cntscl <- scale(object@meta.data[row.names(group.info),"nCount_RNA"])
  
  ## Adding customed covariables
  if(!is.null(cov)){
    group.info <- cbind(group.info, object@meta.data[row.names(group.info),cov])
    colnames(group.info) <- c("group","cntscl",cov)
  }else{
    colnames(group.info) <- c("group","cntscl")
  }
  
  ## Constructing single cell dataset
  fit_data <- MAST::FromMatrix(
    exprsArray = as.matrix(object[[assay]]@data[features,row.names(group.info)]),
    cData = group.info,
    fData = as.data.frame(as.matrix(features,nrow = 1))
  )
  
  data = as.matrix(object[[assay]]@data)
  
  thresh.min = 0
  pct.1 <- round(
    rowSums(data[features, cells.1, drop = FALSE] > thresh.min) /
      length(cells.1),
    digits = 3
  )
  
  pct.2 <- round(
    rowSums(data[features, cells.2, drop = FALSE] > thresh.min) /
      length(cells.2),
    digits = 3
  )
  
  # feature selection (based on average difference)
  # if selecting data slot is scale.data, then calculating the rowMeans as avg gene expression
  # if selecting data slot is data, then using formula log(x = rowMeans(x = expm1(x = x)) + pseudocount.use)
  # if selecting data slot is count, then using formula log(x = rowMeans(x = x) + pseudocount.use)
  mean.fxn <- if (is.null(x = reduction) && slot != "scale.data") {
    switch(
      EXPR = slot,
      'data' = function(x) {
        return(log2(x = rowMeans(x = expm1(x = x)) + pseudocount.use))
      },
      function(x) {
        return(log2(x = rowMeans(x = x) + pseudocount.use))
      }
    )
  } else {
    rowMeans
  }
  
  ## Calculating the avg_logFC between two groups
  data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(data[features, cells.2, drop = FALSE])
  total.diff <- (data.1 - data.2)
  
  ## Constructing zlm model and using Group1 as reference group
  cond<-factor(colData(fit_data)$group)
  cond<-relevel(cond,"Group1")
  colData(fit_data)$group<-cond
  fmla <- as.formula(
    object = paste0(" ~ ", paste(colnames(group.info), collapse = "+"))
  )
  zlmCond <- MAST::zlm(formula = fmla, fit_data)
  
  ## Running a likelihood ratio test here, testing for differences when we drop the group factor
  ## only test the group coefficient.
  summaryCond <- summary(object = zlmCond, doLRT = 'groupGroup2')
  summaryDt <- summaryCond$datatable
  # fcHurdle <- merge(
  #   summaryDt[contrast=='groupGroup2' & component=='H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
  #   summaryDt[contrast=='groupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid'
  # ) #logFC coefficients
  # fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle <- merge(summaryDt[contrast=='groupGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='groupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
  fcHurdle = fcHurdle[features,]
  p_val <- fcHurdle[,`Pr(>Chisq)`]
  p_val_adj = p.adjust(
    p = p_val,
    method = "bonferroni",
    n = nrow(object)
  )
  Coef <- fcHurdle[,`coef`]
  genes.return <- fcHurdle[,`primerid`]
  
  to.return <- data.frame(pct.1 = pct.1, pct.2 = pct.2, Coef = Coef, avg_logFC = total.diff, p_val = p_val, p_val_adj = p_val_adj, row.names = genes.return)
  return(to.return)
}

#load in the goblet cells; Seurat object = "goblet"

#check the idents - ident based on inflammation status 
Idents(goblet) <- goblet$inflam_status

#run MAST 
goblet_mast <- MAST_runfun(goblet, assay = "RNA", ident.1 = "inflamed", ident.2 = "non_inflam", cov = c("sex", "ancestry", "age"))

#save the output
write.table(goblet_mast, "goblet_MAST_run.txt", sep="\t",quote=F,row.names=T)


