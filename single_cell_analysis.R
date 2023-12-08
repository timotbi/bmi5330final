library(Seurat)
library(dplyr)

samplenames <- c("AML1012-D0", "AML210A-D0", "AML419A-D0", "AML916-D0", "AML314-D0", "AML371-D0", "AML475-D0", "AML722B-D0", "AML870-D0", "AML997-D0", "AML329-D0", 
                 "AML420B-D0", "AML556-D0", "AML328-D0", "AML707B-D0")
filepath <- "./data/GSE116256_RAW/"
data_mtx <- list()
md_mtx <- list()

for(i in 1:length(samplenames)){
  print(samplenames[i])
  sample_files <- list.files(path = filepath, pattern = samplenames[i])
  data_mtx[[i]] <- t(read.table(paste0(filepath, sample_files[grep(pattern = "dem", sample_files)]), header = T, row.names = 1))
  md_mtx[[i]] <- read.table(paste0(filepath, sample_files[grep(pattern = "anno", sample_files)]), header = T, row.names = 1, sep = "\t")[, 1:27]
}

data_mtx <- do.call(rbind, data_mtx)
md_mtx <- do.call(rbind, md_mtx)

aml_object <- CreateSeuratObject(counts = t(data_mtx))
aml_object <- AddMetaData(aml_object, md_mtx)

aml_object[["percent.mt"]] <- PercentageFeatureSet(aml_object, pattern = "^MT-")
VlnPlot(aml_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(aml_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(aml_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

aml_object <- SCTransform(aml_object, method = "glmGamPoi", vars.to.regress = "percent.mt") %>% 
  RunPCA(reduction.name = "pca",assay = 'SCT')
aml_object <- RunHarmony(aml_object,group.by.vars = c("orig.ident"), reduction = 'pca', assay.use = 'SCT',
                   project.dim = FALSE,  reduction.save = "harmony_SCT", plot_convergence = T)
aml_object <- RunUMAP(aml_object, reduction = "harmony_SCT", dims = 1:14, assay = 'SCT', reduction.name = 'hsct.umap', reduction.key = 'hsctUMAP_')
aml_object <- RunUMAP(aml_object, reduction = "pca", dims = 1:14, assay = "SCT", reduction.name = "sct.umap", reduction.key = "sctUMAP_")
aml_object <- FindNeighbors(aml_object, reduction = "harmony_SCT")
aml_object <- FindClusters(aml_object, resolution = 1)
DimPlot(aml_object, reduction = "hsct.umap", group.by = "PredictionRefined", cols = c("red", "green", "grey"))
DimPlot(aml_object, reduction = "hsct.umap", label = T)
DimPlot(aml_object, reduction = "hsct.umap", group.by = "CellType")

aml_object <- NormalizeData(aml_object)
aml_object_malignant <- subset(aml_object, PredictionRefined == "malignant")

DimPlot(aml_object_malignant, reduction = "hsct.umap", group.by = "CellType")

Idents(aml_object_malignant) <- "SCT_snn_res.1"

# 556 is NPM1 only, 210A is both, 328 is FLT3 only, and 475 is neither
markers <- list()
markers[[1]] <- FindMarkers(aml_object_malignant, ident.1 = "AML556.D0", ident.2 = "AML210A.D0")
markers[[2]] <- FindMarkers(aml_object_malignant, ident.1 = "AML556.D0", ident.2 = "AML328.D0")
markers[[3]] <- FindMarkers(aml_object_malignant, ident.1 = "AML556.D0", ident.2 = "AML475.D0")
markers[[4]] <- FindMarkers(aml_object_malignant, ident.1 = "AML210A.D0", ident.2 = "AML328.D0")
markers[[5]] <- FindMarkers(aml_object_malignant, ident.1 = "AML210A.D0", ident.2 = "AML475.D0")
markers[[6]] <- FindMarkers(aml_object_malignant, ident.1 = "AML328.D0", ident.2 = "AML475.D0")
var_markers_all <- intersect(intersect(intersect(intersect(intersect(rownames(markers[[1]]), rownames(markers[[2]])),
        rownames(markers[[3]])),
      rownames(markers[[4]])),
    rownames(markers[[5]])),
  rownames(markers[[6]]))
grep("BRD", var_markers_all)

data_clusters <- rbind(prop.table(table(aml_object_malignant$CellType[aml_object_malignant$orig.ident=="AML556.D0"])),
                       prop.table(table(aml_object_malignant$CellType[aml_object_malignant$orig.ident=="AML475.D0"])),
                       prop.table(table(aml_object_malignant$CellType[aml_object_malignant$orig.ident=="AML328.D0"])),
                       prop.table(table(aml_object_malignant$CellType[aml_object_malignant$orig.ident=="AML210A.D0"])))
data_clusters[2, 6] <- 0
rownames(data_clusters) <- c("AML556.D0", "AML475.D0", "AML328.D0", "AML210A.D0")

heatmap(t(data_clusters), cexRow = 1, cexCol = 1)

aml_object_interest <- subset(aml_object_malignant, orig.ident == "AML556.D0" | orig.ident == "AML210A.D0" |
                                orig.ident == "AML328.D0" | orig.ident == "AML475.D0")
Idents(aml_object_interest) <- "orig.ident"
VlnPlot(aml_object_interest, features = c("FLT3", "NPM1", "S100A9", "S100A8", "IFI16", "JUN", "ATF4", "FOS"), pt.size = 0, ncol = 4)
