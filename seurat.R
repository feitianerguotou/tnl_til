library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(MySeuratWrappers)
library(SeuratObject)

#单细胞
TNL_1 <- Read10X("singlecellr/1_Cellranger_result/TNL_1/")
TNL_2 <- Read10X("singlecellr/1_Cellranger_result/TNL_2/")
TNL_3 <- Read10X("singlecellr/1_Cellranger_result/TNL_3/")
TNL_4 <- Read10X("singlecellr/1_Cellranger_result/TNL_4/")
TNL_5 <- Read10X("singlecellr/1_Cellranger_result/TNL_5/")
TIL_1 <- Read10X("singlecellr/1_Cellranger_result/TIL_1/")
TIL_2 <- Read10X("singlecellr/1_Cellranger_result/TIL_2/")
TIL_3 <- Read10X("singlecellr/1_Cellranger_result/TIL_3/")
TIL_4 <- Read10X("singlecellr/1_Cellranger_result/TIL_4/")###
TIL_5 <- Read10X("singlecellr/1_Cellranger_result/TIL_5/")###

TNL_1object <- CreateSeuratObject(counts = TNL_1, project = "TNL_1",min.cells = 3, min.features = 200)
TNL_2object <- CreateSeuratObject(counts = TNL_2, project = "TNL_2",min.cells = 3, min.features = 200)
TNL_3object <- CreateSeuratObject(counts = TNL_3, project = "TNL_3",min.cells = 3, min.features = 200)
TNL_4object <- CreateSeuratObject(counts = TNL_4, project = "TNL_4",min.cells = 3, min.features = 200)
TNL_5object <- CreateSeuratObject(counts = TNL_5, project = "TNL_5",min.cells = 3, min.features = 200)
TIL_1object <- CreateSeuratObject(counts = TIL_1, project = "TIL_1",min.cells = 3, min.features = 200)
TIL_2object <- CreateSeuratObject(counts = TIL_2, project = "TIL_2",min.cells = 3, min.features = 200)
TIL_3object <- CreateSeuratObject(counts = TIL_3, project = "TIL_3",min.cells = 3, min.features = 200)
TIL_4object <- CreateSeuratObject(counts = TIL_4, project = "TIL_4",min.cells = 3, min.features = 200)
TIL_5object <- CreateSeuratObject(counts = TIL_5, project = "TIL_5",min.cells = 3, min.features = 200)

TNL_1object[["percent.mt"]] <- PercentageFeatureSet(TNL_1object, pattern = "^MT-")
TNL_2object[["percent.mt"]] <- PercentageFeatureSet(TNL_2object, pattern = "^MT-")
TNL_3object[["percent.mt"]] <- PercentageFeatureSet(TNL_3object, pattern = "^MT-")
TNL_4object[["percent.mt"]] <- PercentageFeatureSet(TNL_4object, pattern = "^MT-")
TNL_5object[["percent.mt"]] <- PercentageFeatureSet(TNL_5object, pattern = "^MT-")
TIL_1object[["percent.mt"]] <- PercentageFeatureSet(TIL_1object, pattern = "^MT-")
TIL_2object[["percent.mt"]] <- PercentageFeatureSet(TIL_2object, pattern = "^MT-")
TIL_3object[["percent.mt"]] <- PercentageFeatureSet(TIL_3object, pattern = "^MT-")
TIL_4object[["percent.mt"]] <- PercentageFeatureSet(TIL_4object, pattern = "^MT-")
TIL_5object[["percent.mt"]] <- PercentageFeatureSet(TIL_5object, pattern = "^MT-")

VlnPlot(TIL_3object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(TNL_1object, feature1 = "nFeature_RNA", feature2 = "percent.mt")

TNL_1object <- subset(TNL_1object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TNL_2object <- subset(TNL_2object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TNL_3object <- subset(TNL_3object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TNL_4object <- subset(TNL_4object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TNL_5object <- subset(TNL_5object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TIL_1object <- subset(TIL_1object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TIL_2object <- subset(TIL_2object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TIL_3object <- subset(TIL_3object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TIL_4object <- subset(TIL_4object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
TIL_5object <- subset(TIL_5object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)

allobject <- c(TNL_1object,TNL_2object,TNL_3object,TNL_4object,TNL_5object,TIL_1object,TIL_2object,TIL_3object,TIL_4object,TIL_5object)
names(allobject) <- c('TNL_1','TNL_2','TNL_3','TNL_4','TNL_5','TIL_1','TIL_2','TIL_3','TIL_4','TIL_5')

TNL_TIL <- lapply(allobject, function(x) {x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)})

features <- SelectIntegrationFeatures(object.list = TNL_TIL)
anchors <- FindIntegrationAnchors(object.list = TNL_TIL, anchor.features = features)
TNL_TIL <- IntegrateData(anchorset = anchors)

DefaultAssay(TNL_TIL) <- "integrated"
TNL_TIL <- ScaleData(TNL_TIL, verbose = FALSE)
TNL_TIL <- RunPCA(TNL_TIL, npcs = 30, verbose = FALSE)
TNL_TIL <- RunUMAP(TNL_TIL, reduction = "pca", dims = 1:30)
TNL_TIL <- FindNeighbors(TNL_TIL, reduction = "pca", dims = 1:30)
TNL_TIL <- RunTSNE(TNL_TIL, dims = 1:30)

TNL_TIL@meta.data$groups[TNL_TIL@meta.data$orig.ident=="TIL_1"]=c("TIL")

TNL_TIL <- FindClusters(TNL_TIL, resolution = 1)

saveRDS(TNL_TIL, file = "singlecellr/np_tnl/TNL_TIL.rds")
TNL_TIL <- readRDS("singlecellr/til_tnl/til_tnl.rds")#####
saveRDS(TNL_TIL, file = "singlecellr/til_tnl/til_tnl.rds")#####

Idents(TNL_TIL) <- "integrated_snn_res.1"
p1 <- DimPlot(TNL_TIL, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(TNL_TIL, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL, reduction = "umap", split.by = "orig.ident")
DimPlot(TNL_TIL, cells.highlight = CellsByIdentities(object = TNL_TIL, idents = c(40)), reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(TNL_TIL, cells.highlight = CellsByIdentities(object = TNL_TIL, idents = c(33)), reduction = "tsne", label = TRUE, repel = TRUE)
#质检图
VlnPlot(TNL_TIL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
Idents(TNL_TIL) <- "groups"
FeatureScatter(TNL_TIL, feature1 = "nFeature_RNA", feature2 = "percent.mt")

table(Idents(TNL_TIL))
table(TNL_TIL$orig.ident)
table(TNL_TIL$celltype)

table(Idents(TNL_TIL), TNL_TIL$groups)
prop.table(table(Idents(TNL_TIL), TNL_TIL$groups))
cellnumber <- prop.table(table(Idents(TNL_TIL), TNL_TIL$orig.ident))
write.csv(cellnumber,file="singlecellr/til_tnl/cellnumber.csv")

genes_to_check3 = c('HLA-G', 'KRT18', 'KRT8', 'GATA3','PAPPA2', 'KRT7')
genes_to_check3 = c('CD163', 'CD86', 'CD14', 'PTPRC','IL1B', 'MRC1','S100A9','S100A8')#巨噬
genes_to_check3 = c('CD3E', 'CD3D', 'CD7', 'KLRD1','KLRB1','CD4','CD8A','CD8B')#T NK
genes_to_check3 = c('CD79A', 'MS4A1', 'IGHG1', 'MZB1','RCVRN','PTPRC', 'HBB','JCHAIN')#B 红细胞
genes_to_check3 = c('HLA-G', 'KRT18', 'KRT8', 'GATA3','KRT17', 'KRT7')#滋养
genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'DCN','LUM', 'OXTR', 'GJA1')#smc
genes_to_check3 = c('EPCAM','CD24','CCL21', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#上皮 纤维 内皮 淋巴内皮
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('SLC2A1', 'HK2', 'HIF1A', 'LDHA','OXTR', 'GJA1')#代谢
genes_to_check3 <- c('VWF', 'PECAM1','HLA-G','KRT7','DCN', 'LUM', 'TAGLN', 'ACTA2', 'HBB', 'HBA1','TRAC', 'CD3D','CD14', 'KLRD1','NCAM1', 'CD34')#all
genes_to_check3 = c('FCGR2A', 'CD33', 'KIT', 'ENPP3','C3AR1', 'CCR3')#other
genes_to_check3 = c('UBE2C', 'CDC20', 'TOP2A', 'STMN1','CD28', 'UCP2','ABCG2')#other
genes_to_check3 = c('IL7R', 'EBI3', 'SAT1', 'IL2RA','IL12A', 'IL10','IL27')#Th1 th2 th3 treg tr1 tfh th17 th9 th22
genes_to_check3 <- c("ITGAX","ANPEP","FUT4",'MPO',"CD33","CMA1","ITGA2B","TPSAB1","TPSB2")#树突，粒细胞，肥大细胞
genes_to_check3 = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                    'CCR7', 'SELL' , 'TCF7','CXCR6' , 'ITGA1',
                    'FOXP3', 'IL2RA',  'CTLA4','GZMB', 'GZMK','CCL5',
                    'IFNG', 'CCL4', 'CCL3' ,
                    'PRF1' , 'NKG7') 
genes_to_check3 = c("PRL","IGFBP1","WNT4","BMP2", "PRL8A2")#蜕膜
genes_to_check3 = c("ACKR1","ACKR2","ACKR3","IFNG", "IL7R", "IL7","CCR7","CXCR1","CXCR2")
genes_to_check3 <- c('PECAM1','VWF','TAGLN','ACTA2','DCN','LUM','TFF3','CCL21','LYZ','PTPRC','HBB','HLA-G')
genes_to_check3 <- c('CSF3R','CXCR2','FCGR3B','CD27','CD79B','KIT','SOD2','SNX10','SLPI','SPI1')
genes_to_check3 <- c('CCR1','CCR3','CCR4','CCR5')#CCLshouti
VlnPlot(TNL_TIL,features = genes_to_check3, pt.size = 0,ncol = 2)
VlnPlot(TNL_TIL,features = "DDX21", pt.size = 0,idents = c("SMC_1TNL","SMC_2TIL"))
FeaturePlot(TNL_TIL, features = genes_to_check3, reduction = "umap",min.cutoff = 0)
DefaultAssay(TNL_TIL) <- "RNA"
DefaultAssay(TNL_TIL) <- "SCT"
Idents(TNL_TIL) <- "celltype.group.new"
TNL_TIL_markers <- FindAllMarkers(TNL_TIL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TNL_TIL_markers21 <- FindMarkers(TNL_TIL, ident.1 = "21",only.pos = TRUE, logfc.threshold = 0.25)
TNL_TIL_markersSMC <- FindMarkers(TNL_TIL, ident.1 = "SMC_2TIL",ident.2 = "SMC_1TNL",only.pos = F)
write.csv(TNL_TIL_markersSMC,file="singlecellr/TNL_TIL_markers_smc.csv")
markers <- FindMarkers(pbmc_small, ident.1 = "g1", group.by = 'groups', subset.ident = "2")
top10 <- TNL_TIL_markers %>% group_by(cluster) %>% top_n(n = 20)
write.csv(top10,file="singlecellr/TNL_TIL_markers_top20.csv")

library(scCATCH)
clu_markers2 <- findmarkergenes(TNL_TIL,
                                species = "Human",
                                cluster = c('5','30','32','36','10'),
                                match_CellMatch = TRUE,
                                cancer = NULL,
                                tissue = c("Blood","Embryo","Myometrium","Uterus","Amniotic fluid","Placenta",
                                           "Endometrium","Endometrium stroma","Blood vessel","Epithelium","Pluripotent stem cell","Embryonic stem cell"),
                                cell_min_pct = 0.25,
                                logfc = 0.25,
                                pvalue = 0.05)
clu_ann <- scCATCH(clu_markers2$clu_markers,
                   species = "Human",
                   cancer = NULL,
                   tissue = c("Blood","Embryo","Myometrium","Uterus","Amniotic fluid","Placenta",
                              "Endometrium","Endometrium stroma","Blood vessel","Epithelium","Pluripotent stem cell","Embryonic stem cell")
)
saveRDS(clu_ann, file = "singlecellr/til_tnl/clu_ann.rds")
clu_ann <- readRDS("singlecellr/til_tnl/clu_ann.rds")

#标注
SMC=c(8,12,13,17,22,28,35)
Macrophages=c(0,1,3,15,26,29,33,38)
Trophectoderm_cells=c(6,16)
Endothelial_cells=c(2,9,10,20,23,39)
Fibroblasts=c(4,19,25,27)
LEC=c(11,14)
Red_blood_cells=c(18)
Neutrophils=c(5,30)
NK_cells=c(24)
T_cells=c(7,21,37,40)
B_cells=c(32)
Mast_cells=c(36)
Unknown=c(31,34)

current.cluster.ids <- c(SMC,Macrophages,Neutrophils,NK_cells,T_cells,B_cells,Mast_cells,Unknown,
                         Endothelial_cells,
                         Fibroblasts,
                         Trophectoderm_cells,
                         LEC,Red_blood_cells)

new.cluster.ids <- c(rep("SMC",length(SMC)),rep("Macrophages",length(Macrophages)), rep("Neutrophils",length(Neutrophils)),
                     rep("NK_cells",length(NK_cells)),
                     rep("T_cells",length(T_cells)),rep("B_cells",length(B_cells)),
                     rep("Mast_cells",length(Mast_cells)),rep("Unknown",length(Unknown)),
                     rep("Endothelial_cells",length(Endothelial_cells)),
                     rep("Fibroblasts",length(Fibroblasts)),
                     rep("Trophectoderm_cells",length(Trophectoderm_cells)),
                     rep("LEC",length(LEC)),
                     rep("Red_blood_cells",length(Red_blood_cells)))

TNL_TIL@meta.data$celltype.new <- plyr::mapvalues(x = as.integer(as.character(TNL_TIL@meta.data$integrated_snn_res.1)), from = current.cluster.ids, to = new.cluster.ids)
TNL_TIL_macro@meta.data$Macrotype[TNL_TIL_macro@meta.data$integrated_snn_res.0.1=="6"]=c("Redblood-AM")
TNL_TIL@meta.data$celltype.new[TNL_TIL@meta.data$celltype.new == "Myeloid_cells"]  <- c("Monocytic")
TNL_TIL@meta.data$celltype.new[TNL_TIL@meta.data$celltype.new == "Trophectoderm_cells"]  <- c("Trophoblast_cells")
  
DimPlot(TNL_TIL, reduction = "umap",group.by = "celltype.new", split.by = "groups2",order = NULL)
DimPlot(TNL_TIL, reduction = "umap",group.by = "celltype.new")

Idents(TNL_TIL) <- "celltype.new"
TNL_TIL_7cells.markers <- FindAllMarkers(TNL_TIL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_7cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="singlecellr/til_tnl/top10_TNL_TIL_7.csv")
write.csv(TNL_TIL_7cells.markers,file="singlecellr/til_tnl/marker_TNL_TIL_7.csv")
TNL_TIL_7cells.markers <- read.csv("singlecellr/til_tnl/marker_TNL_TIL_7.csv")
top10 <- read.csv("singlecellr/til_tnl/top10_TNL_TIL_7.csv")
DoHeatmap(TNL_TIL,features = top10$gene, cells = sample(1:50000,10000))
set.seed(42)
TNL_TIL_300 <- subset(TNL_TIL, downsample = 300)
DoHeatmap(TNL_TIL_300, features = top10$gene, angle = 15)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))


#new celltype
TNL_TIL_clear <- subset(TNL_TIL, celltype.new == "Unknown", invert = T)
TNL_TIL_clear <- readRDS("singlecellr/til_tnl/TNL_TIL_clear.rds")#####
saveRDS(TNL_TIL_clear, file = "singlecellr/til_tnl/TNL_TIL_clear.rds")#####

table(Idents(TNL_TIL_clear))
table(TNL_TIL$celltype.group.new)
table(TNL_TIL_clear$celltype.new)
table(TNL_TIL_clear$celltype.new, TNL_TIL_clear$orig.ident)
prop.table(table(Idents(TNL_TIL_clear), TNL_TIL$groups))
prop.table(table(TNL_TIL_clear$celltype.new, TNL_TIL_clear$orig.ident))

DimPlot(TNL_TIL_clear, reduction = "umap",group.by = "celltype.new", split.by = "groups",order = NULL)
DimPlot(TNL_TIL_clear, reduction = "umap",group.by = "celltype.new",label = T,repel = T)
DimPlot(TNL_TIL_clear, reduction = "umap",order = NULL)
Idents(TNL_TIL_clear) <- "celltype.new"
DefaultAssay(TNL_TIL_clear) <- "SCT"
TNL_TIL_12cells.markers <- FindAllMarkers(TNL_TIL_clear, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_12cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(top10,file="singlecellr/til_tnl/top10_TNL_TIL_12.csv")
write.csv(TNL_TIL_12cells.markers,file="singlecellr/til_tnl/marker_TNL_TIL_12.csv")

TNL_TIL_13cells.markers <- FindAllMarkers(TNL_TIL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_13cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- top10[order(top10$cluster),]
write.csv(top10,file="singlecellr/til_tnl/top10_TNL_TIL_13.csv")
write.csv(TNL_TIL_12cells.markers,file="singlecellr/til_tnl/marker_TNL_TIL_13.csv")

TNL_TIL_12cells.markers <- read.csv("singlecellr/til_tnl/marker_TNL_TIL_12.csv")
top10 <- read.csv("singlecellr/til_tnl/top10_TNL_TIL_12.csv")
DoHeatmap(TNL_TIL_clear,features = top10$gene, cells = sample(1:50000,10000))
set.seed(42)
TNL_TIL_clear_300 <- subset(TNL_TIL_clear, downsample = 300)
DefaultAssay(TNL_TIL_clear_300) <- "SCT"
Idents(TNL_TIL_clear_300) <- "celltype.new"
DoHeatmap(TNL_TIL_clear_300, features = top10$gene,slot = "scale.data", angle = 20,label = F,group.by = "celltype.new")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))


mapgene <- c('CD79A','VWF','PECAM1','LUM','COL1A1','TFF3','CCL21','TPSB2','CD14','LYZ','FCGR3B','CXCL8','KLRD1','HBB','ACTA2','TAGLN','CD3D','IL7R','KRT7')
mapgene <- c('CD14','LYZ','ACTA2','TAGLN','TFF3','CCL21','VWF','PECAM1','LUM','COL1A1')
mapgene <- c('TPSB2','KLRD1','CD3D','IL7R','HBB','CD79A','KRT7','PAPPA2','FCGR3B','CXCL8')
mapgene <- c('CD79A','VWF','PECAM1','LUM','COL1A1','TFF3','CCL21','TPSB2','CD14','LYZ')
mapgene <- c('FCGR3B','CXCL8','KLRD1','HBB','ACTA2','TAGLN','CD3D','IL7R','KRT7')
mapgene <- c('CD14','LYZ','ACTA2','TAGLN','TFF3','CCL21','VWF','PECAM1','LUM','COL1A1','TPSB2','KLRD1','CD3D','CD2','HBB','CD79A','KRT7','PAPPA2','FCGR3B','CXCL8')
library(MySeuratWrappers)
MySeuratWrappers::VlnPlot(TNL_TIL_clear,group.by = "celltype.new",
                          features = mapgene,stacked=T,sort = F,direction = "vertical",pt.size=0,axis.ticks.size = 1,axis.text.size = 12,axis.title.size=15)+NoLegend()
MySeuratWrappers::MultiFeaturePlot(TNL_TIL_clear, features = mapgene, reduction = "umap",min.cutoff = 0)






#找差异基因
TNL_TIL$celltype.group.new <- paste(Idents(TNL_TIL),TNL_TIL$groups2,sep = "_")
Idents(TNL_TIL) <- "celltype.group.new"
Idents(TNL_TIL) <- "celltype.new"
DefaultAssay(TNL_TIL) <- "SCT"
DefaultAssay(TNL_TIL) <- "RNA"
##SMC,Monocytic,Neutrophils,NK_cells,T_cells,B_cells,Mast_cells,Unknown,Endothelial_cells,Fibroblasts,Trophectoderm_cells,LEC,Red_blood_cells
FindMarkers(TNL_TIL, ident.1 = "Monocytic_1TNL", ident.2 = "Monocytic_2TIL",test.use = "wilcox", verbose = FALSE, min.pct = 0.5) %>% write.csv(file="singlecellr/til_tnl/Red_blood_cells_TNL_TIL_new.csv")
genedeg<-FindMarkers(TNL_TIL, ident.1 = "Neutrophils_1TNL", ident.2 = "Neutrophils_2TIL", verbose = FALSE, assay = 'SCT',slot = 'counts',test.use = "DESeq2") 
des_genes <- read.csv("singlecellr/til_tnl/Immune_cells_TNL_TIL.csv") %>% filter(, p_val_adj < 0.05) %>% filter(, avg_log2FC > 0.584962501 | avg_log2FC < -0.584962501)
des_genes <- arrange(des_genes, avg_log2FC)
DoHeatmap(object = subset(TNL_TIL, celltype == "Immune_cells"),features = des_genes$X, group.by = "groups2")

Idents(TNL_TIL) <- "celltype.new"
mapgene <- c('MT2A', 'JUN','TAGLN', 'RGS5','TPM1', 'TIMP1','HIF1A')
mapgene <- c('MALAT1', 'IFI6', 'VCAM1', 'ENPP2','FTL', 'TIMP3','TAGLN', 'DES','TPM1', 'CLU')
mapgene <- c('SNCG', 'IFI6', 'VCAM1', 'ENPP2','FTL', 'TIMP3','TAGLN', 'DES','TPM1', 'CLU')
mapgene <- c('CD14', 'CST3', 'C1QB', 'CD74','STMN1', 'HLA-DRA','IL32', 'CYP1B1','CTSB', 'SELENOP')
mapgene <- c('S100A8', 'S100A9', 'C1QB', 'C1QA','SELENOP', 'RNASE1')
mapgene <- c('COL1A2', 'COL3A1', 'COL1A1', 'C1QA','SELENOP', 'RNASE1')
mapgene <- c('PLSCR4', 'COL3A1', 'FABP4', 'C1QA','SELENOP', 'RNASE1')
mapgene <- c('CRYAB', 'TFF3', 'EDN1', 'EPCAM','FN1', 'H19')
mapgene <- c('HLA-G', 'FLT1', 'PAPPA2', 'PTGDS','EBI3', 'CDKN1C')
mapgene <- c('CCL3', 'CCL4', 'CXCL10', 'MALAT1','EBI3', 'CDKN1C')
mapgene <- c('IGFBP7', 'NOTCH3', 'DES', 'RGS5','ADIRF', 'EGR1')
mapgene <- c('ACTA2', 'TIMP1', 'NNMT', 'VCAN','PTGDS', 'SAT1')
mapgene = c('MYLK', 'OXTR','GJA1','PTGS2','PTGER3','ACTA2','PTGDS','GJA4')#收缩基因1
mapgene = c('JUN', 'RGS5','ATF3','ACTA2','HIF1A','SRGN','NNMT','TIMP1')#收缩基因1
mapgene = c('S100A8', 'S100A12','VCAN','IL1B','C1QB','CCL4','CD163','SELENOP')
mapgene = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
mapgene = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
mapgene = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2','CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子
mapgene = c('CXCR4','TGFBR2','IL6ST','CSF1R','IL13RA1','IL1R2','CSF2RB','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2','IL2RB','IL1RAP','IL7R')#炎性受体
mapgene = c('CXCR4','TGFBR2','IL6ST','CSF1R','IL13RA1','IL1R2','CSF2RB','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2','IL2RB','IL1RAP','ACKR1','ACKR2','ACKR3','IFNG','IL7','CCR7','CXCR1','CXCR2')#炎性受体
mapgene = c('OXT','OXTR','PTGS2','PTGER1','PTGER2','PTGER3','PTGER4',	'PTGFR')
mapgene = c('PTGS2','PTGER3','OXTR','GJA1')#top
mapgene = c('PTGER3','PTGER2','OXT','PGR')
mapgene = c('CXCL8','CXCR1','CXCR2')#因子受体
mapgene = c('TNF','TNFRSF1A','TNFRSF1B')
mapgene = c('IL1A','IL1B','IL1R1','IL1R2','IL1RAP')
mapgene = c('IL4','IL10','IFNG')
mapgene = c('IL6','IL6R','IL6ST')
mapgene = c('JUN','CALD1','FOS','CALM1')
mapgene = c('CXCL8','CXCR1','CXCR2','TNF','TNFRSF1A','TNFRSF1B','IL1A','IL1B','IL1R1','IL1R2','IL1RAP','IL6','IL6R','IL6ST','CCL2','CCR2','CCL4','CCR5')
mapgene = c('ACKR1','PNP','ADAMTS9','SPRY1','TIMP1')
VlnPlot(TNL_TIL_clear, features = mapgene, pt.size = 0,ncol = 2,idents = c("SMC_1TNL","SMC_2TIL"))
VlnPlot(TNL_TIL_clear, features = mapgene, pt.size = 0,ncol = 4,group.by = "celltype.new")
FeaturePlot(TNL_TIL_clear, features = mapgene, reduction = "umap",min.cutoff = 0,ncol = 4)
FeaturePlot(TNL_TIL_clear, features = mapgene, reduction = "umap",min.cutoff = 0,ncol = 2,split.by = "groups2")
Idents(TNL_TIL_clear) <- "celltype.new"
Idents(TNL_TIL_clear) <- "celltype.group.new"
DefaultAssay(TNL_TIL_clear) <- "SCT"
MySeuratWrappers::VlnPlot(TNL_TIL_clear,features = mapgene,stacked=T,pt.size=0,group.by = "celltype.new")
MySeuratWrappers::MultiFeaturePlot(TNL_TIL_clear,features = mapgene,ncol=4,split.by = "groups2")
VlnPlot(TNL_TIL_clear, features = mapgene, pt.size = 0,ncol = 4,sort = T, split.by = "groups2",group.by = "celltype.new")
FeaturePlot(TNL_TIL_clear, features = mapgene, reduction = "umap",min.cutoff = 0,ncol = 4)
MySeuratWrappers::MultiFeaturePlot(TNL_TIL_clear, features = mapgene, reduction = "umap",min.cutoff = 0,ncol = 8)
DotPlot(TNL_TIL_clear, features = mapgene, group.by = "celltype.group.new",col.min=0)+ theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+coord_flip()

#相关性
TIL_v_ptgs2 <- subset(TIL_v, subset = PTGS2 >1,slot = "counts")
TIL_v_sct_ptgs2 <- GetAssayData(TIL_v_ptgs2[["Spatial"]], slot = "counts")
col2 <- as.matrix(TIL_v_sct_ptgs2)
a <- t(col2)
b <- cor(a)
col1 <- colorRampPalette(c("navy", "white", "firebrick3"))

TNL_TIL_clear_ptgs2 <- subset(TNL_TIL_clear, subset = PTGS2 >1,slot = "counts")
TNL_TIL_clear_ptgs2_nm <- subset(TNL_TIL_clear_ptgs2, idents = c("Neutrophils","Monocytic"))
TNL_TIL_clear_sct_ptgs2 <- GetAssayData(TNL_TIL_clear_ptgs2[["SCT"]], slot = "data")
col2 <- as.matrix(TNL_TIL_clear_sct_ptgs2)
a <- t(col2)
b <- cor(a)
col1 <- colorRampPalette(c("navy", "white", "firebrick3"))


##SMC
TNL_TIL_SMC <- subset(TNL_TIL, celltype.new == "SMC")
TNL_TIL_SMC <- ScaleData(TNL_TIL_SMC, verbose = FALSE)
TNL_TIL_SMC <- RunPCA(TNL_TIL_SMC, npcs = 30, verbose = FALSE)
TNL_TIL_SMC <- RunUMAP(TNL_TIL_SMC, reduction = "pca", dims = 1:30)
TNL_TIL_SMC <- FindNeighbors(TNL_TIL_SMC, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_SMC) <- "celltype.new"
DimPlot(TNL_TIL_SMC, reduction = "umap")

p1 <- DimPlot(TNL_TIL_SMC, reduction = "umap", group.by = "groups2")
p2 <- DimPlot(TNL_TIL_SMC, reduction = "umap", group.by = "SMCtype",label = TRUE, repel = TRUE)
p1 + p2
DefaultAssay(TNL_TIL_SMC) <- "integrated"
TNL_TIL_SMC <- FindClusters(object = TNL_TIL_SMC, resolution = 0.1)
TNL_TIL_SMC@meta.data$SMCtype[TNL_TIL_SMC@meta.data$integrated_snn_res.0.1=="0"]=c("SMC_1")
TNL_TIL_SMC$SMCtype.group.new <- paste(Idents(TNL_TIL_SMC),TNL_TIL_SMC$groups2,sep = "_")
Idents(TNL_TIL_SMC) <- "SMCtype"
Idents(TNL_TIL_SMC) <- "SMCtype.group.new"
DefaultAssay(TNL_TIL_SMC) <- "SCT"
genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'OXTR','GJA1','LUM','DCN')#smc
genes_to_check3 = c('EPCAM', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#上皮 纤维 内皮 淋巴内皮
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('SLC2A1', 'HK2', 'HIF1A', 'LDHA','OXTR', 'GJA1')#代谢
genes_to_check3 = c('PTGS2', 'OXTR','GJA1','ACTA2','PTGER3','PTGER2','PTGDS','GJA4')#收缩基因
genes_to_check3 = c('OXTR','GJA1','GJA4','ACTA2','PTGER3','PTGER2','IGFBP3','IGFBP5','GNG11','RGS2','RGS16','RAMP1','RAMP2','RAMP3','CNN1','MYH11')#收缩基因
genes_to_check3 = c('MYLK', 'OXTR','GJA1','PTGER3','PTGDS','GJA4','ELANE','IFNG','PGR')
FeaturePlot(TNL_TIL_SMC, features = genes_to_check3, reduction = "umap",min.cutoff = 0)
MySeuratWrappers::MultiFeaturePlot(TNL_TIL_SMC,features = genes_to_check3,ncol=4)
DotPlot(TNL_TIL_SMC, features = genes_to_check3, group.by = "SMCtype",col.min=0)+ theme(axis.text.x = element_text(angle= 90,  vjust = 0.5, hjust=0.5))+coord_flip()
MySeuratWrappers::VlnPlot(TNL_TIL_SMC,features = genes_to_check3,stacked=T,pt.size=0,group.by = "SMCtype")
saveRDS(TNL_TIL_SMC, file = "singlecellr/til_tnl/TNL_TIL_SMC.rds")
TNL_TIL_SMC <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_SMC.rds")#######

TNL_TIL_SMC.markers <- FindAllMarkers(TNL_TIL_SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_SMC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_order <- top10[order(top10$cluster),]
write.csv(top10,file="singlecellr/til_tnl/top10_TNL_TIL_13_smc.csv")
write.csv(TNL_TIL_SMC.markers,file="singlecellr/til_tnl/marker_TNL_TIL_13_smc.csv")
top10 <- read.csv("singlecellr/til_tnl/top10_TNL_TIL_13_smc.csv")
table(TNL_TIL_SMC$groups2,TNL_TIL_SMC$integrated_snn_res.0.1)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]
TNL_TIL_SMC_300 <- subset(TNL_TIL_SMC, downsample = 300)
DefaultAssay(TNL_TIL_SMC_300) <- "integrated" #integrated
DefaultAssay(TNL_TIL_SMC_300) <- "SCT" #integrated
DoHeatmap(TNL_TIL_SMC_300, features = top10_order$gene,slot = "scale.data",group.by = "SMCtype", angle = 20,label = F)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))
FeaturePlot(TNL_TIL_SMC, features = top10$gene, reduction = "umap",min.cutoff = 0)

VlnPlot(TNL_TIL_SMC, features = genes_to_check3, pt.size = 0,ncol = 2, split.by = "groups2")
p <- DoHeatmap(TNL_TIL_SMC, features = top10$gene,slot = "scale.data", angle = 20,label = F)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))
ggsave("x.png", p , width = 7, height = 4, dpi = 300)

library(monocle)#######
data <- GetAssayData(TNL_TIL_SMC, assay = 'RNA', slot = 'counts')
cell_metadata <- TNL_TIL_SMC@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

data_matrix <- as(as.matrix(TNL_TIL_SMC@assays$RNA@data), 'sparseMatrix')
feature_ann <- data.frame(gene_id=rownames(data_matrix),gene_short_name=rownames(data_matrix))
rownames(feature_ann)<-rownames(data_matrix)
data_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-TNL_TIL_SMC@meta.data
rownames(sample_ann)<-colnames(data_matrix)
data_pd <- new("AnnotatedDataFrame", data =sample_ann)
data_cds <- newCellDataSet(data_matrix,phenoData =data_pd, featureData =data_fd, expressionFamily=negbinomial.size())
data_cds <- estimateSizeFactors(data_cds)
data_cds <- estimateDispersions(data_cds, parallel = T)
disp_table <- dispersionTable(data_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
data_cds <- setOrderingFilter(data_cds, unsup_clustering_genes$gene_id)
data_cds <- reduceDimension(
  data_cds,
  max_components = 2,
  method = 'DDRTree')
data_cds <- orderCells(data_cds)
p1 <- plot_cell_trajectory(data_cds,cell_size = 1,color_by = "groups2")
p2 <- plot_cell_trajectory(data_cds,cell_size = 1,color_by = "SMCtype")
p3 <- plot_cell_trajectory(data_cds,cell_size = 1,color_by = "Pseudotime")
p1|p2|p3
saveRDS(data_cds, file = "singlecellr/til_tnl/TNL_TIL_SMC_cds.rds")
data_cds <- readRDS(file = "singlecellr/til_tnl/TNL_TIL_SMC_cds.rds")#######

disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
ordering_genes = disp.genes$gene_id
plot_pseudotime = plot_pseudotime_heatmap(data_cds[ordering_genes, ], num_clusters = 3,cores = 2, return_heatmap = T, show_rownames = F)

diff_test <- differentialGeneTest(data_cds[ordering_genes,], cores = 4, fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))
p = plot_pseudotime_heatmap(data_cds[sig_gene_names,], cores = 10,num_clusters=3, show_rownames=T, return_heatmap=T)
p$tree_row
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering,file="singlecellr/til_tnl/TNL_TIL_smc_pseudo_cluster5.csv")
s.genes <- c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'OXTR', 'GJA1')
des.genes <- rownames(clustering)

contraction.genes <- read.table("singlecellr/til_tnl/WP_MYOMETRIAL_RELAXATION_AND_CONTRACTION_PATHWAYS.txt")
contraction.genes <- data.frame(contraction.genes)
contraction.genes <- contraction.genes$V1
overlape.gene <- c('ACTA2', 'ACTB', 'ADM', 'ATF3','CNN1', 'ETS2', 'FOS','IGFBP2')
overlape.gene <- c('IGFBP3', 'IGFBP4', 'IGFBP5', 'IGFBP6','JUN', 'RAMP1', 'RGS16','RGS2','RGS5')
overlape.gene <-intersect(des.genes,contraction.genes)
plot_genes_in_pseudotime(data_cds[overlape.gene,], color_by = "SMCtype")
plot_pseudotime_heatmap(data_cds[overlape.gene, ],cores = 2, return_heatmap = T, show_rownames = T)
table(data_cds@phenoData@data[["State"]],data_cds@phenoData@data[["groups2"]])

clusters <- cutree(plot_pseudotime$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.csv(clustering,file="singlecellr/til_tnl/SMC-pse-cluster.csv")
test_genes <- c('CXCR4', 'CCL3','CCL4','CCL8', 'CCL2','CCL3L1','CCL4L2','CXCL3','CXCL10','TNFSF10','CSF1R','IL7R')#1
test_genes <- c('CCL21')#2
test_genes <- c('CSF3R','CXCL8','IL1B','CXCL2')#3
test_genes = c('MYLK', 'OXTR','GJA1','HIF1A','PTGER3','PTGDS','GJA4','PTGS2')
test_genes <- c('CXCR4', 'CCL3','CCL4','CCL8', 'CCL2','CCL3L1','CCL4L2','CXCL3','CXCL10','TNFSF10','CSF1R','IL7R','CCL21','CSF3R','CXCL8','IL1B','CXCL2')#all
plot_genes_in_pseudotime(data_cds[test_genes,],
                         color_by = "groups2",
                         cell_size=0.5,
                         ncol = 4)

##Macrophages  myeloid
DefaultAssay(TNL_TIL) <- "integrated"
TNL_TIL_macro <- subset(TNL_TIL, celltype.new == "Macrophages")
TNL_TIL_macro <- ScaleData(TNL_TIL_macro, verbose = FALSE)
TNL_TIL_macro <- RunPCA(TNL_TIL_macro, verbose = FALSE)
TNL_TIL_macro <- RunUMAP(TNL_TIL_macro, reduction = "pca", dims = 1:30)
TNL_TIL_macro <- FindNeighbors(TNL_TIL_macro, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_macro) <- "celltype.new"
DimPlot(TNL_TIL_macro, reduction = "umap",split.by = "groups2")

p1 <- DimPlot(TNL_TIL_macro, reduction = "umap", group.by = "groups2",split.by = "groups2")
p2 <- DimPlot(TNL_TIL_macro, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL_macro, reduction = "umap", label = TRUE, repel = TRUE,group.by = "Macrotype",split.by = "groups2")
DefaultAssay(TNL_TIL_macro) <- "integrated"
TNL_TIL_macro <- FindClusters(object = TNL_TIL_macro, resolution = 0.1)
TNL_TIL_macro@meta.data$Macrotype[TNL_TIL_macro@meta.data$integrated_snn_res.0.1=="6"]=c("Redblood-AM")
TNL_TIL_macro$Macrotype.group.new <- paste(Idents(TNL_TIL_macro),TNL_TIL_macro$groups2,sep = "_")
table(Idents(TNL_TIL_macro))
table(TNL_TIL_macro$Macrotype,TNL_TIL_macro$groups)
Idents(TNL_TIL_macro) <- "Macrotype.group.new"
Idents(TNL_TIL_macro) <- "Macrotype"
DefaultAssay(TNL_TIL_macro) <- "SCT"
genes_to_check3 = c('CD163', 'CD86', 'CD14', 'PTPRC','IL1B', 'MRC1','S100A9','S100A8','LYZ','TYROBP')#巨噬
genes_to_check3 = c('FCGR3A', 'NLRP3', 'CD14', 'FCN1','IL1B', 'C1QC','SPP1','FCGR3A')#巨噬
genes_to_check3 = c('CD40', 'CD80', 'CD86','LILRA4', 'CD1C', 'BATF3')#DC
genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'OXTR','GJA1','LUM','DCN')#smc
genes_to_check3 = c('EPCAM', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#上皮 纤维 内皮 淋巴内皮
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL3','CCL3L1','CCL4','CCL4L2','CSF1R','IL1B')#macro炎性
genes_to_check3 = c('CD1C', 'LUM','IL1B','MRC1','MKI67','HBA1','TAGLN')#top
genes_to_check3 = c('CD86','IL1B','CD163', 'MRC1')#top
genes_to_check3 = c('C1QB', 'SPP1','IL1B','FCN1','COL1A1','CCL21','TAGLN','IGFBP5','CD1C','FCER1A','STMN1','MKI67','HBA1','HBA2')
genes_to_check3 = c('C1QB', 'MRC1','IL1B','FCN1','COL1A1','LUM','TAGLN','RGS5','CD1C','FCER1A','MKI67','HBA1')
FeaturePlot(TNL_TIL_macro, features = genes_to_check3, reduction = "umap",min.cutoff = 0,ncol = 4)
saveRDS(TNL_TIL_macro, file = "singlecellr/til_tnl/TNL_TIL_macro.rds")#####
TNL_TIL_macro <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_macro.rds")#####

MySeuratWrappers::MultiFeaturePlot(TNL_TIL_macro,features = genes_to_check3,ncol=4)
VlnPlot(TNL_TIL_macro, features = genes_to_check3, pt.size = 0,ncol = 1, group.by = "Macrotype")
VlnPlot(TNL_TIL_macro, features = genes_to_check3, pt.size = 0,ncol = 2, split.by = "groups2",group.by = "Macrotype")
MySeuratWrappers::VlnPlot(TNL_TIL_macro,features = genes_to_check3,stacked=T,pt.size=0,group.by = "Macrotype.group.new",
                          cols = color7)
color7 <- c("#F8766D","#F8766D", "#C49A00","#C49A00", "#53B400","#53B400", "#00C094","#00C094", "#00B6EB","#00B6EB", "#A58AFF","#A58AFF", "#FB61D7", "#FB61D7")
TNL_TIL_macro.markers <- FindAllMarkers(TNL_TIL_macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_macro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = TNL_TIL_macro, features = top10$gene)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

write.csv(top10,file="singlecellr/til_tnl/top10_TNL_TIL_13_macro.csv")
top10 <- read.csv("singlecellr/til_tnl/top10_TNL_TIL_13_macro.csv")
write.csv(TNL_TIL_macro.markers,file="singlecellr/til_tnl/marker_TNL_TIL_13_macro.csv")
VlnPlot(TNL_TIL_macro, features = genes_to_check3, pt.size = 0,ncol = 2, split.by = "groups2")
table(TNL_TIL_macro$integrated_snn_res.0.1, TNL_TIL_macro$groups)


##T NK
DefaultAssay(TNL_TIL) <- "integrated"
Idents(TNL_TIL) <- TNL_TIL$celltype.new
TNL_TIL_TNK <- subset(TNL_TIL, idents = c("T_cells","NK_cells"))
TNL_TIL_TNK <- ScaleData(TNL_TIL_TNK, verbose = FALSE)
TNL_TIL_TNK <- RunPCA(TNL_TIL_TNK, verbose = FALSE)
TNL_TIL_TNK <- RunUMAP(TNL_TIL_TNK, reduction = "pca", dims = 1:30)
TNL_TIL_TNK <- FindNeighbors(TNL_TIL_TNK, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_TNK) <- "celltype.new"
DimPlot(TNL_TIL_TNK, reduction = "umap",split.by = "groups2",group.by = "Ttype")

p1 <- DimPlot(TNL_TIL_TNK, reduction = "umap", group.by = "groups2",split.by = "groups2")
p2 <- DimPlot(TNL_TIL_TNK, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL_TNK, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.2")
DefaultAssay(TNL_TIL_TNK) <- "integrated"
TNL_TIL_TNK <- FindClusters(object = TNL_TIL_TNK, resolution = 0.2)

TNL_TIL_TNK@meta.data$Ttype[TNL_TIL_TNK@meta.data$integrated_snn_res.0.2=="8"]=c("CD8+ T_cells")
TNL_TIL_TNK$Ttype.group.new <- paste(Idents(TNL_TIL_TNK),TNL_TIL_TNK$groups2,sep = "_")
table(Idents(TNL_TIL_TNK))
table(TNL_TIL_TNK$Macrotype,TNL_TIL_macro$groups)
Idents(TNL_TIL_TNK) <- "Ttype"
Idents(TNL_TIL_TNK) <- "Ttype.group.new"
Idents(TNL_TIL_TNK) <- "integrated_snn_res.0.2"
DefaultAssay(TNL_TIL_TNK) <- "SCT"
genes_to_check3 = c('CD3E', 'CD3D', 'CD2')#T
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1','FOXP3', 'CTLA5')#免疫抑制
genes_to_check3 = c('CD4','GZMA', 'CCL5','CD69', 'TCF7','LEF1','FOXP3','IL2RA', 'ICA1','TOX2', 'IL17A','CTSH')
genes_to_check3 = c('CD4','CD8A', 'BTLA','CD69', 'CCR7','CCL5','KLRD1','C1QB', 'GSN','CXCL14', 'GNLY','FCGR3A', 'CD2','TOP2A')
genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A','MKI67','FTL','IL7R','KIT')
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('CCL4','TNFSF10','IL7R','CCL3','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL4','TNFSF10','IL7R','CCL3','CCL5','XCL1','XCL2','IL2RB')#T NK炎性
genes_to_check3 = c('CD1C', 'LUM','IL1B','MRC1','MKI67','HBA1','TAGLN')#top
genes_to_check3 = c('CCR7', 'IL7R','FOXP3','CXCL13','ZFP36','GZMK','IFIT1','IFNG','LAG3','AREG','TRDC','CD7')#T marker
genes_to_check3 = c('CCR7','IL7R','CD8A','TAGLN','KLRD1','CD4','MKI67','KIT')#top
genes_to_check3 = c('CD4','CD8A','ACTA2','KLRD1','MAFB','MKI67','KIT','FOXP3')#top
genes_to_check3 = c('CD38','IFNG','ICOS','CCR6','ID2','CD2')#top
genes_to_check3 = c('C1QB', 'SPP1','IL1B','FCN1','COL1A1','CCL21','TAGLN','IGFBP5','CD1C','FCER1A','STMN1','MKI67','HBA1','HBA2')
genes_to_check3 = c('C1QB', 'MRC1','IL1B','FCN1','COL1A1','LUM','TAGLN','RGS5','CD1C','FCER1A','MKI67','HBA1')
FeaturePlot(TNL_TIL_TNK, features = genes_to_check3, reduction = "umap",min.cutoff = 0,ncol = 2)
MySeuratWrappers::MultiFeaturePlot(TNL_TIL_TNK,features = genes_to_check3,ncol=4)
saveRDS(TNL_TIL_TNK, file = "singlecellr/til_tnl/TNL_TIL_TNK.rds")#####
TNL_TIL_TNK <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_TNK.rds")#####

VlnPlot(TNL_TIL_TNK, features = genes_to_check3, pt.size = 0,ncol = 4)
VlnPlot(TNL_TIL_TNK, features = genes_to_check3, pt.size = 0,ncol = 4)
MySeuratWrappers::VlnPlot(TNL_TIL_TNK,features = genes_to_check3,stacked=T,pt.size=0,group.by = "Ttype.group.new",
                          cols = color7)
color7 <- c("#F8766D","#F8766D", "#CD9600","#CD9600", "#7CAE00","#7CAE00", "#00BE67","#00BE67", "#00BFC4","#00BFC4", "#00A9FF","#00A9FF", "#C77CFF", "#C77CFF", "#FF61CC", "#FF61CC")
TNL_TIL_TNK.markers <- FindAllMarkers(TNL_TIL_TNK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_TNK.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]
set.seed(42)
TNL_TIL_TNK_80 <- subset(TNL_TIL_TNK, downsample = 80)
DefaultAssay(TNL_TIL_TNK_80) <- "integrated"
DoHeatmap(TNL_TIL_TNK_80, features = top10_order$gene, angle = 15,group.by = "Ttype")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))+theme(axis.text.y = element_text(size = 5))

DoHeatmap(object = TNL_TIL_TNK, features = genes_to_check3)+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

write.csv(top10,file="singlecellr/til_tnl/top10_TNL_TIL_13_TNK.csv")
top10 <- read.csv("singlecellr/til_tnl/top10_TNL_TIL_13_TNK.csv")
write.csv(TNL_TIL_TNK.markers,file="singlecellr/til_tnl/marker_TNL_TIL_13_new_TNK.csv")
VlnPlot(TNL_TIL_macro, features = genes_to_check3, pt.size = 0,ncol = 2, split.by = "groups2")
table(TNL_TIL_macro$integrated_snn_res.0.1, TNL_TIL_macro$groups)



#粒细胞(neutrophils)
DefaultAssay(TNL_TIL_neu) <- "integrated"
TNL_TIL_neu <- subset(TNL_TIL, celltype.new == c("Neutrophils"))
TNL_TIL_neu <- ScaleData(TNL_TIL_neu, verbose = FALSE)
TNL_TIL_neu <- RunPCA(TNL_TIL_neu, verbose = FALSE)
TNL_TIL_neu <- RunUMAP(TNL_TIL_neu, reduction = "pca", dims = 1:30)
TNL_TIL_neu <- FindNeighbors(TNL_TIL_neu, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_neu) <- "celltype.new"
DimPlot(TNL_TIL_neu, reduction = "umap",split.by = "groups2",group.by = "TNL_TIL_neu",label.size = 50)

p1 <- DimPlot(TNL_TIL_neu, reduction = "umap", group.by = "groups2")
p2 <- DimPlot(TNL_TIL_neu, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL_TNK, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.2",split.by = "orig.ident")
DefaultAssay(TNL_TIL_neu) <- "integrated"
TNL_TIL_neu <- FindClusters(object = TNL_TIL_neu, resolution = 0.2)

TNL_TIL_neu@meta.data$TNL_TIL_neu[TNL_TIL_neu@meta.data$integrated_snn_res.0.2=="5"]=c("neutrophils_6")
TNL_TIL_neu$TNL_TIL_neu.group.new2 <- paste(Idents(TNL_TIL_neu),TNL_TIL_neu$groups2,sep = "_")
table(Idents(TNL_TIL_neu))
table(TNL_TIL_neu$Macrotype,TNL_TIL_neu$groups)
Idents(TNL_TIL_neu) <- "TNL_TIL_neu"
Idents(TNL_TIL_neu) <- "TNL_TIL_neu.group.new2"
DefaultAssay(TNL_TIL_neu) <- "SCT"

genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1','FOXP3', 'CTLA5')#免疫抑制
genes_to_check3 = c('CD4','GZMA', 'CCL5','CD69', 'TCF7','LEF1','FOXP3','IL2RA', 'ICA1','TOX2', 'IL17A','CTSH')
genes_to_check3 = c('CD4','CD8A', 'BTLA','CD69', 'CCR7','CCL5','KLRD1','C1QB', 'GSN','CXCL14', 'GNLY','FCGR3A', 'CD2','TOP2A')
genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A','MKI67','FTL','IL7R','KIT')
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL3','CCL4','CXCL8','TNFSF10','CCL2','IL6ST','CSF1R','IL13RA1','IL1R2','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2')#NEU炎性p
FeaturePlot(TNL_TIL_neu, features = genes_to_check3, reduction = "umap",min.cutoff = 0,ncol = 4)
MySeuratWrappers::MultiFeaturePlot(TNL_TIL_neu,features = genes_to_check3,ncol=4)
DotPlot(TNL_TIL_neu, features = genes_to_check3, group.by = "TNL_TIL_neu",split.by="groups",scale = TRUE)+ theme(axis.text.x = element_text(angle= 90,  vjust = 0.5, hjust=0.5))+coord_flip()
saveRDS(TNL_TIL_neu, file = "singlecellr/til_tnl/TNL_TIL_neu.rds")#####
TNL_TIL_neu <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_neu.rds")#####

VlnPlot(TNL_TIL_neu, features = genes_to_check3, pt.size = 0,ncol = 4, group.by = "TNL_TIL_neu.group.new2")
VlnPlot(TNL_TIL_TNK, features = genes_to_check3, pt.size = 0,ncol = 2, split.by = "groups2",group.by = "Macrotype")
MySeuratWrappers::VlnPlot(TNL_TIL_neu,features = genes_to_check3,stacked=T,pt.size=0,group.by = "TNL_TIL_neu.group.new2",
                          cols = color7)
color7 <- c("#F8766D","#F8766D", "#C49A00","#C49A00", "#53B400","#53B400", "#00C094","#00C094", "#00B6EB","#00B6EB", "#A58AFF","#A58AFF", "#FB61D7", "#FB61D7")
TNL_TIL_neu.markers <- FindAllMarkers(TNL_TIL_neu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_neu.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10$cluster <-as.character(top10$cluster)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(object = TNL_TIL_neu, features = top10_order$gene,group.by = "TNL_TIL_neu")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

write.csv(top10,file="singlecellr/til_tnl/top10_TNL_TIL_13_neu.csv")
top10 <- read.csv("singlecellr/til_tnl/top10_TNL_TIL_13_neu.csv")
write.csv(TNL_TIL_neu.markers,file="singlecellr/til_tnl/marker_TNL_TIL_13_neu.csv")
VlnPlot(TNL_TIL_macro, features = genes_to_check3, pt.size = 0,ncol = 2, split.by = "groups2")
table(TNL_TIL_macro$integrated_snn_res.0.1, TNL_TIL_macro$groups)

#B(B_cells)
DefaultAssay(TNL_TIL_B) <- "integrated"
TNL_TIL_B <- subset(TNL_TIL, celltype.new == c("B_cells"))
TNL_TIL_B <- ScaleData(TNL_TIL_B, verbose = FALSE)
TNL_TIL_B <- RunPCA(TNL_TIL_B, verbose = FALSE)
TNL_TIL_B <- RunUMAP(TNL_TIL_B, reduction = "pca", dims = 1:30)
TNL_TIL_B <- FindNeighbors(TNL_TIL_B, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_B) <- "celltype.new"
DimPlot(TNL_TIL_B, reduction = "umap",split.by = "groups2")

p1 <- DimPlot(TNL_TIL_B, reduction = "umap", group.by = "groups2")
p2 <- DimPlot(TNL_TIL_B, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL_B, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.5",split.by = "orig.ident")
DefaultAssay(TNL_TIL_B) <- "integrated"
TNL_TIL_B <- FindClusters(object = TNL_TIL_B, resolution = 0.5)

TNL_TIL_B@meta.data$B_type[TNL_TIL_B@meta.data$integrated_snn_res.0.1=="0"]=c("Memory_B_cells")
TNL_TIL_B$B_type.group.new <- paste(Idents(TNL_TIL_B),TNL_TIL_B$groups2,sep = "_")
table(Idents(TNL_TIL_B))
table(TNL_TIL_B$Macrotype,TNL_TIL_neu$groups)
Idents(TNL_TIL_B) <- "B_type"
Idents(TNL_TIL_neu) <- "integrated_snn_res.0.2"
Idents(TNL_TIL_B) <- "B_type.group.new"
DefaultAssay(TNL_TIL_B) <- "SCT"

genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1','FOXP3', 'CTLA5')#免疫抑制
genes_to_check3 = c('CD4','GZMA', 'CCL5','CD69', 'TCF7','LEF1','FOXP3','IL2RA', 'ICA1','TOX2', 'IL17A','CTSH')
genes_to_check3 = c('CD4','CD8A', 'BTLA','CD69', 'CCR7','CCL5','KLRD1','C1QB', 'GSN','CXCL14', 'GNLY','FCGR3A', 'CD2','TOP2A')
genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A','MKI67','FTL','IL7R','KIT')
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('CD79A','MS4A1','RFTN1', 'ITGB1','IGHA1','TXNIP','IGHD','CD27','SDC1','CD38')#B type
genes_to_check3 = c('CD19','CD38','CD24', 'CD27','BLIMP1','XBP1','IRF4','SDC1')#B plasmablasts and plasma cells.
genes_to_check3 = c('CD79A','MS4A1','CD27','SDC1','CD38','JCHAIN')#B type
genes_to_check3 = c('MS4A1','HLA-DRA', 'JCHAIN','XBP1')
genes_to_check3 = c('CCL21','CXCR4')#b炎性因子2
genes_to_check3 = c('CD79A','MS4A1','CD22','CD19','IGHM','TNFRSF17','IGHG1','CD27','CD24','SDC1','PAX5','PRDM1','FCER2','TCL1A','TUBB','STMN1')
genes_to_check3 = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL3','CCL4','CXCL8','TNFSF10','CCL2','IL6ST','CSF1R','IL13RA1','IL1R2','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2')#NEU炎性p
FeaturePlot(TNL_TIL_B, features = genes_to_check3, reduction = "umap",min.cutoff = 0,ncol = 5)
saveRDS(TNL_TIL_B, file = "singlecellr/til_tnl/TNL_TIL_B.rds")#####
TNL_TIL_B <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_B.rds")#####

MySeuratWrappers::MultiFeaturePlot(TNL_TIL_B,features = genes_to_check3,ncol=4)
VlnPlot(TNL_TIL_B, features = genes_to_check3, pt.size = 0,ncol = 7, group.by = "Macrotype")
VlnPlot(TNL_TIL_B, features = genes_to_check3, pt.size = 0,ncol = 2, split.by = "groups2",group.by = "Macrotype")
MySeuratWrappers::VlnPlot(TNL_TIL_B,features = genes_to_check3,stacked=T,pt.size=0,group.by = "B_type.group.new",cols = color2)

color2 <- c("#F8766D","#F8766D","#00BFC4","#00BFC4")
TNL_TIL_B.markers <- FindAllMarkers(TNL_TIL_B, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(TNL_TIL_B.markers,file="singlecellr/til_tnl/marker_TNL_TIL_13_B.csv")
top10 <- TNL_TIL_B.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(object = TNL_TIL_neu, features = top10$gene,group.by = "TNL_TIL_neu")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

#肥大细胞(Mast_cells)
TNL_TIL_Mast <- subset(TNL_TIL, celltype.new == c("Mast_cells"))
TNL_TIL_Mast <- ScaleData(TNL_TIL_Mast, verbose = FALSE)
TNL_TIL_Mast <- RunPCA(TNL_TIL_Mast, verbose = FALSE)
TNL_TIL_Mast <- RunUMAP(TNL_TIL_Mast, reduction = "pca", dims = 1:30)
TNL_TIL_Mast <- FindNeighbors(TNL_TIL_Mast, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_Mast) <- "celltype.new"
DimPlot(TNL_TIL_Mast, reduction = "umap",split.by = "groups2")

p1 <- DimPlot(TNL_TIL_Mast, reduction = "umap", group.by = "groups2")
p2 <- DimPlot(TNL_TIL_Mast, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL_Mast, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.2",split.by = "orig.ident")
DefaultAssay(TNL_TIL_Mast) <- "integrated"
TNL_TIL_Mast <- FindClusters(object = TNL_TIL_Mast, resolution = 0.1)

TNL_TIL_B@meta.data$B_type[TNL_TIL_B@meta.data$integrated_snn_res.0.1=="0"]=c("B_cells")
TNL_TIL_B$B_type.group.new <- paste(Idents(TNL_TIL_B),TNL_TIL_B$groups2,sep = "_")
table(Idents(TNL_TIL_B))
table(TNL_TIL_B$Macrotype,TNL_TIL_neu$groups)
Idents(TNL_TIL_B) <- "B_type"
Idents(TNL_TIL_neu) <- "integrated_snn_res.0.2"
Idents(TNL_TIL_B) <- "B_type.group.new"
DefaultAssay(TNL_TIL_B) <- "SCT"

genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1','FOXP3', 'CTLA5')#免疫抑制
genes_to_check3 = c('CD4','GZMA', 'CCL5','CD69', 'TCF7','LEF1','FOXP3','IL2RA', 'ICA1','TOX2', 'IL17A','CTSH')
genes_to_check3 = c('CD4','CD8A', 'BTLA','CD69', 'CCR7','CCL5','KLRD1','C1QB', 'GSN','CXCL14', 'GNLY','FCGR3A', 'CD2','TOP2A')
genes_to_check3 = c('CD3E', 'CD3D', 'CD2','CD4','CD8A','MKI67','FTL','IL7R','KIT')
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('CD79A','MS4A1','RFTN1', 'ITGB1','IGHA1','TXNIP','IGHD','CD27')#B type
genes_to_check3 = c('CD19','CD38','CD24', 'CD27','BLIMP1','XBP1','IRF4','SDC1')#B plasmablasts and plasma cells.
genes_to_check3 = c('MS4A1','HLA-DRA', 'JCHAIN','XBP1')
genes_to_check3 = c('CCL3','TGFBR2','IL6ST')#mast炎性因子2
genes_to_check3 = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL3','CCL4','CXCL8','TNFSF10','CCL2','IL6ST','CSF1R','IL13RA1','IL1R2','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2')#NEU炎性p
FeaturePlot(TNL_TIL_Mast, features = genes_to_check3, reduction = "umap",min.cutoff = 0,ncol = 2)
saveRDS(TNL_TIL_Mast, file = "singlecellr/til_tnl/TNL_TIL_B.rds")#####
TNL_TIL_B <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_B.rds")#####

VlnPlot(TNL_TIL_B, features = genes_to_check3, pt.size = 0,ncol = 7, group.by = "Macrotype")
VlnPlot(TNL_TIL_B, features = genes_to_check3, pt.size = 0,ncol = 2, split.by = "groups2",group.by = "Macrotype")
MySeuratWrappers::VlnPlot(TNL_TIL_B,features = genes_to_check3,stacked=T,pt.size=0,group.by = "B_type.group.new",cols = color2)

color2 <- c("#F8766D","#F8766D","#00BFC4","#00BFC4")
TNL_TIL_Mast.markers <- FindAllMarkers(TNL_TIL_Mast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_Mast.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(object = TNL_TIL_neu, features = top10$gene,group.by = "TNL_TIL_neu")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))


#(Endothelial_cells)
TNL_TIL_endo <- subset(TNL_TIL, celltype.new == c("Endothelial_cells"))
TNL_TIL_endo <- ScaleData(TNL_TIL_endo, verbose = FALSE)
TNL_TIL_endo <- RunPCA(TNL_TIL_endo, verbose = FALSE)
TNL_TIL_endo <- RunUMAP(TNL_TIL_endo, reduction = "pca", dims = 1:30)
TNL_TIL_endo <- FindNeighbors(TNL_TIL_endo, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_endo) <- "celltype.new"
DimPlot(TNL_TIL_endo, reduction = "umap",split.by = "groups2")

p1 <- DimPlot(TNL_TIL_endo, reduction = "umap", group.by = "groups2")
p2 <- DimPlot(TNL_TIL_endo, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.3")
DefaultAssay(TNL_TIL_endo) <- "integrated"
TNL_TIL_endo <- FindClusters(object = TNL_TIL_endo, resolution = 0.3)

TNL_TIL_endo@meta.data$endo_type[TNL_TIL_endo@meta.data$integrated_snn_res.0.3=="9"]=c("Venous")
TNL_TIL_endo$endo_type.group.new <- paste(Idents(TNL_TIL_endo),TNL_TIL_endo$groups2,sep = "_")
DimPlot(TNL_TIL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "endo_type.group.new")
table(Idents(TNL_TIL_B))
table(TNL_TIL_endo$endo_type,TNL_TIL_endo$groups)
Idents(TNL_TIL_endo) <- "integrated_snn_res.0.2"
Idents(TNL_TIL_endo) <- "endo_type"
DefaultAssay(TNL_TIL_endo) <- "SCT"
genes_to_check3 = c('PROX1','CCL21', 'PROX1', 'HEY1','CXCL12', 'IGFBP3','CD36', 'CA4','ACKR1')#ymphatic ECs (LECs; CCL21, PROX1).arteries (HEY1, IGFBP3), capillaries (CD36, CA4), veins (ACKR1)
genes_to_check3 = c('ACKR1', 'IGFBP3','CD36','TNFSF10')#ymphatic 
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('APLNR', 'INSR', 'ESM1', 'KDR','VWA1', 'COL4A1')#angiogenic
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('VWF','THBD','EDN1', 'SELP','OLR1')#血管内皮损伤marker
genes_to_check3 = c('CD19','CD38','CD24', 'CD27','BLIMP1','XBP1','IRF4','SDC1')#B plasmablasts and plasma cells.
genes_to_check3 = c('MGP','BMX', 'PRSS23','SLC6A2')
genes_to_check3 = c('TNFSF10')#endo炎性因子2
genes_to_check3 = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL3','CCL4','CXCL8','TNFSF10','CCL2','IL6ST','CSF1R','IL13RA1','IL1R2','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2')#NEU炎性p
FeaturePlot(TNL_TIL_endo, features = genes_to_check3, reduction = "umap",min.cutoff = 0)
MySeuratWrappers::MultiFeaturePlot(TNL_TIL_endo,features = genes_to_check3,ncol=2)
VlnPlot(TNL_TIL_endo, features = genes_to_check3, pt.size = 0, group.by = "endo_type.group.new")
saveRDS(TNL_TIL_endo, file = "singlecellr/til_tnl/TNL_TIL_endo.rds")#####
TNL_TIL_endo <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_endo.rds")#####

TNL_TIL_endo.markers <- FindAllMarkers(TNL_TIL_endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TNL_TIL_endo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_order <- top10[order(top10$cluster),]
DoHeatmap(object = TNL_TIL_endo, features = top10$gene,group.by = "TNL_TIL_neu")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))
#Venous Arterial
Idents(TNL_TIL_endo) <- "endo_type.group.new"
FindMarkers(TNL_TIL_endo, ident.1 = "Arterial_1TNL", ident.2 = "Arterial_2TIL",test.use = "wilcox", verbose = FALSE, min.pct = 0.5) %>% write.csv(file="singlecellr/til_tnl/endo_cells_TNL_TIL_Arterial.csv")



#(fibroblast)
TNL_TIL_fibro <- subset(TNL_TIL, celltype.new == c("Fibroblasts"))
TNL_TIL_fibro <- ScaleData(TNL_TIL_fibro, verbose = FALSE)
TNL_TIL_fibro <- RunPCA(TNL_TIL_fibro, verbose = FALSE)
TNL_TIL_fibro <- RunUMAP(TNL_TIL_fibro, reduction = "pca", dims = 1:30)
TNL_TIL_fibro <- FindNeighbors(TNL_TIL_fibro, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_fibro) <- "celltype.new"
DimPlot(TNL_TIL_fibro, reduction = "umap",split.by = "groups2")

p1 <- DimPlot(TNL_TIL_endo, reduction = "umap", group.by = "groups2")
p2 <- DimPlot(TNL_TIL_endo, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.2")
DefaultAssay(TNL_TIL_endo) <- "integrated"
TNL_TIL_endo <- FindClusters(object = TNL_TIL_endo, resolution = 0.2)

TNL_TIL_B@meta.data$B_type[TNL_TIL_B@meta.data$integrated_snn_res.0.1=="0"]=c("B_cells")
TNL_TIL_B$B_type.group.new <- paste(Idents(TNL_TIL_B),TNL_TIL_B$groups2,sep = "_")
table(Idents(TNL_TIL_B))
table(TNL_TIL_B$Macrotype,TNL_TIL_neu$groups)
Idents(TNL_TIL_endo) <- "integrated_snn_res.0.2"
Idents(TNL_TIL_endo) <- "B_type.group.new"
DefaultAssay(TNL_TIL_endo) <- "SCT"

genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('CD79A','MS4A1','RFTN1', 'ITGB1','IGHA1','TXNIP','IGHD','CD27')#B type
genes_to_check3 = c('CD19','CD38','CD24', 'CD27','BLIMP1','XBP1','IRF4','SDC1')#B plasmablasts and plasma cells.
genes_to_check3 = c('MS4A1','HLA-DRA', 'JCHAIN','XBP1')
genes_to_check3 = c('TGFBR2')#fibro炎性因子2
genes_to_check3 = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL3','CCL4','CXCL8','TNFSF10','CCL2','IL6ST','CSF1R','IL13RA1','IL1R2','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2')#NEU炎性p
FeaturePlot(TNL_TIL_fibro, features = genes_to_check3, reduction = "umap",min.cutoff = 0,split.by = "groups2")
saveRDS(TNL_TIL_fibro, file = "singlecellr/til_tnl/TNL_TIL_endo.rds")#####
TNL_TIL_endo <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_B.rds")#####


#(LEC)
TNL_TIL_LEC <- subset(TNL_TIL, celltype.new == c("LEC"))
TNL_TIL_LEC <- ScaleData(TNL_TIL_LEC, verbose = FALSE)
TNL_TIL_LEC <- RunPCA(TNL_TIL_LEC, verbose = FALSE)
TNL_TIL_LEC <- RunUMAP(TNL_TIL_LEC, reduction = "pca", dims = 1:30)
TNL_TIL_LEC <- FindNeighbors(TNL_TIL_LEC, reduction = "pca", dims = 1:30)
Idents(TNL_TIL_LEC) <- "celltype.new"
DimPlot(TNL_TIL_LEC, reduction = "umap",split.by = "groups2")

p1 <- DimPlot(TNL_TIL_endo, reduction = "umap", group.by = "groups2")
p2 <- DimPlot(TNL_TIL_endo, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(TNL_TIL_endo, reduction = "umap", label = TRUE, repel = TRUE,group.by = "integrated_snn_res.0.2",split.by = "orig.ident")
DefaultAssay(TNL_TIL_endo) <- "integrated"
TNL_TIL_endo <- FindClusters(object = TNL_TIL_endo, resolution = 0.2)

TNL_TIL_B@meta.data$B_type[TNL_TIL_B@meta.data$integrated_snn_res.0.1=="0"]=c("B_cells")
TNL_TIL_B$B_type.group.new <- paste(Idents(TNL_TIL_B),TNL_TIL_B$groups2,sep = "_")
table(Idents(TNL_TIL_B))
table(TNL_TIL_B$Macrotype,TNL_TIL_neu$groups)
Idents(TNL_TIL_endo) <- "integrated_snn_res.0.2"
Idents(TNL_TIL_endo) <- "B_type.group.new"
DefaultAssay(TNL_TIL_endo) <- "SCT"

genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 = c('CD79A','MS4A1','RFTN1', 'ITGB1','IGHA1','TXNIP','IGHD','CD27')#B type
genes_to_check3 = c('CD19','CD38','CD24', 'CD27','BLIMP1','XBP1','IRF4','SDC1')#B plasmablasts and plasma cells.
genes_to_check3 = c('MS4A1','HLA-DRA', 'JCHAIN','XBP1')
genes_to_check3 = c('TGFBR2','IL6ST')#fibro炎性因子2
genes_to_check3 = c('CCL3','CCL21','CCL3L1','CCL4L2','CCL5','XCL1','XCL2','TNFSF10')#炎性因子2
genes_to_check3 = c('CCL3','CCL4','CXCL8','TNFSF10','CCL2','IL6ST','CSF1R','IL13RA1','IL1R2','ACVRL1','TNFRSF1B','TNFRSF1A','IFNGR2','IFNAR1','IL1R1','BMPR2')#NEU炎性p
FeaturePlot(TNL_TIL_LEC, features = genes_to_check3, reduction = "umap",min.cutoff = 0,split.by = "groups2")
saveRDS(TNL_TIL_fibro, file = "singlecellr/til_tnl/TNL_TIL_endo.rds")#####
TNL_TIL_endo <-readRDS(file = "singlecellr/til_tnl/TNL_TIL_B.rds")#####



#KONGJIAN
img_n <- Read10X_Image("singlecellr/spaceranger_result/D1/spatial/", image.name = "tissue_hires_image.png")
TNL_v <- Load10X_Spatial("singlecellr/spaceranger_result/D1/",filename = "filtered_feature_bc_matrix.h5", assay = "Spatial",image = img_n,
                        filter.matrix = TRUE)
TNL_v@images[["slice1"]]@scale.factors[["lowres"]] <- TNL_v@images[["slice1"]]@scale.factors[["hires"]]

img_i <- Read10X_Image("singlecellr/spaceranger_result/B1/spatial/", image.name = "tissue_hires_image.png")
TIL_v <- Load10X_Spatial("singlecellr/spaceranger_result/B1/",filename = "filtered_feature_bc_matrix.h5", assay = "Spatial",image = img_i,
                        filter.matrix = TRUE)
TIL_v@images[["slice1"]]@scale.factors[["lowres"]] <- TIL_v@images[["slice1"]]@scale.factors[["hires"]]

Idents(TNL_v) <- "orig.ident"
Idents(TIL_v) <- "orig.ident"

VlnPlot(TNL_v, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
VlnPlot(TIL_v, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
FeatureScatter(TNL_v,feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", cols = 'blue')
FeatureScatter(TIL_v,feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", cols = 'blue')
SpatialFeaturePlot(TNL_v, features = "nCount_Spatial")+theme(legend.position = "right")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))
SpatialFeaturePlot(TIL_v, features = "nCount_Spatial")+theme(legend.position = "right")+scale_fill_gradientn(colors = c( "blue", "white", "firebrick"))

TNL_v <- SCTransform(TNL_v, assay = "Spatial", verbose = FALSE)
TIL_v <- SCTransform(TIL_v, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(TNL_v, features = c("CD163", "KRT8","TAGLN"))
SpatialFeaturePlot(TIL_v, features = c("KRT7", "KRT8","PTPRC","KRT18"),alpha = c(0.1, 1))

TNL_v <- RunPCA(TNL_v, assay = "SCT", verbose = FALSE)
TNL_v <- FindNeighbors(TNL_v, reduction = "pca", dims = 1:30)
TNL_v <- FindClusters(TNL_v, verbose = FALSE)
TNL_v <- RunUMAP(TNL_v, reduction = "pca", dims = 1:30)

TIL_v <- RunPCA(TIL_v, assay = "SCT", verbose = FALSE)
TIL_v <- FindNeighbors(TIL_v, reduction = "pca", dims = 1:30)
TIL_v <- FindClusters(TIL_v, verbose = FALSE)
TIL_v <- RunUMAP(TIL_v, reduction = "pca", dims = 1:30)

p1 <- DimPlot(TNL_v, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(TNL_v, label = TRUE, label.size = 3)
p1 + p2

LinkedDimPlot(TNL_v)

TNL_v <- FindSpatiallyVariableFeatures(TNL_v, assay = "SCT", features = VariableFeatures(TNL_v)[1:1000], 
                                      selection.method = "markvariogram")
TIL_v <- FindSpatiallyVariableFeatures(TIL_v, assay = "SCT", features = VariableFeatures(TIL_v)[1:1000], 
                                      selection.method = "markvariogram")
VariableFeaturePlot(
  TIL_v,
  cols = c("black", "red"),
  pt.size = 1
)
top.features <- head(SpatiallyVariableFeatures(TIL_v, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(TIL_v, features = top.features, ncol = 3, alpha = c(0.1, 1))
saveRDS(TNL_v, file = "singlecellr/np_tnl/TNL_v_h.rds")
saveRDS(TIL_v, file = "singlecellr/np_tnl/TIL_v_h.rds")

TNL_v <- readRDS("singlecellr/np_tnl/TNL_v.rds")####
TIL_v <- readRDS("singlecellr/np_tnl/TIL_v.rds")####
DefaultAssay(TIL_v) <- "SCT"
DefaultAssay(TNL_v) <- "SCT"
mapgene = c('IL1A', 'IL1B','IL6','CXCL8','TNF','LTA','CCL4','CCL2')#炎性因子
genes_to_check3 <- c('CD79A','VWF','PECAM1','LUM','COL1A1','TFF3','CCL21','TPSB2','CD14','LYZ')
genes_to_check3 <- c('FCGR3B','CXCL8','KLRD1','HBB','ACTA2','TAGLN','CD3D','IL7R','KRT7')
genes_to_check3 <- c('PECAM1','VWF','TAGLN','ACTA2','DCN','LUM','TFF3','CCL21','LYZ','PTPRC','HBB','HLA-G')
genes_to_check3 = c('MS4A1','JCHAIN')
genes_to_check3 <- c('CD79A','VWF','PECAM1','LUM','COL1A1','TFF3','CCL21','TPSB2','CD14','LYZ','FCGR3B','CXCL8','KLRD1','HBB','ACTA2','TAGLN','CD3D','IL7R','KRT7')
SpatialFeaturePlot(TIL_v, features = 'HBB', alpha = c(0.1, 1), min.cutoff = 0)+scale_fill_gradientn(colors = c( "white", "blue"))+theme(legend.title = element_text(size = 20))
SpatialFeaturePlot(TNL_v, features = 'CD79A', alpha = c(0.1, 1), min.cutoff = 0)+scale_fill_gradientn(colors = c( "white", "blue"))+theme(legend.title = element_text(size = 20))

SpatialFeaturePlot(TIL_v, features = 'KIT')+scale_fill_gradientn(colors = c( "white", "blue"))
SpatialFeaturePlot(TNL_v, features = 'KIT')
SpatialPlot(TNL_v, features = genes_to_check3)
FeaturePlot(TIL_v, features = c("IL1B","PTGS2"), blend = T)
VlnPlot(TIL_v,features = c("IL1B","PTGS2"))



TNL_TIL_v <- merge(TNL_v,TIL_v)
SpatialFeaturePlot(TNL_TIL_v, features = 'ACTA2', alpha = c(0.1, 1), min.cutoff = 0)+scale_fill_gradientn(colors = c( "white", "blue"))


SpatialDimPlot(TNL_v,crop = TRUE)
TNL <- subset(TNL_TIL, groups == "TNL")
TIL <- subset(TNL_TIL, groups == "TIL")
TNL <- SCTransform(TNL)
TIL <- SCTransform(TIL)
DefaultAssay(TNL) <- "SCT"
DefaultAssay(TIL) <- "SCT"
Idents(TNL) <- "celltype.new"
Idents(TIL) <- "celltype.new"

anchors_TNL <- FindTransferAnchors(reference = TNL, query = TNL_v, normalization.method = "SCT")
predictions.assay_TNL <- TransferData(anchorset = anchors_TNL, refdata = TNL$celltype.new, prediction.assay = TRUE, 
                                     weight.reduction = TNL_v[["pca"]], dims = 1:30)
anchors_TIL <- FindTransferAnchors(reference = TIL, query = TIL_v, normalization.method = "SCT")
predictions.assay_TIL <- TransferData(anchorset = anchors_TIL, refdata = TIL$celltype.new, prediction.assay = TRUE, 
                                     weight.reduction = TIL_v[["pca"]], dims = 1:30)
##Endothelial_cells Fibroblasts SMC Macrophages Red_blood_cells Trophectoderm_cells LEC B_cells Mast_cells
seleccells <- c("SMC","LEC","Endothelial-cells","Fibroblasts","Monocytic","T-cells")
seleccells <- c("SMC","LEC","Endothelial-cells","Fibroblasts","Monocytic","T-cells","B-cells","Mast-cells","Neutrophils","NK-cells","Trophectoderm-cells","Red-blood-cells")
TNL_v[["predictions"]] <- predictions.assay_TNL
DefaultAssay(TNL_v) <- "predictions"
SpatialFeaturePlot(TNL_v,features = seleccells,ncol = 3)

seleccells <- c("SMC","Endothelial-cells","Fibroblasts")
TIL_v[["predictions"]] <- predictions.assay_TIL
DefaultAssay(TIL_v) <- "predictions"
SpatialFeaturePlot(TIL_v, features = seleccells,ncol = 4)

p <- SpatialFeaturePlot(TIL_v,features = seleccells,ncol = 3)
ggsave("til-v.jpg", p , width = 15, height = 12, dpi = 300)
p <- SpatialFeaturePlot(TNL_v,features = seleccells,ncol = 3)
ggsave("tnl-v.jpg", p , width = 15, height = 12, dpi = 300)

#AddModuleScore
DefaultAssay(TNL_v) <- "SCT"
DefaultAssay(TIL_v) <- "SCT"
cd_features <- list(c(
  'CD14',
  'CD68',
  'CD163'))
TIL_v2 <-AddModuleScore(TIL_v,features = cd_features,name = "SMC")
TIL_v2@meta.data <- merge(TIL_v2@meta.data,m,by = 'cellnames')
Idents(TIL_v2) <- 'SMC1'
SpatialFeaturePlot(TIL_v2)
SpatialDimPlot(TIL_v2)

#配体受体 
TIL_v@reductions$spatial = TIL_v@reductions$umap
TIL_v@reductions$spatial@key = 'spatial_'
TIL_v@reductions$spatial@cell.embeddings = as.matrix(TIL_v@images$slice1@coordinates[,c(3,2)])
TIL_v@reductions$spatial@cell.embeddings[,1] = -TIL_v@reductions$spatial@cell.embeddings[,1]
colnames(TIL_v@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
FeaturePlot(TIL_v, features = c("IL1B","MRC1"), blend = T,cols = c('lightgrey','blue','red'),reduction = 'spatial',combine = T)+coord_flip()

TNL_v@reductions$spatial = TNL_v@reductions$umap
TNL_v@reductions$spatial@key = 'spatial_'
TNL_v@reductions$spatial@cell.embeddings = as.matrix(TNL_v@images$slice1@coordinates[,c(3,2)])
TNL_v@reductions$spatial@cell.embeddings[,1] = -TNL_v@reductions$spatial@cell.embeddings[,1]
colnames(TNL_v@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
FeaturePlot(TNL_v, features = c("IL1B","MRC1"), blend = T,cols = c('lightgrey','blue','red'),reduction = 'spatial',combine = T)+coord_flip()


#kongjian联合
TNL_TIL_v <- merge(TNL_v,TIL_v)

SpatialDimPlot(TNL_TIL_v)
saveRDS(NP_TP_v, file = "singlecellr/np_tnl/NP_TP_v.rds")
NP_TP_v <- readRDS("singlecellr/np_tnl/NP_TP_v.rds")

genes_to_check3 = c('CD163', 'CD86', 'CD14', 'CD209','IL1B', 'MRC1')#巨噬
genes_to_check3 = c('CD3E', 'CD3D', 'CD7', 'KLRD1','KLRB1','CD4')#T NK
genes_to_check3 = c('CD79A', 'MS4A1', 'PTPRC', 'HBB','HBA1')#B 红细胞
genes_to_check3 = c('HLA-G', 'KRT18', 'KRT8', 'GATA3','KRT17', 'KRT7')#滋养
genes_to_check3 = c('TAGLN', 'DES', 'CALD1', 'PLN','ACTA2', 'ACTA1')#smc
genes_to_check3 = c('EPCAM', 'DCN', 'LUM', 'PECAM1','VWF', 'TFF3')#上皮 纤维 内皮 淋巴内皮
genes_to_check3 = c('IL2RA', 'ENTPD1', 'NT5E', 'EBI3','CD274', 'PDCD1')#免疫抑制
genes_to_check3 = c('SLC2A1', 'HK2', 'HIF1A', 'LDHA','OXTR', 'GJA1')#代谢
genes_to_check3 = c('STEAP4')#SMC biaoji, 'C7', 'CCN1', 'ADIRF'
genes_to_check3 = c('IGFBP7', 'ADIRF', 'PLAC9', 'TIMP3','OXTR', 'GJA1')
genes_to_check3 = c('SLC2A1', 'HK2', 'HIF1A', 'LDHA', 'CD79A', 'JCHAIN', 'CD3D', 'TAGLN', 'PTPRC', 'CD44', 'FGF7', 'DES')
genes_to_check3 <- c('LUM','CD14','TAGLN','CD3E','TFF3','VWF','PECAM1','KLRD1','TPSB2','CD79A','FCGR3B','PAPPA2')
genes_to_check3 <- c('C1QA', 'MRC1')
genes_to_check3 <- c('OXTR', 'GJA1','ELANE', 'IFNG','SERPINE1')
genes_to_check3 = c('C1QA', 'MRC1','S100A8','FCN1','CALD1','RGS5','MMRN1','CAVIN2','MKI67','TOP2A')#
genes_to_check3 = c('CD14','CD163', 'CD86')
genes_to_check3 = c('CD3E','CD4','CD8A')
VlnPlot(NP_TP_v,features = genes_to_check3, pt.size = 0,ncol = 2,split.by = "group")
FeaturePlot(NP_TP_v, features = genes_to_check3, reduction = "umap")
SpatialFeaturePlot(TP_v, features = genes_to_check3, min.cutoff = 0,max.cutoff = 5, alpha = c(0.1, 1))


#合并亚群
111 <- left_join(data.frame(TNL_TIL_clear@meta.data),data.frame(TNL_TIL_TNK@meta.data))
TNL_TIL_clear <- AddMetaData(TNL_TIL_clear,TNL_TIL_TNK@meta.data)
DimPlot(TNL_TIL_clear, reduction = "umap", group.by = "Ttype",label = T,repel = T)

TNL_TIL_SMC@meta.data$Ttype <- TNL_TIL_SMC@meta.data$SMCtype
TNL_TIL_macro@meta.data$Ttype <- TNL_TIL_macro@meta.data$Macrotype
TNL_TIL_neu@meta.data$Ttype <- TNL_TIL_neu@meta.data$TNL_TIL_neu
TNL_TIL_B@meta.data$Ttype <- TNL_TIL_B@meta.data$B_type
TNL_TIL_clear_nonimmune <- subset(TNL_TIL_clear,celltype.new == c("Endothelial_cells","Fibroblasts","LEC","Mast_cells","Red_blood_cells","Trophectoderm_cells"))
TNL_TIL_clear_nonimmune <- subset(TNL_TIL_clear, celltype.new == c("B_cells"), invert = T)
TNL_TIL_clear_nonimmune <- subset(TNL_TIL_clear_nonimmune, celltype.new == c("SMC"), invert = T)
saveRDS(TNL_TIL_clear_nonimmune, file = "singlecellr/til_tnl/TNL_TIL_clear_nonimmune.rds")
TNL_TIL_clear_nonimmune<- readRDS(file = "singlecellr/til_tnl/TNL_TIL_clear_nonimmune.rds")
#"Monocytic","Neutrophils","NK_cells","T_cells","SMC"
TNL_TIL_clear_nonimmune@meta.data$Ttype <- TNL_TIL_clear_nonimmune@meta.data$celltype.new
metadataall <- rbind(TNL_TIL_TNK@meta.data[,c("orig.ident","Ttype")],TNL_TIL_SMC@meta.data[,c("orig.ident","Ttype")],TNL_TIL_macro@meta.data[,c("orig.ident","Ttype")],TNL_TIL_neu@meta.data[,c("orig.ident","Ttype")],TNL_TIL_B@meta.data[,c("orig.ident","Ttype")],TNL_TIL_clear_nonimmune@meta.data[,c("orig.ident","Ttype")])
TNL_TIL_clear2 <- AddMetaData(TNL_TIL_clear,metadataall)
DimPlot(TNL_TIL_clear2, reduction = "umap",split.by = "groups2",group.by = "Ttype")
saveRDS(TNL_TIL_clear2, file = "singlecellr/til_tnl/TNL_TIL_clear2.rds")#######
TNL_TIL_clear2 <- readRDS(file = "singlecellr/til_tnl/TNL_TIL_clear2.rds")#######
TNL_TIL_clear2_av <- AverageExpression(TNL_TIL_clear2,group.by = "Ttype",assays = "SCT")
TNL_TIL_clear2_av <-TNL_TIL_clear2_av[[1]]
head(TNL_TIL_clear2_av)
cg <- names(tail(sort(apply(TNL_TIL_clear2_av,1,sd)),1000))
pheatmap::pheatmap(cor(TNL_TIL_clear2_av[cg,],method = "spearman"))



#SPOTLIGHT 要用findmarker中的roc方法
library(SPOTlight)
markers_all <- Seurat::FindAllMarkers(object = TNL_TIL,assay = 'SCT', verbose = TRUE,test.use = "roc",
                                              only.pos = TRUE,
                                              logfc.threshold = 1,
                                              min.pct = 0.8)
write.csv(markers_all,file="singlecellr/til_tnl/marker_TNL_TIL_13_roc.csv")
markers_all <- read.csv("singlecellr/til_tnl/marker_TNL_TIL_13_roc.csv")

sc_p <- DimPlot(TIL,reduction = 'umap',label = T,group.by = 'celltype.new')
st_p <- Seurat::SpatialDimPlot(TIL_v,label = T)
sc_p+st_p

set.seed(123)
TNL_3000 <- subset(TIL, downsample = 3000)
TIL_3000 <- subset(TNL, downsample = 3000)

spotlight_TIL <- SPOTlight(x=TIL_3000,
                           y = TIL_v@assays$Spatial@counts,
                           groups = TIL_3000$celltype.new,
                           mgs = markers_all,
                           weight_id = "myAUC",
                           group_id = "cluster",
                           gene_id = "gene")
saveRDS(spotlight_TIL, file = "singlecellr/til_tnl/spotlight_TIL.rds")
head(mat <- spotlight_TIL$mat)[, seq_len(3)]
mod <- spotlight_TIL$NMF

plotTopicProfiles(
  x = mod,
  y = TIL_300$celltype.new,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

plotTopicProfiles(
  x = mod,
  y = TIL_300$celltype.new,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 6)

library(NMF)
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)
plotCorrelationMatrix(mat)
plotSpatialScatterpie(
  x = TIL_v,
  y = mat)
plotCorrelationMatrix(mat)
plotInteractions(mat, "network")




#CellChat
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(mindr)
library(NMF)
library(circlize)
library(Seurat)
library(ComplexHeatmap)
TNL <- subset(TNL_TIL_clear2, groups == "TNL")
TIL <- subset(TNL_TIL_clear2, groups == "TIL")
DefaultAssay(TNL) <- "SCT"
DefaultAssay(TIL) <- "SCT"

cellchat_tnl <- createCellChat(object = TNL, group.by = "Ttype", assay = "RNA")
cellchat_til <- createCellChat(object = TIL, group.by = "Ttype", assay = "RNA")
groupSize_tnl <- as.numeric(table(cellchat_tnl@idents))
groupSize_til <- as.numeric(table(cellchat_til@idents))
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB

CellChatDB_interaction <-CellChatDB$interaction
write.csv(CellChatDB_interaction, file = "singlecellr/np_tnl/CellChatDB_interaction.csv")
cellchat_tnl@DB <- CellChatDB.use
cellchat_til@DB <- CellChatDB.use

cellchat_tnl <- subsetData(cellchat_tnl) 
cellchat_tnl <- identifyOverExpressedGenes(cellchat_tnl)
cellchat_tnl <- identifyOverExpressedInteractions(cellchat_tnl)
cellchat_tnl <- projectData(cellchat_tnl, PPI.human)  
cellchat_tnl <- computeCommunProb(cellchat_tnl)
cellchat_tnl <- filterCommunication(cellchat_tnl, min.cells = 10)
cellchat_tnl <- computeCommunProbPathway(cellchat_tnl)
cellchat_tnl <- aggregateNet(cellchat_tnl)

cellchat_tnl <- netAnalysis_computeCentrality(cellchat_tnl, slot.name = "netP")
cellchat_tnl@netP$pathways
cellchat_tnl@LR$LRsig$pathway_name
cellchat_tnl@LR$LRsig$antagonist
selectK(cellchat_tnl, pattern = "outgoing")
nPatterns = 3
cellchat_tnl <- identifyCommunicationPatterns(cellchat_tnl, pattern = "outgoing", k = nPatterns)
netVisual_heatmap(cellchat_tnl, signaling = c("BMP"), color.heatmap = "Reds")


cellchat_til <- subsetData(cellchat_til) 
cellchat_til <- identifyOverExpressedGenes(cellchat_til)
cellchat_til <- identifyOverExpressedInteractions(cellchat_til)
cellchat_til <- projectData(cellchat_til, PPI.human)  
cellchat_til <- computeCommunProb(cellchat_til)
cellchat_til <- filterCommunication(cellchat_til, min.cells = 10)
cellchat_til <- computeCommunProbPathway(cellchat_til)
cellchat_til <- aggregateNet(cellchat_til)

cellchat_til <- netAnalysis_computeCentrality(cellchat_til, slot.name = "netP")
cellchat_til@netP$pathways
cellchat_til@LR$LRsig$pathway_name
cellchat_til@LR$LRsig$antagonist
selectK(cellchat_til, pattern = "outgoing")
nPatterns = 3
cellchat_til <- identifyCommunicationPatterns(cellchat_til, pattern = "outgoing", k = nPatterns)
netVisual_heatmap(cellchat_til, signaling = c("CCL"), color.heatmap = "Reds")

saveRDS(cellchat_tnl, file = "singlecellr/til_tnl/cellchat_tnl.rds")
saveRDS(cellchat_til, file = "singlecellr/til_tnl/cellchat_til.rds")
cellchat_tnl <- readRDS("singlecellr/til_tnl/cellchat_tnl.rds")#####
cellchat_til <- readRDS("singlecellr/til_tnl/cellchat_til.rds")#####

#merge
til_tnl_object.list <- list(tnl = cellchat_tnl, til = cellchat_til)

cellchat_til_tnl_merge <- mergeCellChat(til_tnl_object.list, add.names = names(til_tnl_object.list))
cellchat_til_tnl_merge
gg1 <- compareInteractions(cellchat_til_tnl_merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_til_tnl_merge, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


netVisual_diffInteraction(cellchat_til_tnl_merge, comparison = c(1, 2), weight.scale = T)
par(mfrow = c(1,1))
netVisual_heatmap(cellchat_tnl)
netVisual_heatmap(cellchat_til)
netVisual_heatmap(cellchat_til_tnl_merge, measure = "weight")
netVisual_heatmap(cellchat_til_tnl_merge)

#选11个 #红色]（或[蓝色]边表示信号在第二个数据集中增加或[减少]）
cellchat_tnl_11<- subsetCellChat(object = cellchat_tnl, idents.use = c("Red_blood_cells", "Unknown"),invert = T)
cellchat_til_11<- subsetCellChat(object = cellchat_til, idents.use = c("Red_blood_cells", "Unknown"),invert = T)
netVisual_heatmap(cellchat_tnl_11,signaling = "CCL")
netVisual_heatmap(cellchat_til_11)
netVisual_circle(cellchat_tnl_11@net$count, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_tnl_11@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat_til_11@net$count, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_til_11@net$weight, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



til_tnl_11_object.list <- list(TNL = cellchat_tnl_11, TIL = cellchat_til_11)
cellchat_til_tnl_11_merge <- mergeCellChat(til_tnl_11_object.list, add.names = names(til_tnl_11_object.list))
cellchat_til_tnl_11_merge
netVisual_diffInteraction(cellchat_til_tnl_11_merge,weight.scale = T,label.edge = T)
netVisual_diffInteraction(cellchat_til_tnl_11_merge,weight.scale = T,label.edge = T,sources.use ="SMC")
netVisual_diffInteraction(cellchat_til_tnl_11_merge, weight.scale = T, measure = "weight",label.edge = T)
netVisual_circle(cellchat_tnl_11@net$count,sources.use ="SMC",weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_aggregate(cellchat_tnl_11, signaling = "CCL")

saveRDS(cellchat_til_tnl_11_merge, file = "singlecellr/til_tnl/cellchat_til_tnl_11_merge.rds")
cellchat_til_tnl_11_merge <- readRDS("singlecellr/til_tnl/cellchat_til_tnl_11_merge.rds")#####

gg1 <- rankNet(cellchat_til_tnl_11_merge, sources.use = "SMC",mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_til_tnl_11_merge, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
rankNet(cellchat_til_tnl_11_merge, sources.use = "SMC",mode = "comparison",stacked = T, do.stat = F)
netVisual_bubble(cellchat_til_tnl_11_merge,targets.use = "SMC",comparison = c(1,2))
netVisual_bubble(cellchat_til_tnl_11_merge,targets.use = "SMC",comparison = c(1,2),signaling = c("FGF","IGF","CXCL","COLLAGEN","THBS","CCL","THBS"),thresh = 0.01)
num.link <- sapply(til_tnl_11_object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(til_tnl_11_object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(til_tnl_11_object.list[[i]], title = names(til_tnl_11_object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
netAnalysis_signalingRole_scatter(cellchat_tnl_11,signaling = pathway.union)

cellchat_til_tnl_11_merge <- computeNetSimilarityPairwise(cellchat_til_tnl_11_merge, type = "functional")
cellchat_til_tnl_11_merge <- netEmbedding(cellchat_til_tnl_11_merge, type = "functional",umap.method = c("uwot"))
cellchat_til_tnl_11_merge <- netClustering(cellchat_til_tnl_11_merge, type = "functional")
netVisual_embeddingPairwise(cellchat_til_tnl_11_merge, type = "functional", label.size = 3.5)

rankSimilarity(cellchat_til_tnl_11_merge, type = "functional")
netVisual_aggregate(cellchat_til_11, signaling = c("FN1"), layout = "circle",sources.use = "SMC")
##"FGF","IGF","CXCL","COLLAGEN","THBS","CCL","THBS"
netVisual_chord_gene(cellchat_til_11,sources.use = "SMC",slot.name = 'net',signaling = "VEGF")
netVisual_chord_gene(cellchat_til_11,targets.use = "SMC",slot.name = 'net',signaling = "VEGF")
pathway.union <-c("FGF","IGF","CXCL","COLLAGEN","THBS","CCL","THBS","ANXA1","VISFATIN","MIF","ADGRE5","MPZ","TENASCIN","TWEAK")
pathway.union <- union(til_tnl_11_object.list[[1]]@netP$pathways, til_tnl_11_object.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(til_tnl_11_object.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(til_tnl_11_object.list)[1], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(til_tnl_11_object.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(til_tnl_11_object.list)[2], width = 5, height = 6)
ht1 + ht2
netVisual_bubble(cellchat_til_tnl_11_merge,signaling = pathway.union, sources.use = "SMC",comparison = c(1,2),max.dataset = 2, remove.isolate = T)

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
cellchat_tnl_11 <- netAnalysis_computeCentrality(cellchat_np_5, slot.name = "netP")
cellchat_til_11 <- netAnalysis_computeCentrality(cellchat_tp_5, slot.name = "netP")

pathway.union <- union(cellchat_tnl_11@netP$pathways, cellchat_til_11@netP$pathways)
netAnalysis_signalingRole_heatmap(cellchat_tnl_11, width = 5, height = 16,
                                  pattern = "all")
netAnalysis_signalingRole_heatmap(cellchat_til_11, width = 5, height = 16,
                                  pattern = "outgoing", signaling = pathway.union)

gg1 <- rankNet(cellchat_til_tnl_11_merge, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_til_tnl_11_merge, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)
# combining all the identified signaling pathways from different datasets 
netAnalysis_signalingRole_heatmap(cellchat_tnl_11, width = 5, height = 12,
                                  pattern = "incoming")
netAnalysis_signalingRole_heatmap(cellchat_til_11, width = 5, height = 18,
                                  pattern = "incoming")
netAnalysis_signalingRole_heatmap(cellchat_tnl_11, width = 5, height = 12,
                                  pattern = "outgoing",color.heatmap = "GnBu")
netAnalysis_signalingRole_heatmap(cellchat_til_11, width = 5, height = 18,
                                  pattern = "outgoing",color.heatmap = "GnBu")
netVisual_bubble(cellchat_tnl_11, sources.use = "SMC", angle.x = 45, thresh = 0.01)
netVisual_bubble(cellchat_til_11, sources.use = "SMC",, angle.x = 45, remove.isolate = T)

cellchat_til_tnl_11_merge <- identifyOverExpressedGenes(cellchat_til_tnl_11_merge, 
                                                        pos.dataset = "TIL", features.name = "TIL",
                                                  only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat_til_tnl_11_merge, features.name = "TIL")
net.up <- subsetCommunication(cellchat_til_tnl_11_merge, net = net, datasets = "TIL",ligand.logFC = 0.2, receptor.logFC = 0.2)
net.down <- subsetCommunication(cellchat_til_tnl_11_merge, net = net, datasets = "TIL",ligand.logFC = -0.2)
####TNP = cellchat_tnl_11, TIP = cellchat_til_11
netVisual_chord_gene(cellchat_til_tnl_11_merge, targets.use = "SMC",
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5 )
netVisual_chord_gene(cellchat_til_tnl_11_merge, sources.use = "SMC",
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5 )
netVisual_chord_gene(cellchat_til_tnl_11_merge, targets.use = "SMC",
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5 )
netVisual_chord_gene(cellchat_til_tnl_11_merge, sources.use = "SMC",
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5 )

#选3个
cellchat_tnl_3<- subsetCellChat(object = cellchat_tnl, idents.use = c("SMC", "Endothelial_cells","Monocytic"))
cellchat_til_3<- subsetCellChat(object = cellchat_til, idents.use = c("SMC", "Endothelial_cells","Monocytic"))

netVisual_circle(cellchat_tnl_3@net$count, weight.scale = T, label.edge= T, title.name = "Number of interactions",vertex.size.max=20,arrow.size=2,edge.label.cex=1)
netVisual_circle(cellchat_til_3@net$count, weight.scale = T, label.edge= T, title.name = "Number of interactions",vertex.size.max=20,arrow.size=2,edge.label.cex=1)

til_tnl_3_object.list <- list(TNP = cellchat_tnl_3, TIL = cellchat_til_3)
cellchat_til_tnl_3_merge <- mergeCellChat(til_tnl_3_object.list, add.names = names(til_tnl_3_object.list))
cellchat_til_tnl_3_merge
netVisual_diffInteraction(cellchat_til_tnl_3_merge,weight.scale = T,label.edge = T)
netVisual_circle(cellchat_til_tnl_3_merge@net$count,sources.use ="SMC",weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_aggregate(cellchat_tnl_11, signaling = pathway.union,show.legend=T)
netVisual_aggregate(cellchat_til_11, signaling = pathway.union)
pathway.union <-c("CXCL","CCL","MHC-I","MHC-II","IFN-I","IL1","IL4","IL10","IL2","COMPLEMENT","IL17","IL12","IL6","TNF","IFN-II","IL16","PD-L1","XCR")
pathway.union2 <-c("CXCL","CCL","IL1","IL6","TNF","XCR")
rankNet(cellchat_til_tnl_11_merge, sources.use = "SMC",mode = "comparison",stacked = T, do.stat = F)
netVisual_bubble(cellchat_til_tnl_3_merge,targets.use = "SMC",comparison = c(1,2))
netVisual_bubble(cellchat_til_tnl_3_merge,targets.use = "SMC",comparison = c(1,2),signaling = c("FGF","IGF","CXCL","COLLAGEN","THBS","CCL","THBS"),thresh = 0.01)
netVisual_bubble(cellchat_til_tnl_3_merge, sources.use = celluse, targets.use = celluse, comparison = c(1, 2), angle.x = 45,thresh = 0.01,max.dataset = 2,title.name = "Increased signaling in TIL",remove.isolate = T)
netVisual_bubble(cellchat_til_tnl_3_merge, sources.use = celluse, targets.use = celluse, comparison = c(1, 2), angle.x = 45,thresh = 0.01,max.dataset = 1,title.name = "Increased signaling in TIL",remove.isolate = T)

cellchat_til_tnl_3_merge <- identifyOverExpressedGenes(cellchat_til_tnl_3_merge, group.dataset = "datasets", pos.dataset = "TIL", features.name = "TIL", only.pos = FALSE, thresh.pc = 0.25, thresh.fc = 0.58492501, thresh.p = 0.05)
net <- netMappingDEG(cellchat_til_tnl_3_merge, features.name = "TIL")
net.up <- subsetCommunication(cellchat_til_tnl_3_merge, net = net, datasets = "TIL",ligand.logFC = 0.58492501, receptor.logFC = 0.58492501)
net.down <- subsetCommunication(cellchat_til_tnl_3_merge, net = net, datasets = "TIL",ligand.logFC = -0.58492501, receptor.logFC = -0.58492501)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_til_tnl_3_merge)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_til_tnl_3_merge)
pairLR.use.up = net.up[, "interaction_name", drop = F]
netVisual_bubble(cellchat_til_tnl_3_merge, pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 45, remove.isolate = T)
pairLR.use.down = net.down[, "interaction_name", drop = F]
netVisual_bubble(cellchat_til_tnl_3_merge, pairLR.use = pairLR.use.down, comparison = c(1, 2),  angle.x = 45, remove.isolate = T)

netVisual_chord_gene(cellchat_tnl_3,slot.name = 'net', net = net.down)
netVisual_chord_gene(cellchat_til_3,slot.name = 'net', net = net.up)





#10个
cellchat_tnl_10<- subsetCellChat(object = cellchat_tnl, idents.use = c("Red_blood_cells", "Unknown","Trophectoderm_cells"),invert = T)
cellchat_til_10<- subsetCellChat(object = cellchat_til, idents.use = c("Red_blood_cells", "Unknown","Trophectoderm_cells"),invert = T)
netVisual_aggregate(cellchat_tnl_10, signaling = pathway.union2,arrow.size=0.5,weight.scale = T)
netVisual_aggregate(cellchat_til_10, signaling = pathway.union2,arrow.size=0.5)
netAnalysis_signalingRole_scatter(cellchat_til_10, signaling = "CCL")

tnl_til_10.list <- list(TNL = cellchat_tnl_10, TIL = cellchat_til_10)
tnl_til_10_cellchat <- mergeCellChat(tnl_til_10.list, add.names = names(tnl_til_10.list))
netVisual_diffInteraction(tnl_til_10_cellchat, weight.scale = T)
rankNet(tnl_til_10_cellchat, mode = "comparison",signaling = pathway.union2,stacked = T, do.stat = TRUE)
netAnalysis_signalingRole_heatmap(cellchat_tnl_10,pattern = "outgoing",signaling = pathway.union)
netVisual_bubble(tnl_til_10_cellchat,  comparison = c(1, 2), angle.x = 45,signaling = pathway.union2)
saveRDS(tnl_til_10_cellchat, file = "singlecellr/til_tnl/tnl_til_10_cellchat.rds")
tnl_til_10_cellchat <- readRDS("singlecellr/til_tnl/tnl_til_10_cellchat.rds")#####
celluse <- c("Endothelial_cells","Monocytic","SMC")
netVisual_bubble(tnl_til_10_cellchat, sources.use = celluse, targets.use = celluse, comparison = c(1, 2), angle.x = 45,signaling = pathway.union2)

netAnalysis_river(cellchat_tnl_10, pattern = "outgoing")
netVisual_aggregate

#10个subtype
cellchat_tnl_10s<- subsetCellChat(object = cellchat_tnl, idents.use = c("Red_blood_cells", "Unknown","Trophectoderm_cells"),invert = T)
cellchat_til_10s<- subsetCellChat(object = cellchat_til, idents.use = c("Red_blood_cells", "Unknown","Trophectoderm_cells"),invert = T)
netVisual_heatmap(cellchat_tnl_10s, signaling = "MHC-I")
netVisual_heatmap(cellchat_til_10s, signaling = "MHC-I")
netVisual_aggregate(cellchat_tnl_10s, signaling = pathway.union2,arrow.size=0.5)
netVisual_aggregate(cellchat_til_10s, signaling = pathway.union2,arrow.size=0.5)
pathway.union2 <-c("CXCL","CCL","IL1","IL6","TNF","XCR")

tnl_til_10s.list <- list(TNL = cellchat_tnl_10s, TIL = cellchat_til_10s)
tnl_til_10s_cellchat <- mergeCellChat(tnl_til_10s.list, add.names = names(tnl_til_10s.list))
netVisual_diffInteraction(tnl_til_10s_cellchat, weight.scale = T)
rankNet(tnl_til_10s_cellchat, mode = "comparison",signaling = pathway.union2,stacked = T, do.stat = TRUE)
netAnalysis_signalingRole_heatmap(cellchat_tnl_10,pattern = "outgoing",signaling = pathway.union)
netVisual_bubble(tnl_til_10s_cellchat,  comparison = c(1, 2), angle.x = 45,signaling = pathway.union2)
saveRDS(tnl_til_10s_cellchat, file = "singlecellr/til_tnl/tnl_til_10s_cellchat.rds")
tnl_til_10s_cellchat <- readRDS("singlecellr/til_tnl/tnl_til_10s_cellchat.rds")#####
celluse <- c("Endothelial_cells","Monocytic","SMC")
netVisual_bubble(tnl_til_10s_cellchat, sources.use = celluse, targets.use = celluse, comparison = c(1, 2), angle.x = 45,signaling = pathway.union2)
netVisual_heatmap(tnl_til_10s_cellchat, comparison = c(1, 2),signaling = "IL6")
rankNet(tnl_til_10s_cellchat, mode = "comparison", stacked = T, do.stat = TRUE,signaling = pathway.union2)
netVisual_chord_gene(cellchat_tnl_10s,signaling = pathway.union2)

tnl_til_10s_cellchat <- identifyOverExpressedGenes(tnl_til_10s_cellchat, group.dataset = "datasets", pos.dataset = "TIL", features.name = "TIL", only.pos = FALSE, thresh.pc = 0.25, thresh.fc = 0.58492501, thresh.p = 0.05)
net <- netMappingDEG(tnl_til_10s_cellchat, features.name = "TIL")
net.up <- subsetCommunication(tnl_til_10s_cellchat, net = net, datasets = "TIL",ligand.logFC = 0.58492501, receptor.logFC = 0.58492501)
net.down <- subsetCommunication(tnl_til_10s_cellchat, net = net, datasets = "TIL",ligand.logFC = -0.58492501, receptor.logFC = -0.58492501)
gene.up <- extractGeneSubsetFromPair(net.up, tnl_til_10s_cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, tnl_til_10s_cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
netVisual_bubble(tnl_til_10s_cellchat, pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 45, remove.isolate = T)
pairLR.use.down = net.down[, "interaction_name", drop = F]
netVisual_bubble(tnl_til_10s_cellchat, pairLR.use = pairLR.use.down, comparison = c(1, 2),  angle.x = 45, remove.isolate = T)

netVisual_chord_gene(cellchat_tnl_10s,slot.name = 'net', net = net.down,signaling = pathway.union2)
netVisual_chord_gene(cellchat_til_10s,slot.name = 'net', net = net.up,signaling = pathway.union2)
