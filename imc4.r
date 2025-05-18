
##imc analysis FOR HCC

library(CATALYST)
library(SpatialExperiment)
library(SingleCellExperiment)
library(scuttle)
library(scater)
library(imcRtools)
library(cytomapper)
library(dittoSeq)
library(tidyverse)
library(bluster)
library(scran)
library(lisaClust)
library(caret)
###### Processed multiplexed imaging data
##1.read data #### read IMC Segmentation Pipeline generated data
spe <- read_steinbock('./')
counts(spe)[1:5,1:5]
head(colData(spe))
head(rowData(spe))
assay(spe, "exprs") <- asinh(counts(spe)/1)
assay(spe, "scaled") <- t(scale(t(assay(spe, "exprs"))))

## Single-cell processing {#cell-processing}
colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)

read.csv("../../panel1.csv")->panel
spe[panel$name,]->spe
 
 readRDS("../../meta.rds")->meta
spe$PID <- meta$PID[match(spe$sample_id,meta$SID)]
spe$Response <- meta$response[match(spe$sample_id,meta$SID)]
spe$ROI <- as.factor(meta$ROI[match(spe$sample_id,meta$SID)])

library(dittoSeq)
dittoRidgePlot(spe, var = "CD3", group.by = "sample_id", assay = "counts") +
    ggtitle("CD3 - before transformation")
assay(spe, "exprs") <- asinh(counts(spe)/1)
dittoRidgePlot(spe, var = "CD3", group.by = "sample_id", assay = "exprs") +
    ggtitle("CD3 - inverse hyperbolic sine function")
assay(spe, "logcounts") <- log(counts(spe) + 0.01)
dittoRidgePlot(spe, var = "CD3", group.by = "sample_id", assay = "logcounts") +
    ggtitle("CD3 - logarithm")


images <- loadImages("./img/")
masks <- loadImages("./masks/", as.is = TRUE)
read.csv("./panel.csv")$name->channelNames(images)

all.equal(names(images), names(masks))

mcols(images) <- mcols(masks) <-meta

spe$SID=spe$sample_id

df <- data.frame(matrix(ncol = 3, nrow = length(images)))
colnames(df) <- c("sample_id", "width_px", "height_px")

# Fill in the dataframe with object names and feature dimensions
for (i in seq_along(images)) {
  image <- images[[i]]
  df[i, "sample_id"] <- names(images)[i]
  df[i, "width_px"] <- image@dim[1]
  df[i, "height_px"] <- image@dim[2]
}

spe$width_px <- df$width_px[match(spe$sample_id,df$sample_id)]
spe$height_px <- df$height_px[match(spe$sample_id,df$sample_id)]


dittoDimPlot(Tcell, var = "g_type", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3) + 
    ggtitle("Patient ID on UMAP")

dittoDimPlot(spe, var = "Response", reduction.use = "UMAP", size = 0.1, legend.size = 3) + 
    ggtitle("Patient ID on UMAP")
dittoDimPlot(spe, var = "GATA_3", reduction.use = "UMAP", size = 0.1, legend.size = 3) + 
    ggtitle("CD3 on UMAP")
###
rowData(spe)$use_channel <- !grepl("DNA", rownames(spe))

features_oi <- c('CD45','CD45RA','CD3','CD4','CD8','Collagen_I','CD31','CD138','T_bet','CD23','CD68','CD11c','CD20','GATA_3','CD117','CD208','CD45RO','CD57','PNAd','Pan_cytokeratin','aSMA','LYVE_1','PD_L1','LTA')

rowData(spe)$features_oi <- rownames(spe) %in% features_oi
set.seed(12345)
library(BiocParallel)
multicore=MulticoreParam(workers=16)
spe <- runUMAP(spe, subset_row = rowData(spe)$features_oi, exprs_values = "exprs", BPPARAM=multicore) 

library(dplyr)

colData(spe) %>%
    as.data.frame() %>%
    group_by(sample_id) %>%
    summarize(cell_area = sum(area),
           no_pixels = mean(width_px) * mean(height_px)) %>%
    mutate(covered_area = cell_area / no_pixels) %>%
    ggplot() +
        geom_point(aes(reorder(sample_id,covered_area), covered_area)) + 
        theme_minimal(base_size = 15) +
        ylim(c(0, 1)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
        ylab("% covered area") + xlab("")

colData(spe) %>%
    as.data.frame() %>%
    group_by(sample_id) %>%
    ggplot() +
        geom_boxplot(aes(sample_id, area)) +
        theme_minimal(base_size = 15) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
        ylab("Cell area") + xlab("")

summary(spe$area)
sum(spe$area < 5)
spe <- spe[,spe$area >= 5]
saveRDS(spe,'./spe_denoised.rds')
saveRDS(images, "./images.rds")
saveRDS(masks, "./masks.rds")


library(batchelor)
set.seed(12345)
out <- fastMNN(spe, batch = spe$PID,
               auto.merge = TRUE,
               subset.row = rowData(spe)$features_oi,
               assay.type = "exprs",
               BPPARAM = multicore)
reducedDim(spe, "fastMNN") <- reducedDim(out, "corrected")
set.seed(12345)
spe <- runUMAP(spe, dimred= "fastMNN", name = "UMAP_mnnCorrected") 


library(harmony)
library(BiocSingular)
spe <- runPCA(spe, 
              subset_row = rowData(spe)$features_oi, 
              exprs_values = "exprs", 
              ncomponents = 30,
              BSPARAM =  ExactParam())
set.seed(230616)
out <- RunHarmony(spe, group.by.vars = "SID")
reducedDim(spe, "harmony") <- reducedDim(out, "HARMONY")
set.seed(220228)
spe <- runUMAP(spe, dimred = "harmony", name = "UMAP_harmony") 

library(scran)
set.seed(220620)
clusters <- clusterCells(spe, 
                         use.dimred = "harmony", 
                         BLUSPARAM = SNNGraphParam(k = 15, 
                                        cluster.fun = "louvain",
                                        type = "rank"))
spe$h_clusters_corrected <- clusters


library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
mat <- reducedDim(spe, "fastMNN")
set.seed(230619)
out <- Rphenograph(mat, k = 45)
clusters <- factor(membership(out[[2]]))
spe$p_clusters_corrected <- clusters


library(CATALYST)

# Run FlowSOM and ConsensusClusterPlus clustering
set.seed(220410)
som.out <- clusterRows(mat, SomParam(100), full = TRUE)
# Cluster the 100 SOM codes into larger clusters
ccp <- ConsensusClusterPlus(t(som.out$objects$som$codes[[1]]),
                            maxK = 30,
                            reps = 100, 
                            distance = "euclidean", 
                            seed = 220410, 
                            plot = NULL)
CATALYST:::.plot_delta_area(ccp)
som.cluster <- ccp[[12]][["consensusClass"]][som.out$clusters]
spe$som_clusters_corrected <- as.factor(som.cluster)

dittoDimPlot(spe, var = "g_type", reduction.use = "UMAP_mnnCorrected", size = 0.2, legend.show = T,do.label=T) + 
    ggtitle("Patient ID on UMAP")


dittoDimPlot(spe, var = "g_type", reduction.use = "UMAP_harmony", size = 0.2, legend.show = T,do.label=T) + 
    ggtitle("Patient ID on UMAP")
spe <- runTSNE(spe, dimred= "fastMNN", name = "TSNE_mnnCorrected") 

p1<-dittoDimPlot(spe, var = "GATA_3", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3)
p2<-dittoDimPlot(spe, var = "CD68", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3)
p3<-dittoDimPlot(spe, var = "CD3", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3)
p4<-dittoDimPlot(spe, var = "PD_L1", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3)
 p1<-plotReducedDim(spe,dimred='UMAP_mnnCorrected', colour_by = c('CD45RA'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
 p2<-plotReducedDim(spe,dimred='UMAP_mnnCorrected', colour_by = c('LTA'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")

p3<-plotReducedDim(spe,dimred='UMAP_mnnCorrected', colour_by = c('CD68'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
p4<-plotReducedDim(spe,dimred='UMAP_mnnCorrected', colour_by = c('Pan_cytokeratin'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")


write.csv(t(exprs(spe)), './exprs.csv')
write.csv(t(counts(spe)), './counts.csv')
write.csv(t(spe@assays@data$logcounts), './logcounts.csv')
write.csv(rowData(spe), './features_metadata.csv')
write.csv(colData(spe), './observations_metadata.csv')

for (x in names(reducedDims(spe))) {
  write.csv(reducedDim(spe,x), paste0('./',x,'.csv'))
}

assay(spe, "scaled") <- t(scale(t(assay(spe, "exprs"))))

 p1<-plotReducedDim(spe,dimred='UMAP_harmony', colour_by = c('CD20'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
 p2<-plotReducedDim(spe,dimred='UMAP_harmony', colour_by = c('GATA_3'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")

p3<-plotReducedDim(spe,dimred='UMAP_harmony', colour_by = c('CD68'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
p4<-plotReducedDim(spe,dimred='UMAP_harmony', colour_by = c('CD57'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")

plotReducedDim(spe,dimred='UMAP_harmony', colour_by = c('CD57'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
mat <- reducedDim(spe, "fastMNN")
set.seed(12345)
mBIC <- mclust::mclustBIC(mat, G=2:30)
set.seed(12345)
clusters = mclust::Mclust(mat, x=mBIC)



plotReducedDim(Tcell,dimred='UMAP_mnnCorrected', colour_by = c('TIM3'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")

library(pheatmap)
spe$mclust_cluster =clusters$classification
spe$mclust_cluster=factor(spe$mclust_cluster)
dittoDimPlot(spe, var = "subyp", reduction.use = "UMAP_harmony", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")

mean_sce <- aggregateAcrossCells(as(Tcell, "SingleCellExperiment"), ids = DataFrame(sample = Tcell$som_clusters_corrected),average = TRUE)
colnames(mean_sce) <- mean_sce$sample
assay(mean_sce, "arcsinh") <- asinh(assay(mean_sce, "counts")) 
dittoHeatmap(mean_sce[marker1,], assay = "arcsinh",
            annot.by = c('sample'), 
            cluster_cols = F, heatmap.colors = colorRampPalette(c("blue", "white", "red"))(100),
            breaks = seq(-3, 3, length.out = 101))




spe$mclust_label <- factor(spe$mclust_cluster)
levels(spe$mclust_label)=aa$V2
dittoDimPlot(spe, var = "subtype", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")


mean_sce <- aggregateAcrossCells(as(Tcell, "SingleCellExperiment"), ids = DataFrame(sample = Tcell$subtype),average = TRUE)
colnames(mean_sce) <- mean_sce$sample
assay(mean_sce, "arcsinh") <- asinh(assay(mean_sce, "counts")) 
dittoHeatmap(mean_sce[marker1,], assay = "arcsinh",
            annot.by = c('sample'), 
            cluster_cols = F, heatmap.colors = colorRampPalette(c("dark blue", "white", "dark red"))(100),
            breaks = seq(-3, 3, length.out = 101))
dittoHeatmap(Tcell, 
             genes = marker1,
             assay = "exprs", 
             cluster_cols = F, cluster_rows=F,
             scale = "none",
             heatmap.colors = inferno(100), 
             annot.by = c("subtype"))



cur_cells <- sample(seq_len(ncol(spe)), 4000)

# Heatmap visualization - DittoHeatmap
dittoHeatmap(spe_p, 
             genes = marker,
             assay = "exprs", 
             cluster_cols = F, cluster_rows=F,
             scale = "none",
             heatmap.colors = inferno(100), 
             annot.by = c("g_type"))


celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                     ids = spe$g_type, 
                     statistics = "mean",
                     use.assay.type = "exprs")
dittoHeatmap(celltype_mean[marker,],
             assay = "exprs",   annot.by = c("g_type", "ncells"),
             cluster_cols = TRUE,  scale = "none",
             heatmap.colors = viridis(100)
             )

library(ggplot2)
library(cytomapper)
library(SingleCellExperiment)
library(EBImage)
# 假设您已经有以下对象：
# masks: CytoImageList 对象，包含分割掩膜
# spe: SingleCellExperiment 对象，包含细胞元数据
# color: 命名颜色向量，用于不同的细胞亚型
# 选择特定的图像 ID
img_id <- "2119322_1_003"

# 提取对应的掩膜图像
mask_img <- masks[[img_id]]

# 获取掩膜图像的尺寸
dims <- dim(mask_img)

# 将掩膜图像转换为数据框
mask_df <- expand.grid(x = 1:dims[1], y = 1:dims[2])
mask_df$cell_id <- as.vector(t(mask_img))
# 过滤掉背景（cell_id 为 0 的像素）
mask_df <- mask_df[mask_df$cell_id != 0, ]

# 提取与当前图像对应的细胞元数据
spe_df <- as.data.frame(colData(spe))
spe_df <- spe_df[spe_df$SID == img_id, ]

# 合并掩膜数据和细胞元数据
plot_df <- merge(mask_df, spe_df, by.x = "cell_id", by.y = "ObjectNumber")

# 创建 ggplot 图像
ggplot(plot_df, aes(x = x, y = y, fill = subtype)) +
  geom_tile() +
  scale_fill_manual(values = color) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "right", panel.background = element_rect(fill = "white"))













cur_images <- images['1913912_1_003']
cur_images <- normalize(cur_images, separateImages = TRUE)
cur_images <- normalize(cur_images, inputRange = c(0, 0.2))
plotPixels(cur_images,
          mask = masks[c('1913912_1_003')],
           img_id = "SID",
           missing_colour = "white",
           colour_by = c("DNA1",'Pan_cytokeratin'),
           colour = list(DNA1 = c("black", "blue"),'Pan_cytokeratin'=c('black','white')),
           image_title = NULL,
           legend = list(colour_by.title.cex = 0.7,
                         colour_by.labels.cex = 0.7))

cur_images <- images['2004778_1_005']
plotPixels(cur_images, 
           #mask = masks[c('2004778_1_005')],
           img_id = "SID",
           colour_by = c('DNA1',"Pan_cytokeratin",'CD20','CD3','CD68','LTA'),
           bcg = list(DNA1=c(0,5,1),
           Pan_cytokeratin = c(0, 10, 2),
           CD20=c(0,5,1),CD3=c(0,5,1),
           CD68=c(0,5,1),
           LTA=c(0,5,1)),
            colour = list(
                         DNA1=c('black','blue'),Pan_cytokeratin= c("black", "#E21A1C"),CD20=c('black','#EB008B'),CD3=c("black",'#6CCDE3'),CD68=c("black",'#01A14B'),LTA=c("black",'#FFF100')
                         ))


pbMDS(spe, fun='sum',
      by = "sample_id", 
      color_by = "Response",dims=c(1,2),
      features =rownames(spe),label=NULL)+
      scale_color_manual(values=c('#6dbf5a','#0097e9'))+
      ylim(-1,1)+
      xlim(-1.5,1.5)

##############T cell cluster
library(batchelor)
set.seed(12345)
out <- fastMNN(Tcell, batch = Tcell$PID,
               auto.merge = TRUE,
               subset.row = rowData(Tcell)$use_channel,
               assay.type = "exprs",
               BPPARAM = multicore)
 reducedDim(out, "corrected")->corr

set.seed(12345)
mBIC <- mclust::mclustBIC(corr, G=2:10)
set.seed(12345)
clusters = mclust::Mclust(corr, x=mBIC)

spe_Tcell$T_cluster=clusters$classification
spe_Tcell$T_cluster=factor(spe_Tcell$T_cluster)
dittoDimPlot(spe_Tcell, var = "T_cluster", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")
dittoDimPlot(spe_Tcell, var = "Response", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")
mean_sce <- aggregateAcrossCells(as(spe_Tcell, "SingleCellExperiment"), ids = DataFrame(sample = spe_Tcell$T_cluster),average = TRUE)
colnames(mean_sce) <- mean_sce$sample
assay(mean_sce, "arcsinh") <- asinh(assay(mean_sce, "counts")) 
dittoHeatmap(mean_sce[rowData(spe_Tcell)$use_channel,mean_sce$sample!='1'], assay = "arcsinh",
            annot.by = c('sample'), 
            cluster_cols = TRUE, heatmap.colors = colorRampPalette(c("dark blue", "white", "dark red"))(100),
            breaks = seq(-3, 3, length.out = 101))
spe_Tcell$T_label=factor(spe_Tcell$T_cluster)
levels(spe_Tcell$T_label)

dittoDimPlot(spe, var = "mclust_label", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")

spe[,spe$mclust_label %in% Mac]->spe_Mac

out <- fastMNN(spe_Mac, batch = spe_Mac$PID,
               auto.merge = TRUE,
               subset.row = rowData(spe_Mac)$use_channel,
               assay.type = "exprs",
               BPPARAM = multicore)
 reducedDim(out, "corrected")->corr
set.seed(12345)
mBIC <- mclust::mclustBIC(corr, G=2:5)
set.seed(12345)
clusters = mclust::Mclust(corr, x=mBIC)
spe_Mac$Mac_cluster=clusters$classification
spe_Mac$Mac_cluster=factor(spe_Mac$Mac_cluster)
dittoDimPlot(spe_Mac, var = "Mac_cluster", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")

mean_sce <- aggregateAcrossCells(as(spe_Mac, "SingleCellExperiment"), ids = DataFrame(sample = spe_Mac$Mac_cluster),average = TRUE)
colnames(mean_sce) <- mean_sce$sample
assay(mean_sce, "arcsinh") <- asinh(assay(mean_sce, "counts")) 
dittoHeatmap(mean_sce[rowData(spe_Mac)$use_channel,], assay = "arcsinh",
            annot.by = c('sample'), 
            cluster_cols = TRUE, heatmap.colors = colorRampPalette(c("dark blue", "white", "dark red"))(100),
            breaks = seq(-3, 3, length.out = 101))


ggplot(umap,aes(x=umap_1,y=umap_2,color=type))+
geom_point(size=0.05)+theme_classic()+
scale_color_manual(values=co)


dittoDimPlot(spe, var = "Response", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3, do.label = F)+
scale_color_manual(values=c("R"='#6ab559','NR'='#218bcd')) 

cur_masks <- masks[names(masks) %in% '1913912_1_002']
plotCells(cur_masks,
          object = spe, 
          cell_id = "ObjectNumber", 
          img_id = "SID",
          colour_by = "SID",colour=list(mclust_label=co))



###########T cell reclustering



 colnames(spe)[spe$clust_label %in% c("CD4T",'CD8T')]->Tcell
spe[,Tcell]->Tcell


out <- fastMNN(spe_Mac, batch = spe_Mac$PID,
               auto.merge = TRUE,
               subset.row = rowData(spe_Mac)$use_channel,
               assay.type = "exprs",
               BPPARAM = multicore)



library(harmony)
library(BiocSingular)
Tcell <- runPCA(Tcell, 
              subset_row = rowData(Tcell)$features_t, 
              exprs_values = "exprs", 
              ncomponents = 30,
              BSPARAM =  ExactParam())
set.seed(230616)
out <- RunHarmony(Tcell, group.by.vars = "SID")
reducedDim(Tcell, "harmony") <- reducedDim(out, "HARMONY")
set.seed(220228)
Tcell <- runUMAP(Tcell, dimred = "harmony", name = "UMAP_harmony") 

library(scran)
set.seed(220620)
clusters <- clusterCells(Tcell, 
                         use.dimred = "harmony", 
                         BLUSPARAM = SNNGraphParam(k =10, 
                                        cluster.fun = "louvain",
                                        type = "rank"))
Tcell$T_cluster <- clusters
dittoDimPlot(Tcell, var = "T_cluster", reduction.use = "UMAP_harmony", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")
dittoDimPlot(Tcell, var = "clust_label1", reduction.use = "UMAP_harmony", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")
mean_sce <- aggregateAcrossCells(as(Tcell, "SingleCellExperiment"), ids = DataFrame(sample = Tcell$T_cluster),average = TRUE)
colnames(mean_sce) <- mean_sce$sample
assay(mean_sce, "arcsinh") <- asinh(assay(mean_sce, "counts")) 
dittoHeatmap(mean_sce[rowData(Tcell)$features_t,], assay = "arcsinh",
            annot.by = c('sample'), 
            cluster_cols = TRUE, heatmap.colors = colorRampPalette(c("dark blue", "white", "dark red"))(100),
            breaks = seq(-3, 3, length.out = 101))

plotReducedDim(Tcell,dimred='UMAP_harmony', colour_by = c('BCA1'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")

dittoDimPlot(spe, var = "subtype", reduction.use = "UMAP_harmony", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")

spe[,Myeloid]->Myeloid

Myeloid <- runPCA(Myeloid, 
              subset_row = rowData(Myeloid)$features_t, 
              exprs_values = "exprs", 
              ncomponents = 30,
              BSPARAM =  ExactParam())
set.seed(230616)
out <- RunHarmony(Myeloid, group.by.vars = "SID")
reducedDim(Myeloid, "harmony") <- reducedDim(out, "HARMONY")
set.seed(220228)
Myeloid <- runUMAP(Myeloid, dimred = "harmony", name = "UMAP_harmony") 

clusters <- clusterCells(Myeloid, 
                         use.dimred = "harmony", 
                         BLUSPARAM = SNNGraphParam(k =10, 
                                        cluster.fun = "louvain",
                                        type = "rank"))
Myeloid$M_cluster <- clusters
dittoDimPlot(Myeloid, var = "M_cluster", reduction.use = "UMAP_harmony", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")
dittoDimPlot(Myeloid, var = "clust_label", reduction.use = "UMAP_harmony", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")

mean_sce <- aggregateAcrossCells(as(Myeloid, "SingleCellExperiment"), ids = DataFrame(sample = Myeloid$M_cluster),average = TRUE)
colnames(mean_sce) <- mean_sce$sample
assay(mean_sce, "arcsinh") <- asinh(assay(mean_sce, "counts")) 
dittoHeatmap(mean_sce[rowData(Myeloid)$features_m,], assay = "arcsinh",
            annot.by = c('sample'), 
            cluster_cols = TRUE, heatmap.colors = colorRampPalette(c("dark blue", "white", "dark red"))(100),
            breaks = seq(-3, 3, length.out = 101))
p1<-plotReducedDim(Myeloid,dimred='UMAP_harmony', colour_by = c('CD68'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
p2<-plotReducedDim(Myeloid,dimred='UMAP_harmony', colour_by = c('LTA'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
p3<-plotReducedDim(Myeloid,dimred='UMAP_harmony', colour_by = c('LYVE_1'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
p4<-plotReducedDim(Myeloid,dimred='UMAP_harmony', colour_by = c('GATA_3'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
p5<-plotReducedDim(Myeloid,dimred='UMAP_harmony', colour_by = c('CD3'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")
p6<-plotReducedDim(Myeloid,dimred='UMAP_harmony', colour_by = c('IL_4'), by_exprs_values = "scaled", point_size = 0.1, point_alpha = 1)  + 
 scale_colour_gradient2( low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-3, 3), na.value = "#b2182b")


mean_sce <- aggregateAcrossCells(as(spe, "SingleCellExperiment"), ids = DataFrame(sample = spe$g_type),average = TRUE)
colnames(mean_sce) <- mean_sce$sample
assay(mean_sce, "arcsinh") <- asinh(assay(mean_sce, "counts")) 
dittoHeatmap(mean_sce[a,], assay = "arcsinh",
            annot.by = c('sample'), 
            cluster_cols = TRUE, heatmap.colors = colorRampPalette(c("dark blue", "white", "dark red"))(100),
            breaks = seq(-3, 3, length.out = 101))

color=c('B cell'='#355BA7','PD1+TCF1+ Tstem'='#EB2533','#5EB35B',
        'CD4+Tbet+ Th1-like'='#09113E','CD23+ DC'='#E14794','CD4+PD1+CXCL13+ Tfh'='#0C4F29',
        'CD4+CD57+ T'='#FCD31B','CD4+ Treg'='#3B8FC9','CD8+GZMB+ T'='#984C42',
        'Mast Cell'='#56B195','NK'='#665296','Macrophage'='#019092',
        'DC-LAMP'='#D9A5C2','CAF'='#9ABF75','high endothelial venules'='#EA2865',
        'Tumor'='#EC8C52','unknown'='#D8DADC','CD8+PD1+TIM3+ Tex'='#4292C6','Other CD8+ T'='#8DA0CB','Other CD4+ Th'='#1E90FF',
        'Endothelial'='#98FB98','ILC2'='#8A2BE2','Normal Epi'='#666666','CD8+CD57+ T'='#FDB462')






 color=c("Tumor"='#EC8C52',"Macrophage"='#019092',"Endothelial"='#98FB98',"CD8+CD69+ T"='#8DA0CB',"B cell"='#355BA7',"CD4+CD45RO+ T"='#09113E',"Normal Epi"='#666666',"CD8+ Tex"='#4292C6',"other CD8+ T"='#FDB462',
  "CAF"='#9ABF75', "CD4+CCR7+ T"='#FCD31B',"HEVs"='#EA2865',"CD4+HLA-DR+ T"='#1E90FF',"ILC2"='#8A2BE2',"GZMB+ T"='#984C42',"CD8+CXCL13+Ki67+ T"='#D8DADC',"CD4+ Treg"='#3B8FC9',"NK"='#665296',"Mast cell"='#56B195',"CD8+ Tstem"='#EB2533',
  "CD4+ Tfh"='#0C4F29',"LAMP3+ DC"='#D9A5C2')





spe <- runPCA(spe, 
              subset_row = rowData(spe)$features_oi1, 
              exprs_values = "exprs", 
              ncomponents = 50,
              BSPARAM =  ExactParam())
set.seed(230616)
out <- RunHarmony(spe, group.by.vars = "SID")
reducedDim(spe, "harmony1") <- reducedDim(out, "HARMONY")
set.seed(220228)
spe <- runUMAP(spe, dimred = "harmony1", name = "UMAP_harmony1") 
dittoDimPlot(spe, var = "subtype", reduction.use = "UMAP_harmony1", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")


library(batchelor)
set.seed(220228)
out <- fastMNN(spe, batch = spe$SID,
               auto.merge = TRUE,
               subset.row = rowData(spe)$features_oi,
               assay.type = "exprs")
reducedDim(spe, "fastMNN") <- reducedDim(out, "corrected")
library(scater)

set.seed(220228)
spe <- runUMAP(spe, dimred= "fastMNN", name = "UMAP_mnnCorrected") 
dittoDimPlot(spe, var = "subtype", reduction.use = "UMAP_mnnCorrected", size = 0.1, legend.size = 3, do.label = TRUE) + 
    ggtitle("Patient ID on UMAP")

ggplot(UMAP,aes(x=UMAP_1,y=UMAP_2))+
geom_point(aes(color=type),size=2)+
scale_color_manual(values=color)+
theme_classic()



library(SpatialExperiment)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "knn", k = 20)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "expansion", threshold = 20)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "delaunay", max_dist = 20)
colPairNames(spe)

library(ggplot2)
library(viridis)
plotSpatial(spe[,spe$sample_id == "2004778_1_004"], 
            node_color_by = "subtype1", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "neighborhood", 
            nodes_first = FALSE, 
            edge_color_fix = "grey")+
scale_color_manual(values=color)


set.seed(230621)
spe <- detectCommunity(spe, 
                       colPairName = "neighborhood", 
                       size_threshold = 10)

plotSpatial(spe[,spe$sample_id == "2004778_1_001"], 
            node_color_by = "spatial_community", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
    theme(legend.position = "none") +
    ggtitle("Spatial tumor communities") +
    scale_color_manual(values = rev(colors()))


spe <- aggregateNeighbors(spe, 
                          colPairName = "knn_interaction_graph", 
                          aggregate_by = "metadata", 
                          count_by = "subtype")

spe <- aggregateNeighbors(spe, 
                          colPairName = "delaunay_interaction_graph", 
                          aggregate_by = "metadata", 
                          count_by = "subtype")
set.seed(220705)
cn_1 <- kmeans(spe$aggregatedNeighbors, centers =8)
spe$cn_celltypes <- as.factor(cn_1$cluster)


plotSpatial(spe, 
            node_color_by = "tmp", 
            img_id = "sample_id", 
            node_size_fix = 0.1) +
            scale_color_manual(values=c("#EBEBEB",'#4292C6','#EB2533','#984C42','#8A2BE2'))
    ggsci::scale_color_d3('category20')

plotSpatial(spe, 
            node_color_by = "subtype", 
            img_id = "sample_id", 
            node_size_fix = 0.01)+
        scale_color_manual(values=color)+
  theme(
    panel.background = element_rect(fill = "black"),  # 设置面板背景为黑色
    plot.background = element_rect(fill = "black")    # 可选：设置整个绘图区域背景为黑色
  )

for_plot <- prop.table(table(as.character(spe$cn_celltypes), spe$subtype), margin = 1)
pheatmap(for_plot, 
         cluster_rows=F,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100), 
         scale = "column")


# By expression
spe <- aggregateNeighbors(spe, 
                          colPairName = "knn_interaction_graph", 
                          aggregate_by = "expression", 
                          assay_type = "exprs",
                          subset_row = rowData(spe)$use_channel)


set.seed(220705)

cn_2 <- kmeans(spe$mean_aggregatedExpression, centers = 10)
spe$cn_expression <- as.factor(cn_2$cluster)

plotSpatial(spe, 
            node_color_by = "cn_expression", 
            img_id = "sample_id", 
            node_size_fix = 0.1) +
    ggsci::scale_color_d3('category20')

for_plot <- prop.table(table(as.character(spe$cn_expression), spe$subtype), margin = 1)
pheatmap(for_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")

spe <- buildSpatialGraph(spe, 
                         img_id = "sample_id", 
                         type = "knn", 
                         name = "knn_spatialcontext_graph", 
                         k = 40)

# Compute the fraction of cellular neighborhoods around each cell
spe <- aggregateNeighbors(spe, 
                          colPairName = "knn_spatialcontext_graph",
                          aggregate_by = "metadata",
                          count_by = "cn_celltypes",
                          name = "aggregatedNeighborhood")

# Detect spatial contexts
spe <- detectSpatialContext(spe, 
                            entry = "aggregatedNeighborhood",
                            threshold = 0.90,
                            name = "spatial_context")

library(patchwork)

spe <- filterSpatialContext(spe, 
                            entry = "spatial_context",
                            group_by = "PID", 
                            group_threshold = 3,cells_threshold = 100,
                            name = "spatial_context_filtered")

spe <- patchDetection(spe, 
                      patch_cells = spe$subtype == "Tumor",
                      img_id = "sample_id",
                      expand_by = 1,
                      min_patch_size = 10,
                      colPairName = "neighborhood",
                      BPPARAM = MulticoreParam())
spe <- minDistToCells(spe, 
                      x_cells = !is.na(spe$patch_id), 
                      img_id = "sample_id")

plotSpatial(spe, 
            node_color_by = "patch_id", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
    theme(legend.position = "none") 

plotSpatial(spe, 
            node_color_by = "distToCells", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
    scale_color_gradient2(low = "dark blue", mid = "white", high = "dark red")


spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "delaunay", max_dist = 50)

out <- testInteractions(spe, 
                        group_by = "sample_id",
                        label = "subtype", 
                        colPairName = "neighborhood",method='class',
                        iter = 500,
                        BPPARAM = SerialParam(RNGseed = 221029))
cell_type_num <- 22
interaction_num <- cell_type_num * cell_type_num
num_roi <- out@nrows / interaction_num

for (nn in 1:num_roi) {
    start_ind <- (nn - 1) * interaction_num + 1
    end_ind <- nn * interaction_num
    cur_sigvals <- out$sigval[start_ind:end_ind]
    cell_exists <- rep(0, cell_type_num)
    for (mm in 1:cell_type_num) {
        cur_cell_vals <- cur_sigvals[(mm-1)*cell_type_num+1:mm*cell_type_num]
        if (all(is.na(cur_cell_vals)))
            cell_exists[mm] <- 1
    }
    cur_sigvals <- replace(cur_sigvals, is.na(cur_sigvals), -1)
    missing_inds <- which(cell_exists == 1)
    if (length(missing_inds) > 0) {
        ind_combs <- crossing(missing_inds, missing_inds)
        for (jj in 1:dim(ind_combs)[1]) {
            ele_ind <- (ind_combs[[jj,1]]-1)* cell_type_num + ind_combs[[jj,2]]
            cur_sigvals[ele_ind] <- 0
        }
    }
    out$sigval[start_ind:end_ind] <- cur_sigvals
}


readRDS("../../meta.rds")->meta
subset_roi_R=meta$SID[meta$response=='R']
subset_roi_NR=meta$SID[meta$response=='NR']
R_subset_out <- out[out$group_by %in% subset_roi_R, ]
NR_subset_out <- out[out$group_by %in% subset_roi_NR, ]


from_order <-unique(dat$from_label)[c(1:19,22)]
to_order <- from_order[20:1]

max_per_val <- 1
min_per_val <- -0.7

R_subset <- R_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(subset_roi_R)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
R_subset$to_label <- paste0(R_subset$to_label, "-R")

NR_subset <- NR_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(subset_roi_NR)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
NR_subset$to_label <- paste0(NR_subset$to_label, "-NR")


merge_subset <- rbind(R_subset,NR_subset)
merge_to_order <- c()


merge_to_order <- c()
for (cell_type in to_order){
    for (R_type in c("NR", "R")){
        merge_to_order <- append(merge_to_order, paste(cell_type, R_type, sep="-"))
    }
}

order_merge_set <- merge_subset %>% as_tibble() %>% 
    mutate(from_label=factor(from_label, levels=from_order)) %>%
    mutate(to_label=factor(to_label, levels=merge_to_order))
order_merge_set[!is.na(order_merge_set$from_label),]->order_merge_set
order_merge_set[!is.na(order_merge_set$to_label),]->order_merge_set

order_merge_set %>%
    ggplot() +
    geom_tile(aes(from_label, to_label, fill=sum_sigval)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) 







library(scales)
out <- testInteractions(spe, 
                        group_by = "sample_id",
                        label = "subtype", 
                        colPairName = "neighborhood",method='class',
                        iter = 200,
                        BPPARAM = SerialParam(RNGseed = 221029))

readRDS("../../meta.rds")->meta
subset_roi_R=meta$SID[meta$response=='R']
subset_roi_NR=meta$SID[meta$response=='NR']
R_subset_out <- out[out$group_by %in% subset_roi_R, ]
NR_subset_out <- out[out$group_by %in% subset_roi_NR, ]



from_order <-unique(out$from_label)[c(1:19,22)]
to_order <- from_order[20:1]

max_per_val <- 1.00
min_per_val <- -0.8

# update R
R_subset <- R_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(subset_roi_R)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
R_subset$to_label <- paste0(R_subset$to_label, "-R")

NR_subset <- NR_subset_out %>% as_tibble() %>% group_by(from_label, to_label) %>%
    summarize(sum_sigval = sum(sigval, na.rm = TRUE) / length(subset_roi_NR)) %>%
    mutate(across(starts_with("sum"), ~case_when(.x >= 0 ~ .x / max_per_val, TRUE ~ - .x / min_per_val)))
NR_subset$to_label <- paste0(NR_subset$to_label, "-NR")


merge_subset <- rbind(R_subset,NR_subset)
merge_to_order <- c()


merge_to_order <- c()
for (cell_type in to_order){
    for (R_type in c("NR", "R")){
        merge_to_order <- append(merge_to_order, paste(cell_type, R_type, sep="-"))
    }
}

order_merge_set <- merge_subset %>% as_tibble() %>% 
    mutate(from_label=factor(from_label, levels=from_order)) %>%
    mutate(to_label=factor(to_label, levels=merge_to_order))
order_merge_set[!is.na(order_merge_set$from_label),]->order_merge_set
order_merge_set[!is.na(order_merge_set$to_label),]->order_merge_set

order_merge_set %>%
    ggplot() +
    geom_tile(aes(from_label, to_label, fill=sum_sigval)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) 



 plotSpatial(spe[, spe$sample_id == '2119322_1_003'],
            node_color_by = "subtype",
            img_id = "sample_id",
            node_size_fix = 1,
            draw_edges = TRUE,
            colPairName = "expansion_interaction_graph",
            nodes_first = FALSE,
            directed = FALSE,
            edge_color_fix = "grey")+
scale_color_manual(values=color)

heat_color=c("#71C8C3",'#9BD4C7','#C0E5DB','#DAEEEA','#F3FAF8','#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF','#FBE9D8','#F7CFAE','#F4B17E','#F39B4F','#ED822B')




cur_masks=masks[c('2001616_1_001','2200377_1_001')]
Endo <- spe[,spe$subtype == "Endothelial"]
Endo$IL33=assays(Endo)[['scaled']]['IL_33',]
Endo$IL33 <- pmin(Endo$IL33, 3)
plotCells(cur_masks,
          object = Endo, 
          cell_id = "ObjectNumber", thick = TRUE,missing_colour = "black", 
          img_id = "SID",
          colour_by = "IL33",
          exprs_values = "scaled",colour = list('IL33'=c('blue','white','red')))


as.data.frame(rowData(R4))->var_R4
as.data.frame(colData(R4))->obs_R4
as.data.frame(spatialCoords(R4))->spatial_R4
 write.csv(spatial_R4,'spatial_R4.csv')
write.csv(obs_R4,file = "obs_R4.csv")
write.csv(var_R4,file = "var_R4.csv")
write.csv(t(assays(R4)[['exprs']]),'expr_R4.csv')
