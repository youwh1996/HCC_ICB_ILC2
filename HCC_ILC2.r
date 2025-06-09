HCC3T<-SCTransform(HCC3T,assay='Spatial')
HCC3T <- RunPCA(HCC3T, assay = "SCT", verbose = FALSE)
HCC3T <- FindNeighbors(HCC3T, reduction = "pca", dims = 1:30)
HCC3T <- FindClusters(HCC3T, resolution = 1.5)
HCC3T <- RunUMAP(HCC3T, reduction = "pca", dims = 1:30)

sc_meta=data.frame(cellID=colnames(dat),cellType=dat$type1,sampleInfo=dat$PID,row.names=colnames(dat))

sc_meta=data.frame(cellID=colnames(dat),cellType=dat$sub_cluster,sampleInfo=dat$patient,row.names=colnames(dat))
sc_meta=data.frame(cellID=colnames(dat),cellType=dat$DefineTypes,sampleInfo=dat$Pat_Tissues,row.names=colnames(dat))

CARD_obj = createCARDObject(
	sc_count = dat@assays$RNA@counts,
	sc_meta = sc_meta,
	spatial_count = HCC3T@assays$Spatial@counts[,rownames(cord)],
	spatial_location = cord,
	ct.varname = "cellType",
	ct.select = unique(sc_meta$cellType),
	sample.varname = "sampleInfo",
	minCountGene = 100,
	minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
CARD_obj@Proportion_CARD->prop
saveRDS(CARD_obj,'P1T_spot_celltype.rds')

library(SpaTopic)
HCC3T<-HCC3T[,rownames(CARD_obj@Proportion_CARD)]
result_list <- CellTopic(
  CARD_obj@Proportion_CARD,
  HCC3T@meta.data["seurat_clusters"],
  cluster = "seurat_clusters",
  num_topics =12,
  percent = 0.7
)
my_colors <- c("#66b803", "#E5AA9B", "#FABD3E", "#2B8CBE", "#DE77AE", "#9970AB", "gray", "#D5E4A2", "#71D0F5", "#B1746F", "#ADE2D0", "#20DE8BFF", "#CCDE8BFF", "#FFDE8BFF", "#FFA88BFF", "#FF6A8BFF")
HCC3T <- AddMetaData(HCC3T, result_list[["CellTopic"]])
SpatialDimPlot(HCC3T, group.by = "CellTopic", image.alpha = 0, pt.size.factor = 1.1) + scale_fill_manual(values = my_colors)
CellTopic_plot(HCC3T, 
               celltype_topic = result_list[["celltype_topic"]], 
               celltopic = paste0("CellTopic", c(1:15)),
               cols.highlight = c("#DE2D26", "#FEE0D2"),
               #cols.celltype = my_colors, 
               pt.size.factor = 1.2)
CellTopic_plot(ST1, 
               celltype_topic = result_list[["celltype_topic"]], 
               celltopic = "CellTopic2", 
               cols.highlight = c("#DE2D26", "#FEE0D2"),
               pt.size.factor = 4.2)

result_list[["celltype_topic"]]->celltype_topic
bar_plot_data <- celltype_topic
bar_plot_data$CellType <- rownames(celltype_topic)
bar_plot_data$LA=apply(bar_plot_data[,c('CellTopic3','CellTopic7','CellTopic12','CellTopic16','CellTopic10')],1,mean)
rownames(bar_plot_data)[order(-bar_plot_data[,'LA'])][1:15]->aa
ggplot2::ggplot(bar_plot_data[aa,], aes(x = stats::reorder(.data$CellType,
            .data[['LA']], decreasing = TRUE), y = .data[['LA']])) +
            ggplot2::geom_bar(aes(fill = .data$CellType), stat = "identity",
                width = 0.7) + ggplot2::xlab("CellType") + ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none", axis.text.x = element_text(angle = 90,
                hjust = 1, vjust = 0, size = 8), axis.title.y = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        coord_flip()


AddModuleScore(HCC3T,features=list(sig))->HCC3T
SpatialFeaturePlot(P8T,features='OTUD5')+
scale_fill_gradientn(colours = colorRampPalette(c("#20AFD9","#1B96BA",'#167997','#15667D','#135A6E','#114050','#133C4B','#16343E','#1F262A','#2D2A22','#343024','#414023','#545325','#71722A','#8B8C31','#9B9D32','#B6B632','#BBBD2C','#CFCD28','#E9E320'))(100))

###magic diff
library(ggplot2)
Rcells=rownames(meta)[meta$Res=='R']
NRcells=rownames(meta)[meta$Res=='NR']

ds <- min(c(length(Rcells),length(NRcells)))
g1 <- sample(Rcells,ds)
g3 <- sample(NRcells,ds)
##umis magic
umis=spot@assays$Magic
logFC = as.vector(apply(umis, 1, function(x) log2((mean(x[g1]))/(mean(x[g3])))))
logFCx = as.vector(apply(umis, 1, function(x) log2(mean(x[g1]))))
logFCy = as.vector(apply(umis, 1, function(x) log2(mean(x[g3]))))

wilc = p.adjust(as.vector(apply(umis, 1, function(x) wilcox.test(x[g1], x[g3])$p.value)), method = "fdr")
pb <- data.frame(rownames(umis), logFC, logFCx, logFCy, wilc)
colnames(pb) <- c("gene", "fc", "x", "y", "pval")
pb$pval[pb$pval=="NaN"] <- NA
pb <- na.omit(pb)
pb[pb$x!=(-Inf) & pb$y!=(-Inf),]->pb

pb$x2<-c()
pb$y2<-c()
for(i in 1:nrow(pb)){
    pb$x2[i]=mean( spot@assays$Spatial@data[rownames(pb)[i],g1])
    pb$y2[i]=mean( spot@assays$Spatial@data[rownames(pb)[i],g3])
}


ggplot(data=pb,aes(x=rank,y=log2(fc+1),label=show))+
geom_point(size=0.5,color='#BEBEBE')+
geom_hline(yintercept = 0, linetype = "dashed")+
geom_text_repel(data = pb[pb$show != "",],col='#3551AA', size=3.5,max.overlaps=15)+
geom_point(data=pb[pb$show!='',],aes(x=rank,y=log2(fc+1)),size=3,color='red')+
ylim(-5,5)+
theme_classic()

gl<-list()
for(i in as.character(unique(dat$subtype))){
tmp<-subset(dat,subtype==i)
FindMarkers(tmp,group.by='patient',ident.1='R',min.pct=0.2,logfc.threshold=0.25)->gl[[i]]
}


ggplot(meta,aes(x=Confirmed.Response_IRF,y=IL33))+geom_boxplot(outlier.color=NA)+geom_jitter(aes(color=Confirmed.Response_IRF))+ggpubr::stat_compare_means()+theme_classic()+
scale_color_manual(values=c("#6DBF5A",'#0097E9'))

ggplot(meta,aes(x=group,y=TP))+geom_boxplot(outlier.color=NA)+geom_jitter(aes(color=group))+ggpubr::stat_compare_means()+theme_classic()+
scale_color_manual(values=c("#6DBF5A",'#0097E9'))
#########GSE235863
sceasy::convertFormat("~/data/JH_HCC/GSE235863/GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad",'anndata','seurat')->scRNA

FindAllMarkers(dat,only.pos=T,logfc.threshold=0.25,min.pct=0.25)->deg
library(dplyr)
data.frame(deg %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC))->top50
unique(top50$gene)->gene

intersect(rownames(pb)[pb$pval<0.05 & pb$fc>0.1],gene)->sig

###########
dataimpu <- magic(P9T, genes = c("IL33"),
                  npca = 50, n.jobs = 30, assay = "Spatial")

data <- GetAssayData(dataimpu, slot = "data",
                     assay = "MAGIC_Spatial")
data <- data.frame(t(data), check.names = F)
data <- rownames_to_column(data, "id")
 P9T@images$image@coordinates$row->data$row
P9T@images$image@coordinates$col->data$col
ggplot(data,aes(y=-row,x=col,fill=IL33))+geom_point(size=2.5,shape = 21,color = "black",stroke =0.1)+
scale_fill_gradientn(colours = colorRampPalette(c('#FDFDFD','#FDFDFD','#FDFDFD','#DBECEC','#A6CFE3','#7AA1E2','#7850B3','#622163','#581952','#370F24','#360F24','#380C23'))(100))+theme_classic()

p1<-ggplot(data,aes(x=row,y=col,fill=LTA))+geom_point(size=2.5,shape = 21,color = "black",stroke =0.1)+
scale_fill_gradientn(colours = colorRampPalette(c('#FDFDFD','#FDFDFD','#FDFDFD','#DBECEC','#A6CFE3','#7AA1E2','#7850B3','#622163','#581952','#370F24','#360F24','#380C23'))(100))+theme_classic()
p2<-ggplot(data,aes(x=row,y=col,fill=LTB))+geom_point(size=2.5,shape = 21,color = "black",stroke =0.1)+
scale_fill_gradientn(colours = colorRampPalette(c('#FDFDFD','#FDFDFD','#FDFDFD','#DBECEC','#A6CFE3','#7AA1E2','#7850B3','#622163','#581952','#370F24','#360F24','#380C23'))(100))+theme_classic()
p3<-ggplot(data,aes(x=row,y=col,fill=CD70))+geom_point(size=2.5,shape = 21,color = "black",stroke =0.1)+
scale_fill_gradientn(colours = colorRampPalette(c('#FDFDFD','#FDFDFD','#FDFDFD','#DBECEC','#A6CFE3','#7AA1E2','#7850B3','#622163','#581952','#370F24','#360F24','#380C23'))(100))+theme_classic()
SpatialFeaturePlot(P9T,features='ILC2a1')+
scale_fill_gradientn(colours = colorRampPalette(c('#FDFDFD','#FDFDFD','#FDFDFD','#DBECEC','#A6CFE3','#7AA1E2','#7850B3','#622163','#581952','#370F24','#360F24','#380C23'))(100))+theme_classic()

ifelse(dat$ILC2a1>quantile(dat$ILC2a1,0.9),'pos','neg')->dat$ILC2



ggplot()+
geom_point(data=meta[meta$ILC2=='neg' &meta$Stem=='neg',],aes(x=col,y=-row),shape = 21,color = "black",stroke =0.1,size=2.5)+
geom_point(data=meta[meta$ILC2=='pos',],aes(x=col,y=-row),shape = 22,color = "green",stroke =0.1,size=2.5)+
geom_point(data=meta[meta$Stem=='pos',],aes(x=col,y=-row),shape = 24,color = "red",stroke =0.1,size=2.5)+
theme_classic()


p1<-ggplot(data=HCC3L@meta.data)+geom_point(aes(x=col,y=-row,color=sig1),stroke =0.1,size=2.5)+
scale_color_gradientn(colours = colorRampPalette(c("#20AFD9","#1B96BA",'#167997','#15667D','#135A6E','#114050','#133C4B','#16343E','#1F262A','#2D2A22','#343024','#414023','#545325','#71722A','#8B8C31','#9B9D32','#B6B632','#BBBD2C','#CFCD28','#E9E320'))(100),limits = c(-0.45, 1.15))+
theme_classic()
p2<-ggplot(data=HCC3T@meta.data)+geom_point(aes(x=col,y=-row,color=sig1),stroke =0.1,size=2.5)+
scale_color_gradientn(colours = colorRampPalette(c("#20AFD9","#1B96BA",'#167997','#15667D','#135A6E','#114050','#133C4B','#16343E','#1F262A','#2D2A22','#343024','#414023','#545325','#71722A','#8B8C31','#9B9D32','#B6B632','#BBBD2C','#CFCD28','#E9E320'))(100),limits = c(-0.45, 1.15))+
theme_classic()

meta1=P9T@meta.data
meta1$col= P9T@images$image@coordinates$col
meta1$row=P9T@images$image@coordinates$row
meta1$LA=dat@meta.data[colnames(P9T),'LA1']
p2<-ggplot(data=meta1)+geom_point(aes(x=col,y=-row,color=LA),stroke =0.1,size=2.5)+
scale_color_gradientn(colours = colorRampPalette(c('#FDFDFD','#FDFDFD','#FDFDFD','#DBECEC','#A6CFE3','#7AA1E2','#7850B3','#622163','#581952','#370F24','#360F24','#380C23'))(100),limits = c(-0.15, 0.4))+
theme_classic()

 SpatialDimPlot(dat,group.by='LAS',cols=c('neg'="#FDFDFD",'pos'='#5E4FA2'))

##########CD70+ ILC2 signature
ILCimpu <- magic(ILC, genes = c("LTA",'LTB','CD70'),
                  npca = 50, n.jobs = 30, assay = "RNA")

data <- GetAssayData(ILCimpu, slot = "data",
                     assay = "MAGIC_RNA")
data <- data.frame(t(data), check.names = F)
data <- rownames_to_column(data, "id")

###
P9Timpu <- magic(P9T, 
                  npca = 50, n.jobs = 30, assay = "Spatial")
#读取文件名
exp.file = "TPM.merge.txt"
in.gct.file = "ESTIMATE_input.gct"

#将表达量矩阵转换为gct格式，匹配ssGSEA分析要求的格式
outputGCT(exp.file, in.gct.file)
#基因过滤，将输入的GCT格式的表达谱文件中基因与内置common gene进行比对过滤
filterCommonGenes(input.f = exp.file,output.f = in.gct.file,id = "NAME")
#输出过滤后的gct文件
out.score.file = "ESTIMATE_score.gct"

#ESTIMATE分析
estimateScore(in.gct.file,out.score.file,platform = "illumina")

#platform是数据平台，illumina、agilent平台不输出TumorPurity；affymetrix会输出TumorPurity
#绘制样本的肿瘤细胞纯度以及其与公共数据库中样本比较的散点图
#plotPurity(scores=out.score.file)

library(estimate)
estimate<-function(dat,pro){
input.f=paste0(pro,'_estimate_input.txt')
output.f=paste0(pro,'_estimate_gene.gct')
output.ds=paste0(pro,'_estimate_score.gct')
write.table(dat,file=input.f,sep='\t',quote=F)
library(estimate)
filterCommonGenes(input.f=input.f,
output.f=output.f,
id="GeneSymbol")
estimateScore(input.ds=output.f,
output.ds=output.ds,
platform="illumina")## 注意platform
scores=read.table(output.ds,skip=2,header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
return(scores)
}
estimate(tpm,'HCC')->scores

TumorPurity=cos(0.6049872018+0.0001467884*scores[,3])
head(TumorPurity)




library(survival)
library(survminer)
fit<-survfit(Surv(OS.time,OS)~group,data=meta)

geneList<-data.frame(Gene=diff$Gene.Name,FC=as.numeric(diff$WT.vs..KO.Log2.Fold.Change))
glist <- geneList[,2]
names(glist) <- as.character(geneList[,1])
 glist <- sort(glist,decreasing = T)
gsea<-GSEA(glist, TERM2GENE=all, verbose=T, pvalueCutoff = 1)##


fgseaRes<-fgsea(pathways=gl,
stats=rankData,
minSize=5,
maxSize=500, nperm=1000)




group=factor(c(rep("KO",4),rep("WT",4)))
dge<-DGEList(counts=exp,group=group)

design.matrix<-model.matrix(~0+group)
colnames(design.matrix)<-levels(group)
dge<-estimateDisp(dge,design=design.matrix)
fit<-glmQLFit(y=dge,design=design.matrix)
con<-makeContrasts(WT-KO,levels=design.matrix)
qlf<-glmQLFTest(fit,contrast=con)
dge_results<-as.data.frame(topTags(qlf,n=nrow(dge$counts),sort.by="logFC"))
rankData <- dge_results$logFC
names(rankData) <- rownames(dge_results)
head(rankData)

data.frame(fgseaRes)->res

plotEnrichment(gl[["GOBP_POSITIVE_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_POTENTIAL"]],
               rankData)

c("GOMF_HISTONE_H3_METHYL_LYSINE_4_DEMETHYLASE_ACTIVITY",'GOMF_HISTONE_METHYLTRANSFERASE_ACTIVITY_H3_K4_SPECIFIC','GOBP_HISTONE_H3_K4_MONOMETHYLATION','GOBP_CHROMATIN_REMODELING_AT_CENTROMERE','GOBP_NUCLEOLAR_CHROMATIN_ORGANIZATION','REACTOME_METHYLATION')->name
res[res$pathway %in% name,]->tmp
library(ggplot2)
tmp$name=c('Chromatin remodeling','Histone H3 K4 Monomethlytation','Nucleolar Chromatin organization','Histone Methyltransferase  Activity H3 K4 Specific','Methylation')
ggplot(tmp,aes(x=name,y=NES))+geom_col(fill='#FFA040')+coord_flip()+theme_bw()

c('Lta','Ltb','Cd70','mt-Co1','mt-Co2','mt-Co3','Atp6ap2','Atp6ap1l','Ndufa1','Ndufa10','Ndufa2','Ndufa3','Ndufa4','Ndufb2','Uqcr10','Uqcr11','Uqcrb','Cox6a1','Cox7a2','Cox6c','Atp5d')


pheatmap(aa[,c(2:4,5,6,8)],scale='row',border=NA,cluster_cols=F,cluster_rows=F,color=colorRampPalette(c('white','white','white','white','#F5A5A7','#EF5757','#FF0000','#FF0000'))(1024))

gsea<-GSEA(glist, TERM2GENE=data.frame(term='HALLMARK_OXIDATIVE_PHOSPHORYLATION',genes=gl[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']]), verbose=T, pvalueCutoff = 1)##
enrichplot::gseaplot2(gsea,1)

mainPathways=c('GOBP_POSITIVE_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_POTENTIAL',"WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",'WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY','GOCC_MITOCHONDRIAL_TRICARBOXYLIC_ACID_CYCLE_ENZYME_COMPLEX','GOBP_GLUTAMINE_FAMILY_AMINO_ACID_BIOSYNTHETIC_PROCESS')
topPathways <- fgseaRes[fgseaRes$pathway %in% mainPathways,]
plotGseaTable(gl[topPathways$pathway], rankData, topPathways, 
              gseaParam=0.5)
            

c('ILC1a'='#BDD3EC','ILC1b'='#EB8FA9','ILC1c'='#FBB83B','ILC1d'='#5C9DCA','ILC2a'='#1C5BA5','ILC2b'='#736BB0','ILC2c'='#6FBD47','ILC2d'='#212154','ILC3a'='#CB2028','ILC3b'='#F08122','Cycling ILC'='#2B877C')

ggplot(emb,aes(x=UMAP_1,y=UMAP_2))+
geom_point(aes(color=type))+
scale_color_manual(values=c('ILC1a'='#BDD3EC','ILC1b'='#EB8FA9','ILC1c'='#FBB83B','ILC1d'='#5C9DCA','ILC2a'='#1C5BA5','ILC2b'='#736BB0','ILC2c'='#6FBD47','ILC2d'='#212154','ILC3a'='#CB2028','ILC3b'='#F08122','Cycling ILC'='#2B877C'))+
theme_classic()

ggscplot(object=ILC)+
stat_density2d(geom="raster",aes(fill=..density..),
contour=F,show.legend=F)+
geom_scPoint(color="white",size=0.2)+
facet_wrap(~group,ncol=3)+
theme_bw()+
theme(panel.grid=element_blank(),
axis.ticks=element_blank(),
strip.background=element_blank(),
strip.text=element_blank(),
axis.text=element_blank())+
scale_fill_viridis_c(option="magma",direction=1)+coord_cartesian(expand=F)


library(ggrepel)

ggplot()+
geom_point(data=diff,aes(x=logFC,y=log10P),size=0.5,color='#BEBEBE')+
geom_point(data=diff[diff$gene!='',],aes(x=logFC,y=log10P,color=gene),size=2)+
scale_color_manual(values=c('#810F7C','#F781BF','#FDB462','#FF0000','#0000FF','#00FF00','#6FBD47','#8C6BB1','#E0ECF4','#FFD92F'))+
theme_classic()



 P9T@meta.data->meta_P9

rTLS=c("IL7R","CCL19","MALAT1","GSTP1","COTL1","ACTA2","ZFP36L2","TYMP","RGS1","CAPG","HERPUD1","TCF4","IFI27L2","VCAN","CXCL10","IFI16","TCF7")

ILC2_CD70=c("IL7R","GATA3","CD70","LTA","XCL1","IL13","NFKBIZ","PHLDA1","BIRC3","DUSP6","BHLHE40","ICOS","TNFAIP3","SAMSN1","CXCL8","PTGS2","LIF","ZNF331","SLC2A3","FOSL2","CSF1",
"NR4A2","PFKFB3","PTGER4","CSF2","TNFSF10","VIM","CREM","RGS1","PDE4D","FCER1G","SOCS2","CARD16","PRNP","MAP3K8", "CAPG","IRF2BP2","MYC","IL10RA","NAMPT","NDFIP2","XCL2","ZFP36L1",
"GADD45A","H3F3B","CD82","KLRB1","TAB3","CD69","GALNT4")

AddModuleScore(P1T,features=list(ILC2=ILC2_CD70,rTLS=rTLS))->P1T
P9T@meta.data=meta_P9
select1 <- rownames(meta_P1)[meta_P1$Cluster1 >  quantile(meta_P1$Cluster1, probs = c(0.9))]
select2 <- rownames(meta_P1)[meta_P1$region != "NA"]



P1T@images$image@coordinates->spatial_coords
P1T@images$image@spot.radius->radio
ilc2_spots=spatial_coords[select1,c('row','col')]
ilc2_spots$id=select1
 colnames(ilc2_spots)[1:2]=c("x",'y')

TLS_spots=spatial_coords[select2,c('row','col')]
TLS_spots$id=select2
 colnames(TLS_spots)[1:2]=c("x",'y')

calculate_distance_matrix <- function(ilc2_spots, tls_spots) {
  n_ilc2 <- nrow(ilc2_spots)
  n_tls <- nrow(tls_spots)
  
  # 初始化距离矩阵
  distance_matrix <- matrix(0, nrow = n_ilc2, ncol = n_tls)

  for (i in 1:n_ilc2) {
    for (j in 1:n_tls) {
      # 计算欧式距离
      dist <- sqrt((ilc2_spots$x[i] - tls_spots$x[j])^2 + (ilc2_spots$y[i] - tls_spots$y[j])^2)
      distance_matrix[i, j] <- dist
    }
  }
  
  return(distance_matrix)
}

distance_matrix <- calculate_distance_matrix(ilc2_spots, TLS_spots)
rownames(distance_matrix) <- ilc2_spots$id
colnames(distance_matrix) <- TLS_spots$id

input <- do.call(rbind, lapply(colnames(distance_matrix), function(i){
  data.frame(id = i, dist = min(distance_matrix[,i]))
}))
input[order(input$dist),]->input
input$number <- 1:nrow(input)
ggplot(input)+
  geom_smooth(aes(y=TCF7,x=number),span = 1,  method = "loess", se = TRUE,color='red') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1))


tmp.magic <- t(Mac@assays$Spatial@counts)
tmp.magic <- library.size.normalize(tmp.magic)
tmp.magic <- sqrt(tmp.magic)
sce.all.MAGIC <- magic(tmp.magic, genes="all_genes") 
 
Macimmp <- magic(Mac,  npca = 50, n.jobs = 5, assay = "RNA")


conda create -n rna_metagenomics -c bioconda fastp kneaddata bowtie2 samtools
conda activate rna_metagenomics

gunzip RNA_P01_T1_Pre-treatment_R1.fastq.gz
gunzip RNA_P01_T1_Pre-treatment_R2.fastq.gz
kneaddata \
  -i1 RNA_P01_T1_Pre-treatment_R1.fastq -i2 RNA_P01_T1_Pre-treatment_R2.fastq   -o ./kd/ \
  -db  /home/youwh/reference/hg38_bt2/GRCh38_noalt_as/GRCh38_noalt_as -p 10
gzip RNA_P01_T1_Pre-treatment_R1.fastq
gzip RNA_P01_T1_Pre-treatment_R2.fastq

mkdir rRNA_work
sortmerna \
  --ref /home/youwh/Bioapp/silva-arc-16s-id95.fasta \
  --reads ./kd/RNA_P01_T1_Pre-treatment_R1_kneaddata_paired_1.fastq \
  --reads ./kd/RNA_P01_T1_Pre-treatment_R1_kneaddata_paired_2.fastq \
  --workdir ./rRNA_work --paired \
  --aligned rRNA --other non_rRNA --threads 10


