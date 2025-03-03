library(Seurat)
library(SingleCellExperiment)
library(patchwork)


#consensus abundance
tcga_obj_sub = readRDS('./tcga_obj_sub.rds')
tcga_obj_sub = tcga_obj_sub %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()  %>% RunUMAP( return.model = TRUE, dims = 1:50)
tcga_obj_sub <- tcga_obj_sub %>% FindNeighbors() %>% FindClusters(resolution = 1)
p1=DimPlot(tcga_obj_sub,group.by = 'cancer',reduction = "tsne", pt.size = .1,cols=c(
  "#F79159","#C176B0","#EED785","#D4B595","#EEBED6","#EDB869","#CD7560","#7FC28E","#F48F93",
  "#F9CABE","#C76DA8","#76A1EC","#41ae76","#B6D890","#a6cee3","pink","#33a02c","#F8A895",
  "#b2df8a","#e31a1c","#fb9a99","#fdbf6f","#ff7f00","#6a3d9a","#cab2d6","#ffff99","#F79159","#8C9FD1","#8887B4","#8C9FD1","#BAA1CC","#E6674C","#E87A8C")) 

tcga_obj_sub$MAST = "Heterogeneous"
tcga_obj_sub$MAST[tcga_obj_sub$MAST_xcell > quantile(tcga_obj_sub$MAST_xcell, 0.3333) & tcga_obj_sub$MAST_ssgsea > quantile(tcga_obj_sub$MAST_ssgsea, 0.3333) & tcga_obj_sub$NEU_quantiseq > quantile(tcga_obj_sub$mast_NAKAJIMA_MAST_CELL, 0.3333)] = "High"
tcga_obj_sub$MAST[tcga_obj_sub$MAST_xcell < quantile(tcga_obj_sub$MAST_xcell, 0.6666) & tcga_obj_sub$MAST_ssgsea < quantile(tcga_obj_sub$MAST_ssgsea, 0.6666) & tcga_obj_sub$NEU_quantiseq < quantile(tcga_obj_sub$mast_NAKAJIMA_MAST_CELL, 0.6666)] = "Sparse"
p2=DimPlot(tcga_obj_sub,group.by = 'MAST',pt.size = .1, reduction = "tsne",cols = c("grey80","grey60","grey30"))
p=p1+p2



### scRNA-seq Processing

library(Seurat)
library(dplyr)
library(harmony)

bl_list=read.delim("./bl.list")

sc_obj=CreateSeuratObject(counts = data@assays$RNA@counts, min.cells = 10, min.features = 200)

sc_obj=subset(sc_obj,nCount_RNA >= 600)
sc_obj=subset(sc_obj,nCount_RNA <= 120000)
sc_obj=subset(sc_obj,nFeature_RNA >= 400)
sc_obj=subset(sc_obj,nFeature_RNA <= 8000)
sc_obj=subset(sc_obj,percent.mt <= 10)


sc_obj <- NormalizeData(sc_obj) %>% FindVariableFeatures(nfeatures=2000 )
VariableFeatures(sc_obj)=setdiff(VariableFeatures(sc_obj),bl[,1])
sc_obj <- sc_obj %>% ScaleData() %>% RunPCA(verbose=FALSE)
sc_obj <- RunHarmony(sc_obj, group.by.vars = "class")
sc_obj <-  sc_obj %>% RunUMAP(reduction = "harmony", dims = 1:30)
sc_obj <- sc_obj %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.6)



### sample info
sc_obj = readRDS('./sc_obj.rds')
sc_obj_meta=sc_obj@meta.data

d1=data.frame(table(sc_obj_meta$cancer_type))
d1=d1[order(d1$Freq,decreasing = T),]
d1$Var1=factor(d1$Var1,levels = d1$Var1)
p1=ggplot(d1, aes(x = Var1, y = Freq))+
  geom_bar(stat = 'identity',width = 0.8,fill='#a6cee3')+scale_y_log10()+theme_classic()+
  theme( axis.text.x=element_text(angle=90,hjust = 0,colour="black",size=10),legend.position = 'right')+
  scale_fill_manual(values = c("#F08F26"))+ylab("cell")
d2=unique(data.frame(sc_obj_meta$cancer_type,sc_obj_meta$Sample_id))
d2=data.frame(table(d2$cancer_type))
d2$Var1=factor(d2$Var1,levels = d1$Var1)
p2=ggplot(d2, aes(x = Var1, y = Freq))+
  geom_bar(stat = 'identity',width = 0.8,fill='#a1d99b')+scale_y_log10()+theme_classic()+
  theme( axis.text.x=element_text(angle=90,hjust = 0,colour="black",size=10),legend.position = 'right')+
  scale_fill_manual(values = c("#F08F26"))+ylab("Sample")

d3=unique(data.frame(sc_obj_meta$cancer_type,sc_obj_meta$Patient))
d3=data.frame(table(d3$cancer_type))
d3$Var1=factor(d3$Var1,levels = d1$Var1)
p3=ggplot(d3, aes(x = Var1, y = Freq))+
  geom_bar(stat = 'identity',width = 0.8,fill='#fcbba1')+scale_y_log10()+theme_classic()+
  theme( axis.text.x=element_text(angle=90,hjust = 0,colour="black",size=10),legend.position = 'right')+
  scale_fill_manual(values = c("#F08F26"))+ylab("Patient")

d4=data.frame(da$cancer_type,da$site)
num = d4  %>% group_by(da.cancer_type,da.site)
num <- dplyr::summarise(num,count = n())
num=data.frame(num)
num$da.cancer_type=factor(num$da.cancer_type,levels = d1$Var1)
num$da.site=factor(num$da.site,levels = rev(c("Tumor","Metastasis","Normal","Blood","lymph node")))

p4=ggplot(num,aes(cancer_type,count,fill=da.site))+
  geom_bar(stat="identity",position="fill",width = 0.8)+
  scale_fill_manual(values = rev(c("#F6CE39","#8D9AC6","#A2CA5E","#F68865","#6baed6")
  ))+theme_classic()+ylab("Cell fraction")+
  theme( axis.text.x=element_text(angle=90,hjust = 0,colour="black",size=10),legend.position = 'right')

p=p1+p2+p3+p4+plot_layout(ncol=1)



### marker
c0=c("RPL31","RPL27A","RPL21")
c1=c("HSPA1A","HSPA1B","DNAJB1")
c2=c("NFKB1","NR4A3","VEGFA" )
c3=c("CD69", "RGS1","CMA1" ) 
c4=c("PDIA3","ENO1","MYDGF"  ) 
c5=c("CREM","FTH1","ANXA1" )
c6=c("MS4A2","LTC4S","TXNIP")
c7=c( "CD74","HLA-DRA","HLA-DRB1")
c8=c("UBB","FOSB" ,"DNAJA1") 
c9=c( "CCL5","IL32","GZMB") 

p=DotPlot(sc,features = c(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9),group.by = 'cell')+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+
  scale_colour_gradient2(low = "white", mid = "#c6dbef", high = "#084081")
