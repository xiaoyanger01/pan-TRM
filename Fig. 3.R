library(ggplot2)
library(dplyr)
library(viridis) 
library(ggpointdensity) 

### density plot
umap_data <- sc_obj@reductions$umap@cell.embeddings
umap_df <- as.data.frame(umap_data)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$tissue=sc_obj$site
tumor=umap_df %>%filter(tissue == 'Tumor')
normal=umap_df %>%filter(tissue == 'Normal')
blood=umap_df %>%filter(tissue == 'Blood')
Metastasis=umap_df %>%filter(tissue == 'Metastasis')

p1=ggplot(data = tumor, mapping = aes(x = UMAP_1,y = UMAP_2)) + 
  geom_pointdensity(size=0.1) + ggtitle("Tumor")+
  scale_color_viridis_c(option = "F")+theme_classic()

p2=ggplot(data = blood, mapping = aes(x = UMAP_1,y = UMAP_2)) + 
  geom_pointdensity(size=0.1) + ggtitle("Blood")+
  scale_color_viridis_c(option = "F")+theme_classic()

p3=ggplot(data = normal, mapping = aes(x = UMAP_1,y = UMAP_2)) + 
  geom_pointdensity(size=0.1) + ggtitle("Normal")+
  scale_color_viridis_c(option = "F")+theme_classic()

p4=ggplot(data = Metastasis, mapping = aes(x = UMAP_1,y = UMAP_2)) + 
  geom_pointdensity(size=0.1) + ggtitle("Metastasis")+
  scale_color_viridis_c(option = "F")+theme_classic()

p=p1+p2+p3+p4+plot_layout(ncol=2)


###  cell percentage
inflammation=subset(sc_obj,site == 'inflammation')
counts=rbind(count(inflammation$cancer,inflammation$cell),cancer=table(tumor$cell))
data=apply(counts, 1, function(x){x/sum(x)}*100)
data=data.frame(data)

per=melt(data)

p=ggplot(per,aes(value,Var2,fill=Var1))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values = c('#fb8072',"#1f78b4",'#b2df8a','#ff7f00','#bc80bd','#fdbf6f','#a6cee3','#fb9a99','#cab2d6','#b3de69','#bebada','#fb8072','#80b1d3','#fdb462'
  ))+theme_classic()+ylab("Cell fraction")+theme( axis.text.x=element_text(angle=90,hjust = 0,colour="black",size=10),
                                                  legend.position = 'right')


### Ro/e
out.prefix= '/picb/sysgenomics/huangyingying/01_pan_cancer_singlecell/all_combine_mast/00_hg38_gene/'
Ro.e <- do.tissueDist(cellInfo.tb = meta.tb,
                      out.prefix=sprintf("test",out.prefix),
                      pdf.width=4,pdf.height=6,verbose=1)
Ro.e_matrix = Ro.e$OR.dist.mtx
Ro.e_matrix[Ro.e_matrix>2]=2
Ro.e_matrix=data.frame(Ro.e_matrix)

col = colorRampPalette((brewer.pal(n = 3,  name = "OrRd")))

p=pheatmap(Ro.e_matrix,
           show_colnames = T,
           scale = "none",
           cluster_rows = F, 
           cluster_cols =F, 
           color = colorRampPalette(c("white", "#fee0d2","#fb6a4a"))(50),
           display_numbers = TRUE,gaps_col = 14,
           border_color = "grey70",clustering_distance_cols = "correlation")