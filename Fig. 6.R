library(monocle3)

sc_obj_sub=subset(sc_obj,site == 'Tumor')
data <- GetAssayData(sc_obj_sub, assay = 'RNA', slot = 'counts')
cell_metadata <- sc_obj_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)

cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sc_obj_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="cell") 
cds <- cluster_cells(cds)

cds <- learn_graph(cds)


cds <- order_cells(cds)

p=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, alpha=0.5,cell_size=0.4,graph_label_size = 0,
             label_leaves = FALSE,  label_branch_points = FALSE) +
  scale_color_gradient2(low="#3288bd", mid="#ffffbf", high="#f46d43")


cds_sub <- choose_graph_segments(cds,clear_cds=FALSE, return_list = T)
cds_sub <- cds[,cds_sub$cells]
lin1=pseudotime(cds_sub, reduction_method = "UMAP")

cds_sub2 <- choose_graph_segments(cds,clear_cds=FALSE, return_list = T)
cds_sub2 <- cds[,cds_sub2$cells]
lin2=pseudotime(cds_sub2, reduction_method = "UMAP")


p1=ggplot(data,aes(Pseudotime,Antigen.processing.and.presentation,color=lin ))+stat_smooth(method = "gam",formula = y ~ x + I(x ^ 2),se=T,span=1,fullrange =T)+
  scale_color_manual(values = c("#fad99d","#a98abe"))+theme_classic()
p2=ggplot(data,aes(Pseudotime,Stress,color=lin ))+stat_smooth(method = "gam",formula = y ~ x + I(x ^ 2),se=T,span=1,fullrange =T)+
  scale_color_manual(values = c("#fad99d","#a98abe"))+theme_classic()

###MEBOCOST
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from mebocost import mebocost
adata = sc.read_h5ad('./sc_obj.h5ad')
mebo_obj = mebocost.create_obj(
  adata = adata,
  group_col = ['cell'],
  met_est = 'mebocost',
  config_path = './mebocost.conf',
  exp_mat=None,
  cell_ann=None,
  species='human',
  met_pred=None,
  met_enzyme=None,
  met_sensor=None,
  met_ann=None,
  scFEA_ann=None,
  compass_met_ann=None,
  compass_rxn_ann=None,
  cutoff_exp='auto',
  cutoff_met='auto',
  cutoff_prop=0.25,
  sensor_type=['Receptor', 'Transporter', 'Nuclear Receptor'],
  thread=8
)

mebo_obj._load_config_()
mebo_obj.estimator()


met_mat = pd.DataFrame(mebo_obj.met_mat.toarray(),
                       index = mebo_obj.met_mat_indexer,
                       columns = mebo_obj.met_mat_columns)

met_mat.head()

p=ggdotchart(plot_da, x = "variable", y = "value",
             color = "class",                      
             palette = c('#fb8072',"#1f78b4",'#b2df8a','#ff7f00','#bc80bd','#fdbf6f','#a6cee3','#fb9a99','#cab2d6','#b3de69','#bebada','#fb8072','#80b1d3','#fdb462'),
             add = "segments",  
             sorting ='none',
             add.params = list(color = "class", size = 0.4), 
             dot.size ="value",                                 
             ggtheme = theme_pubr())+ylab("Relative metabolite presence")+xlab("")+
  scale_size_continuous(range = c(0.5,4))

### boxplot

p1=ggplot(data,aes(Dihydrolanosterol,Antigen.processing.and.presentation))+geom_point(color='#6baed6',size=1,alpha=0.5)+
  geom_smooth(method = 'lm', se = F, color = 'black')+
  stat_cor(data=da, method = "spearman")+theme_bw()+xlim(-0.2,0.2)
p2=ggplot(data,aes(Dihydrolanosterol,Antigen.processing.and.presentation,color=cell))+geom_point(size=0.5)+
  theme_bw()+scale_color_manual(values = c('#fb8072',"#1f78b4",'#b2df8a','#ff7f00','#bc80bd','#fdbf6f','#a6cee3','#fb9a99','#cab2d6','#b3de69','#bebada','#fb8072','#80b1d3','#fdb462'))+
  theme(legend.position = 'none')