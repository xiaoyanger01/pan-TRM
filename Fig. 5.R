library(CellChat)
library(patchwork)
library(Seurat)
library(NMF)
library(ggalluvial)
library(openxlsx)
library(tidyr)
library(ggrepel)



### cellchat
options(stringsAsFactors = FALSE)
args=commandArgs(T)
rds=readRDS(args[1])


DefaultAssay(rds)<-"RNA"
rds1<- CreateSeuratObject(counts=rds@assays$RNA,object=args[2])
rds1@meta.data$type<-rds@meta.data$cell
cellchat <- createCellChat(object = rds1, meta = rds1@meta.data, group.by = "type")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)

future::plan("multiprocess", workers = 16)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

data=cellchat@net$weight[2:11,c(1,12:17)][,c(6,1,2,3,4,5,7)]
p=pheatmap(data,cluster_rows = F,border_color = 'white' ,cluster_cols = F,color = colorRampPalette(c("#fff5f0", "#fee0d2","#fc9272","#cb181d"))(50))


groupSize <- table(cellchat@idents) %>% as.numeric()
weight_mat <- cellchat@net$weight
cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
cir_mat=cir_mat[c(9,1,12:17),c(9,1,12:17)]
cir_mat['C7', ] <- weight_mat["C7", c(9,1,12:17)]
p=netVisual_circle( cir_mat, vertex.weight= groupSize[c(9,1,12:17)], weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.4)



## ST
st_traint <- CellTrek::traint(st_data=st_obj, 
                              sc_data=sc_obj, 
                              sc_assay='RNA', 
                              cell_names='cell_copykat')


st_celltrek <- CellTrek::celltrek(st_sc_int=brain_tst_traintraint, int_assay='traint', sc_data=sc_obj, sc_assay = 'RNA', 
                                  reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                  dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek


st_celltrek$cell <- factor(st_celltrek$cell, 
                           levels=sort(unique(st_celltrek$cell)))


CellTrek::celltrek_vis(st_celltrek@meta.data %>% 
                         dplyr::select(coord_x, coord_y, cell_type:id_new),
                       brain_celltrek@images$ST.colon2@image, 
                       brain_celltrek@images$ST.colon2@scale.factors$lowres)


## distance
cell_distance=function(brain_celltrek,ref,que,new_name){
  inp_df <- brain_celltrek@meta.data %>% dplyr::select(cell_names = dplyr::one_of('cell_copykat'), coord_x, coord_y)
  inp_df$coord_x = 270-inp_df$coord_x
  output <- kdist(inp_df = inp_df, 
                  ref = ref, #目标细胞类型
                  ref_type = 'all', 
                  que = que,  #其余的细胞
                  k = 5, 
                  new_name = new_name,
                  keep_nn = F)
  
  res = output$kdist_df
  return(res)
  
}



## density

p3=cell_density(obj,starmap_meta,img_temp,'C7')
p4=cell_density(obj,starmap_meta,img_temp,'T')