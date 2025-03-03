library(UCell)
library(tidyr)
library(clusterProfiler)



#### annotation
marker = FindAllMarkers(object = sc_obj)
marker=marker[which(marker$pct.1 > 0.1 & marker$p_val < 0.05 & marker$avg_log2FC > 0.25),]
top_markers  = marker %>%
  group_by(cluster) %>%
  top_n(100, avg_log2FC) %>%
  mutate(row_num = row_number())


FeatureMatrix  = top_markers %>%
  pivot_wider(names_from = cluster, values_from = gene, id_cols = row_num)
FeatureMatrix$row_num = NULL

c5 = read.gmt("kegg_MSigDB.gmt")
ego_comb = character(0)
cluster.names = colnames(FeatureMatrix)
for (i in 1:(length(cluster.names))) {
  # print(cluster.names[i])
  gene_list = FeatureMatrix[,cluster.names[i]]
  ego = enricher(as.matrix(gene_list), TERM2GENE=c5)
  ego_result = ego@result
  ego_result$GeneRatio.num = ego_result$Count/length(as.matrix(gene_list))
  ego_result = ego_result[order(ego_result$GeneRatio.num, decreasing = T),]
  ego_result = subset(ego_result, ego_result$pvalue < 0.05)
  ego_result = cbind(ego_result, cluster.names[i])
  colnames(ego_result)[ncol(ego0_result)] = "Group"
  ego_comb = rbind(ego_comb, ego_result)
}




c0=c("Ribosome","Mitophagy")
c1=c("UV Response Up","Hypoxia")
c2=c("IL-2/STAT5 Signaling","TNF signaling pathway")  # 诱导其他细胞
c3=c("Glutathione metabolism","Renin-angiotensin system")
c4=c("PI3K/AKT/mTOR  Signaling","Fc gamma R-mediated phagocytosis")
c5=c("Reactive Oxygen Species Pathway","Inflammatory Response")
c6=c("Fc epsilon RI signaling pathway","Cholesterol Homeostasis")
c7=c("Antigen processing and presentation","Interferon Alpha Response")
c8=c("Osteoclast differentiation","p53 Pathway")
c9=c("T cell receptor signaling pathway","PD-L1 expression and PD-1 checkpoint pathway in cancer")
key=c(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
ego_comb_sub = ego_comb %>% filter(Description %in% key)
ego_comb_sub$logp = -log10(ego_comb_sub$pvalue)
ego_comb_sub$Description = factor(ego_comb_sub$Description, levels=key)

data=data.frame(id=ego_comb_sub$Description,value=ego_comb_sub$logp,group=ego_comb_sub$group)
data=dcast(data,id~group)
rownames(data)=data$id


p=pheatmap(data,cluster_rows = F,cluster_cols = F,
           color = colorRampPalette(colors = c("white","#ff7f00","#800026"))(100),border_color = 'grey90')



#### published mast cell 
signature.matrix <- (data.frame(ScoreSignatures_UCell(sc_obj@assays$RNA@data, features = signature, name = "")))

p=DotPlot(sub,features = c("JEM.MC1","JEM.MC2","JEM.MC3","JEM.MC4","JEM.MC5","JEM.MC6","MC3","MC4",'Resting.MC1','activated.MC1'),group.by = 'cell',cols = c("#fff5f090","#a50f15") )+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+scale_size_continuous(range = c(2,6))

