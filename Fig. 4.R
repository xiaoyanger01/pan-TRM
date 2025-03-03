library(dplyr)
library(ggplot2)
library(pheatmap)
library(data.table)
library(circlize)
library(grDevices)
library(survminer)
library(survival)
library(RegParallel)
library(ggpubr)
library(corrplot)
library(patchwork)



### HR and survival curve
geneSets = list(c0=c("TPSAB1","RPL27A","RPS27","RPL41","RPL35"),c1=c("TPSAB1","DNAJB1","HSPH1","JUN","BAG3"),
                c2=c("TPSAB1","NFKB1","NR4A3","VEGFA","KDM6B"),c3=c("TPSAB1","RGS1","CMA1","CTSG","FOS"),
                c4=c("TPSAB1","PKM","MYDGF","PDIA3","PPIB"),c5=c("TPSAB1","CREM","BCL2A1","BIRC3","EIF1"),
                c6=c("TPSAB1","TXNIP","TPSB2","CTSD","S100A4"),c7=c("TPSAB1","HLA-DPB1","HLA-DPA1","HLA-DQB1","HLA-DMB"),
                c8=c("TPSAB1","HSPA1A","EGR1","HSPE1","FOS"),c9=c("TPSAB1","CCL5","IL32","NKG7","IL7R"))

out <- data.frame()
for(i in 1:33){
  print(i)
  expr_data =cancer_fpkm[[i]]
  sparse_expr <- as(as.matrix(expr_data), "sparseMatrix")
  cells_rankings <- AUCell_buildRankings(sparse_expr )
  auc_scores <- AUCell_calcAUC(geneSets, cells_rankings)
  auc_df <- data.frame(id = colnames(expr_data), score = t(auc_scores))
  
  clin_data <- cancer_clinical[[names(cancer_fpkm)[i]]]
  os_time <- ifelse(is.na(clin_data$days_to_death), clin_data$days_to_last_follow_up, clin_data$days_to_death)
  
  clinical_df <- data.frame(
    OS = os_time,
    vital_status = clin_data$vital_status,
    id = clin_data$submitter_id,
    age_at_diagnosis = clin_data$age_at_diagnosis,
    gender = clin_data$gender,
    tumor_stage = clin_data$tumor_grade,
    stringsAsFactors = FALSE
  )
  merged_data <- left_join(auc_df, clinical_df, by = "id")
  merged_data$cancer <- names(cancer_fpkm)[i]
  
  out <- rbind(out, merged_data)
}

do.cox=function(data,cancer,cluster){
  pdat = out[out$cancer == cancer,]
  res.cut <- surv_cutpoint(pdat, time = "OS", 
                           event = "vital_status", 
                           variables = cluster 
  )
  res.cat <- surv_categorize(res.cut)
  pdat$group <- res.cat[,cluster] 
  cox = coxph(Surv(OS,vital_status) ~ group, data=pdat)
  cox.summary = summary(cox)
  nvar = length(unique(pdat$group)) - 1
  HR = round(cox.summary$conf.int[1:nvar,1], 2)
  HR.lower = round(cox.summary$conf.int[1:nvar,3], 2)
  HR.upper = round(cox.summary$conf.int[1:nvar,4], 2)
  HR.range = sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef = cox.summary$coefficients[1:nvar,1]
  coef.se = cox.summary$coefficients[1:nvar,3]
  Pval = round(cox.summary$coefficients[1:nvar,5], 4)
  group = gsub("group","",rownames(cox.summary$conf.int)[1:nvar])
  res=data.frame(cancerType = cancer,cluster=cluster, group=group, HR=HR, HR.range=HR.range, coef=coef, coef.se=coef.se, Pval=Pval)
  return(res)
}

p1=pheatmap(data,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#2171b5","#c6dbef","white", "#fee0d2","#fb6a4a"))(50),border_color = NA)

p2=ggplot(data, aes(HR,cluster,color=HR)) +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(xmin = HR_low, xmax = HR_high), width = 0.15) +
  geom_vline(xintercept = 1, linetype = 3) +
  scale_colour_gradientn(colours = color,limits=c(0,2))+
  theme_classic() +xlab("Harzard Ratio (95% CI)")

fit <- survfit(Surv(OS, vital_status) ~score_c1, data = res.cat)
p3=ggsurvplot(fit, pval=TRUE, 
              legend.labs=c( "High","Low"),
              size = 1,
              risk.table = F,
              palette = c("#fc9272","#6baed6"))+xlab("Time(Days) ")

p=p1+p2+p3


### ICB boxplot
p=ggplot(data)+geom_boxplot(aes(class,score1,fill=class,alpha =0.2))+
  geom_jitter(aes(class,score,fill=class,color=class),width =0.2,shape = 21,size=2.5)+
  scale_fill_manual(values = c("#fec44f","#8D9AC6","#F84C56"))+theme_classic()+ylab("c1 signature")+ggtitle("GSE202069 p = 0.8064")+
  scale_color_manual(values = c("#fec44f","#8D9AC6","#F84C56"))

