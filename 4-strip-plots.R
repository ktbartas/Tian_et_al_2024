setwd('/media/emem/Lab/kbartas/gpe_cocsal')
library(Seurat) 
library(magrittr)
library(ggplot2)
library(data.table)
library(dplyr)
combined.sct<-readRDS('R_outs/neuron-outs/gpe_combined-rna-annotated.rds')
DefaultAssay(combined.sct) <- "RNA"
combined.sct<- NormalizeData(object = combined.sct)
Idents(combined.sct) <- "treatment"
Idents(combined.sct) <- "cluster.broad"
print('set idents')
clusts<-unique(Idents(combined.sct))
degs.markers<-read.csv('R_outs/neuron-outs/DGE/ALL-cells_degs.csv')  
degs.markers$cluster <- 'All nuclei'
print('got all cells');
for (clustname in clusts) {
  plotname1 = paste0(clustname,'_degs.csv')
  plotname2 = paste0('R_outs/neuron-outs/DGE/',plotname1)
  temp2<-read.csv(plotname2)  
  if (dim(temp2)[1]>0) {temp2$cluster <- clustname}
  degs.markers<-bind_rows(temp2,degs.markers)  
  label2<-c('Kcnq5','Kcnq3')
  plot_as <- ggplot(temp2, aes(x=avg_log2FC, y=-log(p_val_adj), col=diffexpressed,label=X)) + geom_point() +geom_text(aes(label=ifelse(X %in% label2,as.character(X),'')),hjust=0.5,vjust=0.5)
  plotname1 = paste0(clustname,'_degs.svg')
  plotname2 = paste0('neuron-figs/Kcnqs/',plotname1)
  ggsave(file=plotname2, plot=plot_as)
}
print('got marker and deg list') 
# make strip plot..

require(scales)
library(ggrepel)
# Create vector with levels of object@ident
cluster <- unique(degs.markers$cluster)
# Create vector of default ggplot2 colors

###### P-val corrected ###########
genesinterest<-c('Kcnq3','Kcnq5')#,'Nrxn3','Kcnq1','Kcnq2','Kcnq4','Drd1','Cry2')'Nxph1',
#                 'Drd2','Nr4a2','Bdnf','Jun','Per1','Per2','Per3','Cry1','Egr1','Cdnk1a', 'Nfkbia' , 'Fosb')
cols <- hue_pal()(length(cluster))
col_dat<-as_tibble(cbind(cols, cluster))
col_dat
cell_DE_res_merge<-full_join(degs.markers, col_dat, by = "cluster") 
cell_DE_res_merge<-mutate(cell_DE_res_merge, FC=  ifelse(avg_log2FC <= -.5 | avg_log2FC >= .5 , TRUE, FALSE))
cell_DE_res_merge<-mutate(cell_DE_res_merge, pval_yes=  ifelse(p_val_adj < 0.05 , TRUE, FALSE))
# cell_DE_res_merge<-mutate(cell_DE_res_merge, fillyn=  ifelse(X %in% genesinterest, TRUE, FALSE))
cell_DE_res_merge<-mutate(cell_DE_res_merge, fillyn=  ifelse(pval_yes & FC , TRUE, FALSE))
cell_DE_res_merge
cell_DE_res_merge<-mutate(cell_DE_res_merge, col_use = ifelse(fillyn== TRUE, cols, "dark gray"))
cell_DE_res_merge<-mutate(cell_DE_res_merge, updown = ifelse(p_val_adj < 0.05 & avg_log2FC >=0.5, "Upregulated in cocaine-treated mice", ifelse(p_val_adj < 0.05 & avg_log2FC <= -0.5, "Downregulated in cocaine-treated mice", "NS")))
de_res_table<-table(cell_DE_res_merge$updown, cell_DE_res_merge$cluster)
de_res_table<-as.data.frame.matrix(de_res_table)
de_res_table
col_for_plot<-as.character(cell_DE_res_merge$col_use)
table(col_for_plot)
x<-unique(cell_DE_res_merge$cluster)
cell_DE_res_merge$cluster<- factor(cell_DE_res_merge$cluster, levels = x)
table(cell_DE_res_merge$fillyn,cell_DE_res_merge$cluster )
strip<-ggplot(cell_DE_res_merge, aes(x = cluster, y = avg_log2FC)) +
   geom_jitter(position = position_jitter(0.3), color = col_for_plot, alpha = 0.5)+
    theme_bw() + ylim(-2,2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text_repel(aes(label=ifelse(X %in% genesinterest,paste("(",X,",",p_val_adj,")"),'')),max.overlaps = 3000)
ggsave(file="neuron-figs/DGE/stripplot-kcnqs.svg", plot=strip,height=6,width=6);
print('complete - p-val corrected')
