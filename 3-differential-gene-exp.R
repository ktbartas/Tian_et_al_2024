setwd('/media/emem/Lab/kbartas/gpe_cocsal')
library(Seurat) 
library(magrittr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggrepel)
combined.sct<-readRDS('R_outs/neuron-outs/gpe_combined-rna-annotated.rds')
DefaultAssay(combined.sct) <- "RNA"
print('data loaded 1')
#define negate func
`%!in%` <- Negate(`%in%`)
print('defined func 2')
##DEGs for each cluster 
#"R_outs/diff_exp/DEGs_each_cluster.csv" ?
Idents(combined.sct) <- "cluster.specific"
clusts<-unique(Idents(combined.sct))
combined.sct<- NormalizeData(object = combined.sct)
combined.sct<- ScaleData(object = combined.sct)
Idents(combined.sct) <- "cluster.broad"

levels(combined.sct) <- rev(c("Oligodendrocytes","OPCs", "Microglia","Astrocytes","GABAergic"))
colors = rev(c("#F8766D",         "#DE8C00","#F564E3",  "#C77CFF", "#00BA38"))
plot2 <- VlnPlot(
  object = combined.sct,
  features = c("Kcnq1", "Kcnq2","Kcnq3","Kcnq4","Kcnq5"),
  pt.size = 1,cols=colors
 ) 
ggsave(file='neuron-figs/Kcnqs/by-types-Kcnqs-vln.svg', plot=plot2,height=7,width=15); 
plot2 <- VlnPlot(
  object = combined.sct,
  features = c("Kcnq1", "Kcnq2","Kcnq3","Kcnq4","Kcnq5"),
  pt.size = 0,cols=colors
 ) 
ggsave(file='neuron-figs/Kcnqs/by-types-Kcnqs-vln-nopts.svg', plot=plot2,height=7,width=15)

g1<-VlnPlot(combined.sct,cols=colors, features = c("nFeature_RNA", "nCount_RNA","percent.mt"))
ggsave(file="neuron-figs/QC/combined_QC_cells.svg", plot=g1,width = 8, height = 7)
print('made init violin plots')


combined.sct<-readRDS('R_outs/neuron-outs/gpe_combined-rna-annotated.rds')
DefaultAssay(combined.sct) <- "RNA"
Idents(combined.sct) <- "cluster.specific"
clusts<-unique(Idents(combined.sct))
combined.sct<- NormalizeData(object = combined.sct);
combined.sct@meta.data$celltype.stim <- paste0(combined.sct@meta.data$cluster.specific, "_", 
    combined.sct@meta.data$treatment)
Idents(combined.sct) <- "celltype.stim"
print('ready for for loop')
for (clustname in clusts) {
  # Volcano plot  
  ident1<-paste0(clustname,'_Control')
  ident2<-paste0(clustname,'_Cocaine')
  as_markers<- FindMarkers(combined.sct, ident.1 = ident2, ident.2 = ident1, only.pos = FALSE, test.use = "MAST",min.pct = 0.01,logfc.threshold = 0.1)
  as_markers$diffexpressed <- "NO"
  as_markers$diffexpressed[as_markers$avg_log2FC > .10 & as_markers$p_val_adj< 0.05] <- "UP"
  as_markers$diffexpressed[as_markers$avg_log2FC < -.10 & as_markers$p_val_adj< 0.05]  <- "DOWN"
  as_markers$X<-row.names(as_markers)
  df<-as_markers[as_markers$diffexpressed =='NO', ] 
  label2<-df$X
  plot_as <- ggplot(as_markers, aes(x=avg_log2FC, y=-log(p_val_adj), col=diffexpressed,label=X)) + geom_point() +geom_text_repel(aes(label=ifelse(X %!in% label2,as.character(X),'')),hjust=0.5,vjust=0.5)
  plotname1 = paste0(clustname,'_degs.svg')
  plotname2 = paste0('neuron-figs/DGE/',plotname1)
  ggsave(file=plotname2, plot=plot_as)
  plotname1 = paste0(clustname,'_degs.csv')
  plotname2 = paste0('R_outs/neuron-outs/DGE/',plotname1)
  write.csv(as_markers,plotname2) 
};
print('check - DGE individual clusters done');

Idents(combined.sct) <- "treatment"; 
as_markers<- FindMarkers(combined.sct, ident.1 = 'Cocaine', ident.2 = 'Control', only.pos = FALSE, test.use = "MAST",min.pct = 0.01,logfc.threshold = 0.1);
print('found markers for all cells')
as_markers$diffexpressed <- "NO";
as_markers$diffexpressed[as_markers$avg_log2FC > .10 & as_markers$p_val_adj< 0.05] <- "UP";
as_markers$diffexpressed[as_markers$avg_log2FC < -.10 & as_markers$p_val_adj< 0.05]  <- "DOWN";
as_markers$X<-row.names(as_markers)
df<-as_markers[as_markers$diffexpressed =='NO', ] 
label2<-df$X
plot_as <- ggplot(as_markers, aes(x=avg_log2FC, y=-log(p_val_adj), col=diffexpressed,label=X)) + geom_point() +geom_text_repel(aes(label=ifelse(X %!in% label2,as.character(X),'')))
ggsave(file='neuron-figs/DGE/ALL-cells_degs.svg', plot=plot_as);
write.csv(as_markers,'R_outs/neuron-outs/DGE/ALL-cells_degs.csv');
print('complete all  cells - now run neurons');

#subset to just GABAergic neurons
Idents(combined.sct) <- "cluster.broad"
combined.sct<-subset(x =combined.sct, idents = c("GABAergic"))
Idents(combined.sct) <- "treatment"; 
combined.sct<- NormalizeData(object = combined.sct);
as_markers<- FindMarkers(combined.sct,  ident.1 = 'Cocaine', ident.2 = 'Control', only.pos = FALSE, test.use = "MAST",min.pct = 0.01,logfc.threshold = 0.1);
as_markers$diffexpressed <- "NO";
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
as_markers$diffexpressed[as_markers$avg_log2FC > .10 & as_markers$p_val_adj< 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
as_markers$diffexpressed[as_markers$avg_log2FC < -.10 & as_markers$p_val_adj< 0.05]  <- "DOWN"
# label= row.names(as_markers[1:10,]
as_markers$X<-row.names(as_markers)
df<-as_markers[as_markers$diffexpressed =='NO', ] 
#label2<-df$X
write.csv(as_markers,'R_outs/neuron-outs/DGE/GABAergic_degs.csv');

###########make plots for neuron
as_markers<-read.csv('R_outs/neuron-outs/DGE/GABAergic_degs.csv');
label2<-list('Celf2',
'Phactr1',
'Rarb',
'Rbfox1',
'Cacnb2',
'Gm10754',
'Cacna2d3',
'Kcnd2',
'Dgkb',
'Dlgap2',
'Kcnq3',
'Kcnq5',
'Pde10a',
'Rgs9',
'Itpr1',
'Auts2',
'Nrg1',
'Pdzd2',
'Dlg2',
'Slit3',
'Mctp1',
'Foxp1',
'Ank3',
'Ano3',
'Adcy5',
'Neto1',
'Robo2',
'Cdh8', 
'Fbxl7')
as_markers$plotcolor <- "Not significant";
as_markers$plotcolor[as_markers$p_val_adj< 0.05] <- "Significant, |log2FC|<0.1"
as_markers$plotcolor[as_markers$avg_log2FC > .10 & as_markers$p_val_adj< 0.05] <- "Significant, |log2FC|>0.1"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
as_markers$plotcolor[as_markers$avg_log2FC < -.10 & as_markers$p_val_adj< 0.05]  <- "Significant, |log2FC|>0.1"
plot_as <- ggplot(as_markers, aes(x=avg_log2FC, y=-log(p_val_adj), col=plotcolor,label=X)) + theme_minimal() + geom_point() +geom_text_repel(aes(label=ifelse(X %in% label2,as.character(X),'')),max.overlaps =200)+ geom_vline(xintercept=-0.1, linetype="dashed", color = "red")+ geom_vline(xintercept=0.1, linetype="dashed", color = "red") + geom_hline(yintercept=-log(0.05), linetype="dashed", color = "red")+ xlim(-0.8, 0.8)+ ylim(-1, 35)
ggsave(file='neuron-figs/DGE/GABAergic_degs.svg', plot=plot_as,height=7,width=8)

label2<-list("Cacna2d3",
"Rbfox1",
"Dgkb",
"Celf2",
"Rarb",
"Nrg1",
"Phactr1",
"Lrrc7",
"Ptprn",
"Rgs9",
"Dlgap2",
"Slit3",
"Sgcz",
"Kcnd2",
"Robo2",
"Syt1",
"Ano3",
"Kcnip4",
"Cacnb2",
"Fgf14",
"Adgrb3",
"Kcnq5",
"Gabrg3",
"Erc2",
"Mctp1",
"Sema3c",
"Me3",
"Arglu1",
"Kcnma1",
"Dcc",
"Grm5",
"Homer1",
"Dlg2",
"Gabrb3",
"Plekha5",
"Camkv",
"Lrfn5",
"Ryr3",
"Cit",
"Pcsk2",
"Dach1",
"Nebl",
"Unc13c",
"Asic2",
"Lrrtm3",
"Cdh7",
"Zfp385b",
"Epha6",
"Sv2c",
"Elmo1",
"Akap6",
"Tenm4",
"Caln1",
"Frmpd4",
"Phlpp1",
"Mast4",
"Grm1",
"Chn1",
"Thrb",
"Rogdi",
"Mirg",
"Gabra4",
"Ptprq",
"Rnf112",
"Rprd1a",
"Matk",
"Dab1",
"C2cd5",
"Ptpn5",
"Epha7",
"Kcnab1",
"Lrrk2")
as_markers<-read.csv('R_outs/neuron-outs/DGE/ALL-cells_degs.csv');
as_markers$plotcolor <- "Not significant";
as_markers$plotcolor[as_markers$p_val_adj< 0.05] <- "Significant, |log2FC|<0.1"
as_markers$plotcolor[as_markers$avg_log2FC > .10 & as_markers$p_val_adj< 0.05] <- "Significant, |log2FC|>0.1"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
as_markers$plotcolor[as_markers$avg_log2FC < -.10 & as_markers$p_val_adj< 0.05]  <- "Significant, |log2FC|>0.1"
plot_as <- ggplot(as_markers, aes(x=avg_log2FC, y=-log(p_val_adj), col=plotcolor,label=X)) + theme_minimal() + geom_point() +geom_text_repel(aes(label=ifelse(X %in% label2,as.character(X),'')))+ geom_vline(xintercept=-0.1, linetype="dashed", color = "red")+ geom_vline(xintercept=0.1, linetype="dashed", color = "red") + geom_hline(yintercept=-log(0.05), linetype="dashed", color = "red")+ xlim(-0.8, 0.8)+ ylim(-1, 35)
ggsave(file='neuron-figs/DGE/all-cells_degs-clearbg.svg', plot=plot_as,height=7,width=8)
print('complete (DGE)');
#

#subset to just GABAergic neurons
Idents(combined.sct) <- "cluster.specific"
combined.sct<-subset(x =combined.sct, idents = c("GABAergic 1","GABAergic 2"))
Idents(combined.sct) <- "treatment"; 
combined.sct<- NormalizeData(object = combined.sct);
as_markers<- FindMarkers(combined.sct,  ident.1 = 'Cocaine', ident.2 = 'Control', only.pos = FALSE, test.use = "MAST",min.pct = 0.01,logfc.threshold = 0.1);
as_markers$diffexpressed <- "NO";
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
as_markers$diffexpressed[as_markers$avg_log2FC > .10 & as_markers$p_val_adj< 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
as_markers$diffexpressed[as_markers$avg_log2FC < -.10 & as_markers$p_val_adj< 0.05]  <- "DOWN"
# label= row.names(as_markers[1:10,]
as_markers$X<-row.names(as_markers)
df<-as_markers[as_markers$diffexpressed =='NO', ] 
#label2<-df$X
write.csv(as_markers,'R_outs/neuron-outs/DGE/GABAergic12_degs.csv');

###########make plots for neuron
as_markers<-read.csv('R_outs/neuron-outs/DGE/GABAergic12_degs.csv');
label2<-list('Kcnq3',
'Kcnq5')
as_markers$plotcolor <- "Not significant";
as_markers$plotcolor[as_markers$p_val_adj< 0.05] <- "Significant, |log2FC|<0.1"
as_markers$plotcolor[as_markers$avg_log2FC > .10 & as_markers$p_val_adj< 0.05] <- "Significant, |log2FC|>0.1"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
as_markers$plotcolor[as_markers$avg_log2FC < -.10 & as_markers$p_val_adj< 0.05]  <- "Significant, |log2FC|>0.1"
plot_as <- ggplot(as_markers, aes(x=avg_log2FC, y=-log(p_val_adj), col=plotcolor,label=X)) + theme_minimal() + geom_point() +geom_text_repel(aes(label=ifelse(X %in% label2,as.character(X),'')),max.overlaps =2000)+ geom_vline(xintercept=-0.1, linetype="dashed", color = "red")+ geom_vline(xintercept=0.1, linetype="dashed", color = "red") + geom_hline(yintercept=-log(0.05), linetype="dashed", color = "red")+ xlim(-0.8, 0.8)+ ylim(-1, 35)
ggsave(file='neuron-figs/DGE/GABAergic12_degs.svg', plot=plot_as);
print('GABA 1 2 done')
