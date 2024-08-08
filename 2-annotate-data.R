# conda activate signac 
setwd('/media/emem/Lab/kbartas/gpe_cocsal')
library(Seurat)   
library(magrittr)
library(ggplot2)
library(data.table)
combined.sct<-readRDS('R_outs/neuron-outs/gpe_combined-rna.rds') 
print('2 - data loaded')

DefaultAssay(combined.sct) <- "RNA";
combined.sct<- NormalizeData(object = combined.sct);
p1<-FeaturePlot(combined.sct,features = c("Cnr1","Cnr2")) ;
ggsave(file="neuron-figs/annotations/cannabinoid_markers.svg", plot=p1,width=14) 
#opc markers c("Pdgfra","Cspg4")
p1<-FeaturePlot(combined.sct,features = c("Pdgfra","Cspg4", "Olig2")) ;
ggsave(file="neuron-figs/annotations/opc_markers.svg", plot=p1);
#astro markers c("Gfap","Aqp4","Slc1a2")) #astrocyte markers
p1<-FeaturePlot(combined.sct ,features =c("Gfap","Aqp4","Slc1a2")); #astrocyte markers
ggsave(file="neuron-figs/annotations/astro_markers.svg", plot=p1);

#GABAergic neurons markers gad1 gad2
p1<-FeaturePlot(combined.sct,features =c("Gad2","Gad1","Slc32a1")) ;#GABA markers
ggsave(file="neuron-figs/annotations/GABA_markers.svg", plot=p1);

#glutamatergic neurons markers Slc17a7 Foxp2	Fezf2	Slc17a7
p1<-FeaturePlot(combined.sct,features =c("Slc17a7","Rorb","Foxp2",	"Fezf2"))  ;
ggsave(file="neuron-figs/annotations/Gluta_markers.svg", plot=p1);

#oligodendrocyte markers  c("Mobp","Mbp","Mog")) #oligodendrocyte markers
p1<-FeaturePlot(combined.sct,features =  c("Mobp","Mbp","Mog")) ;#oligodendrocyte markers
ggsave(file="neuron-figs/annotations/oligodendro_markers.svg", plot=p1);

#micro glia markers c("Csf1r","P2ry12")) #microglia markers ctss 	tyrobp 	cd83	c1qa
p1<-FeaturePlot(combined.sct,features =  c("Csf1r","P2ry12","Ctss","Cd83")); #microglia markers
ggsave(file="neuron-figs/annotations/microglia_markers.svg", plot=p1); 

#endothelial markers  
p1<-FeaturePlot(combined.sct,features =  c("Itm2a","Bsg","Zic3")); # 
ggsave(file="neuron-figs/annotations/endothelial_markers.svg", plot=p1); 
print('plots made - general types')


#neuron markers 1
p1<-FeaturePlot(combined.sct,features = c("Chat","Gad1","Gad2","Slc17a8")) ;
ggsave(file="neuron-figs/annotations/neurons_markers1.svg", plot=p1); 

#neuron markers 2
p1<-FeaturePlot(combined.sct ,features =c("Ngfr","Nefm","Tac1","Cadm2"));  
ggsave(file="neuron-figs/annotations/neurons_markers2.svg", plot=p1); 

p1<-FeaturePlot(combined.sct ,features =c("Neurod2","C1ql3","Rorb","Slc17a7"));  
ggsave(file="neuron-figs/annotations/neurons_markers3.svg", plot=p1); 

#neuron markers 4
p1<-FeaturePlot(combined.sct ,features =c("Grem","Pvalb","Fibcd1","Trdn"));  
ggsave(file="neuron-figs/annotations/neurons_markers4.svg", plot=p1); 

#neuron markers 5
p1<-FeaturePlot(combined.sct ,features =c("Elfn1","Grik3","Pappa","Cyp26b1"));  
ggsave(file="neuron-figs/annotations/neurons_markers5.svg", plot=p1); 
print('plots made - neuron figs 1-4')



#neuron markers 6 Slc17a6	Cck Six3 Vip
p1<-FeaturePlot(combined.sct ,features =c("Slc17a6","Cck","Six3","Vip"));  
ggsave(file="neuron-figs/annotations/neurons_markers6.svg", plot=p1); 
#neuron markers 7
p1<-FeaturePlot(combined.sct ,features =c("Cplx3","Trh","Igfbp4","Rgs4"));  
ggsave(file="neuron-figs/annotations/neurons_markers7.svg", plot=p1); 
print('neuron markes 7 done')
#neuron markers 8
p1<-FeaturePlot(combined.sct ,features =c("Hpcal4","Sst","Tac2","Adora1"));  
ggsave(file="neuron-figs/annotations/neurons_markers8.svg", plot=p1); 

#neuron markers 9
p1<-FeaturePlot(combined.sct ,features =c("Satb1","Cbln1","Tcf7l2","Rspo3"));  
ggsave(file="neuron-figs/annotations/neurons_markers9.svg", plot=p1); 
print('neuron markes 9 done')
#neuron markers 10
p1<-FeaturePlot(combined.sct ,features =c("Nefm","Vipr2","Drd1","Pdyn"));  
ggsave(file="neuron-figs/annotations/neurons_markers10.svg", plot=p1); 

#neuron markers 11
p1<-FeaturePlot(combined.sct ,features =c("Adora2a","Cadm2","Pde1c","Sphkap"));  
ggsave(file="neuron-figs/annotations/neurons_markers11.svg", plot=p1); 
#neuron markers 12
p1<-FeaturePlot(combined.sct ,features =c("Tmem255a","Cpne4","Th"));  
ggsave(file="neuron-figs/annotations/neurons_markers12.svg", plot=p1); 

p1<-FeaturePlot(combined.sct,features = c("Penk")) ;
ggsave(file="neuron-figs/annotations/neurons_markers_PENK.svg", plot=p1); 

p1<-FeaturePlot(combined.sct,features = c("Npas1")) ;
ggsave(file="neuron-figs/annotations/neurons_markers_Npas1.svg", plot=p1); 
print('neuron plots done')
DefaultAssay(combined.sct) <- "integrated"
combined.sct <- FindClusters(combined.sct, resolution = 0.4)
#view data 
p2 <- DimPlot(combined.sct, reduction = "umap", label = TRUE, repel = TRUE);
ggsave(file="neuron-figs/integrated_clusters_umap.svg", plot=p2);
DefaultAssay(combined.sct) <- "RNA"
md <- combined.sct@meta.data %>% as.data.table;
md[, .N, by = c("seurat_clusters")]; 
#########################################################################################
#0-10
# 0,6  = 'Oligodendrocytes' Mbp, Mobp
#  7,10  = 'Microglia' 'Ctss','Csf1r'
#  5   = 'OPCs' Pdgfra, Tnr 
# 1,9   = 'Astrocytes' Gja1, Aqp4 , Slc1a2
# NA = 'Glutamatergic' (not expecting any in GPe but we confirmed there are none) 'Slc17a6' Tcf7l2
#   = 'GABAergic ' Gad1, Gad2
##'GABAergic 1' = 3  .....'GABAergic 2' = 4  ...'GABAergic 3' = 8  ....'GABAergic 4' =2
###################use below line only 4238
#'Oligodendrocytes', 'Astrocytes','GABAergic 4','GABAergic 1','GABAergic 2','OPCs','Oligodendrocytes' , 'Microglia' ,'GABAergic 3', 'Astrocytes','Microglia'
new.cluster.ids <- c('Oligodendrocytes', 'Astrocytes','GABAergic 4','GABAergic 1','GABAergic 2','OPCs','Oligodendrocytes' , 'Microglia' ,'GABAergic 3', 'Astrocytes','Microglia');
#merge gaba 4 and 5
names(new.cluster.ids) <- levels(combined.sct)
combined.sct <- RenameIdents(combined.sct, new.cluster.ids)
combined.sct <- StashIdent(combined.sct, save.name = "cluster.specific")
p1<-DimPlot(combined.sct, reduction = "umap")
ggsave(file="neuron-figs/annotations/separateGABA_labelled_umap.svg", plot=p1);  
Idents(combined.sct) <- "cluster.specific"
levels(Idents(combined.sct))
# "Oligodendrocytes" ,"Astrocytes","GABAergic" ,"GABAergic" ,"GABAergic","OPCs","Microglia" ,"GABAergic"
new.cluster.ids <- c( "Oligodendrocytes" ,"Astrocytes","GABAergic" ,"GABAergic" ,"GABAergic","OPCs","Microglia" ,"GABAergic");
names(new.cluster.ids) <- levels(combined.sct)
combined.sct <- RenameIdents(combined.sct, new.cluster.ids)
combined.sct <- StashIdent(combined.sct, save.name = "cluster.broad")
p1<-DimPlot(combined.sct, reduction = "umap")
ggsave(file="neuron-figs/annotations/groupedGABA_labelled_umap.svg", plot=p1)
print('check plots')
## ONLY UNCOMMENT BELOW ONCE YOU ARE SURE THIS IS THE LAST ANNOTATION
Idents(combined.sct) <- "cluster.specific"
saveRDS(combined.sct,'R_outs/neuron-outs/gpe_combined-rna-annotated.rds')
print('saved RDS')
##########
p1<-DimPlot(combined.sct, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(file="neuron-figs/annotations/labelled_umap_neurons.svg", plot=p1); 
all<-FindAllMarkers(combined.sct,only.pos = TRUE,logfc.threshold = 0.3)
write.csv(all,'R_outs/neuron-outs/marker_genes_each_cluster.csv')
print('found marker genes')
DefaultAssay(combined.sct) <- "RNA"
combined.sct<- NormalizeData(object = combined.sct)
Idents(combined.sct) <- "cluster.broad"
#dotplot vis 
levels(combined.sct) <- c("Oligodendrocytes","OPCs", "Microglia","Astrocytes","GABAergic"); 
feats<-rev(c('Mog','Mobp','Pdgfra', 'Cspg4','Ctss','Csf1r','Gja1','Aqp4', 
        'Gad1','Gad2'));
p1<-DotPlot(combined.sct, features = feats)+ theme(axis.text.x = element_text(angle = 45, hjust=1));
ggsave(file="neuron-figs/annotations/Dotplot_markergenes-general.svg", plot=p1,width = 8, height = 4.5); 
p1<-DimPlot(combined.sct, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend(); 
ggsave(file="neuron-figs/annotations/labelled_umap_broad.svg", plot=p1); 
print('made general dot plot')
#see our transgene
combined.sct2 <-subset(combined.sct,idents = c("GABAergic"));
combined.sct<- NormalizeData(object = combined.sct)
Idents(combined.sct2) <- "cluster.specific"; 
levels(combined.sct2) <- rev(c('GABAergic 1','GABAergic 2', 'GABAergic 3', 'GABAergic 4'));
feats<-rev(c("Col8a1","Slc44a5", "Foxp2",  "Pde1c", "Syt1",
  "Tafa2","Dpp10","Asic2","Robo2",
        "Robo1","Gabrb2","Gad2","Gad1"));  
##4 = Sst, 
p1<-DotPlot(combined.sct2, features = feats)+ theme(axis.text.x = element_text(angle = 45, hjust=1));
ggsave(file="neuron-figs/annotations/neuron-Dotplot_markergenes.svg", plot=p1,width = 9, height = 4); 
p1<-DimPlot(combined.sct2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend(); 
ggsave(file="neuron-figs/annotations/labelled_umap_neurons-only.svg", plot=p1); 
print('complete');
# make heatmap 

library(svglite)
combined.sct<-readRDS('R_outs/neuron-outs/gpe_combined-rna-annotated.rds')
DefaultAssay(combined.sct) <- "RNA";
combined.sct<- NormalizeData(combined.sct)
# npas1 between pde1c and foxp2
combined.sct<- ScaleData(combined.sct)
Idents(combined.sct) <- "cluster.specific";  
levels(combined.sct) <- c("Oligodendrocytes","OPCs", "Microglia","Astrocytes","GABAergic 4","GABAergic 3","GABAergic 2","GABAergic 1"); 
feats<-rev(c('Mog','Mobp','Pdgfra', 'Cspg4','Tmem119', 'C1qb','Gja1','Aqp4', 
        'Col8a1','Slc44a5', 'Foxp2', 'Npas1',  'Pde1c', 'Chat',
        'Dpp10','Asic2','Lhx6',
        'Robo1','Gabrb2','Pvalb','Gad1','Gad2')); 
svglite('neuron-figs/annotations/heatmap.svg',width=10, height=8)
DoHeatmap(
  combined.sct,
  features = feats, 
  group.bar = TRUE, 
  slot = "scale.data", 
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = FALSE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)+ scale_fill_gradientn(colors = c("grey4", "grey70", "turquoise3"))
dev.off() 
svglite('neuron-figs/annotations/dotplot-all.svg',width=10, height=5)
DotPlot(combined.sct, features = feats)+ theme(axis.text.x = element_text(angle = 45, hjust=1));
dev.off() 

##subset of genes
combined.sct<-readRDS('R_outs/neuron-outs/gpe_combined-rna-annotated.rds')
DefaultAssay(combined.sct) <- "RNA";
combined.sct<- NormalizeData(combined.sct)
# npas1 between pde1c and foxp2
combined.sct<- ScaleData(combined.sct)
Idents(combined.sct) <- "cluster.specific";  
levels(combined.sct) <- rev(c("Oligodendrocytes","OPCs", "Microglia","Astrocytes","GABAergic 4","GABAergic 3","GABAergic 2","GABAergic 1"))
colors = rev(c("#F8766D",         "#DE8C00","#F564E3",  "#C77CFF",   "#00B4F0",    "#00BFC4",    "#00BA38",    "#7CAE00"))
p1<-DimPlot(combined.sct, reduction = "umap",cols=colors)
ggsave(file="neuron-figs/annotations/separateGABA_labelled_umap.svg", plot=p1,height=7,width=8);  

Idents(combined.sct) <- "cluster.specific";  
feats<-rev(c('Mog','Mobp','Pdgfra', 'Cspg4','Tmem119','Csf1r','Gja1','Aqp4', 
        'Meis2','Foxp2', 'Npas1',
        'Etv1',  'Lhx6','Pvalb',
             'Gad1','Gad2')); 
levels(combined.sct) <- rev(c("Oligodendrocytes","OPCs", "Microglia","Astrocytes","GABAergic 4","GABAergic 3","GABAergic 2","GABAergic 1"))
colors = rev(c("#F8766D",         "#DE8C00","#F564E3",  "#C77CFF",   "#00B4F0",    "#00BFC4",    "#00BA38",    "#7CAE00"))

### all genes are known marker genes except:
## Tafa1 for  GABA3
## Kcnmb2 for GABA1
## Cpne4 for  GABA2
## Mctp1 for GABA4
feats<-rev(c('Mog','Mobp','Pdgfra', 'Cspg4',"Tmem119",
	     'Csf1r','Gja1','Aqp4', 'Gad1','Gad2',
	     'Meis2','Mctp1','Foxp2','Tafa1',"Npas1",
             'Pvalb',"Cpne4","Etv1","Lhx6","Kcnmb2")); 
p1<-DoHeatmap(
  combined.sct,
  features = feats, 
  group.bar = TRUE, group.colors=colors,
  slot = "scale.data", 
  label = TRUE,
  hjust = 0,
  angle = 45,raster=FALSE,
  draw.lines = TRUE, 
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)+ scale_fill_gradientn(colors = c("grey100","lightsteelblue2",  "red4"))
ggsave(file="neuron-figs/annotations/heatmap-dpi100.svg", plot=p1,dpi=100,height=7,width=9);  
ggsave(file="neuron-figs/annotations/heatmap-dpi80.svg", plot=p1,dpi=80,height=7,width=8);   
ggsave(file="neuron-figs/annotations/heatmap-dpi50.svg", plot=p1,dpi=50,height=7,width=7);  
ggsave(file="neuron-figs/annotations/heatmap-dpi40.svg", plot=p1,dpi=40,height=7,width=6);  
ggsave(file="neuron-figs/annotations/heatmap-dpi30.svg", plot=p1,dpi=30,height=7,width=5);   
p1<-DoHeatmap(
  combined.sct,
  features = feats, 
  group.bar = TRUE, group.colors=colors,
  slot = "scale.data", 
  label = TRUE,
  hjust = 0,
  angle = 45, 
  draw.lines = TRUE, 
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)+ scale_fill_gradientn(colors = c("grey100","lightsteelblue2",  "red4"))
ggsave(file="neuron-figs/annotations/heatmap.svg", plot=p1,height=7,width=9);   
ggsave(file="neuron-figs/annotations/heatmap-auto.svg", plot=p1,height=7,width=8);  
ggsave(file="neuron-figs/annotations/heatmap-subset.svg", plot=p1,height=7,width=7);
ggsave(file="neuron-figs/annotations/heatmap1.svg", plot=p1,height=7,width=7);   
ggsave(file="neuron-figs/annotations/heatmap2.svg", plot=p1,height=7,width=6);  
ggsave(file="neuron-figs/annotations/heatmap3.svg", plot=p1,height=7,width=5);
#################FIND MARKER GENES
all<-FindAllMarkers(combined.sct,logfc.threshold = 0.5,only.pos = TRUE)
write.csv(all,"neuron-figs/annotations/marker_genes-each-clust.csv") 
#################FIND MARKER GENES done
levels(combined.sct) <- c("Oligodendrocytes","OPCs", "Microglia","Astrocytes","GABAergic 4","GABAergic 3","GABAergic 2","GABAergic 1"); 
svglite('neuron-figs/annotations/dotplot-reset.svg',width=10, height=7)
DotPlot(combined.sct, features = feats)+ theme(axis.text.x = element_text(angle = 45, hjust=1));
dev.off() 
## make Drd1-5 pltos
feats<-c("Drd1","Drd2","Drd3","Drd4","Drd5" ,
        'Etv1',  'Lhx6','Pvalb',
             'Gad1','Gad2'); 
p1<-FeaturePlot(combined.sct,features =c("Drd1","Drd2","Drd3","Drd4","Drd5")) ;
ggsave(file="neuron-figs/annotations/Drd1-5_markers-umap.svg", plot=p1,width=10, height=15);
svglite('neuron-figs/annotations/dotplot-Drd1-5_markers.svg',width=10, height=6)
DotPlot(combined.sct, features = feats)+ theme(axis.text.x = element_text(angle = 45, hjust=1));
dev.off()
plot_as<-VlnPlot(
  combined.sct,cols=colors,
  features="Drd1",pt.size = 1)
ggsave(file='neuron-figs/annotations/Drd1-vln.svg', plot=plot_as,height=3.5,width=5)
plot_as<-VlnPlot(
  combined.sct,cols=colors,
  features="Drd2",pt.size = 1)
ggsave(file='neuron-figs/annotations/Drd2-vln.svg', plot=plot_as,height=3.5,width=5)
plot_as<-VlnPlot(
  combined.sct,cols=colors,
  features="Drd3",pt.size = 1)
ggsave(file='neuron-figs/annotations/Drd3-vln.svg', plot=plot_as,height=3.5,width=5)
plot_as<-VlnPlot(
  combined.sct,cols=colors,
  features="Drd5",pt.size = 1)
ggsave(file='neuron-figs/annotations/Drd5-vln.svg', plot=plot_as,height=3.5,width=5)
#pt.size = 0

Idents(combined.sct) <- "cluster.broad"
neurs.sct<-subset(x =combined.sct, idents = c("GABAergic"))
Idents(neurs.sct) <-"cluster.specific"
neurs.sct<- NormalizeData(object =neurs.sct);
levels(neurs.sct) <- rev(c("GABAergic 4","GABAergic 3","GABAergic 2","GABAergic 1")) 
neurcolors<-rev(c("#00B4F0",    "#00BFC4",    "#00BA38",    "#7CAE00"))
plot_as<-VlnPlot(
  neurs.sct,cols=neurcolors,
  features="Pvalb",pt.size = 1)
ggsave(file='neuron-figs/annotations/Pvalb-exp-GABAs.svg', plot=plot_as,height=3.5,width=5)
#pt.size = 0
plot_as<-VlnPlot(
  neurs.sct,cols=neurcolors,
  features="Pvalb",pt.size = 0)
ggsave(file='neuron-figs/annotations/Pvalb-exp-GABAs-noDots.svg', plot=plot_as,height=3.5,width=5)
plot_as<-VlnPlot(
  neurs.sct,cols=neurcolors,
  features=c("Kcnq3","Kcnq5"),pt.size = 1)
ggsave(file='neuron-figs/annotations/Kcnqs-exp-GABAs.svg', plot=plot_as,height=3.5,width=5)
#pt.size = 0
#F8766D = Oligodendrocytes
#DE8C00 = OPCs

#7CAE00 = GABAergic 1 
#00BA38 = GABAergic 2 
#00BFC4 = GABAergic 3
#00B4F0 = GABAergic 4

#F564E3 = Microglia
#C77CFF = Astrocytes

# make vln plots
combined.gab12<-subset(x =neurs.sct, idents = c("GABAergic 1","GABAergic 2"))
colors3 = rev(c( "#00BA38",    "#7CAE00"))

levels(combined.gab12) <-rev( c("GABAergic 2","GABAergic 1"))
plot_as<-VlnPlot(
  combined.gab12,cols=colors3,
  features=c("Drd1","Drd2","Drd3","Drd5"),pt.size = 1,ncol=4)
ggsave(file='neuron-figs/annotations/Drd1-5-vln-gaba1-22.svg', plot=plot_as,height=4,width=8)

plot_as<-VlnPlot(
  combined.gab12,cols=colors3,
  features=c("Kcnq3","Kcnq5"),pt.size = 1,ncol=4)
ggsave(file='neuron-figs/annotations/Kcnq3-5-vln-gaba1-2.svg', plot=plot_as,height=4,width=8)

# save metadata file with raw values for kevin
GENES<-rownames(combined.gab12) 
indMyGene<-which(GENES=="Kcnq3") 
CountMyGeneMyCell<-GetAssayData(object = combined.gab12, slot = 'counts')[indMyGene, ]
combined.gab12[['Kcnq3.raw']]<-CountMyGeneMyCell
indMyGene<-which(GENES=="Kcnq5") 
CountMyGeneMyCell<-GetAssayData(object = combined.gab12, slot = 'counts')[indMyGene, ]
combined.gab12[['Kcnq5.raw']]<-CountMyGeneMyCell
kcnqs35 <- combined.gab12@meta.data   
kcnqs35<- subset(kcnqs35, select=c("cluster.specific","treatment","Kcnq3.raw", "Kcnq5.raw"))
head(kcnqs35)
write.csv(kcnqs35,"neuron-figs/annotations/GABA1-2_Kcnq3-5.csv")


# save metadata file with raw values for kevin
GENES<-rownames(neurs.sct) 
indMyGene<-which(GENES=="Kcnq3") 
CountMyGeneMyCell<-GetAssayData(object = neurs.sct, slot = 'counts')[indMyGene, ]
neurs.sct[['Kcnq3.raw']]<-CountMyGeneMyCell
indMyGene<-which(GENES=="Kcnq5") 
CountMyGeneMyCell<-GetAssayData(object = neurs.sct, slot = 'counts')[indMyGene, ]
neurs.sct[['Kcnq5.raw']]<-CountMyGeneMyCell
kcnqs35 <- neurs.sct@meta.data   
kcnqs35<- subset(kcnqs35, select=c("cluster.specific","treatment","Kcnq3.raw", "Kcnq5.raw"))
head(kcnqs35)
write.csv(kcnqs35,"neuron-figs/annotations/GABA1-4_Kcnq3-5.csv")

library(dplyr) 
GABA1<-subset(kcnqs35, cluster.specific=="GABAergic 1")
GABA2<-subset(kcnqs35, cluster.specific=="GABAergic 2")
GABA3<-subset(kcnqs35, cluster.specific=="GABAergic 3")
GABA4<-subset(kcnqs35, cluster.specific=="GABAergic 4")
head(GABA1)
write.csv(GABA1,"neuron-figs/annotations/GABA1_Kcnq3-5.csv")
write.csv(GABA2,"neuron-figs/annotations/GABA2_Kcnq3-5.csv")
write.csv(GABA3,"neuron-figs/annotations/GABA3_Kcnq3-5.csv")
write.csv(GABA4,"neuron-figs/annotations/GABA4_Kcnq3-5.csv")
