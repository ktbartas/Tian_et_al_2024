library(Seurat)
library(Matrix) 
library(tidyverse)
library(cowplot)
library(patchwork)
setwd('/media/emem/Lab/kbartas/gpe_cocsal');  
library(ggplot2); 
library(WGCNA)
library(hdWGCNA)
library(svglite)
set.seed(1)
print('loaded packages')
combined.sct<-readRDS('R_outs/neuron-outs/gpe_combined-rna-annotated.rds')
DefaultAssay(combined.sct) <- "RNA";
combined.sct<- NormalizeData(combined.sct)
print('loaded data')
Idents(combined.sct) <- "cluster.broad"
combined.sct<- SetupForWGCNA(
  combined.sct,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "WGCNA.broad" # the name of the hdWGCNA experiment
)
print('set up for WGCNA')
combined.sct<- MetacellsByGroups(
  seurat_obj = combined.sct,
  group.by = c("cluster.broad", "treatment"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cluster.broad' # set the Idents of the metacell seurat object
)
print('meta cells by groups')
# Co-expression network analysis
# set up the expression matrix
###########################################GABAergic first to test
combined.sct <- SetDatExpr(
  combined.sct,
  group_name = "GABAergic", # the name of the group of interest in the group.by column
  group.by='cluster.broad' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)
print('test gaba first')
# Test different soft powers:
combined.sct<- TestSoftPowers(
  combined.sct,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
# plot the results:
svglite('neuron-figs/hdWGCNA/GABAergic-soft-power-thresh.svg',width=10, height=8)
plot_list <- PlotSoftPowers(combined.sct)
wrap_plots(plot_list, ncol=2)
dev.off() 
power_table <- GetPowerTable(combined.sct)
head(power_table)
print('check soft power val')
###########
########################################################################################################################## EDITED TO HERE#######################################################

# construct co-expression network: 
combined.sct<- ConstructNetwork(
  combined.sct, soft_power=8,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'GABAergic' # name of the topoligical overlap matrix written to disk
)
print('network constructed')
svglite('neuron-figs/hdWGCNA/GABAergic-dendrogram.svg',width=10, height=8)
PlotDendrogram(combined.sct, main='GABAergic hdWGCNA Dendrogram')
dev.off()

## compute module eigengenes with batch correction 
combined.sct <- ScaleData(combined.sct, features=VariableFeatures(combined.sct))

# compute all MEs in the full single-cell dataset
combined.sct <- ModuleEigengenes(
 combined.sct
)
print('compute all MEs')
# harmonized module eigengenes:
hMEs <- GetMEs(combined.sct)
# compute eigengene-based connectivity (kME):
combined.sct <- ModuleConnectivity(
  combined.sct,
  group.by = 'cluster.broad', group_name = 'GABAergic'
)
combined.sct <- ResetModuleNames(
  combined.sct,
  new_name = "GABA-M"
)
print('reset gaba  module names')

svglite('neuron-figs/hdWGCNA/GABAergic-kME-plots.svg',width=9, height=12)
p <- PlotKMEs(combined.sct , ncol=3)
p
dev.off()
print('maade kme plot')

modules <- GetModules(combined.sct)
plot_list <- ModuleFeaturePlot(
  combined.sct,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
svglite('neuron-figs/hdWGCNA/GABAergic-module-feat-plots.svg',width=12, height=5)
wrap_plots(plot_list, ncol=5)
dev.off()
print('made module feature plots')
svglite('neuron-figs/hdWGCNA/GABAergic-module-correlogram.svg',width=10, height=10)
ModuleCorrelogram(combined.sct)
dev.off()
print('made module correlogram plot')

Idents(combined.sct) <- "cluster.broad"; 
gabs <-subset(combined.sct,idents = c('GABAergic'));#
targets <- gabs@meta.data 
print('made targets')
PCvalues=hMEs[rownames(targets),]
PCvalues=PCvalues[ , -which(names(PCvalues) %in% c("grey"))]
print('made pcvals')
plot_df <- cbind(select(targets, c(cluster.specific, treatment)), PCvalues)
plot_df <- reshape2::melt(plot_df, id.vars = c('cluster.specific', 'treatment'))
plot_df$cluster.specific <- factor(plot_df$cluster.specific, levels=c('GABAergic 1', 'GABAergic 2','GABAergic 3', 'GABAergic 4'))
plot_df$metacell <- paste0(plot_df$cluster.specific,plot_df$treatment)

nmodules<-4
nclusters<-5
print('defined nmods and nclusts')
library(ggsignif)
colors <- sub('GABA-M', '', as.character(levels(plot_df$variable)))
p <- ggplot(plot_df, aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot(notch=FALSE) +
  RotatedAxis() + ylab('GABA Module Eigengene') + xlab('') +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
  )
# plot width and height:
w=2*nmodules; h=2.2*nclusters;
svglite('neuron-figs/hdWGCNA/GABAergic-ME_trajectory.svg',width=w,height=h )
p + facet_wrap(cluster.specific~variable, scales='free', ncol=nmodules)+
geom_signif(aes(x=treatment,y=value),comparisons = list(c("Cocaine", "Control")),
           vjust=2,test = "wilcox.test", color="black", map_signif_level = TRUE)
dev.off()

print('made ME traject plots')
##############!! need to test& run below
library(igraph)
library(enrichR) 
ModuleNetworkPlot(combined.sct,outdir = "neuron-figs/hdWGCNA/GABAergic-ModuleNetworks")
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
combined.sct <- RunEnrichr(
  combined.sct,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)
# retrieve the output table
enrich_df <- GetEnrichrTable(combined.sct)
EnrichrBarPlot(
  combined.sct,
  outdir = "neuron-figs/hdWGCNA/GABAergic-enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)
print('made enrichedr plots')
### MAKE SVG version
wgcna_name = "WGCNA.broad"
modules <- GetModules(combined.sct, wgcna_name)
mods <- levels(modules$module)
mods <- mods[mods != "grey"]
enrichr_df <- GetEnrichrTable(combined.sct, wgcna_name)
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), 
        USE.NAMES = FALSE)
}
plot_size<-list(5,5)
outdir<- "neuron-figs/hdWGCNA/GABAergic-enrichr_plots"
n_terms<-10
logscale<-TRUE
if (!dir.exists(outdir)) {
    dir.create(outdir)
}
for (i in 1:length(mods)) {
    cur_mod <- mods[i]
    cur_terms <- subset(enrichr_df, module == cur_mod)
    print(cur_mod)
    cur_color <- modules %>% subset(module == cur_mod) %>% 
        .$color %>% unique %>% as.character
    if (nrow(cur_terms) == 0) {
        next
    }
    cur_terms$wrap <- wrapText(cur_terms$Term, 45)
    plot_list <- list()
    for (cur_db in dbs) {
        plot_df <- subset(cur_terms, db == cur_db) %>% top_n(n_terms, 
            wt = Combined.Score)
        if (cur_mod == "black") {
            text_color = "grey"
        }
        else {
            text_color = "black"
        }
        if (logscale) {
            plot_df$Combined.Score <- log(plot_df$Combined.Score)
            lab <- "Enrichment log(combined score)"
            x <- 0.2
        }
        else {
            lab <- "Enrichment (combined score)"
            x <- 5
        }
        plot_list[[cur_db]] <- ggplot(plot_df, aes(x = Combined.Score, 
            y = reorder(wrap, Combined.Score))) + geom_bar(stat = "identity", 
            position = "identity", color = "white", fill = cur_color) + theme_light() + 
            geom_text(aes(label = wrap), x = x, color = text_color, 
              size = 3.5, hjust = "left") + ylab("Term") + 
            xlab(lab) + ggtitle(cur_db) + theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), legend.title = element_blank(), 
            axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            plot.title = element_text(hjust = 0.5))
    }
    svglite(paste0(outdir, "/", cur_mod, "-bio-proc.svg"), width = plot_size[[1]], 
        height = 7) 
        print(plot_list[1]) 
    dev.off()
    svglite(paste0(outdir, "/", cur_mod, "-cell-com.svg"), width = plot_size[[1]], 
        height = 7) 
        print(plot_list[[2]]) 
    dev.off()
    svglite(paste0(outdir, "/", cur_mod, "-mol-func.svg"), width = plot_size[[1]], 
        height = 7) 
        print(plot_list[[3]]) 
    dev.off()}
combined.sct <- ModuleExprScore(
  combined.sct,
  n_genes = 25,
  method='Seurat'
)
print('calculate expr score')
MEs <- GetMEs(combined.sct, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
combined.sct@meta.data <- cbind(combined.sct@meta.data, MEs)
p <- DotPlot(combined.sct, features=mods, group.by = 'cluster.specific')
svglite("neuron-figs/hdWGCNA/GABAergic-Dotplot-specific.svg")
# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
p
# plot output
dev.off()


svglite("neuron-figs/hdWGCNA/GABAergic-Dotplot-broad.svg")
p <- DotPlot(combined.sct, features=mods, group.by = 'cluster.broad')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
p
# plot output
dev.off()
svglite("neuron-figs/hdWGCNA/GABAergic-all-mod-network.svg")
HubGeneNetworkPlot(
  combined.sct,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()
combined.sct <- RunModuleUMAP(
  combined.sct,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
umap_df <- GetModuleUMAP(combined.sct)

# plot with ggplot
svglite("neuron-figs/hdWGCNA/GABAergic-UMAP-network.svg")
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color, # color each point by WGCNA module
   size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()
svglite("neuron-figs/hdWGCNA/GABAergic-UMAP-labeled-network.svg")
ModuleUMAPPlot(
  combined.sct,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=3 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)
dev.off()

### lollipop plots

######################################################
# cocaine vs saline all cells combined.sct
######################################################
# get the modules table 
mods <- levels(modules$module)
mods <- mods[mods!='grey']
write.csv(modules,"neuron-figs/hdWGCNA/GABA-modules.csv")
# get current genotypes to compare
gt1 <- 'Cocaine'
gt2 <- 'Control'

# get list of cells with these genotypes
g1 <- combined.sct@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- combined.sct@meta.data %>% subset(treatment == gt2) %>% rownames

# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    combined.sct,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox'
)
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2)
DMEs


################################################################################
# Plot DMEs as a lollipop
################################################################################

# plot title
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)

# set plotting attributes for shape, X shape for not significant results
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)

# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))

# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])

# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module

# plot the results
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + 
  ggtitle(cur_title) + NoLegend()

# optionally add the alternating grey stripes
p <- p +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=DMEs) 

fig_dir<- "neuron-figs/hdWGCNA/GABAergic-all-cells-"
# individual plots
pdf(paste0(fig_dir, 'DMEs_lollipop.pdf'), width=4, height=3.5)
p
dev.off()
print('all cells done')

######################################################
# cocaine vs saline all cells combined.sct
######################################################
# get the modules table 
Idents(combined.sct) <- "cluster.broad"; 
gabs <-subset(combined.sct,idents = c('GABAergic'));#
# get current genotypes to compare
gt1 <- 'Cocaine'
gt2 <- 'Control'

# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames

# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    combined.sct,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox'
)
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2)
DMEs
################################################################################
# Plot DMEs as a lollipop
################################################################################

# plot title
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)

# set plotting attributes for shape, X shape for not significant results
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)

# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))

# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])

# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module

# plot the results
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + 
  ggtitle(cur_title) + NoLegend()

# optionally add the alternating grey stripes
p <- p +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=DMEs) 

fig_dir<- "neuron-figs/hdWGCNA/GABAergic-all-neurons-"
# individual plots
pdf(paste0(fig_dir, 'DMEs_lollipop.pdf'), width=4, height=3.5)
p
dev.off()
## MAKE SVG
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + theme_light() + 
  ggtitle(cur_title) + NoLegend()
#
svglite(paste0(fig_dir, 'DMEs_lollipop.svg'), width=5, height=3)
p
dev.off()
print('all neurons done')
######################################################
# cocaine vs saline all GABA
######################################################
Idents(combined.sct) <- "cluster.broad"; 
gabs <-subset(combined.sct,idents = c('GABAergic'));#
# get the modules table 
mods <- levels(modules$module)
mods <- mods[mods!='grey']

# get current genotypes to compare
gt1 <- 'Cocaine'
gt2 <- 'Control'

# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames

# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox'
)
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2)
DMEs
################################################################################
# Plot DMEs as a lollipop
################################################################################

# plot title
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)

# set plotting attributes for shape, X shape for not significant results
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)

# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))

# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])

# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module

# plot the results
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + 
  ggtitle(cur_title) + NoLegend()

# optionally add the alternating grey stripes
p <- p +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=DMEs) 

## plots by cluster ########
clustlist<-c('GABAergic 1','GABAergic 2','GABAergic 3','GABAergic 4')
print('begin for loop')

##########################################################################DELETE LOOP AS IT WRITES BAD PDFS
gabclust<-'GABAergic 1'
Idents(combined.sct) <- "cluster.specific";  
gabs <-subset(combined.sct,idents = gabclust);#
# get the modules table  

######################################################
# cocaine v saline by cluster
######################################################

# get current genotypes to compare
gt1 <- 'Cocaine'
gt2 <- 'Control'

# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames

# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox'
)
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2) 
################################################################################
# Plot DMEs as a lollipop
################################################################################

# plot title
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)

# set plotting attributes for shape, X shape for not significant results
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)

# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))

# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])

# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module

# plot the results
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + 
  ggtitle(cur_title) + NoLegend()

# optionally add the alternating grey stripes
p <- p +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=DMEs) 

fig_dir<-paste0("neuron-figs/hdWGCNA/",gabclust) 
# individual plots
pdf(paste0(fig_dir, 'DMEs_lollipop.pdf'), width=4, height=3.5)
p
dev.off()
## MAKE SVG
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + theme_light() + 
  ggtitle(cur_title) + NoLegend()
#
svglite(paste0(fig_dir, 'DMEs_lollipop.svg'), width=5, height=3)
p
dev.off()
print('all neurons done 1')

#
gabclust<-'GABAergic 2'
Idents(combined.sct) <- "cluster.specific";  
gabs <-subset(combined.sct,idents = gabclust);#
# get the modules table  

######################################################
# cocaine v saline by cluster
######################################################

# get current genotypes to compare
gt1 <- 'Cocaine'
gt2 <- 'Control'

# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames

# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox'
)
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2)
DMEs
################################################################################
# Plot DMEs as a lollipop
################################################################################

# plot title
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)

# set plotting attributes for shape, X shape for not significant results
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)

# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))

# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])

# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module

# plot the results
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + 
  ggtitle(cur_title) + NoLegend()

# optionally add the alternating grey stripes
p <- p +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=DMEs) 

fig_dir<-paste0("neuron-figs/hdWGCNA/",gabclust) 
# individual plots
pdf(paste0(fig_dir, 'DMEs_lollipop.pdf'), width=4, height=3.5)
p
dev.off()
## MAKE SVG
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + theme_light() + 
  ggtitle(cur_title) + NoLegend()
#
svglite(paste0(fig_dir, 'DMEs_lollipop.svg'), width=5, height=3)
p
dev.off()
print('all neurons done 2')
##
#
gabclust<-'GABAergic 3'
Idents(combined.sct) <- "cluster.specific";  
gabs <-subset(combined.sct,idents = gabclust);#
# get the modules table  

######################################################
# cocaine v saline by cluster
######################################################

# get current genotypes to compare
gt1 <- 'Cocaine'
gt2 <- 'Control'

# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames

# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox'
)
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2)
DMEs
################################################################################
# Plot DMEs as a lollipop
################################################################################

# plot title
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)

# set plotting attributes for shape, X shape for not significant results
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)

# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))

# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])

# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module

# plot the results
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + 
  ggtitle(cur_title) + NoLegend()

# optionally add the alternating grey stripes
p <- p +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=DMEs) 

fig_dir<-paste0("neuron-figs/hdWGCNA/",gabclust) 
# individual plots
pdf(paste0(fig_dir, 'DMEs_lollipop.pdf'), width=4, height=3.5)
p
dev.off()
## MAKE SVG
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + theme_light() + 
  ggtitle(cur_title) + NoLegend()
#
svglite(paste0(fig_dir, 'DMEs_lollipop.svg'), width=5, height=3)
p
dev.off()
print('all neurons done 3')
gabclust<-'GABAergic 4'
Idents(combined.sct) <- "cluster.specific";  
gabs <-subset(combined.sct,idents = gabclust);#
# get the modules table  

######################################################
# cocaine v saline by cluster
######################################################

# get current genotypes to compare
gt1 <- 'Cocaine'
gt2 <- 'Control'

# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames

# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox'
)
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2)
DMEs
################################################################################
# Plot DMEs as a lollipop
################################################################################

# plot title
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)

# set plotting attributes for shape, X shape for not significant results
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)

# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))

# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])

# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module

# plot the results
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + 
  ggtitle(cur_title) + NoLegend()

# optionally add the alternating grey stripes
p <- p +  ggforestplot::geom_stripes(aes(y=module), inherit.aes=FALSE, data=DMEs) 

fig_dir<-paste0("neuron-figs/hdWGCNA/",gabclust) 
# individual plots
pdf(paste0(fig_dir, 'DMEs_lollipop.pdf'), width=4, height=3.5)
p
dev.off()
## MAKE SVG
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + theme_light() + 
  ggtitle(cur_title) + NoLegend()
#
svglite(paste0(fig_dir, 'DMEs_lollipop.svg'), width=5, height=3)
p
dev.off()
print('all neurons done 4')
##########################################################################GABA1
gabclust<-'GABAergic 1'
Idents(combined.sct) <- "cluster.specific";  
gabs <-subset(combined.sct,idents = gabclust);
# get current genotypes to compare
gt1 <- 'Cocaine'
gt2 <- 'Control'
# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames
# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox')
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2) 
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)
# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))
# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])
# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module
DMEs1<-DMEs
##########################################################################GABA2
gabclust<-'GABAergic 2'
Idents(combined.sct) <- "cluster.specific";  
gabs <-subset(combined.sct,idents = gabclust);
# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames
# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox')
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2) 
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)
# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))
# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])
# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module
DMEs2<-DMEs
##########################################################################GABA3
gabclust<-'GABAergic 3'
Idents(combined.sct) <- "cluster.specific";  
gabs <-subset(combined.sct,idents = gabclust);
# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames
# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox')
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2) 
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)
# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))
# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])
# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module
DMEs3<-DMEs
##########################################################################GABA4
gabclust<-'GABAergic 4'
Idents(combined.sct) <- "cluster.specific";  
gabs <-subset(combined.sct,idents = gabclust);
# get list of cells with these genotypes
g1 <- gabs@meta.data %>% subset(treatment == gt1) %>% rownames
g2 <- gabs@meta.data %>% subset(treatment == gt2) %>% rownames
# run the DME comparison with a wilcox test
DMEs <- FindDMEs(
    gabs,
    barcodes1 = g1,
    barcodes2 = g2,
    test.use='wilcox')
DMEs$ident.1 <- gt1
DMEs$ident.2 <- gt2
DMEs$comparison <- paste0(gt1, '_vs_', gt2) 
library(ggforestplot)
cur_title <- paste0(gt1, ' vs. ', gt2)
DMEs$shape <- ifelse(DMEs$p_val_adj < 0.05, 21,4)
# order the modules by logFC of the DME test
DMEs <- DMEs %>% arrange(avg_log2FC, descending=TRUE)
DMEs$module <- factor(as.character(DMEs$module), levels=as.character(DMEs$module))
# add number of genes per module 
n_genes <- table(modules$module)
DMEs$n_genes <- as.numeric(n_genes[as.character(DMEs$module)])
# get the module colors
mod_colors <- dplyr::select(modules, c(module, color)) %>%
distinct
cp <- mod_colors$color; names(cp) <- mod_colors$module
DMEs4<-DMEs
print('4dmes made')
DMEs1$specie <- "GABAergic 1"
DMEs2$specie <- "GABAergic 2"
DMEs3$specie <- "GABAergic 3"
DMEs4$specie <- "GABAergic 4"

DMEs = rbind(DMEs1,DMEs2)
# plot the results
fig_dir<-paste0("neuron-figs/hdWGCNA/","GABA1-2") 
## MAKE SVG
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') +facet_wrap(~specie)+ 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + theme_light() + 
  ggtitle(cur_title) + NoLegend()
#
svglite(paste0(fig_dir, 'DMEs_lollipop.svg'), width=5, height=3)
p
dev.off()
print('all neurons done 1')


DMEs = rbind(DMEs1,DMEs2,DMEs3,DMEs4)
# plot the results
fig_dir<-paste0("neuron-figs/hdWGCNA/","GABA1-4") 
## MAKE SVG
p <- DMEs %>%
ggplot(aes(y=module, x=avg_log2FC, size=log(n_genes), color=module)) +
  geom_vline(xintercept=0, color='black') +
  geom_segment(aes(y=module, yend=module, x=0, xend=avg_log2FC), linewidth=0.5, alpha=0.3) +
  geom_point() +
  geom_point(shape=DMEs$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') +facet_wrap(~specie)+ 
  xlab(bquote("Avg. log"[2]~"(Fold Change)")) +
  theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust=0.5, face='plain', size=10)
  ) + theme_light() + 
  ggtitle(cur_title) + NoLegend()
#
svglite(paste0(fig_dir, 'DMEs_lollipop.svg'), width=5, height=3)
p
dev.off()
print('all neurons done 2')
print('complete - all')
