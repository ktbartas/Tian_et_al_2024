setwd('/media/emem/Lab/kbartas/gpe_cocsal');
### REMEMBER TO BE IN CORRECT ENVIRONMENT
# conda activate signac
##load the libraries needed
library(Seurat) ;
library(Signac); 
library(magrittr);
library(ggplot2);
library(data.table);
library(ggplot2);
library(GenomicRanges);
library(future);
library(EnsDb.Mmusculus.v79) ;
library(stringr);
#load 2021 data
cont.data21 <- Read10X(data.dir = "data/raw_feature_bc_matrix_saline2021/raw_feature_bc_matrix_saline2021");
coca.data21 <- Read10X(data.dir = "data/raw_feature_bc_matrix_cocaine2021/raw_feature_bc_matrix_cocaine2021");
#control1 <- CreateSeuratObject(counts = cont.data, min.cells = 1, min.features = 20);
#cocaine1 <- CreateSeuratObject(counts = coca.data, min.cells = 1, min.features = 20);
#load 2023 data

cont.data23 <- Read10X(data.dir = "data/raw_feature_bc_matrix_saline2023/raw_feature_bc_matrix");
coca.data23 <- Read10X(data.dir = "data/raw_feature_bc_matrix_cocaine2023/raw_feature_bc_matrix");
#control2 <- CreateSeuratObject(counts = cont.data, min.cells = 2, min.features = 50);
#cocaine2 <- CreateSeuratObject(counts = coca.data, min.cells = 2, min.features = 50);

List1<-rownames(cont.data23)
List2<-rownames(cont.data21) # get rownames 
keepit<-Reduce(intersect,list(List1,List2)) # keep rownames in both
cont_subsetdata23 <- cont.data23[rownames(cont.data23) %in% keepit, ]  # Extract rows from data 
cont_subsetdata21 <- cont.data21[rownames(cont.data21) %in% keepit, ]  # Extract rows from data 
List1<-rownames(coca.data23)
List2<-rownames(coca.data21)
keepit<-Reduce(intersect,list(List1,List2))
coca_subsetdata23 <- coca.data23[rownames(coca.data23) %in% keepit, ]  # Extract rows from data 
coca_subsetdata21 <- coca.data21[rownames(coca.data21) %in% keepit, ]  # Extract rows from data 
List1<-rownames(coca_subsetdata23)
List2<-rownames(coca_subsetdata21)
identical(List1, List2) #check that they're the same
List1<-rownames(cont_subsetdata23)
List2<-rownames(cont_subsetdata21)
identical(List1, List2)

#######################################^^ reduced to only common rownames. now reduce to only common col names -> next
List1<-colnames(cont.data23)
List2<-colnames(cont.data21)
keepit<-Reduce(intersect,list(List1,List2))
cont_subsetdata23a <- cont_subsetdata23[,colnames(cont_subsetdata23) %in% keepit ]  # Extract cols from data 
cont_subsetdata21a <- cont_subsetdata21[,colnames(cont_subsetdata21) %in% keepit ]  # Extract cols from data 
List1<-colnames(cont_subsetdata23a)
List2<-colnames(cont_subsetdata21a)
identical(List1, List2)
#control reduced to only common colnames

List1<-colnames(coca.data23)
List2<-colnames(coca.data21)
keepit<-Reduce(intersect,list(List1,List2))
coca_subsetdata23a <- coca_subsetdata23[,colnames(coca_subsetdata23) %in% keepit ]  # Extract cols from data 
coca_subsetdata21a <- coca_subsetdata21[,colnames(coca_subsetdata21) %in% keepit ]  # Extract cols from data 
List1<-colnames(coca_subsetdata23a)
List2<-colnames(coca_subsetdata21a)
identical(List1, List2)
#cocaine reduced to only common colnames

coca_use<-coca_subsetdata23a+coca_subsetdata21a
cont_use<-cont_subsetdata23a+cont_subsetdata21a 
control1 <- CreateSeuratObject(counts = cont_use, min.cells = 2, min.features = 50);
cocaine1 <- CreateSeuratObject(counts = coca_use, min.cells = 2, min.features = 50);

control1 <- PercentageFeatureSet(control1, pattern = "^mt-", col.name = "percent.mt");
cocaine1 <- PercentageFeatureSet(cocaine1, pattern = "^mt-", col.name = "percent.mt"); 
##

g1<-VlnPlot(control1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3);
g2<-VlnPlot(cocaine1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3);
ggsave(file="neuron-figs/QC/control_QC_before1.svg", plot=g1,width = 8, height = 12);
ggsave(file="neuron-figs/QC/cocaine_QC_before1.svg", plot=g2,width = 8, height = 12);
print('made init violin plots')
#subset
control1 <- subset(control1, subset = nFeature_RNA < 2000 & nFeature_RNA > 150 & nCount_RNA < 4000 & nCount_RNA > 250 & percent.mt < 10);
cocaine1 <- subset(cocaine1, subset = nFeature_RNA < 2000 & nFeature_RNA > 150 & nCount_RNA < 4000 & nCount_RNA > 250 & percent.mt < 10); 

g1<-VlnPlot(control1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3);
g2<-VlnPlot(cocaine1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3);
ggsave(file="neuron-figs/QC/control_QC_after1.svg", plot=g1,width = 8, height = 12);
ggsave(file="neuron-figs/QC/cocaine_QC_after1.svg", plot=g2,width = 8, height = 12);
print('RNA QC complete');

cocaine1[["treatment"]]<- "Cocaine" ;
control1[["treatment"]]<- "Control" ;
md <- control1@meta.data %>% as.data.table;
md[, .N, by = c("treatment")];
md <- cocaine1@meta.data %>% as.data.table;
md[, .N, by = c("treatment")];
# old : 6,781  saline-treated group  6,203 nuclei from  cocaine-treated group...total of 12,984 nuclei.
saveRDS(control1,'R_outs/control_2023.rds');
saveRDS(cocaine1,'R_outs/cocaine_2023.rds');
print('objects saved')
