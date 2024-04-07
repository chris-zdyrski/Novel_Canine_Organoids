


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ##      ###     ##       ##        ##        ##      ##   ########    ##     #######                               
##      ####     ####    ##      ####       ##         ##    ##   ##           ##    ##                                               
##     ##  ##    ## ##   ##     ##  ##      ##          ##  ##    ##           ##    ##           
##    ##    ##   ##  ##  ##    ##    ##     ##           ####      #######     ##     #######                                           
##    ########   ##   ## ##   ##########    ##            ##             ##    ##           ##                                      
##    ##    ##   ##    ####   ##      ##    ##            ##             ##    ##           ##                                          
##    ##    ##   ##     ###   ##      ##    ########      ##       #######     ##     #######                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#DATE:   2-28-2028
#GOAL:   Unsupervised analysis of organoid vs tissue samples across 4 tissues
#PLAN:

#       1.  UNSUPERVISED ANALYSIS:   
#                           -> CORRELATION MATRIX:   are replicates correlated?   what groups are there?
#                           -> PRINCIPAL COMPONENT ANALYSIS:   what are the 2 major "dimensions" of varaibility in samples
#                           -> PATHWAY ANALYSIS:  transcriptional programs
#
#################################################################################################################################################################################################################################################################### 



load("Tissue-Organoids_parsed.RData")

#REquested to remove ovaries, kidneys and Pancreas_Tissue_B118_1 and Pancreas_Tissue_B118_B118_2
samples2remove<-c("B818_Ovary_Organoid_6_2_2021","b818_ovary_tissue",
                  "B816_Kidney_Organoid","B816_Kidney_Tissue","B818_Kidney_Organoid","B818_Kidney_Tissue",
                  "B818_Pancreas_Tissue_1","B818_Pancreas_Tissue_2")

organoid_zscore<-organoid_zscore[,-which(colnames(organoid_zscore) %in% samples2remove)]

separated_organ_tissue_zscore<-separated_organ_tissue_zscore[,-which(colnames(separated_organ_tissue_zscore) %in% samples2remove)]


library("gplots")
library("RColorBrewer")

###############################################################################################################################################################
#PART 1:   ORGANOIDS VS TISSUE
###############################################################################################################################################################


#####################################################
#  CORR-MATRIX.    Published analysis
#####################################################

organoid_rpm_corr<-cor(organoid_logRPM)

heatmap.2(organoid_rpm_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "rpm corr",key.xlab = "rpm corr",cexRow = 0.6,cexCol = 0.6)

#####################################################
#  CORR-MATRIX.   Refined analysis
#####################################################

organoid_zscore_corr<-cor(organoid_zscore)

heatmap.2(organoid_zscore_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "z-score corr",key.xlab = "z-score corr",cexRow = 0.6,cexCol = 0.6)





#####################################################
#  PCA   published analysis
#####################################################

table(sapply(strsplit(colnames(organoid_zscore),"_"), function(x) x[2]))
#Bladder Endometrium      Kidney       Liver        Lung       ovary       Ovary    Pancreas 
#4       4                4             4           4           1           1           6 

#RUN PCA
organoid_pca<-prcomp(t(organoid_zscore))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)


sample_names_vector[grep("Bladder_Organoid",sample_names_vector)]<-"blue"
sample_names_vector[grep("Endometrium_Organoid",sample_names_vector)]<-"grey50"

sample_names_vector[grep("Liver_Organoid",sample_names_vector)]<-"orange"
sample_names_vector[grep("Lung_Organoid",sample_names_vector)]<-"green3"
sample_names_vector[grep("Pancreas_Organoid",sample_names_vector)]<-"red"


sample_names_vector[grep("Bladder_Tissue",sample_names_vector)]<-"darkblue"
sample_names_vector[grep("Endometrium_Tissue",sample_names_vector)]<-"black"
sample_names_vector[grep("Liver_Tissue",sample_names_vector)]<-"orange3"
sample_names_vector[grep("Lung_Tissue",sample_names_vector)]<-"darkgreen"
sample_names_vector[grep("Pancreas_Tissue",sample_names_vector)]<-"red4"



plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16)
text(organoid_pca_data[,1],organoid_pca_data[,2],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)

plot(organoid_pca_data[,3],organoid_pca_data[,4],col=sample_names_vector,pch=16)
text(organoid_pca_data[,3],organoid_pca_data[,4],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)


#####################################################
#   Heatmap:  top 1000 variance genes
#####################################################

max_var_genes<-names(sort(apply(organoid_logRPM,1,var),decreasing=T))

heatmap.2(t(organoid_zscore[max_var_genes[1:1000],]),col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "log(rpm)",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.2)




#####################################################
#   PATHWAY ANALYSIS
#####################################################

#####################
#LOAD PACKAGES
#####################


library(viper)
library(msigdbr)



#####################
#2 MAKE HALLMARK-OBJECT FOR PATHWAY ANALYSIS
#####################

#EXTRACT MOUSE-HALLMARKS:
HALLMARK_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
HALLMARK_gene_sets$gs_name<-gsub("HALLMARK_","",HALLMARK_gene_sets$gs_name)


#CONVERT MOUSE-HALLMARK DATA-FRAME TO A "REGULON"
msigdbr_list = split(x = HALLMARK_gene_sets$gene_symbol, f = HALLMARK_gene_sets$gs_name)
hallmark_regulon<-list()

for (tf in names(msigdbr_list)){
  sig_genes<-msigdbr_list[[tf]]
  
  hallmark_regulon[[tf]]<-list()
  tmp_tfmode<-rep(1,length(sig_genes))
  names(tmp_tfmode)<-sig_genes
  
  hallmark_regulon[[tf]]$'tfmode'<-tmp_tfmode
  hallmark_regulon[[tf]]$'likelihood'<-as.numeric(tmp_tfmode)
}

#####################
#RUN VIPER-Pathway analysis
#####################
organoid_Hallmarks <- viper(organoid_zscore,hallmark_regulon,method="scale")

#####################
#VISUALIZE PATHWAY ANALYSIS:
#####################

#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(organoid_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:30],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.5,mar=c(5,10))

#####################################################
#FIGURE 3:   PATHWAY BARPLOT:
#####################################################



organ_vs_tissue<-sort(rowMeans(organoid_Hallmarks[,grep("ganoid",colnames(organoid_Hallmarks))])-rowMeans(organoid_Hallmarks[,grep("issue",colnames(organoid_Hallmarks))]))

barplot(organ_vs_tissue[c(1:15,40:50)],horiz=T,las=2,cex.names = 0.4,col="black",mar=c())


###############################################################################################################################################################
#PART 2:   TISSUE VS TISSUE (visualized separatley)
###############################################################################################################################################################



#####################################################
#FIGURE 2B:   PCA   ORGANOID only
#####################################################


#RUN PCA
organoid_pca<-prcomp(t(separated_organ_tissue_zscore[,grep("ganoid",colnames(separated_organ_tissue_zscore))]))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)


sample_names_vector[grep("Bladder_Organoid",sample_names_vector)]<-"blue"
sample_names_vector[grep("Endometrium_Organoid",sample_names_vector)]<-"grey50"
sample_names_vector[grep("Liver_Organoid",sample_names_vector)]<-"orange"
sample_names_vector[grep("Lung_Organoid",sample_names_vector)]<-"green3"
sample_names_vector[grep("Pancreas_Organoid",sample_names_vector)]<-"red"



plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16)
text(organoid_pca_data[,1],organoid_pca_data[,2],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)


#####################################################
#FIGURE 2B:   PCA   TISSUE only
#####################################################


#RUN PCA
organoid_pca<-prcomp(t(separated_organ_tissue_zscore[,grep("issue",colnames(separated_organ_tissue_zscore))]))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)


sample_names_vector[grep("Bladder_Tissue",sample_names_vector)]<-"darkblue"
sample_names_vector[grep("Endometrium_Tissue",sample_names_vector)]<-"black"
sample_names_vector[grep("Liver_Tissue",sample_names_vector)]<-"orange3"
sample_names_vector[grep("Lung_Tissue",sample_names_vector)]<-"darkgreen"
sample_names_vector[grep("Pancreas_Tissue",sample_names_vector)]<-"red4"



plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16)
text(organoid_pca_data[,1],organoid_pca_data[,2],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)




#####################################################
#FIGURE 3:   PATHWAY ANALYSIS - organoid
#####################################################

#####################
#RUN VIPER-Pathway analysis
#####################
organoid_Hallmarks <- viper(separated_organ_tissue_zscore[,grep("ganoid",colnames(separated_organ_tissue_zscore))],hallmark_regulon,method="scale")

#####################
#VISUALIZE PATHWAY ANALYSIS:
#####################

#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(organoid_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:25],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.5,mar=c(5,10))

write.csv(organoid_Hallmarks,file="organoid_pathways.csv")

#####################################################
#FIGURE 3:   PATHWAY ANALYSIS - tissue
#####################################################

#####################
#RUN VIPER-Pathway analysis
#####################
organoid_Hallmarks <- viper(separated_organ_tissue_zscore[,grep("issue",colnames(separated_organ_tissue_zscore))],hallmark_regulon,method="scale")

#####################
#VISUALIZE PATHWAY ANALYSIS:
#####################

#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(organoid_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:25],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.5,mar=c(5,10))


write.csv(organoid_Hallmarks,file="tissue_pathways.csv")




