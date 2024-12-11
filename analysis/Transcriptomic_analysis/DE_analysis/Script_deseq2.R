#### DE Seq analysis using DESeq2
## Developped by Marion Verdenaud and Gabrièle Adam for IPS2

library("DESeq2")
library(ggplot2)
library(tidyr)
library(tibble)
library(plyr)
library(dplyr)
library(reshape2)
library(gplots)
library("ggrepel")
library("pheatmap")
library(stringr)

# define parameters
min_cts<-5
output_dir<-"deseq2_analysis"
lim_pval<-0.05
lim_fc<-0
interaction<-FALSE
nbgenes_clustering<-50
nbgenes_profiles<-20
replicate<-TRUE
targets<-"stage"

cts<-read.csv("input/Cts_COUNTS.csv")
metadata<-read.csv("input/Cts_TARGET.csv")
contrasts<-readRDS("../Quality_Control/DiffAnalysis/Contrasts.rds)
sessioninfopath<-"SessionInfo.txt"

##Filter condition handeling
Project_Name <- "Cts"

##Filter according to filter condition (or not if not defined)
Filter_c=NULL
Filter_v<-NULL
Filter_condition=''
# print (paste("filterdondition",Filter_condition))
##If defined in configuration file
if (Filter_condition != Project_Name)
{	
  fc0<-strsplit(Filter_condition,"_")
  fc<-strsplit(fc0[[1]][2],"-")
  Filter_c<-fc[[1]][1]
  Filter_v<-fc[[1]][2]
}


## MODEL ##

# Define model
# if length list of condition minus Replicates is sup to 1, then create model as addition of term, else, only one factor 


#Here we are using target file
list_meta<-c(targets)

if (! is.null(Filter_c)){
	list_meta<-list_meta[-(which(list_meta == Filter_c))]
}

if (length(list_meta) == 1){
	model<-paste0("~ ",list_meta[1])
}else{
	model<-paste0("~ ",paste0(list_meta, collapse=" + "))
}

#Subset metadata to the target we defined in configuration file + Replicate ?
metadata <- subset(metadata,select=c("labels",list_meta,"Replicate"))

#Subset again condition according to filter value if it is not NULL
if (! is.null(Filter_v)){
  metadata<-metadata %>% filter(grepl(Filter_v,labels))
  cts<-cts %>% select(contains(c("Gene_ID",paste0("_",Filter_v,"_"))))
}

print("#### MODEL #####")
print(model)

# construct DESeqData Object
dds<- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design=formula(model), tidy=TRUE)
dds

# Filter all data sup to a value (default = 0)
KEEP <- rowSums(counts(dds)) > min_cts
dds<- dds[KEEP,]
dds

# run DESeq2 analysis
print("#### Run DE DEseq2 #####")
dds<-DESeq(dds)

# get counts normalized
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts,file=paste0(output_dir,"/Normalised_counts.csv"),row.names=TRUE, quote=FALSE, col.names=NA, sep=",")


## DE ##

# Diff analysis contrats get from DicoExpress 

#Get a vector with the correspondance between value in the contrast and factor name.
which_factor<-c()
for (f in list_meta){
	for (v in unique(metadata[f])){
		wf<-c()
		wf[v]<-f
		which_factor<-c(which_factor,wf)
	}
}

#Deal with interaction - first orientation
#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
#Create a factor with the combinaison of the 2 factors

#Do do only if there is more than one factor
if (length(targets) > 1){
	dds_interact1<-dds
	dds_interact1$group<-factor(paste0(dds[[list_meta[1]]],dds[[list_meta[2]]]))
	design(dds_interact1)<- ~ group
	dds_interact1<-DESeq(dds_interact1)

	#Deal with interaction - second orientation
	dds_interact2<-dds
	dds_interact2$group<-factor(paste0(dds[[list_meta[2]]],dds[[list_meta[1]]]))
	design(dds_interact2)<- ~ group
	dds_interact2<-DESeq(dds_interact2)
}

#for all contrasts (from DicoExpress)
for (cond in contrasts$Contrasts){

	print(paste0("#### Begin work on comparison ",cond," #####"))
	
	# Contrast type [CSP136_FemaleFlower-CSP136_Pulp]-[Mono_FemaleFlower-Mono_Pulp] --> see https://support.bioconductor.org/p/9143932/ to integrate it For now not supported
	if (grepl("\\[.+\\]-\\[.+\\]",cond))
	{
		print ("contrast not done")
		next
	}else if (grepl("_",cond)){
	#It is an interaction contrast

		print ("contrast with interaction")
		#tranform dicoexpress contrast into DEseq2 format
		cond<-str_replace(cond,"\\[","")
		cond<-str_replace(cond,"\\]","")
		contrast<-as.vector(strsplit(cond,"-")[[1]])
		
		first_factor<-strsplit(contrast[1],"_")[[1]][1]
		fact<-c("group")
		contrast<-gsub("_","",contrast)
		contrast<-c(fact,contrast)
		
		#If contrast in first orientation
		if (which_factor[first_factor] == list_meta[1])
		{
			res <-results(dds_interact1,tidy = TRUE,contrast=contrast)
			res_notidy<-results(dds_interact1,tidy = FALSE,contrast=contrast)
		}
		else
		{
			res <-results(dds_interact2,tidy = TRUE,contrast=contrast)
			res_notidy<-results(dds_interact2,tidy = FALSE,contrast=contrast)
		}
	}else{
	
		print ("contrast with no interaction")
	
		#tranform dicoexpress contrast into DEseq2 format
		cond<-str_replace(cond,"\\[","")
		cond<-str_replace(cond,"\\]","")
		contrast<-as.vector(strsplit(cond,"-")[[1]])
		

		#case with no interaction
		fact<-as.vector(which_factor[contrast[1]])
		contrast<-c(fact,contrast)
		res<-results(dds,tidy = TRUE,contrast=contrast)
		res_notidy<-results(dds,tidy = FALSE,contrast=contrast)
    }
	
	# if dir doesn't exist, create it
	if(! file.exists(paste0(output_dir,"/",cond))){
		print("#### Create output directory #####")
		dir.create(paste0(output_dir,"/",cond))
		print("Output directory created")
	}
    
	# reorder results table by p-val and logFC
	# sort results table by pvalue and absolute value of logFC
	res<-res[order(res$padj,abs(res$log2FoldChange)),]
    
	#write result file
	file_name<-paste0(output_dir,"/",cond,"/","result_",cond,".csv")
    
	print("#### Write results in file #####")
	write.table(res,file_name,quote = FALSE,row.names = FALSE,col.names=TRUE,sep = ",")
    
    #end log to print summary
    sink()
    
	# print summary
    summary<-file(paste0(output_dir,"/",cond,"/Summary.txt"), open="wt")
	sink(summary)
# 	sink(summary, type="message")
	print(summary(res))
	sink()
	
	# now create figures
	print("#### Plot figures #####")
	name_pdf<-paste0(output_dir,"/",cond,"/Figures_",cond,".pdf")
	pdf(name_pdf,paper="a4")
	print("#### Plot MA-plot #####")
    
    #plot MA plot
	name_maplot<-paste("MAPlot for ",cond)
	pma<-plotMA(res_notidy, main=name_maplot) # MAplot for the comparison
	print(pma)
    
	#plot volcano plot with ggrepel for top 10 fc
	print("#### Plot volcano plot #####")
	res2<-res_notidy %>%  as.data.frame() %>% add_column(gene=rownames(res_notidy)) %>% mutate(gene_diff_expr= padj < lim_pval & abs(log2FoldChange) >=lim_fc)
	res2<-res2[order(res2$padj,abs(res2$log2FoldChange)),]
	res3<-res2[1:10,]
	p<-ggplot(res2)+geom_point(aes(x = log2FoldChange, y = -log10(padj),label=gene, colour = gene_diff_expr))+
      geom_text_repel(data=res3,aes(x=log2FoldChange,y=-log10(padj),label=gene))+
      xlab("log2 fold change")+ylab("-log10 adjusted pvalue")+
      theme_bw() +ggtitle("Volcano Plot")
	print(p)
    
    #plot barplot des p-values
	print("#### Plot p-values barplot #####")
    p<-ggplot(res2, aes(x=padj))+geom_histogram()+
      xlab("adjusted p-value")+ylab("Counts")+
      theme_bw() +ggtitle("Distribution of adjusted p-value")
	print(p)
	dev.off()
	
	# plot counts for topnbgenes_profiles profile
	print("#### Plot Profile #####")
	name_pdf<-paste0(output_dir,"/",cond,"/Top",nbgenes_profiles,"_profile_",cond,".pdf")
	pdf(name_pdf,paper="a4")
    #par(mar = rep(2, 4))
	for(i in 1:nbgenes_profiles){
      for( meta in list_meta[-length(list_meta)]){
        plotCounts(dds, gene=res2[i,7], intgroup=meta)
      }
	}
	dev.off()

    # Only do this part if there is some DE genes
    if(length(res2[res2$gene_diff_expr=="TRUE","gene"]) > 2){
    # list of genes DEG
		print("#### Plot heatmap #####")
		name_pdf<-paste0(output_dir,"/",cond,"/Top",nbgenes_clustering,"_clustering_",cond,".pdf")
# 		pdf(name_pdf,paper="a4")
		res2<-res2[order(res2$padj,abs(res2$log2FoldChange)),]
		list_DEG<-res2[res2$gene_diff_expr=="TRUE","gene"][1:nbgenes_clustering]
		# plot clustering
		# add 1 to that log is always positive
		DEG_log2NormCounts<-log2(normalized_counts[rownames(normalized_counts) %in% list_DEG,]+1)
		df <- as.data.frame(colData(dds)[,list_meta])
		colnames(df)<-list_meta
		rownames(df)<-rownames(colData(dds))
		
		tryCatch({
		  p<-pheatmap(as.matrix(DEG_log2NormCounts), filename= name_pdf,show_rownames=TRUE, annotation_col=df, cluster_rows=TRUE,cluster_cols=TRUE, cellwidth = 8, fontsize_row = 5)
		},
		error=function(cond){
		  message(paste("error catched"))
		  message(cond)
		}
		)
	} #end of if there is some DEG
}


#print sessionInfo
sink(sessioninfopath)
sessionInfo()
sink()
