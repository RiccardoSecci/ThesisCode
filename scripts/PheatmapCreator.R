
###This File is used to create the final double heatmap ordered in various ways, if type order == Sum_Active then the order will be by how much the drugs
###enriches the desired TFs. in the case of AML leukemia differentiation, then from the way it enriches with its upregulated genes the TF that promote
###Myeloid Differentiation

library(pheatmap)
library(stringr)
source("./script/ClusterMappingCMP_CELL.R")
args = commandArgs(trailingOnly=TRUE)
TypeData = args[1]
file_subset = args[2]
type_order = args[3]
type_matrix = args[4]
clustering8 = args[5]

Pheatmap_create <- function(TypeData,file_subset,type_order,type_matrix = '0.05',clustering8){


	subset_name<-sub("GeneNames","",sub(".txt","",file_subset))
	clustering_methods_initial = read.csv('./output/ClusteringTable/AllTreatmentsClustering.txt',stringsAsFactors = FALSE, sep='\t')
	final_matrix_double_list = readRDS(paste0('./output/',TypeData,'/MeanRank/matrices/DoubleMatrices_',subset_name,'_SignMean_',type_matrix,'.rds'))
	mat_tf_deg = as.data.frame(final_matrix_double_list$final_matrix_DEG_simple)
	print(type_matrix)

	clustering_methods<-CMPClusterMaker(clustering_methods_initial,final_matrix_double_list)
	####creo i vari cluster
	k_means5 <- data.frame('kmeans5'= clustering_methods[,2],row.names =clustering_methods[,1])
	k_means8 <- data.frame('kmeans8'= clustering_methods[,3],row.names =clustering_methods[,1])
	h_clust5 <- data.frame('hclust5'= clustering_methods[,4],row.names =clustering_methods[,1])
	h_clust8 <- data.frame('hclust8' = clustering_methods[,5],row.names =clustering_methods[,1])
	ATC_code <- data.frame('ATC'=clustering_methods[,6],row.names =clustering_methods[,1])
	if(clustering8){
		clustering_methods[,3] <- as.character(clustering_methods[,3])
		clustering_methods[,5] <- as.character(clustering_methods[,5])
		clustering_methods <- clustering_methods[,c(-2,-4),drop = FALSE]
		
	}
	rownames(clustering_methods)= clustering_methods$Drug
	clustering_methods$Drug= NULL

	file_subset = sub("Names","Symbols",file_subset)
	Gene_name_df <- read.csv(paste0('./data/Subsets/',file_subset,".txt"),sep='\t', row.names = 1,header = T)
	Gene_name = rownames(Gene_name_df)

	gene_type = gsub('.txt', '', file_subset)
	colorslist <- c("#FF0000FF" ,"#CCFF00FF", "#00FF66FF" ,"#0066FFFF" ,"#CC00FFFF")

	if(grepl("Apoptosis",file_subset)){
		tfs_active= Gene_name[1:5]
		other_tf = Gene_name[6:30]
		mat_subset = mat_tf_deg[,c(tfs_active, other_tf)]
		aka3 = list(Class = c(Promotes =  colorslist[1], Inhibits=colorslist[4],Both =colorslist[3], Unknown = "grey46"))
		
	}else if(grepl("HematopoiesisGran_Mono", file_subset)){
		tfs_active= Gene_name[1:7]
		other_tf = Gene_name[8:12]
		mat_subset = mat_tf_deg[,c(tfs_active, other_tf)]
		
		aka3 = list(Class = c('Active D' = colorslist[1], 'Active GD' =colorslist[2],'Active MD' = colorslist[3], 'Inhibit D' = colorslist[4], 'Upstream TF' = colorslist[5]))
		
	}else{
		mat_subset = mat_tf_deg
	}
	
	library(RColorBrewer)
	my_palette2 <- c('blue', "green", 'red4',"grey51", 'grey21')
	breaks_heatmap = c(-2.5,-1, 1, 2.5, 3, 4)	


	sum_active <- function(row_number){

	  vector_pos = unlist(row_number[which(row_number>0 & row_number<3 )]%%1)
	  vector_neg = -unlist(row_number[which(row_number<0 )]%%1)

	  VectorSum = sum(c(vector_pos, vector_neg))
	 
	 

	  return(VectorSum)
	}

	###########nel caso vada utilizata solo sulle prime usare col_sum_active
	if(type_order == "Sum_Active"){ ####This variable is to indicate the Order of the drugs in the pHeatmap. If sumactive true then the order is calculate
	#by the way they enrich the most important TFs for the given Phenotype.
		df_active = mat_subset[,tfs_active]
		tfs_sorted = c(tfs_active, other_tf)
		col_sum_active = as.data.frame(sort(apply(df_active,1, sum_active),decreasing = TRUE))
		mat_subset_sort = mat_subset[rownames(col_sum_active),tfs_sorted]
		minim = apply(mat_subset_sort, 1, min)
		index_down = which(minim>=3)
		if(length(index_down)>=1){
			names_sorted =c(rownames(col_sum_active[-index_down,,drop = FALSE]), rownames(col_sum_active[index_down,,drop = FALSE]))
		}else{
			names_sorted = rownames(col_sum_active)
		}
		mat_subset_sort = mat_subset_sort[names_sorted,tfs_sorted ]
		
	}

	if(type_order == "Alphabet" ){
		mat_subset_sort = mat_subset[rownames(col_sum_active),tfs_sorted]
		minim = apply(mat_subset_sort, 1, min)
		index_down = which(minim>=3)
		names_alphabetical<-order(rownames(col_sum_active)[-index_down])
		names_sorted =c(rownames(col_sum_active)[names_alphabetical],rownames(col_sum_active)[index_down])
		mat_subset_sort = mat_subset_sort[names_sorted,tfs_sorted ]
	}

	if(grepl("Limited_to",type_order)){

	  cellLine <- sub("Limited_to_","",type_order)
	  correctindex <- which(grepl(cellLine,rownames(mat_subset)))	  
	  mat_subset_sort <-mat_subset[correctindex,]
	}




	pdf(paste0('./output/',TypeData,'/MeanRank/matrices/Pheatmaps/NewDouble_pheatmap_',subset_name,"_",TypeData,'_Drugs&TF_new', type_matrix, '.pdf'),height = 12)

	

	pheatmap(mat_subset_sort, annotation_row =clustering_methods, annotation_colors=aka3,annotation_col = Gene_name_df, cluster_rows=FALSE, cluster_cols = FALSE, fontsize = 4, color = my_palette2,breaks = breaks_heatmap)

	dev.off()

	return("PheatmapCompletata")
}

Pheatmap_create(TypeData,file_subset,type_order,type_matrix,clustering8)


