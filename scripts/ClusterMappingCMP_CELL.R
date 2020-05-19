
#With this function I calculate the values for the ATC_code, K_means and H_clust of each drug.
#the K_means and H_clust values are the results of a clustering algorithm that clusters drug closer together by calculating various properties.
#5 and 8 refers to the number of clusters created.

library(stringr)

CMPClusterMaker <- function(clustering_data,final_matrix){
    clustering_methods = clustering_data
    final_matrix_double_list = final_matrix
    mat_tf_deg = as.data.frame(final_matrix_double_list$final_matrix_DEG_simple)
    
    cluster_names<-clustering_methods$Drug
    matrix_names<- row.names(mat_tf_deg)
    
    groups <- expand.grid(cluster_names, matrix_names, stringsAsFactors = FALSE)
    colnames(groups)<-c("cluster","matrix")
    groups_mapped <-groups[groups$cluster == sub("_HL60|_PC3|_MCF7","",groups$matrix), ]
    
    new_clustering_methods<-merge(groups_mapped,clustering_methods, by.x="cluster",by.y="Drug",all = TRUE)
    final_clustering_methods<- new_clustering_methods[,-1]
    colnames(final_clustering_methods)<-colnames(clustering_methods)
    final_clustering_methods<-final_clustering_methods[!is.na(final_clustering_methods$Drug),]
  return(final_clustering_methods)
}


ClusterMaker <- function(clustering_data,final_matrix){
    clustering_methods = clustering_data
    final_matrix_double_list = final_matrix
    mat_tf_deg = as.data.frame(final_matrix_double_list)
    
    cluster_names<-clustering_methods$Drug
    matrix_names<- row.names(mat_tf_deg)
    
    groups <- expand.grid(cluster_names, matrix_names, stringsAsFactors = FALSE)
    colnames(groups)<-c("cluster","matrix")
    groups_mapped <-groups[groups$cluster == sub("(_HL60_|_PC3_|_MCF7_)*(UP|DOWN)*","",groups$matrix), ]
    
    new_clustering_methods<-merge(groups_mapped,clustering_methods, by.x="cluster",by.y="Drug",all = TRUE)
    final_clustering_methods<- new_clustering_methods[,-1]
    colnames(final_clustering_methods)<-colnames(clustering_methods)
    final_clustering_methods<-final_clustering_methods[!is.na(final_clustering_methods$Drug),]
	
	k_means5 <- data.frame('kmeans5'= final_clustering_methods[,2],row.names =final_clustering_methods[,1])
	k_means8 <- data.frame('kmeans8'= final_clustering_methods[,3],row.names =final_clustering_methods[,1])
	h_clust5 <- data.frame('hclust5'= final_clustering_methods[,4],row.names =final_clustering_methods[,1])
	h_clust8 <- data.frame('hclust8' = final_clustering_methods[,5],row.names =final_clustering_methods[,1])
	ATC_code <- data.frame('ATC'=final_clustering_methods[,6],row.names =final_clustering_methods[,1])
 	rownames(final_clustering_methods)= final_clustering_methods$Drug
	final_clustering_methods$Drug= NULL
	###########
 
 
	
	
  return(final_clustering_methods)
}
