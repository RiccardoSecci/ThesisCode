library(jsonlite)
library(httr)
base_url <- "https://amp.pharm.mssm.edu/chea3/api/enrich/"
args = commandArgs(trailingOnly=TRUE) #cos? posso fornire gli argomenti
###################
if(FALSE){
setwd("C:/Users/ricca/Desktop/ProgettoDiLaurea/LavoriPython/TF_CMAP")
  filetrattamento = './output/Leukemia/DOWN_logFC/mafenide_HL60_logFC_geni_down.txt'
file_output <- "output/provaMarzo.txt"
}
###############

filetrattamento = args[1]
print(filetrattamento)
file_output = args[2]


if(length(args)==3){
  typeofdatabase = args[3]
}else{
  typeofdatabase = "Literature--ChIP-seq"  # ARCHS4--Coexpression  Integrated--meanRank" "Integrated--topRank"  "GTEx--Coexpression
}


gene_set_almost<-as.matrix(read.csv(filetrattamento, sep='\t'))
gene_set <- paste(gene_set_almost[,1],collapse=",")
call = paste(base_url,gene_set,sep='')
response = GET(call)
json = content(response,"text")
#results as R list of dataframes
results = fromJSON(json)

df_database = results[[typeofdatabase]]
if (typeofdatabase == 'Integrated--meanRank'){
  integrated_method = c('Integrated--topRank', 'Integrated--meanRank')
 
  
  get_sign_value = function(df_database,results, p_value_th){
      ind = which(names(results) %in% integrated_method) 
      p_value_rank = unlist(lapply(names(results)[-ind], function(dataset) length(which(as.numeric(results[[dataset]][,'FET p-value']) < p_value_th))))
      na_rank = unlist(lapply(names(results)[-ind], function(dataset) { data_rank = as.numeric(results[[dataset]][(which(as.numeric(results[[dataset]][,'FET p-value']) == 1 & as.numeric(results[[dataset]][,'Intersect'])==0 )), 'Rank']); if(length(data_rank)>0 ){ min(data_rank)}}))
      na_rank_mean = mean(na_rank)
      p_value_rank_mean = mean(p_value_rank)
      df_database[,paste0('SignMean_', as.character(p_value_th))] = rep(0, nrow(df_database))
      df_database[,paste0('SignMean_', as.character(p_value_th))][which(as.numeric(df_database$Score) < as.integer(p_value_rank_mean))] = 1
      df_database[,paste0('SignMean_', as.character(p_value_th))][which(as.numeric(df_database$Score) > as.integer(na_rank_mean))] = NA
      return(df_database)
      }
  
  df_database = get_sign_value(df_database,results, 0.05)
  df_database = get_sign_value(df_database,results, 0.1)
}
write.table(df_database,file_output, sep='\t', row.names=FALSE, quote=FALSE)


