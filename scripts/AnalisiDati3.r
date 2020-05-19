###########From this script the basic functions for data Analysis (contained in the file "AnalisiDatiFunctions.R") are launched.
###### The matrices cointaining the sets of Enriched TFs for every drug are calculated.
######The matrices are split: one for the TF enriched by the upregulated genes in the signature of each drug, one for the downregulated ones.
#####In the end, one heatmap from the UPMatrix and one for the DOWNMatrix is created
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
##### DICHIARO FUNZIONI UTILI
isEmpty <- function(x) {
  return(length(x)==0)
}


if(FALSE){
 setwd("C:/Users/ricca/Desktop/Cartella28gennaio")
  ngeni = 100
  TFtool = "Chea3"
  TypeData = "Leukemia"
}

source('./script/AnalisiDatiFunctions.R')


args = commandArgs(trailingOnly=TRUE) 
if (length(args) < 2) {
  stop("devi fornire il numero di geni che vuoi che vengano presi in considerazione e se usi AllCells, Leukemia, CellSpecific, opzionale il parametro TFtool .n", call.=FALSE)
} 
if (length(args) < 3) {
  TFtool = "ALL"
} else{
	TFtool = args[3]
}
arg1 = args[1]
if ( arg1!='logFC' & arg1 !='DEG'){
	ngeni = as.numeric(args[1]) # number of up and downregulated genes taken into consideration
} else {
	ngeni=arg1}
#
TypeData = args[2]


library(gplots)# This R library is the one I chose to plot my heatmaps.
  

isEmpty <- function(x) {
  return(length(x)==0)
}
#this function is useful when NA are present in the matrix from which the heatmap is calculated.
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}

#Here various cases are treated, based on the chosen Library of the transcription Factor enrichment tool, Chea3
if(TFtool == "ARCHS4" || TFtool == "ALL" ){
			TFtool ='ARCHS4'
			print(TFtool)
			pathtrattamenti <- paste0('output/',TypeData, "/", TFtool, "/", ngeni, "geni/") 
			elencofile <- list.files(pathtrattamenti) 
			print(head(elencofile))
			  
			matrices_list = chea3_analysis(elencofile,pathtrattamenti)  
			percorsosalvataggio = paste0('output/',TypeData, '/', TFtool,"/matrices/")

			################# 
			dir.create(percorsosalvataggio, showWarnings = FALSE)
			filesalvataggioup <- paste0(percorsosalvataggio,'matriceup_',ngeni, '.txt')
			filesalvataggiosimpleup <-paste0(percorsosalvataggio,'matricesimpleup_',ngeni, '.txt')
			filesalvataggiodown <- paste0(percorsosalvataggio,'matricedown_',ngeni, '.txt')
			filesalvataggiosimpledown <-paste0(percorsosalvataggio,'matricesimpledown_',ngeni, '.txt')

			write.table(matrices_list[[1]],filesalvataggioup,sep = '\t',  quote = FALSE)
			write.table(matrices_list[[3]],filesalvataggiosimpleup,sep = '\t',  quote = FALSE)
			write.table(matrices_list[[2]],filesalvataggiodown, sep = '\t',  quote = FALSE)
			write.table(matrices_list[[4]],filesalvataggiosimpledown,sep = '\t',  quote = FALSE)

			plots_tfs = TF_counts_plots(percorsosalvataggio, TFtool, n_geni= ngeni)
			
			}

			
if (length(args) < 3) {
  TFtool = "ALL"
} 
####################################################################################################################
####################################################################################################################
####################################################################################################################
if (TFtool == "Chea3"  || TFtool == "ALL" ){

			TFtool ='Chea3'
			print(TFtool)
			pathtrattamenti <- paste0('output/',TypeData, "/", TFtool, "/", ngeni, "geni/") 	
			print(pathtrattamenti)
			elencofile <- list.files(pathtrattamenti) 
			print(elencofile)
			
	
			
			matrices_list = chea3_analysis(elencofile, pathtrattamenti)  
			percorsosalvataggio = paste0('output/',TypeData, '/', TFtool,"/matrices/")


			################# 
			dir.create(percorsosalvataggio, showWarnings = FALSE)
			filesalvataggioup <- paste0(percorsosalvataggio,'matriceup_',ngeni, '.txt')
			filesalvataggiosimpleup <-paste0(percorsosalvataggio,'matricesimpleup_',ngeni, '.txt')
			filesalvataggiodown <- paste0(percorsosalvataggio,'matricedown_',ngeni, '.txt')
			filesalvataggiosimpledown <-paste0(percorsosalvataggio,'matricesimpledown_',ngeni, '.txt')

			write.table(matrices_list[[1]],filesalvataggioup,sep = '\t',  quote = FALSE)
			write.table(matrices_list[[3]],filesalvataggiosimpleup,sep = '\t',  quote = FALSE)
			write.table(matrices_list[[2]],filesalvataggiodown, sep = '\t',  quote = FALSE)
			write.table(matrices_list[[4]],filesalvataggiosimpledown,sep = '\t',  quote = FALSE)
			plots_tfs = TF_counts_plots(percorsosalvataggio,TFtool,ngeni) 

} 


if (length(args) < 3) {
  TFtool = "ALL"
} 

sign_column = 'SignMean_0.05' # The value of the DEG gene defined as significant. you could also use 0.1

if (TFtool == "MeanRank"  || TFtool == "ALL" ){
  
      TFtool ='MeanRank'
      print(TFtool)
      pathtrattamenti <- paste0('output/',TypeData, "/", TFtool, "/", ngeni, "geni/") 		
      print(pathtrattamenti)
      elencofile <- list.files(pathtrattamenti) 
      print(elencofile)
      
      
      
      matrices_list = chea3_MeanRank_analysis(elencofile, pathtrattamenti, sign_column)  
      percorsosalvataggio = paste0('output/',TypeData, '/', TFtool,"/matrices/")
      
      
      ################# 
      dir.create(percorsosalvataggio, showWarnings = FALSE)
      filesalvataggiosimpleup <-paste0(percorsosalvataggio,'matricesimpleup_',ngeni,'_',sign_column, '.txt')
      filesalvataggiosimpledown <-paste0(percorsosalvataggio,'matricesimpledown_',ngeni,'_',sign_column, '.txt')
      
      write.table(matrices_list[[3]],filesalvataggiosimpleup,sep = '\t',  quote = FALSE)
      write.table(matrices_list[[4]],filesalvataggiosimpledown,sep = '\t',  quote = FALSE)
      plots_tfs = TF_counts_plots(percorsosalvataggio,TFtool,sign_column, ngeni) 
  
} 


