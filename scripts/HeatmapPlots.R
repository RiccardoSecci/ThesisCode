#This script has to be launched in order to plot the DOUBLE Heatmaps.
#A double heatmap is a heatmap that integrates the enrichment for the UP and DOWN enriched TFs of each drug.

source('./script/data_visualization.R')
args = commandArgs(trailingOnly=TRUE)
input_folder_matrices = args[1]
type_range = 'sign' 
output_folder = input_folder_matrices 

column_sign = 'SignMean_0.05' 
if (grepl('MeanRank' ,input_folder_matrices )){
  
  matrices_list = list.files(input_folder_matrices, pattern=paste0(column_sign, '.txt'))
} else{
  matrices_list = list.files(input_folder_matrices, pattern='.txt')
 
}

if (column_sign=='SignMean_0.05'){
  matrices_list = matrices_list[grepl('SignMean_0.1',matrices_list)==FALSE]
}

matrices_Ngene = list()
for (ngene in c(seq(100,250, 50),'logFC', 'DEG')){
	index_ngene = grep(paste0(ngene),matrices_list )
	#print(ngene)
	if (length(index_ngene)>0){
  	matrices_Ngene[[as.character(ngene)]] = list('simple' =matrices_list[index_ngene][ grep('simple', matrices_list[index_ngene])], 
	  	'tot' = matrices_list[index_ngene][- grep('simple', matrices_list[index_ngene])])
}}

listadendro <- list()
listadendro_ridotta = list()
listamatrici <- list()
listamatrici_ridotte <- list()

th_p_val = 0.05

file_subset = 'Hematopoiesis.txt'
Gene_name <- readLines(paste0('./data/Subsets/',file_subset ))

gene_type = gsub('.txt', '', file_subset)



removeRowCol = function(matrix_input, SIGN_VALUE= 0.5, na_value=2, non_sign_value=1){

        # ROW
        # removing NA value in all tf 
        somma = rowSums(abs(matrix_input))
        index_na = which(somma ==na_value*ncol(matrix_input))
        if (length(index_na)>0){
          cat( 'Number of row with all NA tfs', length(index_na ), '\n')
          matrix_input = matrix_input[-index_na,]
        }
        # COL
        somma = colSums(abs(matrix_input))
        index_na = which(somma ==na_value*nrow(matrix_input))
        if (length(index_na)>0){
          cat( 'Number of col with all NA tfs', length(index_na ), '\n')
          matrix_input = matrix_input[,-index_na]
        }
        
        # ROW
        # removing non sign tf in all tfs per row
        somma = rowSums(abs(matrix_input))
        check_sign = which(apply(matrix_input,1,'min')<=SIGN_VALUE)
        index_non_sign = which(somma ==non_sign_value*ncol(matrix_input))
		common_index = which(index_non_sign %in% check_sign)
		 index_non_sign = index_non_sign[-common_index]
        if (length(index_non_sign)>0){
          cat( 'Number of row with all non sign tfs', length(index_non_sign),'\n')
          matrix_input = matrix_input[-index_non_sign,]
        }
        
        ## COL
        somma = colSums(abs(matrix_input))
        check_sign = which(apply(matrix_input,2,'min')<=SIGN_VALUE)
        index_non_sign = which(somma ==non_sign_value*nrow(matrix_input))
        common_index = which(index_non_sign %in% check_sign)
		 index_non_sign = index_non_sign[-common_index]
        if(length(index_non_sign)>0){
          matrix_input = matrix_input[,-index_non_sign]
          cat( 'Number of col with all non sign tfs', length(index_non_sign),'\n')
        }
        
        
        ### NA and NON_SIGN
  
        ########### removing rows and columns with all non significant tfs 
        # removing non sign tf in all tfs per row
        
        maxRow = apply(abs(matrix_input),1,'max')
        minRow = apply(abs(matrix_input),1,'min')
        index_non_sign_NA = intersect(which(maxRow==na_value), which(minRow==non_sign_value))
        if (length(index_non_sign_NA)>0){
          cat( 'Number of row with all non sign tfs in reduced matrix', length(index_non_sign_NA),'\n')
          matrix_input = matrix_input[-index_non_sign_NA,]
        }
        
        ### COL
        maxCol= apply(abs(matrix_input),2,'max')
        minCol = apply(abs(matrix_input),2,'min')
        index_non_sign_NA = intersect(which(maxCol==na_value), which(minCol==non_sign_value))
        if (length(index_non_sign_NA)>0){
          cat( 'Number of col with all non sign tfs in reduced matrix', length(index_non_sign_NA),'\n')
          matrix_input = matrix_input[,-index_non_sign_NA]
        }
        
		# COL
		all_the_same = unlist(lapply(colnames(matrix_input), function(colNAME) length(as.numeric(unique(matrix_input[,colNAME])))))
		min_value = unlist(lapply(colnames(matrix_input), function(colNAME) min(as.numeric(unique(matrix_input[,colNAME])))))
		index_all_the_same = which(all_the_same==1 & min_value!=SIGN_VALUE)
		 if (length(index_all_the_same)>0){
          cat( 'Number of col with all the same value', length(index_all_the_same),'\n')
          matrix_input = matrix_input[,-index_all_the_same]
        }
		
		
		# ROW
		all_the_same = unlist(lapply(rownames(matrix_input), function(rowNAME) length(as.numeric(unique(matrix_input[rowNAME,])))))
		min_value = unlist(lapply(rownames(matrix_input), function(rowNAME) min(as.numeric(unique(matrix_input[rowNAME,])))))
		index_all_the_same = which(all_the_same==1 & min_value!=SIGN_VALUE)
		 if (length(index_all_the_same)>0){
          cat( 'Number of row with all the same value', length(index_all_the_same),'\n')
          matrix_input = matrix_input[-index_all_the_same,]
        }
		
        return(matrix_input)
  
}




for ( filename in matrices_list){
    	print(filename)
    	filename_split = strsplit(filename, '.txt')[[1]]
    	matrice <- read.csv(paste0(input_folder_matrices, filename), sep = "\t")
    	
    
    
    	if ( type_range=='sign'){
    		if (grepl('MeanRank', input_folder_matrices)){
    			matrice[is.na(matrice)] <- 2
    			matrice[matrice == 1] <- 0.5
    			matrice[ matrice==0 ] <- 1
    		} else{
    			
    			# in the case type_range wont' be sign, a Binary matrix is created.
    			matrice[is.na(matrice)] <- 2
    			matrice[matrice<=th_p_val] <- 0.5
    			matrice[ (matrice> th_p_val) & (matrice!=0.5) & (matrice< 2) ] <- 1
    	} 
    	}
    	
    	matrice_full = removeRowCol(matrice)
    	listamatrici[[filename_split]] = matrice_full

    	colnames(matrice_full) = gsub('_HUMAN', '' , colnames(matrice_full))
  
    	titolo <- paste0("TFs ",gsub('matrice', '', filename_split))
    	titolo = gsub('_',' ',  titolo)
    	print(titolo)
    	threshold=0.05
    	
    	output_folder = input_folder_matrices 
  
    	if ( length(grep('Chea2', output_folder))>=1){
    		new_col = strsplit(colnames(matrice_full), '[.]')
    	} else {
    		new_col = strsplit(colnames(matrice_full), '_')}
    
    	new_cols = unlist(lapply(new_col, function(item) item[[1]][1]))
    	
    	colnames(matrice_full) = new_cols
    
    	filename_split = paste0(filename_split, '_', type_range)
    		listadendro[[filename_split]]<-heatmap_create(as.matrix(matrice_full), output_folder,filename_split,titolo, 2500,1800, 300, TRUE, type_range,threshold=th_p_val)
    		
    #}
    
    	######
    	print('Gene specific reduction')
    	
    	### Heatmap with only the TF chosen for the study. For example, Hematopoietic TFs or Apoptotic TFs.
       	
    	TFutili <- intersect(Gene_name,colnames(matrice))
    	
    	filename_split_ridotto = paste0(filename_split, '_',gene_type)
    	
    	titolo_ridotto <- paste0("TFs ",gsub('matrice', '', filename_split_ridotto))
    	titolo_ridotto = gsub('_', ' ',  titolo_ridotto)
    	matrice_ridotta = matrice[,which(colnames(matrice) %in% TFutili )]
    
    	matrice_ridotta = removeRowCol(matrice_ridotta,SIGN_VALUE= 0.5, na_value=2, non_sign_value=1)
    	
    	print(dim(matrice_ridotta))
    	listamatrici_ridotte[[paste0(filename_split,'_',gene_type)]] = matrice_ridotta 
    	
    		listadendro_ridotta[[filename_split]]<-heatmap_create(as.matrix(matrice_ridotta), output_folder,filename_split_ridotto,titolo_ridotto, 2500,1800, 300, TRUE,type_range, threshold=th_p_val)
		
}

source('./script/MergedMatrices.R')
final_matrix_double_list = list()
listadendro_double_ridotta = list()
for ( ngene in names(matrices_Ngene )){
	print(ngene)
	list_files = matrices_Ngene[[ngene]]
	
	for ( type_tf in names(list_files)){ 
	  if ( length(list_files[[type_tf]])>0){
  		print(type_tf)
  		th_p_val = 0.05
  		list_files_tfs = matrices_Ngene[[ngene]][[ type_tf]]
  		matrixup = list_files_tfs[grep('up', list_files_tfs)]
  		matrixdown = list_files_tfs[grep('down',list_files_tfs)]
  		percorsomatriceup <- paste0(input_folder_matrices ,matrixup)
  		percorsomatricedown <- paste0(input_folder_matrices , matrixdown)
  	
  		input_matrixup <- read.csv(percorsomatriceup,sep='\t')
  		input_matrixdown <- read.csv(percorsomatricedown,sep='\t')
  
  		if (grepl('MeanRank', input_folder_matrices)){
				# UP
      		  input_matrixup[is.na(input_matrixup)] <- 2
      		  input_matrixup[input_matrixup == 1] <- 0.5
      		  input_matrixup[ input_matrixup==0 ] <- 1
        		  # DOWN
      		  input_matrixdown[is.na(input_matrixdown)] <- 2
      		  input_matrixdown[input_matrixdown == 1] <- 0.5
      		  input_matrixdown[ input_matrixdown==0 ] <- 1
  		} else{
			# UP
    		input_matrixup[is.na(input_matrixup)] <- 2
    		input_matrixup[input_matrixup<th_p_val] <- 0.5
    		input_matrixup[(input_matrixup>=th_p_val) & (input_matrixup != 0.5) & (input_matrixup < 2)] <- 1
		# DOWN
    		input_matrixdown[is.na(input_matrixdown)] <- 2
    		input_matrixdown[input_matrixdown<th_p_val] <- 0.5
    		input_matrixdown[(input_matrixdown>=th_p_val) & (input_matrixdown != 0.5) & (input_matrixdown < 2)] <- 1
      }
  		
  
  		### up
  		if ( length(grep('Chea2', input_folder_matrices))>=1){
  				new_col = strsplit(colnames(input_matrixup ), '[.]')
  		} else {
  				new_col = strsplit(colnames(input_matrixup ), '_')
  				}
  
  		new_cols_up = unlist(lapply(new_col, function(item) item[[1]][1]))
  		colnames(input_matrixup) = new_cols_up
  
  		###### down
  
  		if ( length(grep('Chea2', input_folder_matrices))>=1){
  				new_col = strsplit(colnames(input_matrixdown ), '[.]')
  		} else {
  				new_col = strsplit(colnames(input_matrixdown ), '_')
  				}
  
  		new_cols_down = unlist(lapply(new_col, function(item) item[[1]][1]))
  		colnames(input_matrixdown ) = new_cols_down
  
  
  		# double matrix ridotta
  		new_cols = new_cols_up
  		TFutili <- intersect(Gene_name,new_cols)
  
  		# up
  		input_matrixup_ridotta = input_matrixup[,which(new_cols %in% TFutili )]		
  		# down
  		input_matrixdown_ridotta = input_matrixdown[,which(new_cols %in% TFutili )]		
  		
  		final_matrix_double_ridotta = finalmatrixcreate(input_matrixup_ridotta,input_matrixdown_ridotta )
  		
  		final_matrix_double_ridotta = removeRowCol(final_matrix_double_ridotta, SIGN_VALUE = 2.5, na_value=4, non_sign_value = 3)

  		filename_split_ridotto = paste0('final_matrix_ridotta', ngene,  '_', type_tf, column_sign)
  		titolo_ridotto= paste('TFs', ngene, type_tf)
  		#heatmap 
  		listadendro_double_ridotta[[paste0('final_matrix_', ngene,  '_', type_tf)]] <- heatmap_create(as.matrix(final_matrix_double_ridotta), output_folder,filename_split_ridotto,titolo_ridotto, 2500,1800, 300, TRUE,type_range, threshold=th_p_val)
  		
  		#assign(paste0('final_matrix_', ngene,  '_', type_tf), final_matrix_double )
  		final_matrix_double_list[[paste0('final_matrix_', ngene,  '_', type_tf)]] = final_matrix_double_ridotta
  		
  	}
	}
}

saveRDS(final_matrix_double_list, paste0(input_folder_matrices,'DoubleMatrices_', column_sign, '.rds'))


###### Here I plotted the Dendrograms of each Heatmap, I thought it could be useful to understand similarities between drugs.
library(dendextend)

listagruppi = list()
listagruppi_ridotta = list()
output_folder = input_folder_matrices 

for(dendro_nome in names(listadendro)){
	print(dendro_nome)
	dendrogramma <- listadendro[[dendro_nome]]
	listagruppi[[dendro_nome]]<-dendro_create(dendrogramma,output_folder, dendro_nome,dendro_nome,10,10, 3500,1500,300,TRUE)
	titolo_ridotto<- paste0("TFs ",dendro_nome, '_',gene_type, '_', column_sign)
	#print(titolo_ridotto)
	dendrogramma_ridotto <- listadendro_ridotta[[dendro_nome]]
	tf_k = ncol(listamatrici_ridotte[[paste0(dendro_nome,'_',gene_type)]]) -2
	if (tf_k <=1){
		tf_k  = 2	
	}
	listagruppi_ridotta[[dendro_nome]]<-dendro_create(dendrogramma_ridotto,output_folder, titolo_ridotto,titolo_ridotto,10,tf_k , 25,15,500,FALSE)
}

#Here I again plotted the Dendrograms but for the Heatmaps limited to TF subsets.

listagruppi_ridotta = list()
output_folder = input_folder_matrices 

for(dendro_nome in names(listadendro_double_ridotta)){
	print(dendro_nome)
	dendrogramma_ridotta = listadendro_double_ridotta[[dendro_nome]]
	titolo_ridotto<- paste0("TFs ",dendro_nome, '_',gene_type,'_', column_sign)
	#print(titolo_ridotto)
	tf_k = ncol(final_matrix_double_list[[dendro_nome]]) -2
	if (tf_k <=1){
		tf_k  = 2	
	}
	listagruppi_ridotta[[dendro_nome]]<-dendro_create(dendrogramma_ridotta,output_folder,titolo_ridotto, titolo_ridotto,10,tf_k , 25,15,500,FALSE)
	
	
}



