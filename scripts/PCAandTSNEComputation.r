###With this script I create the tSNE and the PCAs of each drug deg genes in the Connectivity Map. I have created it to better understand their behaviour.

Table = readRDS('./data/mylogfc_CMP_CELL_subset_cellLines.rds')

library(readr)
library(Rtsne)
library(ggrepel)

newexecution = FALSE # mettere true se vuoi creare delle nuove distribuzioni

CurrentSubset = "_107_drugs_"    # you could also set it as "alldrugs_all_total_CMAP" if you want to create the tsne and the pca for all the compounds
for(CurrentSubset in c("_107_drugs_","alldrugs","alldrugs_all_total_CMAP")){


		if(CurrentSubset == "_107_drugs_"){
		drugs_useful <-  read.delim("data/Drugs107.txt",header = FALSE,stringsAsFactors = FALSE)
		subsetType = "_107_drugs_"
		repel_size = 1.8
		segment_size = 0.5
		segment_alpha = 0.5
		geom_point_size = 1
		dataframename <- "tsneCurrentSubset.tab"
	
		}else if(CurrentSubset ==  "alldrugs"){drugs_useful <- data.frame(Drugs = sub("_HL60_DEG_geni_down.txt","",list.files(path = "./output/Leukemia_newRun/DOWN_DEG/")),stringsAsFactors=FALSE)
		subsetType = "alldrugs"
		repel_size = 1.4
		segment_size = 0.2
		segment_alpha = 0.2
		geom_point_size = 0.3
		dataframename <- "tsneAllDrugs.tab"
		}else{
		drugs_useful <- read.delim("data/CommonDrugs_all_cell_lines_total_CMAP.txt",header = FALSE,stringsAsFactors = FALSE)
		subsetType = "alldrugs_all_total_CMAP"
		repel_size = 1.4
		segment_size = 0.2
		segment_alpha = 0.2
		geom_point_size = 0.3
		dataframename <- "tsneAll_total_CMAP.tab"
		}
		
		
	

	drugs_useful[,1] <- sort(drugs_useful[,1]) 

	write.table(drugs_useful,file=paste0("output/TSNE_Final/",subsetType,"index.tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA)

	drugs_useful[,1] = gsub('_', ' ' , drugs_useful[,1])


	cols = colnames(Table) 
	correctdrugs = sub("_(HL60|MCF7|PC3)","",colnames(Table)) %in% c(drugs_useful[,1])
	cols_limited = cols[correctdrugs]
	hl60 = grep('_HL60', cols)
	pc3 = grep('_PC3', cols)
	mcf7 =  grep('_MCF7', cols)



	list_subset_drug = list()



	print(subsetType)


	cols_replace = unique(gsub('(_HL60)|(_PC3)|(_MCF7)','',  cols_limited))
	Ndrug = length(cols_replace)

	drug_subset = cols_replace
	list_subset_drug[[subsetType]] = drug_subset

	cols_drug_subset = sort(cols[unlist(lapply(drug_subset, function(name) grep(paste0("^",gsub("+","\\+",gsub(")","\\)",gsub("(","\\(",name,fixed =TRUE),fixed = TRUE),fixed = TRUE),"_"), cols)))])
	data_t = t(Table[,cols_drug_subset])
	my_x <- data.matrix(data_t)

		
	################################################################	
	#### From here I create the tSNE distributions
	################################################################
	set.seed(60)
	if(newexecution){
	my_t<- Rtsne(my_x,PCA=FALSE)}
	
	TSNE_1 = as.numeric(my_t$Y[,1])
	TSNE_2 = as.numeric(my_t$Y[,2])
    }

	samples = unlist(lapply(cols_drug_subset, function(name) { name_split = strsplit(name,'_')[[1]]; name_split[length(name_split)]}))
	drug_names = unlist(lapply(cols_drug_subset, function(name) { name_split = strsplit(name,'_')[[1]]; paste0(name_split[1:(length(name_split)-1)])}))
	drug_indices = match(drug_names,drugs_useful[,1])
    if(newexecution){
	df_TSNE = data.frame('Tsne1'= TSNE_1,'Tsne2' = TSNE_2 ,'CellLine' = samples, 'Drug' = drug_names,'index' = drug_indices)
    }else{
	df_TSNE <- read.table(file=paste0("output/TSNE/TSNETable",subsetType,".tab"),  sep="\t", header = TRUE)
	}
	colors = rep(c('green4','red4','blue4'),Ndrug)
	shapes = rep(c(15,16,17), Ndrug)

	df_TSNE$Color = colors
	df_TSNE$Shape = shapes


	library(RColorBrewer)
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
	col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
	colors_drug = rep(col_vector[1:74],17)
	colors_drug = colors_drug[1:1019]

	set.seed(42)
	library(ggplot2)
	
	pdf(paste0('./output/TSNE_Final/TSNE_Gennaio_withNumbers_',subsetType,'.pdf'), width= 14, height =10)



	if(CurrentSubset == "_107_drugs_"){



	gg1 = ggplot(df_TSNE, aes(x=Tsne1, y=Tsne2,color=CellLine))+
		geom_point( size=geom_point_size)+
		scale_color_manual(values=colors)+	
		geom_text_repel(aes(label = drug_names),force = 0.4,,
						segment.size = segment_size,
						segment.alpha = segment_alpha ,
						size = repel_size,show_guide = F)+
		 xlab(paste0("tSNE dimension 1"))+
		 ylab(paste0("tSNE dimension 2"))+
		theme(legend.position = "right")

	print(gg1)
    variableCHECK = 1
	dev.off()
	}else{
	gg1 = ggplot(df_TSNE, aes(x=Tsne1, y=Tsne2,color=CellLine))+
		geom_point( size=geom_point_size)+
		
		 guides(colour = guide_legend(override.aes = list(size=2)))+
		scale_color_manual(values=colors)+
		geom_text_repel(aes(label = drug_indices),force = 0.4,,
						segment.size = segment_size,
						segment.alpha = segment_alpha ,
						size = repel_size,show_guide = F)+
		 xlab(paste0("tSNE dimension 1"))+
		 ylab(paste0("tSNE dimension 2"))+
		theme(legend.position = "right")
	print(gg1)
    variableCHECK = 0
	dev.off()


	}



	library(ggplot2)
	pdf(paste0('./output/TSNE_Final/TSNE_Gennaio_',subsetType,'.pdf'), width= 14, height =10 )
	gg1 = ggplot(df_TSNE, aes(x=Tsne1, y=Tsne2,color=CellLine))+
		geom_point( size=geom_point_size)+
		guides(colour = guide_legend(override.aes = list(size=2)))+
		scale_color_manual(values=colors)+
		 xlab(paste0("tSNE dimension 1"))+
		 ylab(paste0("tSNE dimension 2"))+
		theme(legend.position = "right")
	print(gg1)

	dev.off()


	library(ggplot2)
	pdf(paste0('./output/TSNE_Final/TSNE_Gennaio_shape_', subsetType,'.pdf'), width= 14, height =10 )
	gg1 = ggplot(df_TSNE, aes(x=Tsne1, y=Tsne2,shape=CellLine, color=drug_names))+
		geom_point( size=1.8)+
		
		scale_color_manual(values=colors_drug )+
		scale_shape_manual(values=shapes)+ 
		xlab(paste0("tSNE dimension 1"))+
		ylab(paste0("tSNE dimension 2"))+
		guides(color = FALSE) +
		theme(legend.position="right")


	print(gg1)

	dev.off()
	if(newexecution){
	write.table(df_TSNE,file=paste0("output/TSNE_Final/TSNETable",subsetType,".tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA)
	}							
								
								
	###################################################					
	#####Here I calculate the PCA distributions			
	###################################################	
	if(newexecution){
	my_p <- prcomp(my_x, center=TRUE)


	PCA_1 = as.numeric(my_p$x[,1])
	PCA_2 = as.numeric(my_p$x[,2])
    }

	samples = unlist(lapply(cols_drug_subset, function(name) { name_split = strsplit(name,'_')[[1]]; name_split[length(name_split)]}))
	drug_names = unlist(lapply(cols_drug_subset, function(name) { name_split = strsplit(name,'_')[[1]]; paste0(name_split[1:(length(name_split)-1)])}))
	drug_indices = match(drug_names,drugs_useful[,1])
	if(newexecution){
	df_pca = data.frame('PC1'= PCA_1,'PC2' = PCA_2 ,'CellLine' = samples, 'Drug' = drug_names,'index' = drug_indices)
	}else{
	df_pca <- read.table(file=paste0("output/PCA/PCATable",subsetType,".tab"),  sep="\t", header = TRUE)
	}
	colors = rep(c('green4','red4','blue4'),Ndrug)
	shapes = rep(c(15,16,17), Ndrug)

	df_pca$Color = colors
	df_pca$Shape = shapes


	if(newexecution){
	eigvs = my_p$sdev^2
	pca1_s = round(( eigvs[1] / sum(eigvs))*100,3)
	pca2_s = round(( eigvs[2] / sum(eigvs))*100,3)
	
	saveRDS(eigvs, file = paste0("output/PCA/PCA",subsetType,"eigenvalues.rds"))
	}else{
	eigvs <- readRDS(paste0("output/PCA/PCA",subsetType,"eigenvalues.rds"))
	pca1_s = round(( eigvs[1] / sum(eigvs))*100,3)
	pca2_s = round(( eigvs[2] / sum(eigvs))*100,3)
	}

	


	library(RColorBrewer)
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
	col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
	colors_drug = rep(col_vector[1:74],17)
	colors_drug = colors_drug[1:1019]
	library(ggplot2)
	pdf(paste0('./output/PCA_Final/PCA_Gennaio_',subsetType,'.pdf'), width= 14, height =10 )
	gg1 = ggplot(df_pca, aes(x=PC1, y=PC2,color=CellLine))+
		geom_point( size=geom_point_size)+
		
	   guides(colour = guide_legend(override.aes = list(size=2))) +
		scale_color_manual(values=colors)+
		xlab(paste('PC1', '(', pca1_s,'%)'))+
		ylab(paste('PC2', '(', pca2_s,'%)'))+
		theme(legend.position = "right")

	print(gg1)

	dev.off()

	set.seed(42)

	if( CurrentSubset == "_107_drugs_" ){



	pdf(paste0('./output/PCA_Final/PCA_Gennaio_withNumbers_',subsetType,'.pdf'), width= 14, height =10 )
	gg1 = ggplot(df_pca, aes(x=PC1, y=PC2,color=CellLine))+
		geom_point( size=geom_point_size)+
		geom_text_repel(aes(label = drug_names),
						size = repel_size,show_guide = F,force = 0.4,
						segment.size = segment_size,
						segment.alpha = segment_alpha )+
		scale_color_manual(values=colors)+
		xlab(paste('PC1', '(', pca1_s,'%)'))+
		ylab(paste('PC2', '(', pca2_s,'%)'))+
		theme(legend.position = "right")
		

	print(gg1)
    variableCHECK = 1
	print(variableCHECK)
	dev.off()
	}else{
	pdf(paste0('./output/PCA_Final/PCA_Gennaio_withNumbers_',subsetType,'.pdf'), width= 14, height =10 )
	gg1 = ggplot(df_pca, aes(x=PC1, y=PC2,color=CellLine))+
		geom_point( size=geom_point_size)+
		geom_text_repel(aes(label = drug_indices),
						size = repel_size,show_guide = F,force = 0.4,
						segment.size = segment_size,
						segment.alpha = segment_alpha )+
		scale_color_manual(values=colors)+
		xlab(paste('PC1', '(', pca1_s,'%)'))+
		ylab(paste('PC2', '(', pca2_s,'%)'))+
		theme(legend.position = "right")
		

	print(gg1)
    variableCHECK= 0
	print(variableCHECK)
	dev.off()




	}

	library(ggplot2)
	pdf(paste0('./output/PCA_Final/PCA_Gennaio_shape_', subsetType,'.pdf'), width= 14, height =10 )
	gg1 = ggplot(df_pca, aes(x=PC1, y=PC2,shape=CellLine, color=drug_names))+
		geom_point( size=1.8)+
		scale_color_manual(values=colors_drug )+
		scale_shape_manual(values=shapes)+ 
		xlab(paste('PC1', '(', pca1_s,'%)'))+
		ylab(paste('PC2', '(', pca2_s,'%)'))+
		guides(color = FALSE) +
		theme(legend.position="right")
		

	print(gg1)

	dev.off()
	if(newexecution){
	write.table(df_pca,file=paste0("output/PCA_Final/PCATable",subsetType,".tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA)
	}
}



