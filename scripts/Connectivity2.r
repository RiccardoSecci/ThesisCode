#With this script I estract from the initial matrix coming from the Connectivity Map Build02 database the various signatures of each drug.


#! /usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE) #così posso fornire gli argomenti

if (length(args) < 3) {
  stop("devi fornire il numero di geni che vuoi che vengano presi in considerazione e la table gene expression di partenza.n", call.=FALSE)
} 
filename = args[2]
typeData = args[3]
print(filename)
pathfiletable = paste0("data/" , filename)
print(pathfiletable)
table_gene_expression = readRDS(pathfiletable) 
trattamento = colnames(table_gene_expression) 

Ngenes = args[1]
if (Ngenes !='logFC' & Ngenes != 'DEG'){
	Ngenes = as.numeric(args[1]) !
} 

completapercorso <- paste0('output/', typeData,'/' )
print(completapercorso)
dir.create(completapercorso)
nomeFileTrattamenti = paste0(completapercorso, 'listaTrattamenti.txt')
write.table(trattamento,nomeFileTrattamenti, sep = '\t' ,quote = FALSE, row.names = FALSE)
output_up_gene = paste0(completapercorso,'/UP_', Ngenes,'/')
output_down_gene =  paste0(completapercorso,'/DOWN_', Ngenes,'/')
dir.create(output_up_gene)
dir.create(output_down_gene)



if ( Ngenes !='DEG'){
	tab_logfc = data.frame()
	i=0         
	for (nome in trattamento){
		i=i+1
		cat('Trattamento ',i,  nome,'\n')
		## ****** sostituzione spazi nome composto con "_" e / con 'tilde' !!! cercare funzione per sostituire!!!
		
		nome2 <- gsub("/", "~", nome)
		nome1 <- gsub(" ","__",nome2)
		nome_corretto = paste0(nome1,"_",Ngenes)
		   #non funziona nome <- readline(prompt="inserisci il nome del composto: ")
		 #tabella_ordinata <- table_gene_expression[order(-nome) , ]
		  ##### WHICH!!!!!!!!! ***************************
		 indice<-match(nome,trattamento) #the index of the treatment of interest

		 oggetto <- table_gene_expression[indice] 

		
		
		geni<-rownames(oggetto)
	
		espressione <- oggetto[,1]

		genedataframe <- data.frame("geni" = geni, "espressione" = espressione) 
		
	
		
		ordinato <- genedataframe[order(-genedataframe$espressione),]  #the data is ordered in a descendent way
		
		if ( Ngenes!='logFC' & Ngenes !='DEG'){
			upregolati <- ordinato[1:Ngenes,]
			size = 11990 # dim(table_gene_expression)[1]
		
			downregolati <- tail(ordinato,Ngenes)
		} else if (Ngenes=='logFC') {
			upregolati <- ordinato[which(ordinato$espressione>1),]
			downregolati <- ordinato[which(ordinato$espressione< -1),]
			nup = nrow(upregolati)
			ndown = nrow(downregolati)
			tab_logfc = rbind(tab_logfc, data.frame('Trattamento'=nome, 'N_up_gene'= nup, 'N_dw_gene' =ndown))
			if ( nup> 500){
				upregolati <- ordinato[1:500,]
			}
			if( ndown > 500){
				downregolati <- tail(ordinato,500)
			} 
		}        
	
		
		filename_UP<-paste0(output_up_gene, nome_corretto, "_geni_up.txt") 
		filename_DOWN<-paste0(output_down_gene,nome_corretto, "_geni_down.txt")
		
		write.table(upregolati,filename_UP, sep = '\t' ,quote = FALSE, row.names = FALSE) 
		write.table(downregolati,filename_DOWN, sep = '\t' ,quote = FALSE, row.names = FALSE)
		


	
	 
	}

	if ( nrow(tab_logfc)>=1 ) {
		cat('The total number of treatments is ', length(trattamento), 'the number of compounds with regulated genes are' , nrow(tab_logfc))
		  write.table(tab_logfc,paste0(completapercorso, 'NumberOfRegulatedGenes.txt'), sep = '\t' ,quote = FALSE, row.names = FALSE)
	

	}
}


tab_deg = data.frame()
if (Ngenes=='DEG'){

	if (typeData=='AllCells'){
		path_data = './data/DEG_genes/CMP/'
		list_files = list.files(path_data)

	} else if (typeData== 'Leukemia'){
		path_data = './data/DEG_genes/Leukemia/'
		list_files = list.files(path_data, pattern='.txt')

	} else if (typeData== 'Leukemia_Combat'){
		path_data = './data/DEG_genes/Leukemia_Combat/'
		list_files = list.files(path_data, pattern='.txt')
		
    } else if (typeData== 'Leukemia_Combat_probe'){
		path_data = './data/DEG_genes/Leukemia_Combat_probe/'
		list_files = list.files(path_data, pattern='.txt')

	} else if (typeData== 'CMP_CELL_Combat_probe'){
		path_data = './data/DEG_genes/CMP_CELL_Combat_probe/'
		list_files = list.files(path_data, pattern='.txt')

	} else if (typeData== 'Leukemia_newRun'){
		path_data = './data/DEG_genes/Leukemia_newRun/'
		list_files = list.files(path_data, pattern='.txt')

	} else if (typeData== 'CellSpecific'){
		path_data = './data/DEG_genes/CMP_CELL/'
		list_files = list.files(path_data)
		
	}
	
	for ( filename in list_files){
		nome = gsub('.txt','', filename)
		print(nome)
		oggetto = read.csv(paste0(path_data , filename),sep='\t', check.names=FALSE)
		ordinato = oggetto[order(-oggetto[,1]),,drop=FALSE]
		upregolati <- ordinato[which(ordinato[,1]>0),,drop=FALSE]
		downregolati <- ordinato[which(ordinato[,1]< 0),,drop=FALSE]
		nup = nrow(upregolati)
		ndown = nrow(downregolati)
		if ( nup> 500){
			upregolati <- ordinato[1:500,, drop=FALSE]
		}
		if( ndown > 500){
			downregolati <-  as.data.frame(tail(ordinato,500))
		} 
		
		nome_corretto = paste0(nome,"_",Ngenes)    
		filename_UP<-paste0(output_up_gene, nome_corretto, "_geni_up.txt") 
		filename_DOWN<-paste0(output_down_gene,nome_corretto, "_geni_down.txt")
		tab_deg = rbind(tab_deg, data.frame('Trattamento'=nome, 'N_up_gene'= nup, 'N_dw_gene' =ndown))
		if ( nup>2 & ndown> 2){
			upregolati$geni = rownames(upregolati)
			downregolati$geni = rownames(downregolati)
			df_up = upregolati[,c(2,1)]
			df_down = downregolati[,c(2,1)]
			write.table(df_up,filename_UP, sep = '\t' ,quote = FALSE, row.names=FALSE) 
			write.table(df_down,filename_DOWN, sep = '\t' ,quote = FALSE, row.names=FALSE)
		}
	}
	if ( nrow(tab_deg)>=1 ) {
		cat('The total number of treatments is ', length(trattamento), 'the number of compounds with sign regulated genes are' , nrow(tab_deg))
		  write.table(tab_deg,paste0(completapercorso, 'NumberOfRegulatedGenes_SIGN.txt'), sep = '\t' ,quote = FALSE, row.names = FALSE)
	
	}
}

