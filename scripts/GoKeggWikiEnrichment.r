#######Script that computes the Enrichment Analysis of the compounds for the Go,Kegg and WikiPathways Terms instead of the TF.
####I thought this kind of enrichment could help me understand better the activity of each drug.


if (!require("enrichR")){

  install.packages("enrichR")
  
}

library(enrichR) 

makeMatrixEnrichr <- function(list_items, db, ApoList, HemaList){
	
	all_terms_db = unique(unlist(lapply(list_items, function(item) as.character(item[[db]][,'Term']))))
	drugs = names(list_items)
	EnrichmentTable <- matrix(nrow = length(all_terms_db),  ncol = length(drugs), dimnames = list(all_terms_db ,drugs))
	for (drug in drugs){
		results_drug = list_items[[drug]][[db]]
		rownames(results_drug) = results_drug$Term
	 	for(term in all_terms_db){
			
			EnrichmentTable[term,drug]<- results_drug[term,8]
		}
	 }
		
	#### Filter table
	# ALL
	EnrichmentTable_NoNAcols <- EnrichmentTable[,colSums(is.na(EnrichmentTable)) != nrow(EnrichmentTable) ,drop=FALSE ]
	EnrichmentTable_NoNArows <- EnrichmentTable_NoNAcols[rowSums(is.na(EnrichmentTable_NoNAcols)) != ncol(EnrichmentTable_NoNAcols),,drop=FALSE  ]
	#Filters out rows and columns that contain only NAs
	write.table(EnrichmentTable, file=paste0('./output/',typedata,'/EnrichmentTables/',db,"_Enrichment_global.tab"), quote=F, sep="\t", row.names=TRUE, col.names=NA)

	# ApoList
	EnrichmentTableApo =  EnrichmentTable[which(rownames(EnrichmentTable) %in% ApoList),]
	write.table(EnrichmentTableApo, file=paste0('./output/',typedata,'/EnrichmentTables/',db,"_APO_Enrichment_total.tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA)
	#Filters out rows and columns that contain only NAs
	EnrichmentTableApo_NoNAcols <- EnrichmentTableApo[,colSums(is.na(EnrichmentTableApo)) != nrow(EnrichmentTableApo),drop=FALSE  ]
	EnrichmentTableApo_NoNArows <- EnrichmentTableApo_NoNAcols[rowSums(is.na(EnrichmentTableApo_NoNAcols)) != ncol(EnrichmentTableApo_NoNAcols) ,,drop=FALSE ]
	write.table(EnrichmentTableApo_NoNArows, file=paste0('./output/',typedata,'/EnrichmentTables/',db,"_APO_Enrichment.tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA)

	# HemaList
	EnrichmentTableHema =  EnrichmentTable[which(rownames(EnrichmentTable) %in% HemaList),]
	write.table(EnrichmentTableHema, file=paste0('./output/',typedata,'/EnrichmentTables/',db,"_Hema_Enrichment_total.tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA)
	#Filters out rows and columns that contain only NAs
	EnrichmentTableHema_NoNAcols <- EnrichmentTableHema[,colSums(is.na(EnrichmentTableHema)) != nrow(EnrichmentTableHema),drop=FALSE ]
	EnrichmentTableHema_NoNArows <- EnrichmentTableHema_NoNAcols[rowSums(is.na(EnrichmentTableHema_NoNAcols)) != ncol(EnrichmentTableHema_NoNAcols) ,,drop=FALSE ]
	write.table(EnrichmentTableHema_NoNArows, file=paste0('./output/',typedata,'/EnrichmentTables/',db,"_Hema_Enrichment.tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA)
	
	return(list('all'=EnrichmentTable, 'Apo'=EnrichmentTableApo_NoNArows, 'Hema'= EnrichmentTableHema_NoNArows))
	
}


split_matrix <- function(list_results_db , db,typedata){
	list_up_matrix = list()
	list_down_matrix = list()
	for (type_subset in names(list_results_db)){
		print(type_subset)
		data_subset = list_results_db[[type_subset]]
		# up and down separation
		data_subset_up = data_subset[,grep('UP',colnames(data_subset)),drop=FALSE]
		data_subset_down = data_subset[,grep('DOWN',colnames(data_subset)),drop=FALSE]
		# split by cell type
		cells = c('HL60', 'PC3','MCF7')
		for(cell_type in cells){
			print(cell_type)
			data_subset_up_cell = data.frame()
			data_subset_down_cell =  data.frame()
			ind_up_cell = grep(cell_type,colnames(data_subset_up))
			if(length(ind_up_cell)>0){
				data_subset_up_cell = data_subset_up[,ind_up_cell,drop=FALSE]
			}
			ind_down_cell = grep(cell_type,colnames(data_subset_down))
			if(length(ind_down_cell)>0){
				data_subset_down_cell = data_subset_down[,ind_down_cell,drop=FALSE]
			}
			list_up_matrix[[paste0(cell_type, '_',type_subset)]]= data_subset_up_cell
			list_down_matrix[[paste0(cell_type, '_',type_subset)]]= data_subset_down_cell
			if(nrow(data_subset_up_cell)>0){
				write.table(data_subset_up_cell,file=paste0('./output/',typedata,'/EnrichmentTables/',db,"_",type_subset,"_",cell_type,"_UP_Enrichment.tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA )
			}
			
			if(nrow(data_subset_down_cell)>0){
				write.table(data_subset_down_cell,file=paste0('./output/',typedata,'/EnrichmentTables/',db,"_",type_subset,"_",cell_type,"_DOWN_Enrichment.tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA) 
			}
		}
	
	}

	return(list('UP'=list_up_matrix, 'DOWN'=list_down_matrix))
}

args = commandArgs(trailingOnly=TRUE)  
if (length(args) < 1) {
  typedata = "CellSpecific"
  enricherlib <- "ALL"
}else if(length(args) == 1) {
 typedata = args[1]
 enricherlib <- "ALL"
}else if(length(args) >= 2){
 typedata = args[1]
 enricherlib <- args[2]
}


filedata <- "data/WikiPathwaysKeggGoTermsList.txt"
TermLists<-(as.list(read.delim(filedata,sep = "\n",header= F,stringsAsFactors = F)))[[1]]

TermIndexes<-grep("%",TermLists) 
TermFamilies<- sub(":","",sub("%","",TermLists[TermIndexes]))

WikiApoTerms<- TermLists[(TermIndexes[1]+1):(TermIndexes[2]-1)]
WikiHemaTerms<-TermLists[(TermIndexes[2]+1):(TermIndexes[3]-1)]
KeggApoTerms<-TermLists[(TermIndexes[3]+1):(TermIndexes[4]-1)]
KeggHemaTerms<-TermLists[(TermIndexes[4]+1):(TermIndexes[5]-1)]
GoApoTerms<-TermLists[(TermIndexes[5]+1):(TermIndexes[6]-1)]
GoHemaTerms<-TermLists[(TermIndexes[6]+1):length(TermLists)]

 

liblist <- c("WikiPathways_2019_Human","KEGG_2019_Human","GO_Biological_Process_2018")

if(enricherlib%in%liblist){
liblist <-enricherlib
}


fileslist <- list.files(path = c(paste0("output/",typedata,"/DOWN_DEG"),paste0("output/",typedata,"/UP_DEG")))

#I need to convert their names with little trick
drugs <-(sub("(DEG_geni)_up.txt","UP",fileslist))
drugs <-(sub("(DEG_geni)_down.txt","DOWN",drugs))

		
if(FALSE){
if(enrichrlib == "WikiApo"){
	libraryterms <- WikiApoTerms
	enrichrlibFull<-"WikiPathways_2019_Human"
}else if(enrichrlib == "WikiHema"){
	libraryterms<-WikiHemaTerms
	enrichrlibFull<-"WikiPathways_2019_Human"
}else if(enrichrlib == "KeggApo"){
	libraryterms<-KeggApoTerms
	enrichrlibFull<-"KEGG_2019_Human"
}else if(enrichrlib == "KeggHema"){
	libraryterms<-KeggHemaTerms
	enrichrlibFull<-"KEGG_2019_Human"
}else if(enrichrlib == "GoApo"){
	libraryterms<-GoApoTerms
	enrichrlibFull<-"GO_Biological_Process_2018"
}else if(enrichrlib == "GoHema"){
	libraryterms<-GoHemaTerms
	enrichrlibFull<-"GO_Biological_Process_2018"
}else{	
	enrichrlibFull <- enrichrlib }

}


drug_enrichR = list()

for(drug in drugs){
	  
	  cat("The drug is: ", drug,'\n')
	  UPDOWN <-grepl("UP",drug)
	  if(UPDOWN){
	  	genefile <-(sub("UP","DEG_geni_up.txt",drug))
	  	genefile <- paste0("output/",typedata,"/UP_DEG/",genefile)
	  }else{
	  	genefile <-(sub("DOWN","DEG_geni_down.txt",drug)) 
	  	genefile <- paste0("output/",typedata,"/DOWN_DEG/",genefile)
	  } 
	  
	  degenes=read.csv(genefile,stringsAsFactors = FALSE, sep='\t', check.names=FALSE)
	 
	  
	  if(UPDOWN == TRUE){
		maxvalue = degenes[1,2]
		degenes[,2] <- degenes[,2] / maxvalue
	  }else{
		maxvalue <- min(degenes[,2])
		degenes[,2] <- degenes[,2] / maxvalue
	  }
				
			results <- try(enrichr(degenes, databases=liblist))
			  if(class(results) == "try-error") { next }   
			
			drug_enrichR[[drug]] = results		
							
			  
			Sys.sleep(3)
		  
	saveRDS(drug_enrichR, paste0('./output/', typedata, '/EnrichmentTables/EnrichR_List_results.rds'))
}






ListLibFiltered = list('GOHema'=GoHemaTerms, 'GOApo'=GoApoTerms, 'KeggHema'=KeggHemaTerms, 'KeggApo'=KeggApoTerms, 'WikiHema'=WikiHemaTerms, 'WikiApo'=WikiApoTerms)


kegg = makeMatrixEnrichr(drug_enrichR, liblist[2], KeggApoTerms, KeggHemaTerms)
wiki = makeMatrixEnrichr(drug_enrichR, liblist[1], WikiApoTerms, WikiHemaTerms)
go = makeMatrixEnrichr(drug_enrichR, liblist[3], GoApoTerms, GoHemaTerms)


kegg_matrices<-split_matrix(kegg , liblist[2],typedata)
wiki_matrices<-split_matrix(wiki , liblist[1],typedata)
go_matrices<-split_matrix(go , liblist[3],typedata)









