###This script has been created to produce Histograms and Density Maps of the results of the Enrichment for the TF.
###I thought it would help me understand how strong the signal of enrichment was.

library(data.table)
typedata <- "CellSpecific"
pathbase <- paste0("output/",typedata,"/EnrichmentTables/")
fileslist <- list.files(path = c(paste0("output/",typedata,"/EnrichmentTables/")))
fileslist <- fileslist[grep("Enrichment", fileslist)]

filename <- fileslist[1]

MatrixDrugs <- read.table("data/Subsets/ElencoDroghe.tab",  sep="\t", header = TRUE, stringsAsFactors = FALSE )
SubsetDrugs <- MatrixDrugs[,1]
SubsetCBPA = list()

matricedimensioni <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(matricedimensioni)<- c("Dataset","Righe","Colonne") 
for(filename in fileslist){
	currentTable<- read.table(paste0(pathbase,filename),  sep="\t", header = TRUE, stringsAsFactors = TRUE,row.names = 1,check.names=FALSE )

	indexSubsets<- sub(("_DOWN|_UP"),"",colnames(currentTable)) %in% SubsetDrugs
	subsetCBPA<-currentTable[,indexSubsets]
	graphname<- sub("_Enrichment.tab","",filename)

	write.table(subsetCBPA,file=paste0('./output/', typedata, '/EnrichmentTables/CEBPALimited/',graphname,'_CEBPA_enrichment.tab'), quote=F, sep="\t", row.names=TRUE,col.names=NA) 
	if(TRUE){
	matricecorrente <- as.matrix(currentTable)
	matricecorrente[is.na(matricecorrente)] <- 0
	vectorForGraph <- as.vector(matricecorrente)
	pdf(paste0(pathbase,'/Densities/',graphname,'.pdf'),height=20,width=10)
    
	print(graphname)
	check<-try(density(vectorForGraph))
	#DENSITY
	segnaleErrore <- ""
	dimensioni <- c(0,0)
	if(class(check) == "try-error"){segnaleErrore <- paste0("there has been an error for this one table: ", graphname," ");dev.off();next}#lineaperilfile <- paste0(segnaleErrore,"Le dimensioni della tabella ",graphname," sono: ", dimensioni[1]," ",dimensioni[2],"\n");filesalvataggio <- (paste0(pathbase,"DimensionEnrichmentTables.txt"));write(lineaperilfile,file=filesalvataggio,append=TRUE);	dev.off();next}
	else{dimensioni <- paste(dim(currentTable))
	plot(check,xlab='Enrichr Score', main = graphname)
	}
		
	dev.off()

	rigadataframe <- data.frame("Dataset" = graphname,"Righe" = dimensioni[1],"Colonne"=dimensioni[2])
	matricedimensioni <- rbind(matricedimensioni,rigadataframe)
	}

}

write.table(matricedimensioni,file=paste0('./output/', typedata, '/EnrichmentTables/enrichmentDimensioni.tab'), quote=F, sep="\t", row.names=TRUE,col.names=NA) 
