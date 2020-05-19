#I created this function to help me visualize the Matrix created just for the DEG genes of Each drug.
#In this case, I chose to visualize the level of expression (not the enrichment) of each TF in the DEG data from the ConnectivityMap
#the argument pvalue refers to the pvalue of the differential expressed gene considered as significant.

BuildDegMatrix <- function(p_value,drugfile,TFs,pathfiles){
pathbase= pathfiles
filestoread<- list.files(path = pathbase)

if(p_value == 0.05){
treatments<- drugfile[which(drugfile[,2]=='SI'),1]
}else{
treatments<- drugfile[which(drugfile[,3]=='SI'),1]
}

namesfiles<-sub("_HL60","",filestoread)
namesfiles<-sub(".txt","",namesfiles)
finaltreatments<- intersect(namesfiles,treatments)


tabletreatments <- matrix(ncol=length(TFs),nrow=length(finaltreatments),dimnames=list(finaltreatments,TFs))
for(currentfile in finaltreatments){
	
	    
		completefile<- read.csv(paste0(pathbase,currentfile,"_HL60.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t", encoding="UTF-8" )
	    TFpresenti<-intersect(row.names(completefile),TFs)
		print(TFpresenti)
		if(length(completefile[TFpresenti,]>0)){
		for(TF in row.names(completefile[TFpresenti,,drop=FALSE])){
		print(TF)
		tabletreatments[currentfile,TF]<-completefile[TF,]
		}
		
		
	}
}

tabletreatments<-tabletreatments[rowSums(is.na(tabletreatments)) != ncol(tabletreatments),]
return(tabletreatments)
}
