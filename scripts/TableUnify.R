###############I've created this script to unify data in a file coming from the Prestwick Library and PubChemFileInfo

library(plyr)
library(dplyr)
library(gtools)
pathbase = "./data/DrugInformation/"

Prestwicktable<-read.csv(paste0(pathbase,"PrestwickModified.tab"),header = TRUE, stringsAsFactors = FALSE, sep = "\t", encoding="UTF-8", )

PubchemInfo<-read.csv(paste0(pathbase,"PubChemFileInfo.tab"),header = TRUE, stringsAsFactors = FALSE, sep = "\t", encoding="UTF-8", )

PubchemInfoReduced<-unique(PubchemInfo[c("CMAP_NAME","catalog_name","CAS")])

mergeddata<- merge(Prestwicktable,PubchemInfoReduced, by.x=c("Chemical.name","CAS.number"),by.y = c("catalog_name","CAS"),all.x = TRUE) #gli NA matchano

mergeddata1<- merge(Prestwicktable,PubchemInfoReduced, by.x="Chemical.name",by.y ="catalog_name",all.x = TRUE) #gli NA matchano
mergeddata2<- merge(Prestwicktable,PubchemInfoReduced, by.x="CAS.number",by.y = "CAS",all.x= TRUE) #gli NA matchano
names(mergeddata2)[1]<-"CAS"
names(mergeddata2)[2]<-"Chemical-name"
totalmergeddata<-smartbind(mergeddata1,mergeddata2)

x<-totalmergeddata

finalmatrix<-x[!is.na(x$Chemical.name),]

write.table(mergeddata2[c("Chemical-name", "CMAP_NAME","catalog_name","CAS","Prestw.number","Plate.Nb....Position.Nb","fmla.structure","Mol.weight..structure","Precautions","Therapeutic.group","Mechanism.of.action", "Side.Effect.s.")], paste0(pathbase,"MergedTables/","MergeonCASTable.tab" ), sep="\t",row.names = FALSE,col.names = TRUE)
write.table(mergeddata1, paste0(pathbase,"MergedTables/","MergeonChemicalnametable.tab" ), sep="\t",row.names = FALSE,col.names = TRUE)
write.table(finalmatrix[c("Chemical.name", "CMAP_NAME","CAS","Prestw.number","Plate.Nb....Position.Nb","fmla.structure","Mol.weight..structure","Precautions","Therapeutic.group","Mechanism.of.action", "Side.Effect.s.")], paste0(pathbase,"MergedTables/","TrueFinalTable.tab" ), sep="\t",row.names = FALSE,col.names = TRUE)

if(FALSE){
PrestReduced <- Prestwicktable[c("Chemical.name","CAS.number","Prestw.number")]

Prestwicklist <-dlply(PrestReduced,.(Chemical.name),c)
names(Prestwicklist)<-Prestwicktable['Chemical.name']
PubchemInfolist <-dlply(PubchemInfoReduced,.(CMAP_NAME),c)
names(PubchemInfolist)<- PubchemInfoReduced["CMAP_NAME"]
PubchemInfolist <- unique(PubchemInfolist)

########I have to change the colnames in order to make the merge operation possible
colnames(PrestReduced)<- c("Compound.name" ,"CAS"  ,  "Prestw.number")
PrestReduced["Compound.name"]<-tolower(PrestReduced[[1]])
PubchemInfoReduced["CMAP_NAME"]<-tolower(PubchemInfoReduced[[1]])
mergedata<-merge(PrestReduced,PubchemInfoReduced, by.x=c("CAS","Compound.name"), by.y=c("CAS","CMAP_NAME"),all = TRUE) #gli NA matchano
mergedataCAS<-merge(PrestReduced,PubchemInfoReduced, by="CAS",all = TRUE) #gli NA matchano

PubchemInfoReduced2<- unique(PubchemInfo[c("CMAP_NAME","CAS","catalog_name")])

mergeddatacompletePrestwick<- merge(Prestwicktable,PubchemInfoReduced2, by.x="CAS.number",by.y = "CAS",all = TRUE) #gli NA matchano
#####The best type of merge is obtained by using catalog_name and Chemical.name
mergeddatacompletePrestwickfinal<- merge(Prestwicktable,PubchemInfoReduced2, by.x=c("Chemical.name","CAS.number"),by.y = c("catalog_name","CAS"),all = TRUE) #gli NA matchano
mergeddatacompletePrestwickfinal<-mergeddatacompletePrestwickfinal[c("Chemical.name", "CMAP_NAME","CAS.number","Prestw.number","Plate.Nb....Position.Nb","fmla.structure","Mol.weight..structure","Precautions","Therapeutic.group","Mechanism.of.action", "Side.Effect.s.")]
mergeddatacompletePrestwick<-mergeddatacompletePrestwick[c("Chemical.name", "catalog_name","CMAP_NAME","CAS.number","Prestw.number","Plate.Nb....Position.Nb","fmla.structure","Mol.weight..structure","Precautions","Therapeutic.group","Mechanism.of.action", "Side.Effect.s.")]


write.table(mergeddatacompletePrestwickfinal, paste0(pathbase,"MergedTables/","completeMergedTable.tab" ), sep="\t",row.names = FALSE,col.names = TRUE)

write.table(mergedata, paste0(pathbase,"MergedTables/","SimpleMergedTable.tab" ), sep="\t",col.names = TRUE, row.names = FALSE )

}









                                                                                