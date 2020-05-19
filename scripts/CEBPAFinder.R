#I have created this script to easily visualize the drugs that enriched the TF CEBPA, that we think could be more promising for our screening

matrice<-readRDS("Leukemia_NewRunDoubleMatrices_SignMean_0.05.rds")[[1]]
compoundsrows<-matrice[which(matrice[,"CEBPA"]>=0 & matrice[,"CEBPA"]<2.6),]

drugs<-as.data.frame(matrice[which(matrice[,"CEBPA"]>=0 & matrice[,"CEBPA"]<2.6),"CEBPA"])
drugs$valori<- drugs[,1]
drugs[,1]<- rownames(drugs)

colnames(drugs)<- c("composti","valori")
rownames(drugs) <- NULL
lineaperilfile <- (drugs)

filesalvataggio <- ("data/Subsets/ElencoDroghe.tab")

write.table(lineaperilfile,file=filesalvataggio, sep="\t", row.names=FALSE)
