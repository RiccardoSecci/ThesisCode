
########In this script I calculate the distance between the points belonging to each drug.
#A given drug can have a maximum of three points, in the case the data is available for all the cell lines (HL60, PC3, MCF7)

distanceCreate <- function(database,ordinazionerandomica){

	distanceTable <- data.frame(averageDistance =double(), Drug = character())

	if(ordinazionerandomica == TRUE){
		set.seed(42)
		tipoordinazione = "random"
		database <-  database[sample(nrow(database)),]
	}else{
	tipoordinazione = "ordinato"
	}

	for(drugrow in seq(1,nrow(database),by=3)){

		x1 <- database[drugrow,1]
		x2 <- database[(drugrow + 1),1]
		x3 <- database[(drugrow + 2),1]
		y1 <- database[drugrow,2]
		y2 <- database[(drugrow + 1),2]
		y3 <- database[(drugrow +2),2]
		distance12 <- sqrt(( x1 - x2)^2 + ( y1 - y2)^2)
		distance23 <- sqrt(( x2 - x3)^2 + ( y2 - y3)^2)
		distance13 <- sqrt(( x1 - x3)^2 + ( y1 - y3 )^2) 


		if(ordinazionerandomica){
		droga = paste0("random", database[drugrow,5])
		}else{droga = database[drugrow, 4 ]
		}
		distanza = (distance12 + distance23 + distance13)/3 
		rowdata = data.frame(  averageDistance =distanza, Drug = droga)
		distanceTable = rbind(distanceTable, rowdata)
	}
	
	write.table(distanceTable,file=paste0("output/TSNE/DistanceTable",tipoordinazione,".tab"), quote=F, sep="\t", row.names=TRUE,col.names=NA) 
	return(distanceTable)

}

subsetType <- "_107_drugs_Novembre"  #  I could also use the "alldrugs" parameter to do it an all the compounds.

database = read.table(file=paste0("output/TSNE/TSNETable",subsetType,".tab"),  sep="\t",header = TRUE,row.names = 1 )
orderedDistance <- distanceCreate(database,FALSE)
randomDistance <- distanceCreate(database,TRUE)

#################I Draw the graph representing the distributions of the distance values

pdf("output/TSNE/DistanceOrdered.pdf")
orderdensity <- density(orderedDistance$averageDistance)
plot(orderdensity, main="Density Plot of drug CellLine Distances")
polygon(orderdensity, col="red", border="blue")

dev.off()

pdf("output/TSNE/DistanceRandom.pdf")
unorderdensity <- density(randomDistance$averageDistance)
plot(unorderdensity, main="Density Plot of random CellLine Distances")
polygon(unorderdensity, col="red", border="blue")

dev.off()


pdf("output/TSNE/BothDistributions.pdf")
plot(unorderdensity,main= "BothDistributions",col = 'red')
lines(orderdensity,col = 'blue')

dev.off()

#########################This test permits me to check if the distributions are significantly different.

x <- orderedDistance[,1]
y <- randomDistance[,1]
risultato <- t.test(x, y, paired = TRUE, alternative = "two.sided")










