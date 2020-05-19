if(FALSE){
#### questo script contiene le funzioni che creano heatmap e dendrogrammi
if(!require("dendextend")){
 install.packages("dendextend",dependencies = TRUE) 

library(dendextend)
}  
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
}

##### DICHIARO FUNZIONI UTILI
isEmpty <- function(x) {
  return(length(x)==0)
}
#funzione per problema heatmap con NA
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1
  return(edist)
}



library(dendextend)
library(gplots)
library(RColorBrewer)
heatmap_create <- function(matx, output_folder, filetitle,main_title,hgt,wdth,reso,flagpng,type_range ,threshold=0.5) {
  
	output_folder_heatmap = paste0( output_folder, 'heatmaps/')
	dir.create(output_folder_heatmap)
  if(flagpng == TRUE){
    png(file = paste0(output_folder_heatmap,filetitle,".png"),width = wdth,height=hgt, res = reso)
  } else{
    pdf(file = paste0(output_folder_heatmap,filetitle,".pdf"),width = wdth,height=hgt)
  }
  #my_palette <- c('gray51', colorRampPalette(c("red", "black", "green"))(n = 200) )
  #breaks=c(seq(min(-log2(matx)), max(-log2(matx)),length=150))
  
  #my_palette2=c(rep('gray51', length(which(breaks<threshold))),
   #        colorRampPalette(c("red", "black", "green"))(n= length(which(breaks>=threshold))-1)) 
  par(cex.main=0.8)
  #breaks = c(-0.5001,-0.5,0.5,1)
  
	if (type_range=='sign'){
		my_palette2 = c('green', 'red', 'grey51')
		breaks_heatmap = c(0,0.5,1,2)
		data = 	as.matrix(matx)
	} else {
		breaks_heatmap = c(0,0.5,1,2)
		matrice[is.na(matrice)] <- -2
		data = 	as.matrix(-log2(matx))
	}
	
	min_data = min(data)
	if (min_data < 0){
		# up sign 0.5  down ( 0.5 sign, 1 not sign , 2 NA)  
		# 1 , 1.5 , 2.5 , 3
		# down sing 0.5 up (0.5 sing , 1 not sign, 2 NA)
		# -1.5, -2.5, -3
		# 1 + 1 (non significativi) =  2
		# 1 + 2 = 3 
		# 2 + 1 -3
		# 2 + 2 = 4 
		#my_palette2 = c('down non sign NA ', 'down sing NA ', 'down sign non sign ', 'NON SIGN' ,'green', 'sign up non sign down', 'sign up e NA', 'non sign up e NA', 'grey51')
		#my_palette2 <-c( colorRampPalette(c("darkblue",'blue', 'white', "green", 'orange', 'orange4'))(8) , "grey51")
		my_palette2 <- c('blue', "green", 'red4',"grey51", 'grey21')
		breaks_heatmap = c(-2.5,-1, 1, 2.5, 3, 4)	

	}

	margin_b = 8
	cexCol = 0.7
	# margin condition
	if (length(grep('Chea3' , output_folder) | grep('ARCHS4' , output_folder) )>0){
		margin_b = 7
		cexCol = 0.6
		if ( ncol(data)>20){
			cexCol = 0.3
	}
	} 
	
  heatobject <- heatmap.2(data,
                          #cellnote = mat_data,  # same data set for cell labels
                          main = main_title, # heat map title
                          notecol="black",      # change font color of cell labels to black
                          density.info="none",  # turns off density plot inside color legend
                          trace="none",         # turns off trace lines inside the heat map
                          #hclustfun = hclust(method = "complete"),
                          cexRow = 0.3,
		          cexCol= cexCol,
                          #breaks=c(seq(min(-log2(matx)), max(-log2(matx)),length=150)), 
                          breaks = breaks_heatmap,
                          #margins =c(12,9),     # widens margins around plot
                          col=my_palette2,       # use on color palette defined earlier
                          #breaks=colors,    # enable color transition at specified limits
                          dendrogram="both",     # only draw a row dendrogram
                          key = TRUE,
                          #breaks=colors,
                         margins= c(margin_b,8),
                          distfun=dist_no_na,
                         # na.color="#9E0142",
                          Colv="Rowv"
                          # symm=F,symkey=F,symbreaks=T, scale="none",
                       
                          #distfun=function(x) dist(x, method="euclidean"), 
                          #hclustfun=function(x) hclust(x, method="ward.D2")
  )            # turn off column clustering
  
  dev.off()
  
  
  
  
  
  
  return(heatobject)
  
}


dendro_create<-function(heatobject,output_folder, filetitle,title,gruppidrug,gruppiTF,wdth,hgt,reso,flagpng,onlycut = FALSE ){
  ########## SALVO I GRUPPI CREATI IN UN FILE APPOSITO
  dendrorows <-heatobject$rowDendrogram
  dendrocolumns <- heatobject$colDendrogram
  
  groupsrows <- as.data.frame(cutree(dendrorows,k = gruppidrug) )
  groupscolumns <- as.data.frame(cutree(dendrocolumns,k = gruppiTF))
  
  colnames(groupsrows) <- "Cluster"
  colnames(groupscolumns) <- "Cluster"
  output_folder_dendro = paste0(output_folder, 'dendrograms/')
dir.create(output_folder_dendro)
  write.table(groupsrows,paste0( output_folder_dendro ,"cutreerowGROUPS_",filetitle, ".tsv"), sep = '\t')
  write.table(groupscolumns,paste0( output_folder_dendro , "cutreecolumnGROUPS_",filetitle,  ".tsv"),sep = '\t')
  ####ho creato questo controllo dato che la funzione non riesce a colorare pi? di 10 gruppi, 
  #### quindi ho reso questa parte di visualizzazione "opzionale"
  if(!onlycut){
    
  dendrocolumns <-color_branches(dendrocolumns,k = gruppiTF) 
  dendrorows <-color_branches(dendrorows,k = gruppidrug)
  
	value_cex_col = 1
	if (nrow(groupsrows)>100){
	value_cex_col = 0.1
	}

    if(flagpng == TRUE){
      png(file = paste0( output_folder_dendro ,"rows_",filetitle,".png"),width = wdth,height=hgt, res = reso)
    } else{
      pdf(file = paste0(output_folder_dendro, "rows_",filetitle,".pdf"),width = wdth)
    }
	par(mar=c(8,5,2,2))
    dendrorows %>% set("labels_cex", value_cex_col) %>% plot(main = title)
   
    dev.off()
	
	value_cex_row = 1
	if (nrow(groupscolumns)>20){
	value_cex_row = 0.3
	}


    if(flagpng == TRUE){
      png(file = paste0(output_folder_dendro, "columns_",filetitle,".png"),width = wdth,height=hgt, res = reso)
    } else{
      pdf(file = paste0(output_folder_dendro, "columns_",filetitle,".pdf"),width = wdth, height = hgt)
    }
	par(mar=c(8,5,2,2))
	 dendrocolumns %>% set("labels_cex",value_cex_row ) %>% plot(main = title)
    #plot(dendrocolumns,main=title)
    dev.off()
  }
  return(c(groupsrows,groupscolumns))
  
}



if(FALSE){
#### RICORDA CHE IL PDF NECESSITA SOLO DELLA WDTH, intorno a 14
####provo la funzione
dendro_create(dendro,"dendrogramfirst100down",wdth = 14,flagpng = FALSE)
heatmap_create(matriceridottaDOWN100,"matriceridottaUP100prova2",'UpRegulatedTfactors', 1500,1000, 300, TRUE)
###########PROVE
dendro<-heatmap_create(matriceridottaDOWN100,"matriceridottaUP100prova2",'UpRegulatedTfactors', 1500,1000, 300, FALSE)
dendrogramma <- dendro$rowDendrogram
albero <- as.hclust(cutree(dendrogramma, k = 30))
plot(dendrogramma)### DEVI FARE COS?
dendrogramma2<- dendro$colDendrogram
plot(dendrogramma2)
dev.off()

###ricorda di usare sempre as.hclust! 

###varie prove
dendrocolonne<- as.hclust(dendro$colDendrogram)
groups=cutree(dendrocolonne, k=4)
x<-cbind(x,groups)
x1<- subset(x, groups==1)
x2<- subset(x, groups==2)
x3<- subset(x, groups==3)
x4<- subset(x, groups==4)
a =cutree(dendrocolonne, k=3)
plot(a)
dev.off()
x<-abline(h=c(30,40)) ###metodo per disegnare una riga
vettore <- c(0,0,0) 
for(TF in TFutili){
  vettore[a[TF]] = vettore[a[TF]] + 1
}
}    
