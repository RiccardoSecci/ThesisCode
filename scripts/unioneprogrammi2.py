#This is the most important script. From here, the TF enriched by the UP/Downregulated genes of each drug are calculated
from ScriptX2K import *
from TF_Tools import *
import argparse
import sys 
import os
import os.path
from os import listdir
from os.path import isfile, join
import re

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--numgeni", help="number of genes")
	parser.add_argument("--tool", help="scrivi Chea2 se vuoi Chea2, Chea3 se vuoi Chea3 ARCHS4")
	parser.add_argument("--GeneData", help="il file con il quale vuoi lavorare. leukemia_matrix.rds se vuoi solo i dati di leukemia, mylogfc_CMP_CELL_subset_cellLines.rds, mylogfc_CMP.rds se vuoi usarli senza le info sulla linea cellulare" )
	parser.add_argument("--TypeData", help='se vuoi i dati riassuntivi senza linea cellulare, scrivi AllCells. Leukemia se interessa soltanto la linea legata alla leucemia, CellSpecific se vuoi gli altri.')
	args = parser.parse_args()

	print("questi sono gli argomenti che hai dato allo script:")
	print(args.numgeni)
	print(args.tool)
	print(args.GeneData)
	print(args.TypeData)
	ngeni = args.numgeni
	TFtool = args.tool 
	pathfiletable = args.GeneData
	typeData = args.TypeData

	completapercorso = 'output/'+ typeData 

	directoryscelta = completapercorso + '/' + TFtool


directoryesistente = completapercorso + '/UP_' + ngeni +'/'
directoryesistenteDOWN = completapercorso + '/DOWN_' + ngeni +'/'


if( not os.path.isdir(directoryesistente)):
	stringacomando = "Rscript  script/Connectivity2.r " + ngeni + " " + pathfiletable + ' ' + typeData # ....nomefile input aggiungerlo dentro connectivity.r 
	os.system(stringacomando) 
	print("ho finito Connectivity2.r")

try:
	os.mkdir(directoryscelta) 	
except Exception as e: 
	pass

onlyfiles = [f for f in listdir(directoryesistente) if isfile(join(directoryesistente, f))]


trattamenti = []
for filename in onlyfiles:
	if typeData=='Leukemia':
		filename = re.sub("_HL60","",filename)
	if ngeni=='logFC':
		trattamenti.append(re.sub("_logFC_geni_(up|down)\.txt","",filename))
	elif ngeni == 'DEG':
		 trattamenti.append(re.sub("_DEG_geni_(up|down)\.txt","",filename))
	else:
		trattamenti.append(re.sub("_[0-9]*_geni_(up|down)\.txt","",filename))
	


directorysalvataggio = directoryscelta + "/" + ngeni +"geni"+"/"
try:
	os.mkdir(directorysalvataggio) 
except Exception as e: 
	print("la directory esisteva gia")


for UP_DOWN in ["UP", "DOWN"]: 
		
	for i, filename in enumerate(onlyfiles):
		trattamento = trattamenti[i]
		daStampare = 'Calculating for ' +TFtool +" " +str(i)+' '+ trattamento +" " + ngeni + " geni"+ " " + UP_DOWN
		trattamento = trattamento.replace("/","~")
		filename = filename.replace("/","~")
		if(UP_DOWN == "DOWN"):
			filename = filename.replace("_up.","_down.")

		
		print(daStampare,  flush=True)

		stringafileinput =   completapercorso + "/" + UP_DOWN + "_" + ngeni +'/' + filename
		pathtofileoutput = directorysalvataggio + "risultatiCheA_" + trattamento +"geni" + UP_DOWN + ngeni +".txt"

		
		try:
			if(TFtool == "Chea2"):
				Chea2_Run(stringafileinput, pathtofileoutput)
			elif(TFtool == "ARCHS4"):
				type_database = 'ARCHS4--Coexpression'
				print("ARCHS4 done")
				Chea3_Run(stringafileinput, pathtofileoutput, type_database)
			elif(TFtool == "MeanRank"):
				type_database = 'Integrated--meanRank'
				print("MeanRank done")
				Chea3_Run(stringafileinput, pathtofileoutput, type_database)
			elif(TFtool == "Chea3"):
				type_database = "Literature--ChIP-seq"
				Chea3_Run(stringafileinput,pathtofileoutput,  type_database)
				print("Chea3 done")
			else:
				print("The tool is either not available or not correctly named by the user. Check your input")
				

		except Exception as e:
			print("C'e stato un errore per questo trattamento, mi spiace")
			print("l'eccezione e' questa:", e)
	
print("finita l'esecuzione")
		
