############ This little tool contains the functions Chea3_Run and Chea2_Run, that enables us to do the enrichment for the TFs. 
#At first in my project I was using Chea2
import argparse
import sys 
import os
import os.path
from os import listdir
from os.path import isfile, join
import re

from ScriptX2K import *


def Chea3_Run(nome_file,nome_output, type_database):
	stringacomando = "Rscript  script/Chea3inR.r '" + nome_file + "' '" + nome_output+"'" + " '"+type_database+"'"
	print(nome_output)
	os.system(stringacomando) #this permits me to launch the script Chea3inR 
	return


def Chea2_Run(stringafileinput,pathtofileoutput):
	stringafileinput = stringafileinput.replace("SSSSSS"," ")
	file_espressione = open(stringafileinput,"r")
 
	contenuto_file =file_espressione.read()
	file_espressione.close()

	raw_file = contenuto_file.split("\n")[1:-1]  
	input_genes = []

	for i in raw_file:
		riga = i.split("\t")
		input_genes.append(riga[0])

	x2k_results = run_X2K(input_genes)
	print("X2K done")
	list_dict = x2k_results['ChEA']
	pathtofile = pathtofileoutput
	out_file = open(pathtofile,"w")
	length_list = len(list_dict) - 1
	table_out = []
	for k, dict in enumerate(list_dict):
		if k==0:
			header = list(dict.keys())
			headers = "\t".join(header)
			out_file.write(headers + "\n")
		row_value = list(dict.values())

		row_value = [ str(element) for element in row_value ]
		string_row = '\t'.join(row_value)
		table_out.append(string_row)
		if k== length_list:
			out_file.write(string_row)

		else:
			out_file.write(string_row + "\n")
	

	out_file.close()
