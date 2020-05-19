#I've created this file for the manipulation and the integration of various tables I had, it's not part of the main project

import numpy as np
import pandas as pd
import os
import csv
pathbase = './data/DrugInformation/'



rows_list = []
with open(pathbase + 'PrestwickLibrary.csv', 'r') as csvFile:
	reader = csv.reader(csvFile)
	for row in reader:
		#row = [x.splitlines() for x in row ]
		rows_list.append(row)

csvFile.close()
tabellacompleta = pd.DataFrame(rows_list[1:],columns = rows_list[0])
tabellacompleta[tabellacompleta=='']=np.nan


tabellacompleta=tabellacompleta.apply(lambda x: x.ffill().bfill())
tabellacompleta.groupby('Chemical name').apply(lambda x: x.ffill().bfill()).drop_duplicates()
tabellacompleta=tabellacompleta.groupby('Chemical name').agg(lambda x: '---'.join(set(x)))

tabellacompleta.iloc[0:10]
tabellacompleta.to_csv(pathbase + 'PrestwickModified.tab', sep='\t')
