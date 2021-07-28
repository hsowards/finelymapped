import numpy as np
import pandas as pd
import sys

file = sys.argv[1]
out =  sys.argv[2]

data = np.load(file)
lst = data.files
print(lst)
print(data['data'][:10,])
print(data['row'][:10,])
print(data['col'][:10,])
print(len(data['row']))
print(len(data['col']))
print(data['shape'])
print(data['format'])

#for item in lst:
	#print(data[item].head())

#np.savetxt(out, data['data'], newline='\n')
