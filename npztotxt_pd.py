import pandas as pd
import numpy as np
import sys

npz = np.load(sys.argv[1])
df = pd.DataFrame.from_dict({item: npz[item] for item in npz.files}, orient='index')
print(df.head())
