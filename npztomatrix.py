import scipy.sparse as sparse
import numpy as np
import pandas as pd
import sys

R = sparse.load_npz(sys.argv[1]).toarray()
R = R + R.T
print(len(R))
print(R[0:10,0:10])

np.savetxt(sys.argv[2], R, delimiter=",")


