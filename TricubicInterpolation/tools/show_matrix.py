import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append('../pyTricubic')
from tricubic_matrix import tricubicMat

tM = tricubicMat[:,:].copy()
tM = 1.*(tM != 0)

def sc(AA,i,j):
    AA[:,[i,j]] = AA[:,[j,i]]
def sr(AA,i,j):
    AA[[i,j],:] = AA[[j,i],:]



rows = []
for j in range(64):
    if sum(tM[j,:]) == 1:
        rows.append(j)
print(rows)

plt.close('all')
plt.matshow(tM)
plt.show()
