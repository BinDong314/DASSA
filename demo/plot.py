import h5py
import matplotlib.pyplot as plt    

import numpy as np

#h5diff -d 0.000000000001 -v  
file_org="/project/projectdirs/m1248/dassa-new/tmatch-verify-data/output/tmach-jun5-1740-2tpt.h5"
file_reg="/project/projectdirs/m1248/dassa-new/tmatch-verify-data/output/tmach-reorg-jun5-1740-2tpt.h5"

f = h5py.File(file_reg, "r")
data = f.get('/dat')[:]
print("Values bigger than 0.5 =", data[data>0.5])
print("Their indices are ", np.nonzero(data > 0.5))
plt.plot(data)
plt.show()



f = h5py.File(file_org, "r")
data = f.get('/dat')[:]
print("Values bigger than 0.5 =", data[data>0.5])
print("Their indices are ", np.nonzero(data > 0.5))
plt.plot(data)
plt.show()
