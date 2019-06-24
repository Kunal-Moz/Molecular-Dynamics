import numpy as np
import matplotlib.pyplot as plt
import csv
Labels = ['Tstep', 'Temp', '<Temp>', 'PE', '<PE>', 'Pressure', '<Pres>', 'Order', 'H-fn', 'Nabors', 'E']
print(Labels)

A = np.loadtxt("run1.dat")#,dtype={'names': ('Tstep', 'Temp', '<Temp>', 'PE', '<PE>', 'Pressure', '<Pres>', 'Order', 'H-fn', 'Nabors', 'E'),
#...                      'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')})
A = np.loadtxt("run_prod1.dat")
Nstep = A[:,0]

i = input("Enter which column to print: ")
ylbl = Labels[i]

plt.plot(Nstep , A[:,i])
plt.xlabel("T")
plt.ylabel(ylbl)
plt.show()





