import numpy as np
import matplotlib.pyplot as plt
import csv
Labels = ["I", "r", "Numerator", "G(r)"]
print(Labels)

#A = np.loadtxt("gofr.dat")#,dtype={'names': ('Tstep', 'Temp', '<Temp>', 'PE', '<PE>', 'Pressure', '<Pres>', 'Order', 'H-fn', 'Nabors', 'E'),
#...                      'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')})
A = np.loadtxt("gofr.dat")
Nstep = A[:,0]

i = input("Enter which column to print: ")
ylbl = Labels[i]
xlbl = Labels[1]
plt.plot(Nstep , A[:,i])
plt.xlabel(xlbl)
plt.ylabel(ylbl)
plt.show()

