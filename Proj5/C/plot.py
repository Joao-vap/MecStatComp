#plot files

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#read data
data = np.loadtxt('c1ene1100.dat', delimiter=',')

print(data)

#plot
plt.plot(data[:,0], data[:,1])
plt.xlabel('Iterações')
plt.ylabel('Energia')
plt.title('step = 0.0001')
plt.show()
