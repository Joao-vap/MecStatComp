# fit exponential to data

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data =  [4,2.0000000000000000,     
        5,4.0000000000000000,   
        6,10.000000000000000,  
        7,25.000000000000000, 
        8,59.000000000000000,
        9,190.00000000000000,
        10,529.00000000000000]   

x = np.array(data[0::2])
y = np.array(data[1::2])

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

popt, pcov = curve_fit(func, x, y, p0=(1, 1e-6, 1))

# print exponential fit parameters
space = np.linspace(0, 10, 1000)
exp = func(space, *popt)
print("a = %s , b = %s, c = %s" % (popt[0], popt[1], popt[2]))

# plot data and fit
plt.plot(x, y, 'ko', label="Original Noised Data")
plt.plot(space, exp, 'g-', label="Fitted Curve")
plt.legend()
plt.show()



