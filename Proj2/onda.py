# plot onda.dat file

import numpy as np
import matplotlib.pyplot as plt
import math

# we will shift the data in the y-axis by 1.01 every line of the data
# so that we can see the lines on top of each other
shift = 0.3
itime = 100

# read the data
data = np.loadtxt('tarefa-6-aux-11820812.out', usecols=range(100))

# create the figure
fig = plt.figure()
ax = fig.add_subplot(111)

# plot the data
for i in range(itime, 0, -1):
    x = np.linspace(0, 1, 100)
    y = data[i,:]
    # black lines with thickness 1
    ax.plot(x, shift*i - y, 'k-', lw=0.2)

# set the limits of the axes

ax.set_xlim(0, 1)

# set the labels of the axes
ax.set_xlabel('Espaço')
ax.set_ylabel('Tempo (n * 3/10)')

# Set the title of the plot
ax.set_title('')

# set aspect ratio so that y is longer than x
ax.set_aspect('auto')

# invert the y-axis
ax.invert_yaxis()

# show the plot
plt.show()