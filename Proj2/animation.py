#plot 'onda.dat'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

# Set up the figure
fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)

# import data
data = np.loadtxt('onda.dat', usecols=range(100))

x = np.linspace(0, 1, 100)

# Initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# Animation function.  This is called sequentially
def animate(i):
    x = np.linspace(0, 1, 100)
    y = data[i,:]
    line.set_data(x, y)
    return line,

# Call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                                 frames=1000, interval=20, blit=True)

# Set title
plt.title('Error Propagation')

# Save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
