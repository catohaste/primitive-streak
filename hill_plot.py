from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data
i = np.arange(0,10,0.01)
p = np.arange(0,1000,1)

alpha_ii = 0.1 # micrometers per second 
k_ii = np.power(2,4) # nano Molar
a_ii = 1 # no units

alpha_pi = 0.3 # micrometers per second 
k_pi = np.power(500,4) # nano Molar
a_pi = 1 # no units

# plot i
p,i = np.meshgrid(p,i)
z = (alpha_pi * np.power(p,4)) / (k_pi + a_pi * np.power(p,4)) + (alpha_ii * np.power(i,4)) / (k_ii + a_ii * np.power(i,4))

# Plot the surface.
surf = ax.plot_surface(p, i, z, cmap=cm.coolwarm,linewidth=0, antialiased=False)
                       


# Customize the z axis.
ax.set_zlim(0.0, 0.5)
ax.set_xlabel('Propagator conc.')
ax.set_ylabel('Inhibitor conc.')
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()