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

alpha_i = 0.3 # micrometers per second
k_i = np.power(350,4) # nano Molar

a_ii = 50 # no units


# plot i
p,i = np.meshgrid(p,i)
z = (alpha_i * (np.power(p,4) + np.power(a_ii*i,4))) / (k_i + np.power(p,4) + np.power(a_ii*i,4))

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