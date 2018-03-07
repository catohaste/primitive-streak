from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
v = np.arange(0,150,0.5)
i = np.arange(0,100,0.5)
v,i = np.meshgrid(v,i)
alpha_v = np.power(30,4)
k_v = 5
a_vv = 1
a_iv = 1
z = (k_v * np.power(v,4)) / (alpha_v + a_vv * np.power(v,4) + a_iv*np.power(i,4))

# Plot the surface.
surf = ax.plot_surface(v, i, z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-0.01, 5.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()