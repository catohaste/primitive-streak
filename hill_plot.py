from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data
v = np.arange(0,150,0.5)
i = np.arange(0,100,0.5)
p = np.arange(0,10,0.02)
alpha_v = np.power(30,4)
alpha_i = np.power(2,4)
alpha_p = np.power(200,4)
k_v = 5
k_i = 0.5
k_p = 0.7
a_vv = 1
a_iv = 1
a_vi = 1
a_pi = 1
a_ii = 1
a_vp = 1
a_pp = 1


# plot v
# v,i = np.meshgrid(v,i)
# z= (k_v * np.power(v,4)) / (alpha_v + a_vv * np.power(v,4) + a_iv*np.power(i,4))

# plot p
# v,p = np.meshgrid(v,p)
# z = (k_p * np.power(v,4)) / (alpha_p + a_vp * np.power(v,4))

# plot i
p,i = np.meshgrid(p,i)
z = (k_i * np.power(p,4)) / (alpha_i + a_pi * np.power(p,4))

# Plot the surface.
surf = ax.plot_surface(p, i, z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
                       


# Customize the z axis.
ax.set_zlim(-0.01, 0.5)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()