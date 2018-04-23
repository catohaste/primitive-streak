from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

# Make data
i = np.arange(0,10,0.01)
p = np.arange(0,1000,1)
v = np.arange(0,10,0.01)

# parameters for dI
alpha_i = 0.005 # micrometers per second
k_i = np.power(400,4) # nano Molar
a_ii = 40 # no units

# parameters for dV
alpha_v = 0.005
k_v = np.power(0.75,4)
a_vv = 0.5
a_iv = 2.5
decay_v = 0.001

c=0.001

# # dI ######################################################
#
# # calculate dI/dt
# p,i = np.meshgrid(p,i)
# dI = (alpha_i * (np.power(p,4) + np.power(a_ii*i,4))) / (k_i + np.power(p,4) + np.power(a_ii*i,4))
#
# fig_dI = plt.figure()
# ax_dI = fig_dI.gca(projection='3d')
#
# surface_dI = ax_dI.plot_surface(p, i, dI, cmap=cm.coolwarm,linewidth=0, antialiased=False)
#
# # Customize the z axis.
# ax_dI.set_zlim(0.0, 0.005)
# ax_dI.set_xlabel('Propagator conc.')
# ax_dI.set_ylabel('Inhibitor conc.')
# ax_dI.zaxis.set_major_locator(LinearLocator(10))
# ax_dI.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# # Add a color bar which maps values to colors.
# fig_dI.colorbar(surface_dI, shrink=0.5, aspect=5)
#
# plt.show()


# dV ###################################################### 

# calculate dV/dt
v,i = np.meshgrid(v,i)
dV = (alpha_v * np.power(a_vv*v + c,4)) / (k_v + np.power(a_iv*i,4) + np.power(a_vv*v + c,4))
 
fig_dV = plt.figure()
ax_dV = fig_dV.gca(projection='3d')

surface_dV = ax_dV.plot_surface(v, i, dV, cmap=cm.coolwarm,linewidth=0, antialiased=False)
                
# Customize the z axis.
ax_dV.set_zlim(0, 0.005)
ax_dV.set_xlabel('Inducer conc.')
ax_dV.set_ylabel('Inhibitor conc.')
ax_dV.zaxis.set_major_locator(LinearLocator(10))
ax_dV.zaxis.set_major_formatter(FormatStrFormatter('%.04f'))
# Add a color bar which maps values to colors.
fig_dV.colorbar(surface_dV, shrink=0.5, aspect=5)

plt.show()
