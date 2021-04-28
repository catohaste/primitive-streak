from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

# Make data
b = np.arange(0,10,0.01)
b_for_v = np.arange(4,9,0.01)
v = np.arange(0,12,0.01)

points_b = b.shape[0]

# parameters for dI
power_b = 4
basal_b = 0
alpha_b = 0.004 # micrometers per second
k_b = np.power(6,power_b) # nano Molar
a_bb = 1 # no units
a_cb = 1 # no units
a_b = 1
calc_rescale = 0.05
# calc_rescale = 1
decay_param_b = alpha_b / 10

# parameters for dV
basal_v = 0.0005
alpha_v = 0.004
k_v = np.power(3,4)
a_vv = 1
a_bv = 1
decay_param_v = alpha_v / 15


# # dI ######################################################
#
# # calculate dI/dt

# Here is the original equation
# The calcium concentration can be rescaled so that is of approximately the same order as the BMP concentration
c_off = 0
c = c_off
production_calcium_off = basal_b + (alpha_b * (np.power(calc_rescale*c,power_b) + np.power(b,power_b))) / (k_b + np.power(a_cb * c,power_b) + np.power(a_bb * b,power_b))

c_on = 126
c = c_on
production_calcium_on = basal_b + (alpha_b * (np.power(calc_rescale*c,power_b) + np.power(b,power_b))) / (k_b + np.power(a_cb * calc_rescale * c,power_b) + np.power(a_bb * b,power_b))

decay_b = - b * decay_param_b

diffusion_b = b * 0.00001 - 0.00001

dB_dt_on = production_calcium_on + decay_b + diffusion_b

dB_dt_off = production_calcium_off + decay_b + diffusion_b

# I also tested an alternative 
# c_off = 20
# c = c_off
# dB_calcium_off = alpha_b * (np.power(c,power_b)+ np.power(b,power_b) + np.power(c * b,power_b)) / (k_b + np.power(a_cb * c,power_b) + np.power(a_bb * b,power_b) + np.power(a_b * b * c,power_b))
#
# c_on = 126
# c = c_on
# dB_calcium_on = alpha_b * (np.power(c,power_b)+ np.power(b,power_b) + np.power(c * b,power_b)) / (k_b + np.power(a_cb * c,power_b) + np.power(a_bb * b,power_b) + np.power(a_b * b * c,power_b))

fig_calcium_on = plt.figure()
ax_calcium_on = fig_calcium_on.gca()

ax_calcium_on.set_title('Calcium on')
ax_calcium_on.set_ylim(-alpha_b*1.1, alpha_b*1.1)
ax_calcium_on.set_xlabel('BMP conc.')
ax_calcium_on.set_ylabel('dB / dt')
ax_calcium_on.axhline(y=0, color='k')
ax_calcium_on.axvline(x=0, color='k')

line_production_on = ax_calcium_on.plot(b,production_calcium_on,label='production')
line_decay_on = ax_calcium_on.plot(b,decay_b,label='decay')
# line_diffusion_on = ax_calcium_on.plot(b,diffusion_b,label='diffusion')
line_net_on = ax_calcium_on.plot(b,dB_dt_on,label='NET')

ax_calcium_on.legend()

plt.show()

fig_calcium_off = plt.figure()
ax_calcium_off = fig_calcium_off.gca()

ax_calcium_off.set_title('Calcium off')
ax_calcium_off.set_ylim(-alpha_b*1.1, alpha_b*1.1)
ax_calcium_off.set_xlabel('BMP conc.')
ax_calcium_off.set_ylabel('dB / dt')
ax_calcium_off.axhline(y=0, color='k')
ax_calcium_off.axvline(x=0, color='k')

line_production_off = ax_calcium_off.plot(b,production_calcium_off,label='production')
line_decay_off = ax_calcium_off.plot(b,decay_b,label='decay')
# line_diffusion_off = ax_calcium_off.plot(b,diffusion_b,label='diffusion')
line_net_off = ax_calcium_off.plot(b,dB_dt_off,label='NET')

ax_calcium_off.legend()

plt.show()


# dV ###################################################### 

# calculate dV/dt
# v,b_for_v = np.meshgrid(v,b_for_v)
# production_v = basal_v + (alpha_v * (np.power(v,4)) / (k_v + np.power(a_bv*b_for_v,4) + np.power(a_vv*v,4)))
#
# decay_v = -v * decay_param_v
#
# net_dV = production_v + decay_v
#
# fig_dV = plt.figure()
# ax_dV = fig_dV.gca(projection='3d')
#
# #surface_production = ax_dV.plot_surface(v, b_for_v, production_v, cmap=cm.coolwarm,linewidth=0.1, antialiased=False)
# #surface_decay = ax_dV.plot_surface(v, b_for_v, decay_v, cmap=cm.viridis,linewidth=0, antialiased=False)
# surface_net = ax_dV.plot_surface(v, b_for_v, net_dV, cmap=cm.magma,linewidth=0, antialiased=False)
#
# # Customize the z axis.
# ax_dV.set_zlim(-alpha_v*1.1, alpha_v*1.1)
# ax_dV.set_xlabel('Inducer conc.')
# ax_dV.set_ylabel('Inhibitor conc.')
# ax_dV.set_zlabel('dV')
# ax_dV.zaxis.set_major_locator(LinearLocator(10))
# ax_dV.zaxis.set_major_formatter(FormatStrFormatter('%.04f'))
# # Add a color bar which maps values to colors
# fig_dV.colorbar(surface_net, shrink=0.5, aspect=5)
#
# plt.show()
