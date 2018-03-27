import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import math

number_of_cells = 1000
time_points = 3000


values = range(number_of_cells)

# set up dummy data
# values = np.zeros( (time_points,number_of_cells) ,dtype=float)
# for cell in range(number_of_cells):
#     for t_point in range(time_points):
#         values[t_point,cell] = (cell + t_point) % number_of_cells


colorValues = np.arange(1024)
colorValues = colorValues.reshape(32,32)

theta = np.linspace(0,2*math.pi,number_of_cells)
r = 3

x1 = [r*math.cos(angle) - 7  for angle in theta]
x2 = [r*math.cos(angle) for angle in theta]
x3 = [r*math.cos(angle) + 7 for angle in theta]
y = [r*math.sin(angle) for angle in theta]

dummy_fig = plt.figure()
dummy_ax = dummy_fig.add_subplot(111)


fig = plt.figure()
ax = fig.add_subplot(111)
cax = fig.add_axes([0.83, 0.05 , 0.05, 0.3])


reds = cm = plt.get_cmap('Reds')
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=reds)
#print(scalarMap.get_clim())

redColorVal = np.zeros((number_of_cells,4),dtype=float)
# for t_point in range(time_points):
for cell in range(number_of_cells):
    redColorVal[cell,:] = scalarMap.to_rgba(values[cell])
        
# print(redColorVal[0,:,:])
        
plotted = ax.scatter(x2,y,c=redColorVal)

im = dummy_ax.imshow(colorValues, cmap='Reds')
fig.colorbar(im, cax=cax, orientation='vertical')

ax.set_aspect('equal')
ax.axis('off')

current_time = 0
time_string = 't = ' + str(current_time) + 's'
chem_string = 'Vg1'

ax.text(-0.35, 0.5, chem_string,fontsize=20)
time_text = ax.text(-0.5, -0.5, time_string,fontsize=16)

def init():
    line_propagator.set_data([], [])
    #line_inducer.set_data([], [])
    #line_inhibitor.set_data([], [])
    return line_propagator,#line_inducer,line_inhibitor,
    
def animate(i):
    x = list(range(number_of_cells))
    y_propagator = propagator[100*i,:]
    #y_inducer = inducer[100*i,:]
    #y_inhibitor = inhibitor[100*i,:]
    line_propagator.set_data(x, y_propagator)
    #line_inducer.set_data(x, y_inducer)
    #line_inhibitor.set_data(x, y_inhibitor)
    return line_propagator,#line_inducer,line_inhibitor,

#anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=30, blit=True)

#anim.save('circle_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

plt.show()


# fig = plt.figure()
# ax = plt.axes(xlim=(-1, 1000),ylim=(-0.5,max_y_value),xlabel='Cell ID',ylabel='Concentration')
# line_propagator, = ax.plot([], [],'o',color='purple',markersize = 0.5,label='Propagator')
# line_inducer, = ax.plot([], [],'go',markersize = 0.5,label='Inducer')
# line_inhibitor, = ax.plot([], [],'ro',markersize = 0.5,label='Inhibitor')
#
# # Now add the legend with some customizations.
# legend = ax.legend(loc='upper center',shadow=False,markerscale=3,numpoints=5)
# for label in legend.get_texts():
#     label.set_fontsize('small')
#
# anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=216, blit=True)
#
# anim.save('model_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])
#
# plt.show()



