import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import animation
import numpy as np
import math

number_of_cells = 1000
time_points = 3000


values = range(number_of_cells)

# set up dummy data
values = np.zeros( (time_points,number_of_cells) ,dtype=float)
for cell in range(number_of_cells):
    for t_point in range(time_points):
        values[t_point,cell] = (cell + t_point) % number_of_cells


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
cNorm  = colors.Normalize(vmin=0, vmax=values[-1,-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=reds)
#print(scalarMap.get_clim())

redColorVal = np.zeros((time_points,number_of_cells,4),dtype=float)
for cell in range(number_of_cells):
    redColorVal[0,cell,:] = scalarMap.to_rgba(values[0,cell])

for t_point in range(1,time_points):
    for cell in range(number_of_cells):
        if cell == 0:
            cell_minus_one = number_of_cells - 1
        else:
            cell_minus_one = cell - 1
        
        redColorVal[t_point,cell,:] = redColorVal[t_point-1,cell_minus_one,:]
        
plotted = ax.scatter(x2,y,s=100)
plotted.set_color(redColorVal[0,:,:])

im = dummy_ax.imshow(colorValues/8, cmap='Reds')
fig.colorbar(im, cax=cax, orientation='vertical',label='units')

ax.set_aspect('equal')
ax.axis('off')
ax.text(-5,3,"Marginal zone\nconcentration", fontsize=20)


current_time = 0
time_string = 't = ' + str(current_time) + 's'
chem_string = 'Vg1'

ax.text(-0.35, 0.5, chem_string,fontsize=20)
time_text = ax.text(-0.7, -0.5,[],fontsize=16)
time_text.set_text(time_string)

def init():
    plotted.set_color(redColorVal[0,:,:])
    current_time = 0
    time_string = 't = ' + str(current_time) + 's'
    time_text.set_text(time_string)
    return plotted, time_text

def animate(i):
    sample_rate = 10
    plotted.set_color(redColorVal[i*sample_rate,:,:])
    current_time = i*sample_rate
    time_string = 't = ' + str(current_time) + 's'
    time_text.set_text(time_string)
    return plotted, time_text

anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=300, blit=True)

anim.save('circle_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

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



