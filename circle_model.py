import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import animation
import numpy as np
import math
import pickle

number_of_cells = 751

inducer = pickle.load( open( "inducer.p", "rb" ) )
inhibitor = pickle.load( open( "inhibitor.p", "rb" ) )
propagator = pickle.load( open( "propagator.p", "rb" ) )

number_of_steps = len(inducer)

inducer_dummy = np.arange(256)
inducer_dummy = inducer_dummy.reshape(16,16)
inhibitor_dummy = np.arange(256)
inhibitor_dummy = inhibitor_dummy.reshape(16,16)
propagator_dummy = np.arange(256)
propagator_dummy = propagator_dummy.reshape(16,16)

theta = np.linspace(0,2*math.pi,number_of_cells) + (3*math.pi)/2
r = 3

x1 = [r*math.cos(angle) - 2*r - 1  for angle in theta]
x2 = [r*math.cos(angle) for angle in theta]
x3 = [r*math.cos(angle) + 2*r + 1 for angle in theta]
y = [r*math.sin(angle) for angle in theta]

dummy_fig = plt.figure()
dummy_ax = dummy_fig.add_subplot(111)

fig = plt.figure()
inducer_ax = fig.add_subplot(311)
inhibitor_ax = fig.add_subplot(312)
propagator_ax = fig.add_subplot(313)

# inducer_cax = fig.add_axes([0.83 - 7, 0.05 , 0.05, 0.3])
# inhibitor_cax = fig.add_axes([0.83, 0.05 , 0.05, 0.3])
# propagator_cax = fig.add_axes([0.83 + 7, 0.05 , 0.05, 0.3])

inducer_cNorm  = colors.Normalize(vmin=np.amin(inducer), vmax=np.amax(inducer))
inducer_scalarMap = cmx.ScalarMappable(norm=inducer_cNorm, cmap ='Blues')
inhibitor_cNorm  = colors.Normalize(vmin=np.amin(inhibitor), vmax=np.amax(inhibitor))
inhibitor_scalarMap = cmx.ScalarMappable(norm=inhibitor_cNorm, cmap ='Reds')
propagator_cNorm  = colors.Normalize(vmin=np.amin(propagator), vmax=np.amax(propagator))
propagator_scalarMap = cmx.ScalarMappable(norm=propagator_cNorm, cmap ='Purples')

# inducer_colors = np.zeros((time_points,number_of_cells,4),dtype=float)
inducer_colors = inducer_scalarMap.to_rgba(inducer)
inhibitor_colors = inhibitor_scalarMap.to_rgba(inhibitor)
propagator_colors = propagator_scalarMap.to_rgba(propagator)
        
inducer_plotted = inducer_ax.scatter(x1,y,s=50)
inducer_plotted.set_color(inducer_colors[0,:,:])
inhibitor_plotted = inhibitor_ax.scatter(x2,y,s=50)
inhibitor_plotted.set_color(inhibitor_colors[0,:,:])
propagator_plotted = propagator_ax.scatter(x3,y,s=50)
propagator_plotted.set_color(propagator_colors[0,:,:])

# im = dummy_ax.imshow(inducer_dummy/8, cmap='Blues')
# fig.colorbar(im, cax=inducer_cax, orientation='vertical',label='units')
# im = dummy_ax.imshow(inducer_dummy/8, cmap='Reds')
# fig.colorbar(im, cax=inhibitor_cax, orientation='vertical',label='units')
# im = dummy_ax.imshow(inducer_dummy/8, cmap='Purples')
# fig.colorbar(im, cax=propagator_cax, orientation='vertical',label='units')

inducer_ax.set_aspect('equal')
inducer_ax.axis('off')
inhibitor_ax.set_aspect('equal')
inhibitor_ax.axis('off')
propagator_ax.set_aspect('equal')
propagator_ax.axis('off')
# ax.text(-5,3,"Marginal zone\nconcentration", fontsize=20)


current_time = 0
time_string = 't = ' + str(current_time) + 's'
inducer_string = 'cVg1'
inhibitor_string = 'BMP'
propagator_string = 'Calcium'

inducer_ax.text(-0.35 - 7, 0.5, inducer_string,fontsize=20)
inhibitor_ax.text(-0.35, 0.5, inhibitor_string,fontsize=20)
propagator_ax.text(-0.35 + 7, 0.5, propagator_string,fontsize=20)
time_text = inhibitor_ax.text(-0.7, -0.5,[],fontsize=16)
time_text.set_text(time_string)

def init():
    inducer_plotted.set_color(inducer_colors[0,:,:])
    inhibitor_plotted.set_color(inhibitor_colors[0,:,:])
    propagator_plotted.set_color(propagator_colors[0,:,:])
    current_time = 0
    time_string = 't = ' + str(current_time) + 's'
    time_text.set_text(time_string)
    return inducer_plotted, inhibitor_plotted, propagator_plotted, time_text

sample_rate = 20
def animate(i):
    sample_rate = 20
    inducer_plotted.set_color(inducer_colors[i*sample_rate,:,:])
    inhibitor_plotted.set_color(inhibitor_colors[i*sample_rate,:,:])
    propagator_plotted.set_color(propagator_colors[i*sample_rate,:,:])
    current_time = i*sample_rate
    time_string = 't = ' + str(current_time) + 's'
    time_text.set_text(time_string)
    return inducer_plotted, inhibitor_plotted, propagator_plotted, time_text
    
number_of_frames = int(np.ceil(number_of_steps / sample_rate))

anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=number_of_frames, blit=True)

anim.save('model_circle_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

plt.show()
    



