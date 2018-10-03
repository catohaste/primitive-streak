
import pickle
import numpy as np
import math

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import animation




file_prefix = 'test' # date_string
pickle_dir = "output/pickles/"+file_prefix
propagator = pickle.load( open( pickle_dir+"/"+file_prefix+"_propagator.p", "rb" ) )
inducer = pickle.load( open( pickle_dir+"/"+file_prefix+"_inducer.p", "rb" ) )
inhibitor = pickle.load( open( pickle_dir+"/"+file_prefix+"_inhibitor.p", "rb" ) )


def create_circle_animation(var, cmap_string, label_string, file_prefix):
    
    number_of_steps = var.shape[0]
    number_of_cells = var.shape[1]
    
    var_dummy = np.arange(256)
    var_dummy = var_dummy.reshape(16,16)
    
    theta = np.linspace(0,2*math.pi,number_of_cells) + (3*math.pi)/2
    radius = 3
    
    x = [radius*math.cos(angle) for angle in theta]
    y = [radius*math.sin(angle) for angle in theta]
    
    dummy_fig = plt.figure()
    dummy_ax = dummy_fig.add_subplot(111)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    cNorm  = colors.Normalize(vmin=np.amin(var), vmax=np.amax(var))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap =cmap_string)
    
    # var_colors = np.zeros((time_points,number_of_cells,4),dtype=float) OLD
    var_colors = scalarMap.to_rgba(var)
    
    plotted = ax.scatter(x,y,s=50)
    plotted.set_color(var_colors[0,:,:])
    
    colormap_ax = fig.add_axes([0.83, 0.05 , 0.05, 0.3])
    im = dummy_ax.imshow(var_dummy/8, cmap=cmap_string)
    fig.colorbar(im, cax=colormap_ax, orientation='vertical',label='units')
    
    ax.set_aspect('equal')
    ax.axis('off')
    
    # ax.text(-5,3,"Marginal zone\nconcentration", fontsize=20)
    
    ax.text(-0.9, 0.5, label_string, fontsize=20)
    
    current_time = 0
    time_string = 't = ' + str(current_time) + 's'
    time_text = ax.text(-0.9, -0.5,[],fontsize=16)
    time_text.set_text(time_string)
    
    def init():
        plotted.set_color(var_colors[0,:,:])
        current_time = 0
        time_string = 't = ' + str(current_time) + 's'
        time_text.set_text(time_string)
        return plotted, time_text

    sample_rate = 20
    def animate(i):
        sample_rate = 20
        plotted.set_color(var_colors[i*sample_rate,:,:])
        current_time = i*sample_rate
        time_string = 't = ' + str(current_time) + 's'
        time_text.set_text(time_string)
        return plotted, time_text

    number_of_frames = int(np.ceil(number_of_steps / sample_rate))

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=number_of_frames, blit=True)

    anim.save("output/videos/"+file_prefix+'_'+label_string+'.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

    #plt.show()


create_circle_animation(propagator,'Purples','Calcium',file_prefix)
create_circle_animation(inhibitor,'Reds','BMP4',file_prefix)
create_circle_animation(inducer,'Greens','Vg1',file_prefix)

