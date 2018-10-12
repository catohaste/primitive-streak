
import pickle
import numpy as np
import math

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
from matplotlib import animation




file_prefix = 'test' # date_string
pickle_dir = "output/pickles/"+file_prefix
propagator = pickle.load( open( pickle_dir+"/"+file_prefix+"_propagator.p", "rb" ) )
inducer = pickle.load( open( pickle_dir+"/"+file_prefix+"_inducer.p", "rb" ) )
inhibitor = pickle.load( open( pickle_dir+"/"+file_prefix+"_inhibitor.p", "rb" ) )

inducer_min=np.amin(inducer)
inducer_max=np.amax(inducer)
inducer_threshold = 2
colors1_proportion = (inducer_threshold - inducer_min) / (inducer_max - inducer_min)

def create_colormap_2colors(proportion_colors1):
    
    number_colors1 = math.floor(proportion_colors1 * 256)
    number_colors2 = 256 - number_colors1
    
    # sample the colormaps that you want to use. Use 128 from each so we get 256
    # colors in total
    colors1 = cmx.Blues(np.linspace(0, 1, number_colors1))
    oranges = cmx.Oranges(np.linspace(0, 1, number_colors2))
    orange = oranges[math.ceil(number_colors2*0.65),:]
    colors2 = np.zeros( (number_colors2,4), dtype=np.float64)
    colors2[:number_colors2] = orange
    
    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
    return mymap
    
def plot_color_gradients(cmap):
    
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(gradient, aspect='auto', cmap=cmap)
    ax.set_axis_off()

    plt.show()
    

def create_circle_animation(var, colormap, label_string, file_prefix):
    
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
    
    cNorm  = mcolors.Normalize(vmin=np.amin(var), vmax=np.amax(var))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap =colormap)
    
    # var_colors = np.zeros((time_points,number_of_cells,4),dtype=float) OLD
    var_colors = scalarMap.to_rgba(var)
    
    plotted = ax.scatter(x,y,s=50)
    plotted.set_color(var_colors[0,:,:])
    
    print(var.max())
    colorbar_max = var.max()
    
    colormap_ax = fig.add_axes([0.83, 0.05 , 0.05, 0.3])
    im = dummy_ax.imshow((var_dummy/256)*colorbar_max, cmap=colormap)
    fig.colorbar(im, cax=colormap_ax, orientation='vertical',label='nM')
    
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
        current_time = int((i*sample_rate) / 2)
        time_string = 't = ' + str(current_time) + 's'
        time_text.set_text(time_string)
        return plotted, time_text

    number_of_frames = int(np.ceil(number_of_steps / sample_rate))

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=number_of_frames, blit=True)

    anim.save("output/videos/"+file_prefix+'_'+label_string+'.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

    #plt.show()


inducer_colormap = create_colormap_2colors(colors1_proportion)

# plot_color_gradients(inducer_colormap)

create_circle_animation(propagator,plt.get_cmap('Oranges'),'Calcium',file_prefix)
create_circle_animation(inhibitor,plt.get_cmap('Reds'),'BMP',file_prefix)
create_circle_animation(inducer,inducer_colormap,'cVg1',file_prefix)

