
import pickle
import numpy as np
from decimal import Decimal

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
from matplotlib import animation
from matplotlib import style


def plot_initial_conditions(variables,params):
    
    number_of_cells = params['number_of_cells']
    inducer_threshold = params['inducer_threshold']

    style.use('fivethirtyeight')
    fig = plt.figure()

    # define colours
    color_inhibitor = '#D50032'
    color_inducer = '#0097A9'
    color_propagator = '#4B384C'

    max_yval_propagator = 6000
    max_yval_protein = 15

    # set up left axes
    ax_protein = plt.axes(xlim=(-1, number_of_cells + 1),ylim=(0.0,max_yval_protein),xlabel='Cell Location')
    ax_protein.set_ylabel('Protein conc. (nM)', color='black')
    ax_protein.tick_params('y', colors='black')
    plt.xticks(np.linspace(0,number_of_cells-1,num=3), ['post.','ant.','post.'])

    # set up right axes
    ax_propagator = ax_protein.twinx()
    ax_propagator.set_ylabel('Calcium conc. (nM)', color=color_propagator)
    ax_propagator.tick_params('y', colors=color_propagator)
    ax_propagator.set_ylim(0,max_yval_propagator)
    ax_propagator.grid(False)

    plt.tight_layout()

    # Define dummy lines for legend
    legend_propagator = ax_protein.plot([],[],'-',color=color_propagator,linewidth = 1.0,label='Calcium')

    # Initiate data to be animated
    line_propagator, = ax_propagator.plot(range(0,number_of_cells), variables.propagator[0,:], '-',color=color_propagator, linewidth = 1.0, label='Calcium')
    line_inhibitor, = ax_protein.plot(range(0,number_of_cells), variables.inhibitor[0,:], linewidth=0,marker='.', color=color_inhibitor, markersize = 1, label='BMP')
    line_inducer, = ax_protein.plot(range(0,number_of_cells), variables.inducer[0,:], linewidth=0,marker='.', color=color_inducer, markersize = 1, label='cVg1')
    time_text = ax_protein.text(0.75*number_of_cells,max_yval_protein*0.95,[],fontsize='medium')

    # draw inducer threshold line
    inducer_threshold_array = inducer_threshold * np.ones((number_of_cells+1,1),dtype=int)
    line_inducer_threshold, = ax_protein.plot(range(0,number_of_cells+1), inducer_threshold_array,'--', color=color_inducer, linewidth = 1,label='cVg1 threshold')
    # threshold_text = ax_protein.text(0.65*number_of_cells,threshold,'approx. threshold',fontsize='x-small')

    # Add the legend
    legend = ax_protein.legend(loc='upper left',shadow=False,markerscale=4,numpoints=5)
    for label in legend.get_texts():
        label.set_fontsize('small')

    plt.show()


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
    frames_per_second = 24
    interval = np.ceil(1000/frames_per_second)

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=interval, frames=number_of_frames, blit=True)

    anim.save("output/videos/"+file_prefix+'_'+label_string+'.mp4', fps=frames_per_second, extra_args=['-vcodec', 'libx264'])

    #plt.show()
    
def show_model_animation(variables,params):
    
    number_of_steps = params['number_of_steps']
    number_of_cells = params['number_of_cells']
    step_length = params['step_length']
    inducer_threshold = params['inducer_threshold']
    file_prefix = params['file_prefix']
    
    style.use('fivethirtyeight')
    fig = plt.figure()

    # define colours
    color_inhibitor = '#D50032'
    color_inducer = '#0097A9'
    color_propagator = '#4B384C'

    max_yval_propagator = 6000
    max_yval_protein = 15

    # set up left axes
    ax_protein = plt.axes(xlim=(-1, number_of_cells + 1),ylim=(0.0,max_yval_protein),xlabel='Cell Location')
    ax_protein.set_ylabel('Protein conc. (nM)', color='black')
    ax_protein.tick_params('y', colors='black')
    plt.xticks(np.linspace(0,number_of_cells-1,num=3), ['post.','ant.','post.'])

    # set up right axes
    ax_propagator = ax_protein.twinx()
    ax_propagator.set_ylabel('Calcium conc. (nM)', color=color_propagator)
    ax_propagator.tick_params('y', colors=color_propagator)
    ax_propagator.set_ylim(0,max_yval_propagator)
    ax_propagator.grid(False)

    plt.tight_layout()

    # Define dummy lines for legend
    legend_propagator = ax_protein.plot([],[],'-',color=color_propagator,linewidth = 1.0,label='Calcium')

    # Initiate data to be animated
    line_propagator, = ax_propagator.plot([], [],'-',color=color_propagator,linewidth = 1.0,label='Calcium')
    line_inhibitor, = ax_protein.plot([], [],linewidth=0,marker='.',color=color_inhibitor,markersize = 1,label='BMP')
    line_inducer, = ax_protein.plot([], [],linewidth=0,marker='.',color=color_inducer,markersize = 1,label='cVg1')
    time_text = ax_protein.text(0.75*number_of_cells,max_yval_protein*0.95,[],fontsize='medium')

    # draw inducer threshold line
    inducer_threshold_array = inducer_threshold * np.ones((number_of_cells+1,1),dtype=int)
    line_inducer_threshold, = ax_protein.plot(range(0,number_of_cells+1), inducer_threshold_array,'--', color=color_inducer, linewidth = 1,label='cVg1 threshold')
    # threshold_text = ax_protein.text(0.65*number_of_cells,threshold,'approx. threshold',fontsize='x-small')

    # Add the legend
    legend = ax_protein.legend(loc='upper left',shadow=False,markerscale=4,numpoints=5)
    for label in legend.get_texts():
        label.set_fontsize('small')


    def init():
        line_propagator.set_data([], [])
        line_inducer.set_data([], [])
        line_inhibitor.set_data([], [])
        time_text.set_text('')
        return time_text,line_propagator,line_inhibitor,line_inducer,

    sample_rate = 20
    def animate(i):
        sample_rate = 20
        x = list(range(number_of_cells))
        y_propagator = variables.propagator[sample_rate*i,:]
        y_inducer = variables.inducer[sample_rate*i,:]
        y_inhibitor = variables.inhibitor[sample_rate*i,:]
        line_propagator.set_data(x, y_propagator)
        line_inducer.set_data(x, y_inducer)
        line_inhibitor.set_data(x, y_inhibitor)
        current_time = Decimal(step_length*sample_rate*i)
        time_string = 't = ' + str(current_time) + 's'
        time_text.set_text(time_string)
        return time_text,line_propagator,line_inhibitor,line_inducer,

    number_of_frames = int(np.ceil(number_of_steps / sample_rate))

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=number_of_frames, blit=True)

    anim.save('output/videos/'+file_prefix+'_model_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

    plt.show()
    
    
# file_prefix = 'test' # date_string
# pickle_dir = "output/pickles/"+file_prefix
# propagator = pickle.load( open( pickle_dir+"/"+file_prefix+"_propagator.p", "rb" ) )
# inducer = pickle.load( open( pickle_dir+"/"+file_prefix+"_inducer.p", "rb" ) )
# inhibitor = pickle.load( open( pickle_dir+"/"+file_prefix+"_inhibitor.p", "rb" ) )


# inducer_min=np.amin(inducer)
# inducer_max=np.amax(inducer)
# colors1_proportion = (inducer_threshold - inducer_min) / (inducer_max - inducer_min)


# inducer_colormap = create_colormap_2colors(colors1_proportion)
#
# # plot_color_gradients(inducer_colormap)
#
# create_circle_animation(propagator,plt.get_cmap('Oranges'),'Calcium',file_prefix)
# create_circle_animation(inhibitor,plt.get_cmap('Reds'),'BMP',file_prefix)
# create_circle_animation(inducer,inducer_colormap,'cVg1',file_prefix)

