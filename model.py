# import tests
from test_model import *

# standard imports
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import style
import time
import datetime
import os
import pickle
from shutil import copyfile
from decimal import Decimal

class Variables:
    def __init__(self, number_of_steps, number_of_cells):
        self.inducer = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
        self.inhibitor = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
        self.propagator = np.zeros( (number_of_steps,number_of_cells), dtype=int)
        self.propagator_pulse = np.zeros( (number_of_steps,number_of_cells), dtype=bool)
        self.propagator_rest = np.zeros( (number_of_steps,number_of_cells), dtype=bool)
        self.propagator_time = np.zeros( (number_of_steps,number_of_cells), dtype=int)
        
    def initialize_inducer(self,inducer_low_conc,inducer_high_conc,high_conc_cells):
        self.inducer[:,:] = inducer_low_conc
        for cell in high_conc_cells:  
            self.inducer[:,cell] = inducer_high_conc
            self.inducer[:,cell] = inducer_high_conc
        
    def remove_portion_of_embryo(self, current_step, portion_to_remove):
        # check current_step is not larger than total number of steps
        
        # check portion_to_remove is a list
        
        # check portion_to_remove is entirely within defined embryo
        
        for cell_to_remove in portion_to_remove:
            self.inducer[current_step - 1,cell_to_remove] = 0
            self.inhibitor[current_step - 1,cell_to_remove] = 0
            self.propagator[current_step - 1,cell_to_remove] = 0
            self.propagator_pulse[current_step - 1,cell_to_remove] = False
            self.propagator_rest[current_step - 1,cell_to_remove] = False
            self.propagator_time[current_step - 1,cell_to_remove] = 0
            
    def pickle_variables(self, pickle_directory, file_prefix):
        # check 'os' has been imported
        
        if not os.path.exists(pickle_directory):
            os.makedirs(pickle_directory)
        pickle.dump( self.inducer, open( pickle_directory+"/"+file_prefix+"_inducer.p", "wb" ) )
        pickle.dump( self.inhibitor, open( pickle_directory+"/"+file_prefix+"_inhibitor.p", "wb" ) )
        pickle.dump( self.propagator, open( pickle_directory+"/"+file_prefix+"_propagator.p", "wb" ) )


# import parameters
from params import *

# initiate variables
var = Variables(number_of_steps,number_of_cells)

inhibitor_space = np.zeros( (number_of_cells), dtype=np.float64)
inhibitor_transform = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)

# initialize inducer
var.initialize_inducer(inducer_low_conc, inducer_high_conc, streak_cells)

# initialize inhibitor
tmp_m1 = (inhibitor_high_conc - inhibitor_low_conc) / (max_cell_id / 2 )
tmp_c1 = inhibitor_low_conc
tmp_m2 = (inhibitor_low_conc - inhibitor_high_conc) / (max_cell_id / 2 )
tmp_c2 = (2 * inhibitor_high_conc) - inhibitor_low_conc
for j in range(1):
    for i in range(math.ceil(number_of_cells/2)):
        var.inhibitor[j,i] = (tmp_m1 * i) + tmp_c1
        inhibitor_space[i] = (tmp_m1 * i) + tmp_c1
    for i in range(math.ceil(number_of_cells/2),number_of_cells):
        var.inhibitor[j,i] = (tmp_m2 * i) + tmp_c2
        inhibitor_space[i] = (tmp_m2 * i) + tmp_c2

# initialize propagator
control_cells = streak_cells
var.propagator_pulse[0,control_cells] = True
var.propagator_time[0,control_cells] = pulse_time - 1
var.propagator[0,control_cells] = pulse_value

print("Number of steps: ",number_of_steps)
#print("Delay steps: ",delay_steps,"\n")

for step in range(1,number_of_steps):
    
    var.propagator_rest[step,:] = var.propagator_rest[step-1,:]
    var.propagator_pulse[step,:] = var.propagator_pulse[step-1,:]
    
    if step > steps_before_cut:
        var.remove_portion_of_embryo(step,cells_to_remove)
        
    for cell in range(number_of_cells):

        # set periodic boundary conditions for strip
        if cell == 0 :
            cell_minus_one = number_of_cells-1
            cell_plus_one = cell+1
        elif cell == number_of_cells-1:
            cell_minus_one = cell-1
            cell_plus_one = 0
        else:
            cell_minus_one = cell-1
            cell_plus_one = cell+1

        var.inducer[step,cell] = var.inducer[step-1,cell] + step_length * ( inducer_basal + ( (inducer_saturation * np.power(var.inducer[step-1,cell],1)) / ( inducer_k + np.power(a_vv * var.inducer[step-1,cell],1) +  np.power(a_iv * var.inhibitor[step-1,cell],1) ) ) - (inducer_decay * var.inducer[step-1,cell]) + (inducer_diffusion * ( var.inducer[step-1,cell_minus_one] - 2 * var.inducer[step-1,cell] + var.inducer[step-1,cell_plus_one] ) ) )
        
        inhibitor_transform[step-1,cell] = a_ii * (var.inhibitor[step-1,cell] + inhibitor_space[cell])

        var.inhibitor[step,cell] = var.inhibitor[step-1,cell] + step_length * ( ( (inhibitor_saturation*( np.power(inhibitor_transform[step-1,cell],4) + np.power(var.propagator[step-1,cell],4) ) ) / ( inhibitor_k + np.power(inhibitor_transform[step-1,cell],4) + np.power(var.propagator[step-1,cell],4) ) ) - (inhibitor_decay * var.inhibitor[step-1,cell]) + (inhibitor_diffusion * ( var.inhibitor[step-1,cell_minus_one] - 2 * var.inhibitor[step-1,cell] + var.inhibitor[step-1,cell_plus_one] ) ) )

        
        if var.propagator_rest[step-1,cell]:
            if var.propagator_time[step-1,cell] > 0:
                var.propagator_time[step,cell] = var.propagator_time[step-1,cell] - 1
            else:
                var.propagator_rest[step,cell] = False
        elif var.propagator_pulse[step-1,cell]:
            if var.propagator_time[step-1,cell] > 0:
                var.propagator_time[step,cell] = var.propagator_time[step-1,cell] - 1
            else:
                var.propagator_pulse[step,cell] = False
                var.propagator_rest[step,cell] = True
                var.propagator_time[step,cell] = rest_time - 1
        elif var.propagator_pulse[step-1,cell_minus_one] or var.propagator_pulse[step-1,cell_plus_one]:
            var.propagator_pulse[step,cell] = True
            var.propagator_time[step,cell] = pulse_time - 1
        elif var.inducer[step-1,cell] > inducer_threshold:
            var.propagator_pulse[step,cell] = True
            var.propagator_time[step,cell] = pulse_time - 1
            
        var.propagator[step,cell] = (
                var.propagator_pulse[step,cell] * pulse_value
            ) + (
                var.propagator_rest[step,cell] * (
                    decay_time - (rest_time - var.propagator_time[step,cell])
                ) * propagator_decay_rate
            )
            

        if var.inducer[step,cell] < 0:
            var.inducer[step,cell] = 0
        if var.inhibitor[step,cell] < 0:
            var.inhibitor[step,cell] = 0
        if var.propagator[step,cell] < 0:
            var.propagator[step,cell] = 0
            

# print("Inducer: ",var.inducer[3000,369])
# print("Inhibitor: ",var.inhibitor[3000,369])
# print("Propagator: ",var.propagator[3000,369])

# pickle data
now = datetime.datetime.now()
date_string = now.strftime("%Y-%m-%d_%H")+now.strftime("%M")
file_prefix = 'changes' # date_string
pickle_dir = "output/pickles/"+file_prefix
var.pickle_variables(pickle_dir, file_prefix)

# # save parameter values
copyfile("params.py", "output/params/"+file_prefix+"_params.py")


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

# # draw inhibitor threshold line
# inhibitor_threshold = 1
# inhibitor_threshold_array = inhibitor_threshold * np.ones((number_of_cells+1,1),dtype=int)
# line_inhibitor_threshold, = ax_protein.plot(range(0,number_of_cells+1), inhibitor_threshold_array,':', color=color_inhibitor, linewidth = 1 ,label='BMP threshold')
# # threshold_text = ax_protein.text(0.65*number_of_cells,threshold,'approx. threshold',fontsize='x-small')

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
    y_propagator = var.propagator[sample_rate*i,:]
    y_inducer = var.inducer[sample_rate*i,:]
    y_inhibitor = var.inhibitor[sample_rate*i,:]
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

