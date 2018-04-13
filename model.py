from test_model import *
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import style
import time
import pickle
    
step_length = 1 # seconds

number_of_cells = 1001
model_runtime = 1 * 60 * 60 # seconds
transcription_delay_time = 15 * 60 # seconds

delay_steps = int(transcription_delay_time / step_length)
number_of_steps = int(model_runtime / step_length)

pulse_time = 3
rest_time = 10

pulse_value = rest_time
pulse_interval = 50

pulse = np.arange(0,number_of_steps/2,pulse_interval)

# model parameters
saturation_inducer = 5
saturation_inhibitor = 0.1
saturation_propagator = 0.7

a_vv = 1
a_iv = 1
a_vi = 1
a_pi = 1
a_ii = 1
a_vp = 1
a_pp = 1

alpha_inducer = np.power(30,1)
alpha_inhibitor = np.power(2,1)
alpha_propagator = np.power(70,4)

decay_inducer = 0.6
decay_inhibitor = 0.001
decay_propagator = 0.6

diffusion_inducer = 0.2
diffusion_inhibitor = 0
diffusion_propagator = 0.4

inducer_source = (list(range(number_of_cells)) - ((np.ones(number_of_cells)-1)*0.5))*0.001

# initiate variables
inducer = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
inhibitor = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
propagator_pulse = np.zeros( (number_of_steps,number_of_cells), dtype=bool)
propagator_rest = np.zeros( (number_of_steps,number_of_cells), dtype=bool)
propagator_time = np.zeros( (number_of_steps,number_of_cells), dtype=int)
propagator = np.zeros( (number_of_steps,number_of_cells), dtype=int)

# inducer initial conditions
inducer_number_of_starting_cells = 20
inducer_high_conc = 1
inducer_low_conc = -1

start_steps = int(np.ceil(number_of_steps*0.5))

inducer[:,:] = inducer_low_conc
inducer[:start_steps,:inducer_number_of_starting_cells] = inducer_high_conc
inducer[:start_steps,-inducer_number_of_starting_cells:] = inducer_high_conc

# inhibitor initial conditions
inhibitor_high_conc = 30
inhibitor_low_conc = 15
max_cell_id = number_of_cells - 1
tmp_m1 = (inhibitor_high_conc - inhibitor_low_conc) / (max_cell_id / 2 )
tmp_c1 = inhibitor_low_conc
tmp_m2 = (inhibitor_low_conc - inhibitor_high_conc) / (max_cell_id / 2 )
tmp_c2 = (2 * inhibitor_high_conc) - inhibitor_low_conc
for j in range(delay_steps):
    for i in range(math.ceil(number_of_cells/2)):
        inhibitor[j,i] = (tmp_m1 * i) + tmp_c1
    for i in range(math.ceil(number_of_cells/2),number_of_cells):
        inhibitor[j,i] = (tmp_m2 * i) + tmp_c2

# propagator initial conditions
control_cell = 0
propagator_pulse[0,control_cell] = True
propagator_time[0,control_cell] = pulse_time - 1
propagator[0,control_cell] = pulse_value - 1

print("Number of steps: ",number_of_steps)
#print("Delay steps: ",delay_steps,"\n")

for step in range(1,number_of_steps):
    
    propagator_rest[step,:] = propagator_rest[step-1,:]
    propagator_pulse[step,:] = propagator_pulse[step-1,:]
    
    if step in pulse:
        propagator_pulse[step-1,control_cell] = True
        propagator_time[step-1,control_cell] = pulse_time - 1
        propagator[step-1,control_cell] = pulse_value - 1
        
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

        #inducer[step,cell] = inducer[step-1,cell] + step_length * ( inducer_source[cell] + ( (saturation_inducer * np.power(inducer[step-delay_steps,cell],1)) / ( alpha_inducer + (a_vv * np.power(inducer[step-delay_steps,cell],1)) + (a_iv * np.power(inhibitor[step-delay_steps,cell],1)) ) ) - (decay_inducer * inducer[step-1,cell]) + (diffusion_inducer * ( inducer[step-1,cell_minus_one] - 2 * inducer[step-1,cell] + inducer[step-1,cell_plus_one] ) ) )

        inhibitor[step,cell] = inhibitor[step-1,cell] + step_length * ( ( (saturation_inhibitor*( np.power(propagator[step-1,cell],1)) ) / ( alpha_inhibitor + (a_pi*np.power(propagator[step-1,cell],1)) ) ) - (decay_inhibitor * inhibitor[step-1,cell]) + (diffusion_inhibitor * ( inhibitor[step-1,cell_minus_one] - 2 * inhibitor[step-1,cell] + inhibitor[step-1,cell_plus_one] ) ) )

        
        if propagator_rest[step-1,cell]:
            if propagator_time[step-1,cell] > 0:
                propagator_time[step,cell] = propagator_time[step-1,cell] - 1
            else:
                propagator_rest[step,cell] = False
        elif propagator_pulse[step-1,cell]:
            if propagator_time[step-1,cell] > 0:
                propagator_time[step,cell] = propagator_time[step-1,cell] - 1
            else:
                propagator_pulse[step,cell] = False
                propagator_rest[step,cell] = True
                propagator_time[step,cell] = rest_time - 1
        elif propagator_pulse[step-1,cell_minus_one] or propagator_pulse[step-1,cell_plus_one]:
            propagator_pulse[step,cell] = True
            propagator_time[step,cell] = pulse_time - 1
            
        propagator[step,cell] = (
                propagator_pulse[step,cell] * pulse_value
            ) + (
                propagator_rest[step,cell] * propagator_time[step,cell]
            )
        
        
        # propagator[step,cell] = propagator[step-1,cell] + step_length * (
        #     (
        #         (
        #             saturation_propagator * (np.power(inducer[step-1,cell],4))
        #         ) / (
        #             alpha_propagator + (a_vp * np.power(inducer[step-1,cell],4))
        #         )
        #     ) - (decay_propagator * propagator[step-1,cell]) + (
        #         diffusion_propagator * (
        #             propagator[step-1,cell_minus_one] - 2 * propagator[step-1,cell] + propagator[step-1,cell_plus_one]
        #         )
        #     )
        # )

        # if inducer[step,cell] < 0:
        #     inducer[step,cell] = 0
        if inhibitor[step,cell] < 0:
            inhibitor[step,cell] = 0
        if propagator[step,cell] < 0:
            propagator[step,cell] = 0
            

#print("Inducer: ",inducer[3000,369])
print("Inhibitor: ",inhibitor[3000,369])
print("Propagator: ",propagator[3000,369])


style.use('fivethirtyeight')
fig = plt.figure()

max_yval_ax1 = 50

# set up left axes
ax1 = plt.axes(xlim=(-1, number_of_cells + 1),ylim=(0.0,max_yval_ax1),xlabel='Cell Location')
ax1.set_ylabel('Propagator conc.', color='purple')
ax1.tick_params('y', colors='purple')
plt.xticks(np.linspace(0,number_of_cells,num=3), ['post.','ant.','post.'])

# set up right axes
ax2 = ax1.twinx()
ax2.set_ylabel('Inhibitor conc.', color='red')
ax2.tick_params('y', colors='red')
ax2.set_ylim(0,150)
ax2.grid('off')

plt.tight_layout()

# draw inhibitor threshold line
threshold = 20
threshold_array = threshold * np.ones((number_of_cells+1,1),dtype=int)
line_threshold, = ax2.plot(range(-1,number_of_cells), threshold_array,'k-',linewidth = 0.75,label='Inhibitor threshold')
threshold_text = ax2.text(0.8*number_of_cells,threshold,'threshold',fontsize='small')

# Initiate data to be animated
line_propagator, = ax1.plot([], [],'-',color='purple',linewidth = 1.0,label='Propagator')
line_inducer, = ax1.plot([], [],'go',markersize = 0.5,label='Inducer')
line_inhibitor, = ax2.plot([], [],'ro',markersize = 0.5,label='Inhibitor')
time_text = ax1.text(0.8*number_of_cells,max_yval_ax1*0.95,[],fontsize='small')

# Define dummy lines for legend
legend_inhibitor = ax1.plot([],[],'ro',markersize = 0.5,label='Inhibitor')

# Add the legend
legend = ax1.legend(loc='upper center',shadow=False,markerscale=3,numpoints=5)
for label in legend.get_texts():
    label.set_fontsize('small')

    
def init():
    line_propagator.set_data([], [])
    line_inducer.set_data([], [])
    line_inhibitor.set_data([], [])
    time_text.set_text('')
    return time_text,line_propagator,line_inhibitor,line_inducer,
    
sample_rate = 10
def animate(i):
    sample_rate = 10
    x = list(range(number_of_cells))
    y_propagator = propagator[sample_rate*i,:]
    y_inducer = inducer[sample_rate*i,:]
    y_inhibitor = inhibitor[sample_rate*i,:]
    line_propagator.set_data(x, y_propagator)
    line_inducer.set_data(x, y_inducer)
    line_inhibitor.set_data(x, y_inhibitor)
    current_time = step_length*sample_rate*i
    time_string = 't = ' + str(current_time) + 's'
    time_text.set_text(time_string)
    return time_text,line_propagator,line_inhibitor,line_inducer
    
number_of_frames = int(np.ceil(model_runtime / sample_rate))

anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=number_of_frames, blit=True)

anim.save('model_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

plt.show()

