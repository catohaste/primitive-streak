import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import style
import time
import pickle

def find_max_y_value(inducer, inhibitor, propagator):
    max_tmp = max(inducer.max(),inhibitor.max())
    max_y_value = math.ceil(max(max_tmp,propagator.max()))
    
    return max_y_value

def init():
    line_propagator.set_data([], [])
    line_inducer.set_data([], [])
    line_inhibitor.set_data([], [])
    return line_propagator,line_inducer,line_inhibitor,
    
def animate(i):
    x = list(range(number_of_cells))
    y_propagator = propagator[100*i,:]
    y_inducer = inducer[100*i,:]
    y_inhibitor = inhibitor[100*i,:]
    line_propagator.set_data(x, y_propagator)
    line_inducer.set_data(x, y_inducer)
    line_inhibitor.set_data(x, y_inhibitor)
    return line_propagator,line_inducer,line_inhibitor,
    
step_length = 1 # seconds

number_of_cells = 1000
model_runtime = 6 * 60 * 60 # seconds
transcription_delay_time = 15 * 60 # seconds

delay_steps = int(transcription_delay_time / step_length)
number_of_steps = int(model_runtime / step_length)

# model parameters
saturation_inducer = 5
saturation_inhibitor = 5
saturation_propagator = 10

a_vv = 1
a_iv = 1
a_vi = 1
a_pi = 1
a_ii = 1
a_vp = 1
a_pp = 1

alpha_inducer = np.power(30,4)
alpha_inhibitor = np.power(30,4)
alpha_propagator = np.power(30,4)

decay_inducer = 0.4
decay_inhibitor = 0.4
decay_propagator = 0.8

diffusion_inducer = 0.2
diffusion_inhibitor = 0.2
diffusion_propagator = 0.8

inducer_source = (list(range(number_of_cells)) - ((np.ones(number_of_cells)-1)*0.5))*0.001

# initiate variables
inducer = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
inhibitor = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
propagator = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)

# inducer initial conditions
inducer_number_of_starting_cells = 10
inducer_high_conc = 125
inducer_low_conc = 35

inducer[:delay_steps,:] = inducer_low_conc
inducer[:delay_steps,:inducer_number_of_starting_cells] = inducer_high_conc
inducer[:delay_steps,-inducer_number_of_starting_cells:] = inducer_high_conc

# set initial gradient of inhibitor as (-x^2) type function
inhibitor_high_conc = 85
inhibitor_low_conc = 52
max_cell_id = number_of_cells - 1
tmp_a = ((max_cell_id/2)**2 - max_cell_id) / (inhibitor_high_conc - inhibitor_low_conc)
tmp_b = (-max_cell_id) * (tmp_a)
for j in range(delay_steps):
    for i in range(number_of_cells):
        inhibitor[j,i] = tmp_a * i**2 + tmp_b*i + inhibitor_low_conc

print("Number of steps: ",number_of_steps)
print("Delay steps: ",delay_steps,"\n")

for step in range(delay_steps,number_of_steps):
    for cell in range(number_of_cells):
    
        # set boundary conditions for strip
        if cell == 0 : 
            cell_minus_one = number_of_cells-1
            cell_plus_one = cell+1
        elif cell == number_of_cells-1:
            cell_minus_one = cell-1
            cell_plus_one = 0
        else:
            cell_minus_one = cell-1
            cell_plus_one = cell+1

        inducer[step,cell] = inducer[step-1,cell] + step_length * ( inducer_source[cell] + ( (saturation_inducer * np.power(inducer[step-delay_steps,cell],4)) / ( alpha_inducer + (a_vv * np.power(inducer[step-delay_steps,cell],4)) + (a_iv * np.power(inhibitor[step-delay_steps,cell],4)) ) ) - (decay_inducer * inducer[step-1,cell]) + (diffusion_inducer * ( inducer[step-1,cell_minus_one] - 2 * inducer[step-1,cell] + inducer[step-1,cell_plus_one] ) ) )
        
        inhibitor[step,cell] = inhibitor[step-1,cell] + step_length * ( ( (saturation_inhibitor*(np.power(inhibitor[step-delay_steps,cell],4) + np.power(propagator[step-delay_steps,cell],4)) ) / ( alpha_inhibitor + (a_vi*np.power(inducer[step-delay_steps,cell],4)) + (a_pi*np.power(propagator[step-delay_steps,cell],4)) + (a_ii*np.power(inhibitor[step-delay_steps,cell],4)) ) ) - (decay_inhibitor * inhibitor[step-1,cell]) + (diffusion_inhibitor * ( inhibitor[step-1,cell_minus_one] - 2 * inhibitor[step-1,cell] + inhibitor[step-1,cell_plus_one] ) ) )
    
        propagator[step,cell] = propagator[step-1,cell] + step_length * ( ( (saturation_propagator*(np.power(inducer[step-delay_steps,cell],4) + np.power(propagator[step-1,cell],4)) ) / ( alpha_propagator + (a_vp*np.power(inducer[step-delay_steps,cell],4)) + (a_pp*np.power(propagator[step-1,cell],4)) ) ) - (decay_propagator * propagator[step-1,cell]) + (diffusion_propagator * ( propagator[step-1,cell_minus_one] - 2 * propagator[step-1,cell] + propagator[step-1,cell_plus_one] ) ) )
        
    
        if inducer[step,cell] < 0:
            inducer[step,cell] = 0
        if inhibitor[step,cell] < 0:
            inhibitor[step,cell] = 0
        if propagator[step,cell] < 0:
            propagator[step,cell] = 0

print("Inducer: ",inducer[3000,369])
print("Inhibitor: ",inhibitor[3000,369])
print("Propagator: ",propagator[3000,369])

max_y_value = find_max_y_value(inducer, inhibitor, propagator)

style.use('fivethirtyeight')

fig = plt.figure()
ax = plt.axes(xlim=(-1, 1000),ylim=(-0.5,max_y_value),xlabel='Cell ID',ylabel='Concentration')
line_propagator, = ax.plot([], [],'o',color='purple',markersize = 0.5,label='Propagator')
line_inducer, = ax.plot([], [],'go',markersize = 0.5,label='Inducer')
line_inhibitor, = ax.plot([], [],'ro',markersize = 0.5,label='Inhibitor')

# Now add the legend with some customizations.
legend = ax.legend(loc='upper center',shadow=False,markerscale=3,numpoints=5)
for label in legend.get_texts():
    label.set_fontsize('small')

anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=216, blit=True)

anim.save('model_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

plt.show()

