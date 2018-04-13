import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import style

number_of_steps = 3000
number_of_cells = 1001

pulse_time = 3
rest_time = 10

pulse_value = rest_time
pulse_interval = 50

propagator_pulse = np.zeros( (number_of_steps,number_of_cells), dtype=bool)
propagator_rest = np.zeros( (number_of_steps,number_of_cells), dtype=bool)
propagator_time = np.zeros( (number_of_steps,number_of_cells), dtype=int)
propagator = np.zeros( (number_of_steps,number_of_cells), dtype=int)

control_cell = 0
propagator_pulse[0,control_cell] = True
propagator_time[0,control_cell] = pulse_time - 1
propagator[0,control_cell] = pulse_value - 1

pulse = np.arange(0,number_of_steps,pulse_interval)

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

max_y_value = 50 + 1

# set up figure and twin axes
fig = plt.figure()
ax1 = plt.axes(xlim=(-1, number_of_cells),ylim=(0.0,max_y_value),xlabel='Cell ID')
ax1.set_ylabel('propagator conc.', color='purple')
ax1.tick_params('y', colors='purple')

plt.tight_layout()

line_propagator, = ax1.plot([], [],'-',color='purple',markersize = 2,label='Propagator')
#line_inducer, = ax1.plot([], [],'go',markersize = 0.5,label='Inducer')
#line_inhibitor, = ax1.plot([], [],'ro',markersize = 0.5,label='Inhibitor')

time_text = ax1.text(0.8*number_of_cells,max_y_value*0.95,[],fontsize='small')

# Now add the legend with some customizations.
legend = ax1.legend(loc='upper center',shadow=False,markerscale=3,numpoints=5)
for label in legend.get_texts():
    label.set_fontsize('small')
    
def init():
    line_propagator.set_data([], [])
    #line_inducer.set_data([], [])
    #line_inhibitor.set_data([], [])
    time_text.set_text('')
    return time_text, line_propagator,#line_inducer,#line_inhibitor,

sample_rate = 10    
def animate(i):
    sample_rate = 10
    x = list(range(number_of_cells))
    y_propagator = propagator[int(np.floor(sample_rate*i)),:]
    #y_inducer = inducer[sample_rate*i,:]
    #y_inhibitor = inhibitor[100*i,:]
    line_propagator.set_data(x, y_propagator)
    #line_inducer.set_data(x, y_inducer)
    #line_inhibitor.set_data(x, y_inhibitor)
    current_time = int(np.floor(sample_rate*i))
    time_string = 't = ' + str(current_time) + 's'
    time_text.set_text(time_string)
    return time_text,line_propagator,#line_inducer,#line_inhibitor,
    
number_of_frames = int(np.floor(number_of_steps / sample_rate))

anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=number_of_frames, blit=True)

anim.save('bool_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

plt.show()


		
			
