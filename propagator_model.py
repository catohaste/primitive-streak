import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import style
import time
import pickle
    
step_length = 0.01 # seconds

number_of_cells = 1
model_runtime = 300 # seconds
# transcription_delay_time = 15 * 60 # seconds

# delay_steps = int(transcription_delay_time / step_length)
number_of_steps = int(model_runtime / step_length)

# model parameters
inducer_source_value = (list(range(number_of_cells)) - ((np.ones(number_of_cells)-1)*0.5))*0.001

# initiate variables
inducer = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
ca_cytoplasm = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
ca_cell = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)

# inducer initial conditions
inducer_number_of_starting_cells = 1
inducer_start_conc = 0.0002
inducer_end_conc = 0.0
inducer_on_time = int(np.ceil(number_of_steps * 0.5))

inducer[:inducer_on_time,:inducer_number_of_starting_cells] = inducer_start_conc
inducer[:inducer_on_time,-inducer_number_of_starting_cells:] = inducer_start_conc

inducer[-inducer_on_time:,:inducer_number_of_starting_cells] = inducer_end_conc
inducer[-inducer_on_time:,-inducer_number_of_starting_cells:] = inducer_end_conc

ca_cytoplasm[0,:] = 0.001
ca_cell[0,:] = 0.02

# parameters
k_1 = 40.0 # s^(-1)
k_2 = 0.02 # s^(-1)
d_1 = 0.3 # M
d_2 = 0.4 # M
d_3 = 0.2 # M
d_p = 0.2 # M 
d_a = 0.4 # M
rho = 0.06 # micro meter^(-1)
alpha = 2
beta = 0.1
gamma = 0 # ********** WARNING ************
v_0 = 0.2 # micro M s^(-1)
v_3 = 9.0 # micro M s^(-1)
v_4 = 3.6 # micro M s^(-1)
v_c = 4.0 # micro M s^(-1)
K_0 = 0.06 # M ********** WARNING ************
K_4 = 0.12 # M
K_3 = 0.12 # M

print("Number of steps: ",number_of_steps)

for step in range(1,number_of_steps):
    for cell in range(number_of_cells):
    
        # set boundary conditions for strip
        # if cell == 0 :
#             cell_minus_one = number_of_cells-1
#             cell_plus_one = cell+1
#         elif cell == number_of_cells-1:
#             cell_minus_one = cell-1
#             cell_plus_one = 0
#         else:
#             cell_minus_one = cell-1
#             cell_plus_one = cell+1

        cell_minus_one = cell
        cell_plus_one = cell
    
        k_r = k_1 * (
            np.power(
                d_2 * inducer[step-1,cell] * ca_cytoplasm[step-1,cell] * (
                    (d_1 + inducer[step-1,cell]) / (d_3 + inducer[step-1,cell])
                )
            ,3) / ( 
                np.power(d_p + inducer[step-1,cell],3) * np.power(d_a + ca_cytoplasm[step-1,cell],3) * np.power(
                    d_2 * ((d_1 + inducer[step-1,cell]) / (d_3 + inducer[step-1,cell])) + ca_cytoplasm[step-1,cell]
                ,3)
            )
        ) + k_2
        
        cell_function = rho * (
            v_0 + v_c * (
                inducer[step-1,cell] / ( inducer[step-1,cell] + K_0 )
            ) - v_4 * (
                ca_cytoplasm[step-1,cell]**2 / ( K_4**2 + ca_cytoplasm[step-1,cell]**2 ) 
            )
        )
        
        
        ca_cell[step,cell] = ca_cell[step-1,cell] + step_length * (
            cell_function + gamma * (
                ca_cytoplasm[step,cell_plus_one] + ca_cytoplasm[step,cell_minus_one] - 2*ca_cytoplasm[step,cell]
            )
        )
        
        
        ca_cytoplasm[step,cell] = ca_cytoplasm[step-1,cell] + step_length * (
            cell_function + rho * alpha * (
                k_r*(1/beta)*(
                    ca_cell[step-1,cell] - (
                        ca_cytoplasm[step-1,cell] * (1-beta)
                    )
                ) - (
                   (v_3 * ca_cytoplasm[step,cell]**2) / (K_3**2 + ca_cytoplasm[step,cell]**2)
                )
            ) + gamma * (
                ca_cytoplasm[step,cell_plus_one] + ca_cytoplasm[step,cell_minus_one] - 2*ca_cytoplasm[step,cell]
            )
        )
        
        if inducer[step,cell] < 0:
            inducer[step,cell] = 0
        if ca_cytoplasm[step,cell] < 0:
            ca_cytoplasm[step,cell] = 0
        if ca_cell[step,cell] < 0:
            ca_cell[step,cell] = 0
            
    # if step % 100 == 0:
    #     print('Step: ',step)
    #     print(np.amin(inducer[step,:]),' < inducer < ',np.amax(inducer[step,:]))
    #     print(np.amin(inhibitor[step,:]),' < inhibitor < ',np.amax(inhibitor[step,:]))
    #     print(np.amin(ca_cytoplasm[step,:]),' < ca_cytoplasm < ',np.amax(ca_cytoplasm[step,:]))

print("Inducer: ",inducer[10,0])
print("Calcium_cytoplasm: ",ca_cytoplasm[10,0])

# def find_max_y_value(inducer, inhibitor, ca_cytoplasm):
#     max_tmp = max(inducer.max(),inhibitor.max())
#     max_y_value = math.ceil(max(max_tmp,ca_cytoplasm.max()))
#   return max_y_value

max_y_value = max(ca_cytoplasm.max(),ca_cell.max())

style.use('fivethirtyeight')

fig = plt.figure()
ax = plt.axes(xlim=(-1, 1),ylim=(0,max_y_value),xlabel='Cell ID',ylabel='Concentration')
line_ca_cytoplasm, = ax.plot([], [],'o',color='purple',markersize = 3,label='ca_cytoplasm')
line_ca_cell, = ax.plot([], [],'o',color='green',markersize = 3,label='ca_cell')

# Now add the legend with some customizations.
legend = ax.legend(loc='upper center',shadow=False,markerscale=1,numpoints=5)
for label in legend.get_texts():
    label.set_fontsize('small')
    
def init():
    time_string = 't = 0s'
    line_ca_cytoplasm.set_data([], [])
    line_ca_cell.set_data([], [])
    ax.text(0.75, 0.05, time_string,fontsize=20)
    return line_ca_cytoplasm, line_ca_cell,

def animate(i):
    time_string = 't = ' + str(i*100) + 's'
    x = list(range(number_of_cells))
    y_ca_cytoplasm = ca_cytoplasm[100*i,:]
    y_ca_cell = ca_cell[100*i,:]
    line_ca_cytoplasm.set_data(x, y_ca_cytoplasm)
    line_ca_cell.set_data(x, y_ca_cell)
    return line_ca_cytoplasm, line_ca_cell,

number_of_frames = int(np.floor(number_of_steps / 100))

anim = animation.FuncAnimation(fig, animate, init_func=init, interval=42, frames=number_of_frames, blit=True)

anim.save('propagator_animation.mp4', fps=24, extra_args=['-vcodec', 'libx264'])

plt.show()

