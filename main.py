# import tests
from test_model import *

# standard imports
import math
import numpy as np
from shutil import copyfile

# import model
from model import *
from params import *
from plots import *


# initiate variables
var = Variables(number_of_steps,number_of_cells)

inhibitor_space = np.zeros( (number_of_cells), dtype=np.float64)
inhibitor_transform = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)

# initialize inducer
var.initialize_inducer_stageXII(inducer_low_conc, inducer_high_conc)

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
control_cells = np.zeros( (number_of_cells), dtype=bool)
for cell in range(number_of_cells):
    if var.inducer[0,cell] > inducer_threshold:
        control_cells[cell] = True
var.propagator_pulse[0,control_cells] = True
var.propagator_time[0,control_cells] = pulse_time - 1
var.propagator[0,control_cells] = pulse_value

# plot_initial_conditions(var,plot_params)

print("Number of steps: ",number_of_steps)

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

        # inhibitor_transform[step-1,cell] = a_ii * (var.inhibitor[step-1,cell] + inhibitor_space[cell])
        inhibitor_transform[step-1,cell] = a_ii * var.inhibitor[step-1,cell]

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


# pickle data
var.pickle_variables(pickle_dir, file_prefix)

# save parameter values
copyfile("params.py", "output/params/"+file_prefix+"_params.py")


show_model_animation(var,plot_params)


