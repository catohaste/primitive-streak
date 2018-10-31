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

# initialize inducer
var.initialize_inducer_stageXII(inducer_low_conc, inducer_high_conc)

# initialize inhibitor
var.initialize_inhibitor(inhibitor_low_conc, inhibitor_high_conc)

# initialize propagator
control_cells = np.zeros( (number_of_cells), dtype=bool)
for cell in range(number_of_cells):
    if var.inducer[0,cell] > inducer_threshold:
        control_cells[cell] = True
var.propagator_pulse[0,control_cells] = True
var.propagator_time[0,control_cells] = pulse_time - 1
var.propagator[0,control_cells] = pulse_value

# print(var.inducer[0,75] - var.inducer[0,76])
# print(var.inhibitor[0,1] - var.inhibitor[0,0])

# plot_initial_conditions(var,plot_params)

print("Number of steps: ",number_of_steps)

cell_dict = {
    'cell':0,
    'cell_minus_one':0,
    'cell_plus_one':0
}

for step in range(1,number_of_steps):

    var.propagator_rest[step,:] = var.propagator_rest[step-1,:]
    var.propagator_pulse[step,:] = var.propagator_pulse[step-1,:]

    # if step > steps_before_cut:
    #     var.remove_portion_of_embryo(step,cells_to_remove)

    for cell in range(number_of_cells):

        cell_dict['cell'] = cell

        # set periodic boundary conditions for strip
        if cell == 0 :
            cell_minus_one = number_of_cells-1
            cell_plus_one = cell+1
            cell_dict['cell_minus_one'] = number_of_cells-1
            cell_dict['cell_plus_one'] = cell+1
        elif cell == number_of_cells-1:
            cell_dict['cell_minus_one'] = cell-1
            cell_dict['cell_plus_one'] = 0
            cell_minus_one = cell-1
            cell_plus_one = 0
        else:
            cell_dict['cell_minus_one'] = cell-1
            cell_dict['cell_plus_one'] = cell+1
            cell_minus_one = cell-1
            cell_plus_one = cell+1

        var.inducer[step,cell] = var.inducer[step-1,cell] + step_length * inducerODE(step-1,cell_dict,var,model_params)

        var.inhibitor[step,cell] = var.inhibitor[step-1,cell] + inhibitorODE(step-1,cell_dict,var,model_params)

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


