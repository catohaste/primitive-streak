
import numpy as np
import matplotlib.pyplot as plt

from model import *
from plots import *
from params import *

inducer_critical = 8
inhibitor_critical = 6
propagator_critical = 130

inducer_threshold = 7.9

# initiate variables
var = Variables(number_of_steps,number_of_cells)

inhibitor_space = np.zeros( (number_of_cells), dtype=np.float64)
inhibitor_transform = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)

# initialize inducer
var.initialize_inducer_stageXII(inducer_critical, inducer_critical)

# initialize inhibitor
var.initialize_inhibitor(inhibitor_critical, inhibitor_critical)

# initialize propagator
control_cells = np.zeros( (number_of_cells), dtype=bool)
for cell in range(number_of_cells):
    if var.inducer[0,cell] > inducer_threshold:
        control_cells[cell] = True
var.propagator_pulse[0,control_cells] = True
var.propagator_time[0,control_cells] = pulse_time - 1
var.propagator[0,control_cells] = pulse_value

# plot_initial_conditions(var,plot_params)

cell_dict = {
    'cell':0,
    'cell_minus_one':0,
    'cell_plus_one':0
}

model_params = {
    'inhibitor_basal': inhibitor_basal,
    'inhibitor_saturation': inhibitor_saturation,
    'inhibitor_k': inhibitor_k,
    'a_ii': a_ii,
    'inhibitor_decay': inhibitor_decay,
    'inhibitor_diffusion': inhibitor_diffusion,

    'inducer_basal': inducer_basal,
    'inducer_saturation': inducer_saturation,
    'inducer_k': inducer_k,
    'a_iv': a_iv,
    'a_vv': a_vv,
    'inducer_decay': inducer_decay,
    'inducer_diffusion': inducer_diffusion
}

dimensions = (100,100,100)

min_axis = np.zeros((len(dimensions)))

inducer_saturation = np.linspace(0.0001,1,dimensions[0])
inducer_decay = np.linspace(0.0001,0.1,dimensions[1])
inducer_diffusion = np.linspace(0.0001,10,dimensions[2])

dV = np.zeros(dimensions)
min_dV = np.zeros(len(dimensions))

for k in range(dimensions[2]):
    for j in range(dimensions[1]):
        for i in range(dimensions[0]):
            model_params['inducer_saturation'] = inducer_saturation[i]
            model_params['inducer_decay'] = inducer_decay[j]
            model_params['inducer_diffusion'] = inducer_diffusion[k]
            dV[i,j,k] = inducerODE(0, cell_dict, var, model_params)
        
flat_dV = dV.flatten()

arg_min_dV = np.argmin(dV)
min_dV = flat_dV[arg_min_dV]

arg_max_dV = np.argmax(dV)
max_dV = flat_dV[arg_max_dV]

print('argmin = ' + str(arg_min_dV))
print('min = ' + str(min_dV))


# plt.plot(range(dimensions[0]),dV[:,50,50])
#
# plt.show()
#
# plt.plot(range(dimensions[1]),dV[50,:,50])
#
# plt.show()
#
# plt.plot(range(dimensions[2]),dV[50,50,:])
#
# plt.show()

