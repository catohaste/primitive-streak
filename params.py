
import numpy as np
import datetime

""" file parameters """
now = datetime.datetime.now()
date_string = now.strftime("%Y-%m-%d_%H")+now.strftime("%M")
file_prefix = 'half_cut' # date_string
pickle_dir = "output/pickles/"+file_prefix

""" time parameters """
step_length = 0.5 # seconds
steps_per_second = 2
model_runtime = 2 * 60 * 60 # seconds

number_of_steps = int(model_runtime / step_length)

""" space parameters """
number_of_cells = 751
max_cell_id = number_of_cells - 1
streak_cells = np.union1d(np.arange(0,30), np.arange(number_of_cells - 30 ,number_of_cells))

""" experiment parameters """
steps_before_cut = int(np.ceil(number_of_steps*0.5))
half_no_of_cells_to_remove = int(np.floor(number_of_cells*0.25))
cells_to_remove = np.union1d(np.arange(0,half_no_of_cells_to_remove), np.arange(number_of_cells - half_no_of_cells_to_remove ,number_of_cells))

""" propagator pulse parameters """
pulse_value = 800 # nano molar
pulse_time = 5 * steps_per_second # steps
rest_time = 180 * steps_per_second # steps
decay_time = 50 * steps_per_second # steps
propagator_decay_rate = pulse_value / decay_time # per step

""" equation parameters : inhibitor """
inhibitor_saturation = 0.011 # micro meters per second
inhibitor_k = np.power(400,4) # nano molar
a_ii = 40 # no units
inhibitor_decay = 0.001
inhibitor_diffusion = 0.0004
inhibitor_m = 0.000022
inhibitor_c = 1

""" equation parameters : inducer """
inducer_saturation = 0.005
inducer_k = np.power(10,4)
a_iv = 6
a_vv = 2
inducer_decay = 0.001
inducer_diffusion = 0.0004
inducer_basal = 0.0008

""" initial conditions : inducer """
inducer_high_conc = 10
inducer_low_conc = 3
inducer_threshold = 8

""" initial conditions : inhibitor """
inhibitor_high_conc = 8
inhibitor_low_conc = 5

plot_params = {
    'number_of_steps': number_of_steps,
    'number_of_cells': number_of_cells,
    'step_length': step_length,
    'inducer_threshold': inducer_threshold,
    'file_prefix': file_prefix
}

