
import numpy as np

""" time parameters """
step_length = 0.5 # seconds
steps_per_second = 2
# model_runtime = 2 * 60 * 60 # seconds
model_runtime = 2 * 60 # seconds

number_of_steps = int(model_runtime / step_length)

steps_before_remove_streak = int(np.ceil(number_of_steps*0.5))

""" space parameters """
number_of_cells = 751
streak_number_of_cells = int(60 / 2)
max_cell_id = number_of_cells - 1

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
inducer_k = np.power(0.55,4)
a_iv = 1
a_vv = 1
inducer_decay = 0.001
inducer_diffusion = 0.0004
inducer_basal = 0.0004

""" initial conditions : inducer """
inducer_high_conc = 2.5
inducer_low_conc = 1
inducer_threshold = 2

""" initial conditions : inhibitor """
inhibitor_high_conc = 8
inhibitor_low_conc = 3.5

