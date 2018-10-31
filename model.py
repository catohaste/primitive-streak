import numpy as np
import math
import os
import pickle

def gradient_and_intercept(x1,y1,x2,y2):
    gradient = (y2 - y1) / (x2 - x1)
    intercept = y1 - gradient*x1
    
    return gradient, intercept

class Variables:
    def __init__(self, number_of_steps, number_of_cells):
        self.inducer = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
        self.inhibitor = np.zeros( (number_of_steps,number_of_cells), dtype=np.float64)
        self.propagator = np.zeros( (number_of_steps,number_of_cells), dtype=int)
        self.propagator_pulse = np.zeros( (number_of_steps,number_of_cells), dtype=bool)
        self.propagator_rest = np.zeros( (number_of_steps,number_of_cells), dtype=bool)
        self.propagator_time = np.zeros( (number_of_steps,number_of_cells), dtype=int)
        
    def initialize_inducer_streak(self,inducer_low_conc,inducer_high_conc,high_conc_cells):
        self.inducer[:,:] = inducer_low_conc
        for cell in high_conc_cells:  
            self.inducer[:,cell] = inducer_high_conc
            
    def initialize_inducer_stageXII(self,inducer_low_conc,inducer_high_conc):
        number_of_cells = self.inducer.shape[1]
        self.inducer[:,:] = inducer_low_conc
        quarter_no_of_cells = int(np.floor(number_of_cells*0.25))
        half_no_of_cells_high_inducer = int(np.floor(number_of_cells*0.1))
        high_conc_cells_inducer = np.union1d(np.arange(0,half_no_of_cells_high_inducer), np.arange(number_of_cells - half_no_of_cells_high_inducer ,number_of_cells))
        left_graded_inducer = np.setdiff1d(np.arange(0,quarter_no_of_cells),np.arange(0,half_no_of_cells_high_inducer))
        right_graded_inducer = np.setdiff1d(np.arange(number_of_cells - quarter_no_of_cells, number_of_cells),np.arange(number_of_cells - half_no_of_cells_high_inducer ,number_of_cells))
        
        for cell in high_conc_cells_inducer:
            self.inducer[0,cell] = inducer_high_conc

        gradient1, intercept1 = gradient_and_intercept(quarter_no_of_cells - 1, inducer_low_conc, half_no_of_cells_high_inducer, inducer_high_conc)
        gradient2, intercept2 = gradient_and_intercept(number_of_cells - quarter_no_of_cells - 1, inducer_low_conc, number_of_cells - half_no_of_cells_high_inducer, inducer_high_conc)
        for cell in left_graded_inducer:
            self.inducer[0,cell] = gradient1 * cell + intercept1
        for cell in right_graded_inducer:
            self.inducer[0,cell] = gradient2 * cell + intercept2
            
            
    def initialize_inhibitor(self,inhibitor_low_conc,inhibitor_high_conc):
        number_of_cells = self.inhibitor.shape[1]
        max_cell_id = number_of_cells - 1
        tmp_m1,tmp_c1 = gradient_and_intercept(0,inhibitor_low_conc,max_cell_id/2,inhibitor_high_conc)
        tmp_m2,tmp_c2 = gradient_and_intercept(max_cell_id/2,inhibitor_high_conc,max_cell_id,inhibitor_low_conc)
        for j in range(1):
            for i in range(math.ceil(number_of_cells/2)):
                self.inhibitor[j,i] = (tmp_m1 * i) + tmp_c1
            for i in range(math.ceil(number_of_cells/2),number_of_cells):
                self.inhibitor[j,i] = (tmp_m2 * i) + tmp_c2
        
        
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
        
        
def inducerODE(step, cell_dict, variables, model_params):
    
    cell = cell_dict['cell']
    cell_minus_one = cell_dict['cell_minus_one']
    cell_plus_one = cell_dict['cell_plus_one']
    
    inducer_basal = model_params['inducer_basal']
    inducer_saturation = model_params['inducer_saturation']
    inducer_k = model_params['inducer_k']
    a_vv = model_params['a_vv']
    a_bv = model_params['a_bv']
    inducer_decay = model_params['inducer_decay']
    inducer_diffusion = model_params['inducer_diffusion']
    
    dV = ( inducer_basal + ( (inducer_saturation * np.power(variables.inducer[step,cell],4)) / ( inducer_k + np.power(a_vv * variables.inducer[step,cell],4) +  np.power(a_bv * variables.inhibitor[step,cell],4) ) ) - (inducer_decay * variables.inducer[step,cell]) + (inducer_diffusion * ( variables.inducer[step,cell_minus_one] - 2 * variables.inducer[step,cell] + variables.inducer[step,cell_plus_one] ) ) )
    
    return dV
    
def inhibitorODE(step, cell_dict, variables, model_params):
    
    cell = cell_dict['cell']
    cell_minus_one = cell_dict['cell_minus_one']
    cell_plus_one = cell_dict['cell_plus_one']
    
    inhibitor_basal = model_params['inhibitor_basal']
    inhibitor_saturation = model_params['inhibitor_saturation']
    inhibitor_k = model_params['inhibitor_k']
    a_bb = model_params['a_bb']
    a_cb = model_params['a_cb']
    inhibitor_decay = model_params['inhibitor_decay']
    inhibitor_diffusion = model_params['inhibitor_diffusion']
    
    dB = ( inhibitor_basal + ( (inhibitor_saturation*( np.power(variables.inhibitor[step,cell],4) + np.power(variables.propagator[step,cell],4) ) ) / ( inhibitor_k + np.power(a_bb * variables.inhibitor[step,cell],4) + np.power(a_cb * variables.propagator[step,cell],4) ) ) - (inhibitor_decay * variables.inhibitor[step,cell]) + (inhibitor_diffusion * ( variables.inhibitor[step,cell_minus_one] - 2 * variables.inhibitor[step,cell] + variables.inhibitor[step,cell_plus_one] ) ) )
    
    return dB
    
