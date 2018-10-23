import numpy as np
import os
import pickle

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

        x1 = quarter_no_of_cells - 1
        y1 = inducer_low_conc
        x2 = half_no_of_cells_high_inducer
        y2 = inducer_high_conc
        gradient = (y2 - y1) / (x2 - x1)
        intercept = y1 - gradient*x1
        for cell in left_graded_inducer:
            self.inducer[0,cell] = gradient * cell + intercept

        x1 = number_of_cells - quarter_no_of_cells - 1
        y1 = inducer_low_conc
        x2 = number_of_cells - half_no_of_cells_high_inducer
        y2 = inducer_high_conc
        gradient = (y2 - y1) / (x2 - x1)
        intercept = y1 - gradient*x1
        for cell in right_graded_inducer:
            self.inducer[0,cell] = gradient * cell + intercept
        
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
        
