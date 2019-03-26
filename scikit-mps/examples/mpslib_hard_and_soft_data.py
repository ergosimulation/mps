'''
mpslib_hard_and_soft_data.py

Example of parsing hard and soft data to MPSLIB algorithms
'''

import mpslib as mps
import numpy as np

if __name__ == '__main__':
        
    
    O1=mps.mpslib(method='mps_snesim_tree')
    
    TI1, TI_filename1 = mps.trainingimages.strebelle(3, coarse3d=1)
    O1.ti=TI1
    
    #%%
    
    O1.par['n_multiple_grids']=2;
    O1.par['n_cond']=49
    O1.par['n_real']=100
    O1.par['debug_level']=-1
    O1.par['simulation_grid_size'][0]=35
    O1.par['simulation_grid_size'][1]=30
    
    # Set hard data
    d_hard = np.array([[ 25.,  10.,   0.,   1.],
           [  5.,  10.,   0.,   1.],
           [  7.,  10.,   0.,   1.],
           [ 14.,  10.,   0.,   1.],
           [ 15.,  10.,   0.,   0.]])
    O1.d_hard = d_hard
    O1.par['hard_data_fnam']='hard.dat'
    
    # Set soft data
    d_soft = np.array([[ 25.,  20.,   0.,   0.99],
           [  5.,  20.,   0.,   0.99],
           [  7.,  20.,   0.,   0.99],
           [ 15.,  20.,   0.,   0.01]])
    O1.d_soft = d_soft
    O1.par['soft_data_fnam']='soft.dat'
    
    #%%O1.remove_gslib_after_simulation=0;
    
    
    O1.par['n_threads']=20#
    O1.run_parallel()
    
    O1.plot_reals()
    O1.plot_etype()
