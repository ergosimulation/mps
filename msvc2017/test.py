# mpslib_example_1.py
# example of using mpslib

#%%
import mpslib as mps
import matplotlib.pyplot as plt
import os
import numpy as np
plt.ion()
plt.set_cmap('hot')

#%% MPS_SNESIM_TREE

O1 = mps.mpslib(method='mps_genesim', n_real = 1)
O1.par['n_max_cpdf_count']=1e+9 # ENESIM
#O1.par['n_max_cpdf_count']=30 # ENESIM
O1.par['simulation_grid_size']=[100,80,1]
O1.par['ti_fnam']='ti.dat'
O1.par['rseed']=0
#O1.mpslib_exe_folder=os.path.join(os.getcwd(),'x64','Release')
O1.mpslib_exe_folder=os.path.join(os.getcwd(),'..')
O1.parameter_filename = 'mps_genesim2.txt'
O1.par['debug_level']=1
#O1.par['n_cond']=36
O1.par['n_cond']=25
#O1.par['n_cond']=81
O1.remove_gslib_after_simulation=0;

O1.run()

#%%

plt.figure(1)
plt.subplot(1,3,1)
plt.imshow(O1.sim[0])
plt.subplot(1,3,2)
plt.imshow(O1.cde[0])
plt.show()
H=np.sum(O1.cde[0])
print("Entropy = %5.2f" % (H) )

