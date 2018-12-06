import numpy as np
import os
import subprocess
from . import eas as eas
from . import trainingimages as trainingimages
# from . import trainingimages
import time


def is_exe(filename):
    return os.path.isfile(filename) and os.access(filename, os.X_OK)


# %% FUNCTIONS FOR PLOTTING OF EAS FILES


# %% THE CLASS
class mpslib:

    def __init__(self, parameter_filename='mps.txt', method='mps_genesim', debug_level=1, n_real=1, rseed=1,
                 out_folder='.', ti_fnam='ti.dat', simulation_grid_size=np.array([80, 40, 1]),
                 origin=np.zeros(3), grid_cell_size=np.array([1, 1, 1]), hard_data_fnam='hard.dat',
                 shuffle_simulation_grid=2, entropyfactor_simulation_grid=4, shuffle_ti_grid=1,
                 hard_data_search_radius=1,
                 soft_data_categories=np.arange(2), soft_data_fnam='soft.dat', n_threads=1, verbose_level=1,
                 template_size=np.array([8, 7, 1]), n_multiple_grids=3, n_cond=36, n_min_node_count=0,
                 n_max_ite=1000000, n_max_cpdf_count=1, distance_measure=1, distance_min=0, distance_pow=1,
                 max_search_radius=10000000, remove_gslib_after_simulation=1, gslib_combine=1, ti=np.empty(0), 
                 colocate_dimension=0):
        '''Initialize variables in Class'''

        mpslib_py_path, fn = os.path.split(__file__)
        self.mpslib_exe_folder = os.path.abspath(os.path.join(mpslib_py_path, '..', '..'))
        # self.mpslib_exe_folder = os.path.join(os.path.dirname('mpslib.py'),'..')
        self.blank_grid = None
        self.blank_val = np.NaN
        self.parameter_filename = parameter_filename.lower()  # change string to lower case
        self.method = method.lower()  # change string to lower case
        self.verbose_level = verbose_level
        self.verbose_level = debug_level
        
        self.remove_gslib_after_simulation = remove_gslib_after_simulation  # remove individual gslib fiels after simulation
        self.gslib_combine = gslib_combine  # combine realzations into one gslib file

        self.sim = None
        self.cde = None

        self.par = {}

        self.par['n_real'] = n_real
        self.par['rseed'] = rseed
        self.par['n_max_cpdf_count'] = n_max_cpdf_count
        self.par['out_folder'] = out_folder
        self.par['ti_fnam'] = ti_fnam.lower()  # change string to lower case
        self.par['simulation_grid_size'] = simulation_grid_size
        self.par['origin'] = origin
        self.par['grid_cell_size'] = grid_cell_size
        self.par['hard_data_fnam'] = hard_data_fnam.lower()  # change string to lower case
        self.par['shuffle_simulation_grid'] = shuffle_simulation_grid
        self.par['entropyfactor_simulation_grid'] = entropyfactor_simulation_grid
        self.par['shuffle_ti_grid'] = shuffle_ti_grid
        self.par['hard_data_search_radius'] = hard_data_search_radius
        self.par['soft_data_categories'] = soft_data_categories
        self.par['soft_data_fnam'] = soft_data_fnam.lower()  # change string to lower case
        self.par['n_threads'] = n_threads
        self.par['debug_level'] = debug_level

        # if the method is GENSIM, add package specific parameters
        if self.method == 'mps_genesim':
            self.par['n_cond'] = n_cond
            self.par['n_max_ite'] = n_max_ite
            self.par['n_max_cpdf_count'] = n_max_cpdf_count
            self.par['distance_measure'] = distance_measure
            self.par['distance_min'] = distance_min
            self.par['distance_pow'] = distance_pow
            self.par['colocate_dimension'] = colocate_dimension
            self.par['max_search_radius'] = max_search_radius
            # self.par['exe'] = 'mps_genesim'

        if self.method[0:10] == 'mps_snesim':
            self.par['template_size'] = template_size
            self.par['n_multiple_grids'] = n_multiple_grids
            self.par['n_min_node_count'] = n_min_node_count
            self.par['n_cond'] = n_cond
            # self.exe = 'mps_mps_snesim_tree'

        # Set verbose_level on eas as well
        eas.debug_level = self.par['debug_level'];
        # Check if on windows
        self.iswin = 0
        if (os.name == 'nt'):
            self.iswin = 1

    def which(self, program):
        '''
        self.which: Locate executable in the following order:
            1) current directotu
            2) MPSlib installation directory
            3) In system path
        '''
        # check if executable is found in the current folder
        if is_exe(program):
            return program
        else:
            # Check of executable is located in MPLSIB folder
            program_path_mpslib = os.path.abspath(os.path.join(self.mpslib_exe_folder, program))
            if is_exe(program_path_mpslib):
                return program_path_mpslib
            else:
                # Check of executable is located in system path
                for path in os.environ["PATH"].split(os.pathsep):
                    path = path.strip('"')
                    exe_file = os.path.join(path, program)
                    if is_exe(exe_file):
                        return exe_file

                print("################################################################")
                print("#")
                print("# mpslib: " + program + " not found !!!!")
                print("# PLEASE ADD THE MPSLIB PROGRAM TO THE SYSTEM PATH")
                print("# OR ADD THE LOCATION OF THE MPSLIB PROGRAMS TO THE SYSTEM PATH")
                print("#")
                print("#################################################################")

                return None

    # % Check parameter file setting using  GENESIM
    def parfile_check_genesim(self):

        self.par.setdefault('n_cond', 25)
        self.par.setdefault('n_max_ite', 10000)
        self.par.setdefault('n_max_cpdf_count', 10)
        self.par.setdefault('distance_measure', 1)
        self.par.setdefault('distance_min', 0)
        self.par.setdefault('distance_pow', 0)
        self.par.setdefault('max_search_radius', 1000000)

    # % Check parameter file setting using  SNESIM
    def parfile_check_snesim(self):

        self.par.setdefault('template_size', np.array([5, 5, 1]))
        self.par.setdefault('n_multiple_grids', 3)
        self.par.setdefault('n_min_node_count', 0)
        self.par.setdefault('n_cond', -1)

    # Change simulation method        
    def change_method(self, method='mps_genesim'):
        print(self.method)
        if method == 'mps_genesim':
            self.parfile_check_genesim()
            self.method = method
        if method[0:10] == 'mps_snesim':
            self.parfile_check_snesim()
            self.method = method

    def par_write(self):
        if self.method == 'mps_genesim':
            self.mps_genesim_par_write()
        elif self.method[0:10] == 'mps_snesim':
            self.mps_snesim_par_write()
        else:
            s = 'unknown method type: {}'.format(self.method)
            raise Exception(s)

    def mps_genesim_par_write(self):
        full_path = os.path.join(self.parameter_filename)

        if (self.verbose_level > 0):
            print('mpslib: Writing GENESIM type parameter file: ' + full_path);

        file = open(full_path, 'w')

        file.write('Number of realizations # %d\n' % self.par['n_real'])
        file.write('Random Seed (0 `random` seed) # %d \n' % self.par['rseed'])
        file.write('Maximum number of counts for condtitional pdf # %d\n' % self.par['n_max_cpdf_count'])
        file.write('Max number of conditional point # %d\n' % self.par['n_cond'])
        file.write('Max number of iterations # %d\n' % self.par['n_max_ite'])
        file.write("Distance measure [1:disc, 2:cont], minimum distance, power # %d %3.1f %3.1f\n" % (
        self.par['distance_measure'], self.par['distance_min'], self.par['distance_pow']))
        file.write('Colocate Dimension # %d \n' % self.par['colocate_dimension'])
        file.write("Max Search Radius for conditional data # %f \n" % self.par['max_search_radius'])
        file.write('Simulation grid size X # %d\n' % self.par['simulation_grid_size'][0])
        file.write('Simulation grid size Y # %d\n' % self.par['simulation_grid_size'][1])
        file.write('Simulation grid size Z # %d\n' % self.par['simulation_grid_size'][2])
        file.write('Simulation grid world/origin X # %g\n' % self.par['origin'][0])
        file.write('Simulation grid world/origin Y # %g\n' % self.par['origin'][1])
        file.write('Simulation grid world/origin Z # %g\n' % self.par['origin'][2])
        file.write('Simulation grid grid cell size X # %g\n' % self.par['grid_cell_size'][0])
        file.write('Simulation grid grid cell size Y # %g\n' % self.par['grid_cell_size'][1])
        file.write('Simulation grid grid cell size Z # %g\n' % self.par['grid_cell_size'][2])
        file.write('Training image file (spaces not allowed) # %s\n' % self.par['ti_fnam'])
        file.write('Output folder (spaces in name not allowed) # %s\n' % self.par['out_folder'])
        file.write(
            'Shuffle Simulation Grid path (1 : random, 0 : sequential) # %d\n' % self.par['shuffle_simulation_grid'])
        file.write('Shuffle Training Image path (1 : random, 0 : sequential) # %d\n' % self.par['shuffle_ti_grid'])
        file.write('HardData filename  (same size as the simulation grid)# %s\n' % self.par['hard_data_fnam'])
        file.write('HardData seach radius (world units) # %g\n' % self.par['hard_data_search_radius'])
        file.write('Softdata categories (separated by ;) # ')
        for c, ii in enumerate(self.par['soft_data_categories']):
            if c == 0:
                file.write('%d' % ii)
            else:
                file.write('; %d' % ii)
        file.write('\n')
        file.write('Soft datafilenames (separated by ; only need (number_categories - 1) grids) # %s\n' %
                   self.par['soft_data_fnam'])
        file.write('Number of threads (not currently used) # %d\n' % self.par['n_threads'])
        file.write(
            'Debug mode(2: write to file, 1: show preview, 0: show counters, -1: no ) # %d\n' % self.par['debug_level'])
        file.close()

    def mps_snesim_par_write(self):
        cwd = os.getcwd();
        full_path = os.path.join(self.parameter_filename)

        if (self.verbose_level > 0):
            print('mpslib: writing SNESIM type parameter file: ' + full_path);

        file = open(full_path, 'w')

        file.write('Number of realizations # %d\n' % self.par['n_real'])
        file.write('Random Seed (0 `random` seed) # %d \n' % self.par['rseed'])
        file.write('Number of mulitple grids (start from 0) # %d\n' % self.par['n_multiple_grids'])
        file.write('Min Node count (0 if not set any limit)# %d\n' % self.par['n_min_node_count'])
        file.write('Maximum number condtitional data (0: all) # %d\n' % self.par['n_cond'])
        if (self.par['template_size'].ndim==2):
            file.write('Search template size X # %d %d\n' % (self.par['template_size'][0, 0], self.par['template_size'][0, 1]))
            file.write('Search template size Y # %d %d\n' % (self.par['template_size'][1, 0], self.par['template_size'][1, 1]))
            file.write('Search template size Z # %d %d\n' % (self.par['template_size'][2, 0], self.par['template_size'][2, 1]))
        else:
            file.write('Search template size X # %d %d\n' % (self.par['template_size'][0], self.par['template_size'][0]))
            file.write('Search template size Y # %d %d\n' % (self.par['template_size'][1], self.par['template_size'][1]))
            file.write('Search template size Z # %d %d\n' % (self.par['template_size'][2], self.par['template_size'][2]))
            #file.write('Search template size X # %d\n' % self.par['template_size'][0])
            #file.write('Search template size Y # %d\n' % self.par['template_size'][1])
            #file.write('Search template size Z # %d\n' % self.par['template_size'][2])

        file.write('Simulation grid size X # %d\n' % self.par['simulation_grid_size'][0])
        file.write('Simulation grid size Y # %d\n' % self.par['simulation_grid_size'][1])
        file.write('Simulation grid size Z # %d\n' % self.par['simulation_grid_size'][2])
        file.write('Simulation grid world/origin X # %g\n' % self.par['origin'][0])
        file.write('Simulation grid world/origin Y # %g\n' % self.par['origin'][1])
        file.write('Simulation grid world/origin Z # %g\n' % self.par['origin'][2])
        file.write('Simulation grid grid cell size X # %g\n' % self.par['grid_cell_size'][0])
        file.write('Simulation grid grid cell size Y # %g\n' % self.par['grid_cell_size'][1])
        file.write('Simulation grid grid cell size Z # %g\n' % self.par['grid_cell_size'][2])
        file.write('Training image file (spaces not allowed) # %s\n' % self.par['ti_fnam'])
        file.write('Output folder (spaces in name not allowed) # %s\n' % self.par['out_folder'])
        file.write(
            'Shuffle Simulation Grid path (1 : random, 0 : sequential) # %d\n' % self.par['shuffle_simulation_grid'])
        file.write('Shuffle Training Image path (1 : random, 0 : sequential) # %d\n' % self.par['shuffle_ti_grid'])
        file.write('HardData filename  (same size as the simulation grid)# %s\n' % self.par['hard_data_fnam'])
        file.write('HardData seach radius (world units) # %g\n' % self.par['hard_data_search_radius'])
        file.write('Softdata categories (separated by ;) # ')
        for c, ii in enumerate(self.par['soft_data_categories']):
            if c == 0:
                file.write('%d' % ii)
            else:
                file.write('; %d' % ii)
        file.write('\n')
        file.write('Soft datafilenames (separated by ; only need (number_categories - 1) grids) # %s\n' %
                   self.par['soft_data_fnam'])
        file.write('Number of threads (not currently used) # %d\n' % self.par['n_threads'])
        file.write(
            'Debug mode(2: write to file, 1: show preview, 0: show counters, -1: no ) # %d\n' % self.par['debug_level'])
        file.close()

    def run(self, normal_msg='Elapsed time (sec)', silent=False):
        """
            *Description:*\n
            This function runs the mpslib executable from python.

            :param normal_msg: last line write to screen upon  successful completion
            :param silent: boolean, determines if modle outputs should be written to the screen

            :return: returns boolean success if the model run has been completed
            """
        import os
        success = False

        # check if TI is set, if not, get the one
        if (not os.path.isfile(self.par['ti_fnam'])) and (not hasattr(self, 'ti')):
            print('mpslib: Training image "%s" not found - USING DEFAULT!' % (self.par['ti_fnam']))
            self.ti = trainingimages.strebelle()[0]

        # write training image if set
        if hasattr(self, 'ti'):
            if self.ti.shape[0] > 0:
                # self.par['ti_fnam']='TrainingImage.dat'
                eas.write_mat(self.ti, self.par['ti_fnam'])

        # write soft data if set
        if hasattr(self, 'd_soft'):
            self.delete_soft_data()
            eas.write(self.d_soft, self.par['soft_data_fnam'])

        # write hard data if set
        if hasattr(self, 'd_hard'):
            self.delete_hard_data()
            eas.write(self.d_hard, self.par['hard_data_fnam'])

        # write parameter file
        self.par_write()


        exe_file = self.method
        if self.iswin:
            exe_file = exe_file + '.exe'

        exe_path = self.which(exe_file)

        if (self.verbose_level > -1):
            print(self.verbose_level)
            print("mpslib: trying to run '%s'  in folder '%s'"%(exe_file,exe_path))

        if exe_path is None:
            s = 'mpslib: The program {} does not exist or is not executable.'.format(exe_file)
            raise Exception(s)
            return -1
        else:
            if (self.verbose_level > -1):
                s = 'mpslib: Using the following executable to run the model: {}'.format(exe_path)
                print(s)


        if not os.path.isfile(os.path.join(self.parameter_filename)):
            s = 'The the mpslib input file does not exists: {}'.format(self.parameter_filename)
            raise Exception(s)

        if (self.verbose_level > 0):
            print("mpslib: trying to run  " + exe_path + " " + self.parameter_filename)

        if self.iswin:
            CREATE_NO_WINDOW = 0x08000000
            # stdout = subprocess.run([cmd,self.parameter_filename], stdout=subprocess.PIPE, creationflags=CREATE_NO_WINDOW)
            proc = subprocess.Popen([exe_path, self.parameter_filename], stdout=subprocess.PIPE,
                                    creationflags=CREATE_NO_WINDOW)
        else:
            proc = subprocess.Popen([exe_path, self.parameter_filename], stdout=subprocess.PIPE)

        while success == False:
            line = proc.stdout.readline()
            c = line.decode('utf-8', errors='ignore')
            if c != '':
                if normal_msg.lower() in c.lower():
                    success = True
                c = c.rstrip('\r\n')
                if not silent:
                    print('{}'.format(c))
            else:
                break

        # read on simulated data
        self.sim = []
        for i in range(0, self.par['n_real']):
            filename = '%s_sg_%d.gslib' % (self.par['ti_fnam'], i)
            OUT = eas.read(filename)
            if (self.verbose_level > 0):
                print('mpslib: Reading: %s' % (filename))
            self.sim.append(OUT['Dmat'])
            success = True
        
        # read conditional entropt simulated data
        self.cde = []
        for i in range(0, self.par['n_real']):
            filename = '%s_cde_%d.gslib' % (self.par['ti_fnam'], i)
            OUT = eas.read(filename)
            if (self.verbose_level > 0):
                print('mpslib: Reading: %s' % (filename))
            self.cde.append(OUT['Dmat'])
            success = True
            
            
            

        # combine gslib output files
        if (self.gslib_combine):
            import os
            n = 0;
            header = []
            for i in range(self.par['n_real']):
                cur_file = '%s_sg_%d.gslib' % (self.par['ti_fnam'], i)
                if os.path.isfile(cur_file):
                    if (self.verbose_level > 1):
                        print('mpslib: Merging %s' % cur_file)
                    d_cur = eas.read(cur_file)
                    if n == 0:
                        n_data = len(d_cur['D'])
                        Dall = np.zeros((n_data, self.par['n_real']))
                    header.append('Real #%d' % (i + 1))
                    Dall[:, i] = d_cur['D']
                    n = n + 1
        Dall = Dall[:, range(n)]

        # remove gslib output files
        if (self.remove_gslib_after_simulation):
            self.delete_gslib()

        filename_out = '%s.gslib' % (self.par['ti_fnam'])
        title = 'Realizations from %s - %s' % (self.method, d_cur['title'])
        eas.write(Dall, filename_out, title=title, header=header)

        return success

    def blank_sim(self):
        if hasattr(self, 'sim') is False:
            s = 'Simulation results not imported'
            raise Exception(s)

        if self.blank_grid is None:
            s = 'No blanking grid defined'
            raise Exception(s)

        nsim = len(self.sim)
        sim_grid_size = np.shape(self.sim[0])
        blank_grid_size = np.shape(self.blank_grid)

        if sim_grid_size != blank_grid_size:
            s = "Blank grid with size {} incompatiable with sim. grid size {} ".format(sim_grid_size, blank_grid_size)
            raise Exception(s)

        blk = self.blank_grid == 0

        self.simblk = np.copy(self.sim)
        for ii in np.arange(nsim):
            self.simblk[ii][blk == True] = self.blank_val

    # Loads existing simulation results into class without running mps algorithm
    def load_sim(self):
        # Check if simulation results are already imported
        if self.sim is not None:
            s = 'Simulation reaults already imported'
            raise Exception(s)

        self.sim = []
        for i in range(0, self.par['n_real']):
            filename = '%s_sg_%d.gslib' % (self.par['ti_fnam'], i)

            if os.path.isfile(filename) is False:
                s = '{} file not found'.format(filename)

            OUT = eas.read(filename)

            if (self.verbose_level > 0):
                print('mpslib: Reading: %s' % (filename))
            self.sim.append(OUT['Dmat'])
    
    # Loads existing conditional entropy
    def load_cde(self):
        # Check if simulation results are already imported
        if self.cde is not None:
            s = 'Conditional Entropy already imported'
            raise Exception(s)

        self.cde = []
        for i in range(0, self.par['n_real']):
            filename = '%s_cde_%d.gslib' % (self.par['ti_fnam'], i)

            if os.path.isfile(filename) is False:
                s = '{} file not found'.format(filename)

            OUT = eas.read(filename)

            if (self.verbose_level > 0):
                print('mpslib: Reading: %s' % (filename))
            self.cest.append(OUT['Dmat'])
            

    # delete gslib files
    def delete_gslib(self, remove_all_gslib=0):
        '''Removes GSLIB output files from folder

        :param remove_all_gslib: (def=0) removes aonly GSLIB files related to the cirrent run
                                 (1) removed ALL GSLIB files.
        :return:
        '''

        import os
        import glob

        if (remove_all_gslib) == 1:
            file_list = glob.glob('*.gslib')
        else:
            file_list = glob.glob('%s/%s*.gslib' % (self.par['out_folder'], self.par['ti_fnam']))
        for file in file_list:
            os.remove(os.path.join(self.par['out_folder'], file))
            if (self.verbose_level > 1):
                print('Removing {0}'.format(os.path.join(self.par['out_folder'], file)))

            # Delete hard data

    def delete_hard_data(self):
        import os
        if (os.path.isfile(self.par['hard_data_fnam'])):
            os.remove(self.par['hard_data_fnam'])

    # Delete soft data
    def delete_soft_data(self):
        import os
        if (os.path.isfile(self.par['soft_data_fnam'])):
            os.remove(self.par['soft_data_fnam'])

    # Delete locard filense
    def delete_local_files(self):
        self.delete_soft_data()
        self.delete_hard_data()
        self.delete_gslib()



    # Select random set of hard data
    def seq_gibbs_set_hard_data(self, step=0.1):
        '''
        Set hard data fro sequential GIbbs resimulation
        Currenyly only works in 2D,
        :param step: Percantage of parameters to resimulate
        :return: d_hard
        '''
        if (hasattr(self, 'x') == 0):
            self.x = np.arange(self.par['simulation_grid_size'][0]) * self.par['grid_cell_size'][0] + self.par['origin'][0]
        if (hasattr(self, 'y') == 0):
            self.y = np.arange(self.par['simulation_grid_size'][1]) * self.par['grid_cell_size'][1] + self.par['origin'][1]
        if (hasattr(self, 'z') == 0):
            self.z = np.arange(self.par['simulation_grid_size'][2]) * self.par['grid_cell_size'][2] + self.par['origin'][2]

        if (hasattr(self, 'sim') == 0):
            d_hard = np.arange(0)
            return d_hard

        D = self.sim[0]
        D = D.ravel()

        if (hasattr(self, 'xx') == 0):
            self.yy, self.xx = np.meshgrid(self.y, self.x, sparse=False, indexing='ij')
            self.xx = self.xx.ravel()
            self.yy = self.yy.ravel()

        N = self.xx.size

        N_hard = np.int(np.ceil((1 - step) * N))

        i_hard = np.random.choice(N, N_hard)

        d_hard = np.zeros((N_hard, 4))
        for i in np.arange(i_hard.size):
            d_hard[i, 0] = self.xx[i_hard[i]]
            d_hard[i, 1] = self.yy[i_hard[i]]
            d_hard[i, 2] = self.z[0]
            d_hard[i, 3] = D[i_hard[i]]
            #print('i=%d, i_hard=%d' % (i, i_hard[i]))

        self.d_hard = d_hard
        return d_hard

    # plot realizations (only in 2D so far)
    def plot_reals(self, nr=25, hardcopy=0, hardcopy_filename='reals'):
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import numpy as np

        nr = np.min((self.par['n_real'], nr))
        nsp = int(np.ceil(np.sqrt(nr)))

        fig = plt.figure(1)
        sp = gridspec.GridSpec(nsp, nsp, wspace=0.1, hspace=0.1)
        plt.set_cmap('hot')
        for i in range(0, nr):
            ax1 = plt.Subplot(fig, sp[i])
            fig.add_subplot(ax1)
            plt.imshow(self.sim[i], interpolation='none')
            plt.title("Real %d" % (i + 1))

        fig.suptitle(self.method + ' - ' + self.parameter_filename, fontsize=16)
        if (hardcopy):
            plt.savefig(hardcopy_filename)

        plt.show(block=False)

    # plot etypes (only in 2D so far)
    def plot_etype(self):

        import matplotlib.pyplot as plt
        import numpy as np
        import os
        from scipy import stats

        # read soft data ('check if it exist')
        use_soft = 0
        if (os.path.isfile(self.par['soft_data_fnam'])):
            d_soft = eas.read(self.par['soft_data_fnam'])
            use_soft = 1

        # read hard data ('check if it exist')
        use_hard = 0
        if (os.path.isfile(self.par['hard_data_fnam'])):
            d_hard = eas.read(self.par['hard_data_fnam'])
            use_hard = 1

        # compute Etype
        emean = np.mean(self.sim, axis=0)
        estd = np.std(self.sim, axis=0)
        emode = stats.mode(self.sim, axis=0)
        emode = emode[0][0]

        # plot the Etypes
        fig = plt.figure(2)
        fig.clf()
        plt.set_cmap('hot')
        plt.subplot(1, 3, 1)
        im = plt.imshow(emean, zorder=-1)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        if (use_hard):
            plt.plot(d_hard['D'][:, 0], d_hard['D'][:, 1], "k.", MarkerSize=25, zorder=0)
            plt.scatter(d_hard['D'][:, 0], d_hard['D'][:, 1], c=d_hard['D'][:, 3], s=25, zorder=1)
        if (use_soft):
            plt.plot(d_soft['D'][:, 0], d_soft['D'][:, 1], ".", color=((.4, .4, .4)), MarkerSize=25, zorder=0)
            plt.scatter(d_soft['D'][:, 0], d_soft['D'][:, 1], c=d_soft['D'][:, 3], s=25, zorder=2)
        plt.title('Etype Mean')

        plt.subplot(1, 3, 2)
        im = plt.imshow(estd)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        plt.title('Etype Std')

        plt.subplot(1, 3, 3)
        im = plt.imshow(emode)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        plt.title('Etype Mode')

        # plt.savefig("soft_ti_example_%s_%s.png" % (O1.method,O1.par['ti_fnam']), dpi=600)
        fig.suptitle(self.method + ' - ' + self.parameter_filename, fontsize=16)
        plt.show(block=False)


'''
    def to_file(self):
        nreal = self.nreal
        outname = self.ti_fnam + '_all.gslib'
        f = open(outname, 'w')
        f.write('%d %d %d\n' % (self.simulation_grid_size[0], self.simulation_grid_size[1], self.simulation_grid_size[2]))
        if nreal == 1:
            f.write('%d\n' % nreal)
            f.write('%s\n' % 'real1')
            for ii in np.arange(self.simulation_grid_size[0] * self.simulation_grid_size[1] * self.simulation_grid_size[2]):
                f.write('%d\n' % self.simulations[ii])
            f.close()
        else:
            f.write('%d\n' % nreal)
            for ii in np.arange(nreal):
                f.write('%s%d\n' % ('real', ii))

            for ii in np.arange(self.simulation_grid_size[0] * self.simulation_grid_size[1] * self.simulation_grid_size[2]):
                for jj in np.arange(nreal):
                    f.write('%d ' % self.simulations[ii, jj])
                f.write('\n')
        f.close()

'''
