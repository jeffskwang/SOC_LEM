import numpy as np
import importlib
import sys
import os
parameters = importlib.import_module(sys.argv[1])
globals().update(parameters.__dict__)

def pad_rasters(input_file,pad_mode):
        s = np.loadtxt(input_file,skiprows=6)
        key = ['' for i in range(0,6)]
        val = ['' for i in range(0,6)]
        newval = ['' for i in range(0,6)]
        
        with open(input_file, 'r') as f:
                for i in range(0,6):
                    key[i], val[i] = f.readline().split()
        
        newval[0] = str(int(val[0]) + 2)
        newval[1] = str(int(val[1]) + 2)
        newval[2] = str(int(float(val[2]) - float(val[4])+0.5))
        newval[3] = str(int(float(val[3]) - float(val[4])+0.5))
##        newval[2] = str(float(val[2]) - float(val[4]))
##        newval[3] = str(float(val[3]) - float(val[4]))
        newval[4] = val[4]
        newval[5] = val[5]

        z = np.zeros((s.shape[0] +2,s.shape[1] +2))

        z[:,:]=-9999.
        z[1:-1,1:-1] = s 
        if pad_mode == 0:
                np.savetxt(results_folder+'\\results\\input\\padded_'+os.path.basename(input_file),s,delimiter='\t',newline='\n'\
                   ,header= key[0]+ '\t' + val[0]+'\n'+\
                   key[1]+ '\t' + val[1]+'\n'+\
                   key[2]+ '\t' + val[2]+'\n'+\
                   key[3]+ '\t' + val[3]+'\n'+\
                   key[4]+ '\t' + val[4]+'\n'+\
                   key[5]+ '\t' + val[5]\
                   , comments='')
        elif pad_mode == 1:
                np.savetxt(results_folder+'\\results\\input\\padded_'+os.path.basename(input_file),z,delimiter='\t',newline='\n'\
                   ,header= key[0]+ '\t' + newval[0]+'\n'+\
                   key[1]+ '\t' + newval[1]+'\n'+\
                   key[2]+ '\t' + newval[2]+'\n'+\
                   key[3]+ '\t' + newval[3]+'\n'+\
                   key[4]+ '\t' + newval[4]+'\n'+\
                   key[5]+ '\t' + newval[5]\
                   , comments='')

