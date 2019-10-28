#importing libaries
import numpy as np
import os
import shutil
import time

#working directory
parent = 'C:\\Users\\jeffs\Desktop\\UMASS-local\\SOC_LEM'

#run name
runname='willis'

#initial condition
ini_file = 'willis_elev.asc'

#physical parameters
U = 0.0            # [m/yr]  uplift
K = 0.0001      # [1/yr] vertical erodibility constant
D = 0.1            #[m^2/yr] hillslope diffusion coefficient
La = 0.1           # [m] active layer
K_SOC = 0.1  # [m] SOC exponent, i.e. SOC[z] = C_SOC * exp(-z/K_SOC)
C_SOC = 0.05# [1/m] SOC coeffcieint, i.e. SOC[z] = C_SOC * exp(-z/K_SOC)

#numerical parameters
T = 150. # [yr] Simulation Time
dt = 1.0 # [yr] model timestep
hole_function = 0# 0 is off and 1 is on

#output parameters
dt_plot = 10. # [yr] plot timestep

#####HOUSEKEEPING#####
#make output folder
if os.path.isdir(parent+'\\results'):
    shutil.rmtree(parent+'\\results')
    time.sleep(1)
os.makedirs(parent+'\\results')
#load landlab components
from landlab.components import FlowAccumulator
from landlab.components import DepressionFinderAndRouter
from landlab.components import FastscapeEroder
from landlab.components import LinearDiffuser
from landlab.components import TaylorNonLinearDiffuser
from landlab import RasterModelGrid, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
from landlab.components import LateralEroder
from landlab.utils import structured_grid
from landlab.io import read_esri_ascii
#load initial condition
(grid, eta) = read_esri_ascii(parent + '\\input\\' + ini_file, name='topographic__elevation')
#define arrays
SOC_La = grid.add_ones('node','SOC_La')
SOC_La *= 0.05
SOC_transfer = grid.add_zeros('node','SOC_transfer')
#grid size and number of cells
dx = grid.dx
dy = grid.dy
nt = int(T / dt)
nt_plot = int(dt_plot/dt)
#output files
data_elevation=[]
data_area=[]
data_soc=[]
#saving parameters and initial condition
parameter_grid = grid.node_vector_to_raster(eta,flip_vertically=True)
np.savetxt(parent+'\\results\\'+ 'parameters_and_IC.asc',parameter_grid,delimiter='\t',newline='\n',header='runname\t' + runname+ '\ncellsx\t'+str(parameter_grid.shape[0])+'\ncellsy\t'+str(parameter_grid.shape[1])+'\ndx\t'+str(grid.dx)+'\ndy\t'+str(grid.dy)+'\nT\t'+str(T)+'\ndt\t'+str(dt)+'\nK\t'+str(K)+'\nD\t'+str(D)+'\nU\t'+str(U), comments='')
#boundary conditions #under construction (it should have mirrored boundaries on the hilltops
for edge in (grid.nodes_at_top_edge,grid.nodes_at_bottom_edge,grid.nodes_at_left_edge, grid.nodes_at_right_edge):
    grid.status_at_node[edge] = FIXED_VALUE_BOUNDARY
###import components###
#Stream Power Incision Model
sp = FastscapeEroder(grid, K_sp = K ,m_sp=0.5, n_sp=1.0)
#Diffusion
ld = LinearDiffuser(grid, linear_diffusivity= D, method = 'simple')
#Flow Accumulation
if hole_function == 0:
    fa = FlowAccumulator(grid,flow_director='FlowDirectorD8')
elif hole_function == 1:
    fa = FlowAccumulator(grid,flow_director='FlowDirectorD8',depression_finder=DepressionFinderAndRouter)
#start time
start_time = time.time()
#########################

#####FUNCTIONS#####
def soil_and_SOC_transport(eta,SOC_La):
    dzdx = grid.calc_grad_at_link(eta)
    SOC_La_uphill = grid.map_value_at_max_node_to_link(eta,SOC_La)
    q = - D * dzdx
    qc = q * SOC_La_uphill
    dqda = grid.calc_flux_div_at_node(q)
    dqcda = grid.calc_flux_div_at_node(qc)

    return dqda,dqcda

def SOC_transfer_function(eta_old,eta_ini,dzdt,SOC_La,SOC_transfer):
    interface = eta_old - eta_ini - La
    SOC_transfer[dzdt>0.0] = SOC_La[dzdt>0.0] 
    SOC_transfer[dzdt<0.0] = C_SOC * np.exp(interface[dzdt<0.0] / K_SOC)

    return SOC_transfer

#####################

eta_ini = eta.copy()
### LOOP START ###
for t in range(0,nt + 1):
    if t%nt_plot == 0:
        print ('Time = ' + str(t * dt) + ' yrs; ' + str(int(float(t)*dt/T*1000.)/10.) + '% done')
        #Append the new data
        data_elevation.append(grid.at_node['topographic__elevation'].copy())
        data_area.append(grid.at_node['drainage_area'].copy())
        data_soc.append(grid.at_node['SOC_La'].copy())

        #Save the files
        np.save(parent+'\\results\\' + 'elevation', data_elevation)
        np.save(parent+'\\results\\' + 'area', data_area)
        np.save(parent+'\\results\\' + 'soc', data_soc)
        
    fa.run_one_step()
    eta_old = eta.copy()
    sp.run_one_step(dt)
    dqda,dqcda = soil_and_SOC_transport(eta,SOC_La)
    eta[grid.core_nodes] += dt *(U - dqda[grid.core_nodes])
    dzdt = (eta - eta_old)/dt
    SOC_transfer_function(eta_old,eta_ini,dzdt,SOC_La,SOC_transfer)
    SOC_La[grid.core_nodes]  += dt/La * (SOC_transfer[grid.core_nodes] * dqda[grid.core_nodes]  - dqcda[grid.core_nodes] )

#end time
stop_time = time.time()
print (str(round((stop_time -start_time )/60.,1))+' mins')
