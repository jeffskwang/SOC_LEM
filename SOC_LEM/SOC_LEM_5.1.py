#importing libaries
import numpy as np
import os
import importlib
import shutil
import time
import sys
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

#####HOUSEKEEPING#####
#import parameters
code_folder = os.getcwd()
sys.path.append(code_folder +'\\drivers')
sys.path.append(code_folder +'\\utilities')
parameters = importlib.import_module(sys.argv[1])
globals().update(parameters.__dict__)

#make output folder
if os.path.isdir(results_folder+'\\results'):
    shutil.rmtree(results_folder+'\\results')
    time.sleep(3)
os.makedirs(results_folder+'\\results')
os.makedirs(results_folder+'\\results\\input')
shutil.copyfile(code_folder +'\\drivers\\'+sys.argv[1]+'.py',results_folder+'\\results\\input\\'+sys.argv[1]+'.py')

#make padded inputs
from pad_rasters import*
pad_rasters(ini_file,pad_mode)
for real_data_file_temp in real_data_file:
    pad_rasters(real_data_file_temp,pad_mode)
for ndvi_data_file_temp in ndvi_data_file:
    pad_rasters(ndvi_data_file_temp,pad_mode)

#load landlab components
from landlab.components import LinearDiffuser
from landlab import RasterModelGrid
from landlab.utils import structured_grid
from landlab.io import read_esri_ascii
from landlab.components import FlowAccumulator, DepressionFinderAndRouter

#load initial conditionf
(grid, eta) = read_esri_ascii(results_folder+'\\results\\input\\padded_'+os.path.basename(ini_file), name='topographic__elevation')
grid.set_nodata_nodes_to_closed(eta, -9999.)
grid.set_fixed_value_boundaries_at_grid_edges(True,True,True,True)
if flow_method == 'D4':
    fa = FlowAccumulator(grid,'topographic__elevation',flow_director='FlowDirectorSteepest',depression_finder=DepressionFinderAndRouter)
elif flow_method == 'D8':
    fa = FlowAccumulator(grid,'topographic__elevation',flow_director='FlowDirectorD8')
elif flow_method == 'Dinf':
    fa = FlowAccumulator(grid,'topographic__elevation',flow_director='FlowDirectorDINF')
nrows = grid.number_of_node_rows
ncols = grid.number_of_node_columns

#define arrays
SOC_La = grid.add_zeros('node','SOC_La')
SOC_transfer = grid.add_zeros('node','SOC_transfer')
dqda = grid.add_zeros('node','dqda')
nz = int(2.0 * Z / dz)
SOC_z = np.zeros((nrows,ncols,nz))
dz_ini = np.zeros((nrows,ncols,nz))
z_coordinates = np.zeros((nrows,ncols,nz))

#grid size and number of cells
dx = grid.dx
dy = grid.dy
nt = int(T / dt)
nt_plot = int(dt_plot/dt)
three_D_nt_plot = int(three_D_dt_plot/dt)

#output files
data_elevation=[]
data_soc=[]
data_area=[]

#timer
start_time = time.time()
#########################

#####FUNCTIONS#####
def curv(eta,curvature,max_curvature):
    for i in range (1,nrows-1):
        for j in range (1,ncols-1):
            curvature[i,j] = (eta[i-1,j] - 2.0 * eta[i,j] + eta[i+1,j]) / dx / dx + (eta[i,j-1] - 2.0 * eta[i,j] + eta[i,j+1]) / dy / dy
            if curvature[i,j] > max_curvature:
                curvature[i,j] = max_curvature
    for i in range (0,nrows):
        curvature[i,0] = 0.0
        curvature[i,ncols-1] = 0.0
    for j in range (0,ncols):
        curvature[0,j] = 0.0
        curvature[nrows-1,j] = 0.0
    return curvature

def IC(eta, dz_ini, SOC_z, SOC_La):
    eta_gauss = np.zeros((nrows,ncols))

    #copy elevation
    eta_temp_long = eta.copy()
    #convert -9999 to numpy NaNs
    eta_temp_long[eta_temp_long==-9999]=np.nan
    #restructure into a 2D array
    eta_temp = eta_temp_long.reshape(nrows,ncols)

    #astropy Gaussian (fills NaNs)
    kernel = Gaussian2DKernel(sigma_value/dx)
    eta_gauss = convolve(eta_temp, kernel,boundary='extend')

    #calculate curvature
    max_curvature = 0.01 #stops negative soci profiles (especially in ditches)    
    curvature = np.zeros((nrows,ncols))
    curvature  = curv(eta_gauss,curvature,max_curvature)
    curvature[np.isnan(eta_temp)]=np.nan
    
    for i in range(0,nrows):
        for j in range(0,ncols):
            dz_ini[i,j,:] = np.linspace(-Z + dz / 2.0, Z - dz / 2.0, nz) 
            B_SOCI = B_SOCI_int + curvature[i,j] * B_SOCI_slope + B_err * np.sqrt(B_MSE* (1./B_n + ((curvature[i,j] -B_average)**2.0)/B_sum_resid))
            C_SOCI = C_SOCI_int + curvature[i,j] * C_SOCI_slope + C_err * np.sqrt(C_MSE* (1./C_n + ((curvature[i,j] -C_average)**2.0)/C_sum_resid))
            SOC_z[i,j,0:int(nz/2)] = A_SOCI + B_SOCI * np.exp(dz_ini[i,j,0:int(nz/2)] * C_SOCI)          
            avg_count = 0
            for k in range(int(nz/2 - La / dz),int(nz/2)):
                avg_count += 1
                SOC_La.reshape(nrows,ncols)[i,j] += SOC_z[i,j,k]
            SOC_La.reshape(nrows,ncols)[i,j] /= float(avg_count)
            
    return dz_ini, SOC_z, SOC_La

def boundary_conditions(eta):
    noflow = grid.map_value_at_min_node_to_link(eta,eta)
    noflow[noflow!=-9999] = 1.0
    noflow[noflow==-9999] = 0.0
    
    return noflow

def soil_and_SOC_transport(eta,SOC_La,dqda,noflow):
    dzdx = grid.calc_grad_at_link(eta)
    SOC_La_uphill = grid.map_value_at_max_node_to_link(eta,SOC_La)    
    q = - D * dzdx * noflow  #q is a link
    qc = q * SOC_La_uphill * noflow
    dqda = grid.calc_flux_div_at_node(q)
    dqcda = grid.calc_flux_div_at_node(qc)

    if flow_method == 'D4' or flow_method == 'D8':
        fa.run_one_step()
        slope = grid.at_node['topographic__steepest_slope']
        area = grid.at_node['drainage_area']
        receiver = grid.at_node['flow__receiver_node']
        
        noflow_sheet = eta[receiver]
        noflow_sheet[noflow_sheet!=-9999] = 1.0
        noflow_sheet[noflow_sheet==-9999] = 0.0
        
        qs = K * slope ** 1.4 * area ** 1.5 * noflow_sheet

        #water
        dqda[receiver] -= qs/dx
        dqda += qs/dx

        dqcda[receiver] -= qs*SOC_La/dx
        dqcda += qs*SOC_La/dx
        
    elif flow_method == 'Dinf':
        fa.run_one_step()
        slope = grid.at_node['topographic__steepest_slope']
        area = grid.at_node['drainage_area']
        area_proportion = grid.at_node['flow__receiver_proportions']
        receiver = grid.at_node['flow__receiver_node']
        
        noflow_sheet = eta[receiver]
        noflow_sheet[noflow_sheet!=-9999] = 1.0
        noflow_sheet[noflow_sheet==-9999] = 0.0

        slope[slope<0.0] = 0.0
        
        qs1 = K * slope[:,0] ** 1.4 * (area * area_proportion[:,0]) ** 1.5 * np.min(noflow_sheet,axis=1)
        qs2 = K * slope[:,1] ** 1.4 * (area * area_proportion[:,1]) ** 1.5 * np.min(noflow_sheet,axis=1)

        #water
        dqda[receiver[:,0]] -= qs1/dx
        dqda[receiver[:,1]] -= qs2/dx
        dqda += (qs1 + qs2)/dx

        dqcda[receiver[:,0]] -= qs1*SOC_La/dx
        dqcda[receiver[:,1]] -= qs2*SOC_La/dx
        dqcda += (qs1 + qs2)*SOC_La/dx
        
    return dqda,dqcda

def find_SOC_cell(interface_z_old,interface_z_new,z_coor,SOC):
    index_locat_old = (np.abs(z_coor - interface_z_old)).argmin()
    index_locat_new = (np.abs(z_coor - interface_z_new)).argmin()
    halfway_point = 0.5 * z_coor[index_locat_old] + 0.5 * z_coor[index_locat_new]    
    SOC_interface = (SOC[index_locat_old] * (interface_z_old - halfway_point)  + SOC[index_locat_new] * (halfway_point - interface_z_new)) / (interface_z_old - interface_z_new)
    
    return SOC_interface  

def find_top_cell_active_layer(interface_z,z_coor,SOC):
    top_cell_active_layer = (np.abs(z_coor - interface_z + dz /2.)).argmin()
    
    return top_cell_active_layer        

def find_bottom_cell_active_layer(interface_z,z_coor,SOC):
    bottom_cell_active_layer = (np.abs(z_coor - interface_z)).argmin()
    
    return bottom_cell_active_layer      

def SOC_transfer_function(eta_old,eta_ini,dzdt,SOC_La,SOC_transfer,SOC_z):
    interface = eta_old - eta_ini - La
    interface_new = eta_old + dzdt * dt - eta_ini - La
    for i in range(0,nrows):
        for j in range(0,ncols):
            if dzdt.reshape(nrows,ncols)[i,j] < 0.0:
                SOC_transfer.reshape(nrows,ncols)[i,j] = find_SOC_cell(interface.reshape(nrows,ncols)[i,j],interface_new.reshape(nrows,ncols)[i,j],dz_ini[i,j,:],SOC_z[i,j,:])
            elif dzdt.reshape(nrows,ncols)[i,j] > 0.0:
                SOC_transfer.reshape(nrows,ncols)[i,j] = SOC_La.reshape(nrows,ncols)[i,j]
                
    return SOC_transfer

def SOC_profile_update(eta,eta_ini,dzdt,SOC_La,SOC_z):
    interface = eta - eta_ini - La
    for i in range(0,nrows):
        for j in range(0,ncols):
            top_cell_active_layer =  find_top_cell_active_layer(interface.reshape(nrows,ncols)[i,j]+La,dz_ini[i,j,:],SOC_z[i,j,:])
            bottom_cell_active_layer =  find_bottom_cell_active_layer(interface.reshape(nrows,ncols)[i,j],dz_ini[i,j,:],SOC_z[i,j,:])
            SOC_z[i,j,top_cell_active_layer+1:] = 0.0
            SOC_z[i,j,bottom_cell_active_layer+1:top_cell_active_layer+1]=SOC_La.reshape(nrows,ncols)[i,j]

            if dzdt.reshape(nrows,ncols)[i,j] > 0.0:
                dz_interface_old = (interface.reshape(nrows,ncols)[i,j] - dt * dzdt.reshape(nrows,ncols)[i,j]) - dz_ini[i,j,bottom_cell_active_layer]  + dz / 2.0
                dz_interface_new = interface.reshape(nrows,ncols)[i,j] - dz_ini[i,j,bottom_cell_active_layer]  + dz / 2.0
                SOC_z[i,j,bottom_cell_active_layer] = (SOC_z[i,j,bottom_cell_active_layer] * dz_interface_old + dt * dzdt.reshape(nrows,ncols)[i,j] *  SOC_La.reshape(nrows,ncols)[i,j]) / (dz_interface_new)
                
    return SOC_z

##### LOOP START #####
noflow = boundary_conditions(eta)
dz_ini, SOC_z, SOC_La = IC(eta, dz_ini, SOC_z, SOC_La)
eta_ini = eta.copy()
for t in range(0,nt + 1):
##    print(t)
##    print (t,np.argmax(SOC_La.reshape(nrows,ncols)),np.max(SOC_La.reshape(nrows,ncols)),(nrows,ncols))
    if t%nt_plot == 0:
        print ('Time = ' + str(t * dt) + '; ' + str(int(float(t)*dt/T*1000.)/10.) + '% done')
        #Append the new data
        data_elevation.append(grid.at_node['topographic__elevation'].copy())
        data_soc.append(grid.at_node['SOC_La'].copy())
        if flow_method != 'noflow':
            data_area.append(grid.at_node['drainage_area'].copy())

        #Save the files
        np.save(results_folder+'\\results\\' + 'elevation', data_elevation)
        np.save(results_folder+'\\results\\' + 'soci', data_soc)
        if flow_method != 'noflow':
            np.save(results_folder+'\\results\\' + 'area', data_area)
            
        if t%three_D_nt_plot == 0:
            np.save(results_folder+'\\results\\' +'3D_SOC_' + '%06d' % + int(t*dt) +'yrs.npy',SOC_z)
            np.save(results_folder+ '\\results\\' +'3D_surface_' + '%06d' % + int(t*dt) +'yrs.npy',eta.reshape(nrows,ncols))
    eta_old = eta.copy()
    dqda,dqcda = soil_and_SOC_transport(eta,SOC_La,dqda,noflow)
    eta[grid.core_nodes] += dt *(- dqda[grid.core_nodes])
    dzdt = (eta - eta_old)/dt
    SOC_transfer = SOC_transfer_function(eta_old,eta_ini,dzdt,SOC_La,SOC_transfer,SOC_z) 
    SOC_La[grid.core_nodes]  += dt/La * (SOC_transfer[grid.core_nodes] * dqda[grid.core_nodes]  - dqcda[grid.core_nodes])
    SOC_z  = SOC_profile_update(eta,eta_ini,dzdt,SOC_La,SOC_z)        

#end time
stop_time = time.time()
print (str(round((stop_time -start_time )/60.,1))+' mins')
