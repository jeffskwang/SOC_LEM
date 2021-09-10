#folders
input_folder = ###LOCATION OF INPUT FOLDER
results_folder = ###LOCATION OF RESULTS FOLDER
#simulation name
runname = 'example_hoffman'
#initial DEM file
ini_file = input_folder +'\\hoffman_DEM.txt'
#observed SOCI file
real_data_file =[input_folder +'\\hoffman_SOCI.txt',]
#ndvi data file
ndvi_data_file = [input_folder +'\\hoffman_NDVI.txt']

###physical parameters###
D = 0.25  #[m^2/yr] diffusion coefficient
La = 0.2 #[m] mixing layer
K = 1.0 * 10. ** (-5.) #[m^0.5/yr] erodibility

###SOCI Profile parameters###
sigma_value = 8.0 #[m] gaussian filter length

#SOCI[depth] = A_SOCI + B_SOCI * exp(-C_SOCI*depth)
#B_SOCI = B_SOCI_int + B_SOCI_slope * curvature
#C_SOCI = C_SOCI_int + C_SOCI_slope * curvature
A_SOCI = 1.5986022611991917
B_SOCI_slope = 190.63
B_SOCI_int = 6.055 - A_SOCI
C_SOCI_slope = -192.29
C_SOCI_int = 1.896

#B_SOCI statistics
B_err = 0.0 #multiple of B_error to add to B_SOCI vs. curvature regression
B_n = 9. #number of B_SOCI observations
B_average = -0.000678733164902306 #average B_SOCI value
B_MSE = 0.6834358539062254 #mean square error of B_SOCI
B_sum_resid = 1.3258417587612205e-05 #sum of residuals of B_SOCI

#C_SOCI statistics
C_err = 0.0#multiple of C_error to add to C_SOCI vs. curvature regression
C_n = 29.#number of C_SOCI observations
C_average = -0.0013086024861925212#average C_SOCI value
C_MSE = 0.593443745455371 #mean square error of C_SOCI
C_sum_resid = 0.0006057960850312982#sum of residuals of C_SOCI

###numerical parameters###
T = 200. # [-] dimensionless simulation time
dt = 1.0 # [-] model timestep
hole_function = 0# 0 is off and 1 is on
dz = 0.01 #[m] soil depth grid step
Z = 2.0 #max deposition, erosion
flow_method = 'noflow' #flow method, either noflow, D4, D8, or Dinf
pad_mode = 1 #0 is a fixed value boundary, 1 is a wall boundary

###output parameters###
dt_plot = 5. # [-] plot timestep
three_D_dt_plot = 40. # [-] 3d data plot timestep, warning: these files are larger



