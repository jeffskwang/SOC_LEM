#python game
#written by jeffrey kwang
#email: jeffskwang@gmail.com
#version = alpha_1

############################
###LIBRARIES###
############################

import pygame
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib
import matplotlib.pyplot as plt
import os
#plt.switch_backend('Agg') fgtf <--- pillsbury's comtribution
folder = 'willis_strat_hole'

############################
###FUNCTIONS###
############################

def plot_setup_color(plot_array,x,y,eta,eta_ini,xlabel,ylabel):
    width_pixel,height_pixel = plot_array.shape[0], plot_array.shape[1]
    fig = Figure(figsize=(float(width_pixel)/100.,float(height_pixel)/100.),dpi=100)
    ax = fig.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    X = np.zeros((nz,pixels+1))
    Y = np.zeros((nz,pixels+1))
    
    for i in xrange(0,x.shape[0]+1):
        if i == 0:
            X[:,i] = x[i]
            Y[:,i] = np.linspace(1.,-1.,200) + eta_ini[i]
        else:
            X[:,i] = x[i-1]
            Y[:,i] = np.linspace(1.,-1.,200) + eta_ini[i-1]
    im = ax.pcolor(X,Y,np.rot90(y),cmap=new_cmap)
    ax.plot(x,eta,color='k')
    ax.plot(x,eta_ini,color='r')
        
    fig.set_tight_layout(True)
    canvas = FigureCanvas(fig)

    canvas.draw()
    
    buf = fig.canvas.tostring_rgb()
    ncols,nrows = fig.canvas.get_width_height()
    buf = np.fromstring(buf, dtype=np.uint8).reshape(nrows, ncols, 3)

    return np.transpose(buf,(1, 0, 2))

def plot_setup(plot_array,x,y,xlabel,ylabel):
    width_pixel,height_pixel = plot_array.shape[0], plot_array.shape[1]
    fig = Figure(figsize=(float(width_pixel)/100.,float(height_pixel)/100.),dpi=100)
    ax = fig.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    ax.plot(x,y,color='r')

    ax.text(x[0],y[0],'A')
    ax.text(x[-1],y[-1],'A\'')
        
    fig.set_tight_layout(True)
    canvas = FigureCanvas(fig)

    canvas.draw()
    
    buf = fig.canvas.tostring_rgb()
    ncols,nrows = fig.canvas.get_width_height()
    buf = np.fromstring(buf, dtype=np.uint8).reshape(nrows, ncols, 3)

    return np.transpose(buf,(1, 0, 2))

############################
###PARAMETERS###
############################
#scale of screen_res / DEM_res
scale = 4

#frame rate
f_rate = 10

#color map
import matplotlib.colors as colors
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

new_cmap = truncate_colormap(matplotlib.cm.BrBG_r, 0.7, 1.0)


############################
###RASTERS###
############################
parent_folder = os.getcwd() +'\\' + folder+'\\'
year_list = []
year_int = 0
for filename in os.listdir(parent_folder):
    if filename.startswith('3D_surface_'):
        year_list.append(int(filename[11:-7]))
    
DEM_INI = np.rot90(np.load(parent_folder+'3D_surface_' + '%06d' % + int(0) +'yrs.npy'))
DEM = np.rot90(np.load(parent_folder+'3D_surface_' + '%06d' % + int(year_list[year_int]) +'yrs.npy'))
SOC = np.rot90(np.load(parent_folder+'3D_SOC_' + '%06d' % + int(year_list[year_int]) +'yrs.npy'))
min_ele = np.min(DEM)
max_ele = np.max(DEM)
SOC[SOC==0]=np.nan
res_width = DEM.shape[0]
res_height = DEM.shape[1]

############################
###PYGAME###
############################

#initialize pygame
zzz = pygame.init()
key_down = 0
change_data = 0

#make game display
gameDisplay = pygame.display.set_mode((int(res_width * scale + res_height * scale),int(res_height * scale)))

#font
font = pygame.font.SysFont('arial',20)
text1 = font.render('A', True, (255,0,0)) 
text2 = font.render('A\'', True, (255,0,0))
text_year = font.render(str(year_list[year_int]) + ' years of erosion', True,(255,0,0))
  
#set title of game
pygame.display.set_caption('dem')

#update screen
pygame.display.update()

#set gameExit to false to enter the main loop
gameExit = False

#set clock for framerate lock
clock = pygame.time.Clock()

############################
###ARRAYS###
############################
game_display_array = np.zeros((int(res_width),int(res_height),3),dtype=int)

plot_1 = np.zeros((int(res_height*scale),int(res_height*scale),3),dtype=int)
plot_2 = np.zeros((int(res_height*scale),int(res_height*scale),3),dtype=int)

x_start,y_start,x_end,y_end = 0.0,0.0,0.0,0.0
pixels = 501
dx = 2.95302823
nz = SOC.shape[2]
SOC_blank = np.empty(nz)
SOC_blank[:] = np.NaN

##############################
#####MAIN LOOP###
##############################
while not gameExit:
    #event handler (temp event) use pygame events, e.g. contains mouse data, keyboard presses
    for event in pygame.event.get():
        #quit and leave the loop
        if event.type == pygame.QUIT:
            gameExit = True

    if event.type == pygame.KEYDOWN:
        if event.key == pygame.K_UP and key_down == 0:
            key_down = 1
            year_int += 1
            if year_int > (len(year_list) - 1):
                year_int = len(year_list) - 1
            else:
                change_data = 1
        if event.key == pygame.K_DOWN and key_down == 0:
            key_down = 1
            year_int -= 1
            if year_int < 0:
                year_int = 0
            else:
                change_data = 1

    if change_data == 1:
        DEM = np.rot90(np.load(parent_folder+'3D_surface_' + '%06d' % + int(year_list[year_int]) +'yrs.npy'))
        SOC = np.rot90(np.load(parent_folder+'3D_SOC_' + '%06d' % + int(year_list[year_int]) +'yrs.npy'))
        min_ele = np.min(DEM)
        max_ele = np.max(DEM)
        SOC[SOC==0]=np.nan
        text_year = font.render(str(year_list[year_int]) + ' years of erosion', True,(255,0,0))
        change_data = 0
        
    if pygame.mouse.get_pressed()[0]:
        (x_mouse,y_mouse) = event.pos
        x_start = int(x_mouse/scale)
        y_start = int(y_mouse/scale)

    if pygame.mouse.get_pressed()[2]:
        (x_mouse,y_mouse) = event.pos
        x_end = int(x_mouse/scale)
        y_end = int(y_mouse/scale)

    if pygame.mouse.get_pressed()[1]:
        x_start_DEM = int(res_width)-x_start-1 
        y_start_DEM = int(res_height)-y_start-1
        x_end_DEM = int(res_width)-x_end-1
        y_end_DEM = int(res_height)-y_end-1
        x_arr = np.linspace(float(x_start_DEM),float(x_end_DEM),pixels)
        y_arr = np.linspace(float(y_start_DEM),float(y_end_DEM),pixels)
        l_arr = np.linspace(0,np.sqrt((float(x_start_DEM )-float(x_end_DEM ))**2.0+(float(y_start_DEM )-float(y_end_DEM ))**2.0),pixels) * dx
        eta_arr = np.zeros(pixels)
        eta_ini_arr = np.zeros(pixels)
        SOC_arr = np.zeros((pixels,nz))
        x_curr,y_curr=int(x_arr[0]),int(y_arr[0])
        for ijk in xrange(0,pixels):
            x_lower = int(x_arr[ijk])
            x_upper = x_lower + 1
            y_lower = int(y_arr[ijk])
            y_upper = y_lower + 1
            
            if x_curr != x_lower or y_curr != y_lower:
                SOC_arr[ijk,:] = SOC_blank
                x_curr,y_curr=x_lower,y_lower
            else:
                SOC_arr[ijk,:] = SOC[int(x_lower + 0.5),int(y_lower + 0.5),:]
            eta_arr[ijk] = DEM[int(x_lower + 0.5),int(y_lower + 0.5)]
            eta_ini_arr[ijk] = DEM_INI[int(x_lower + 0.5),int(y_lower + 0.5)]
        
        plot_1 = plot_setup(plot_1,l_arr,eta_arr,'L [$m$]',r'$\eta$ [$m$]')
        plot_2 = plot_setup_color(plot_2,l_arr,SOC_arr,eta_arr,eta_ini_arr,'L [$m$]',r'z [$m$]')        
        
    if event.type == pygame.KEYUP:
        if event.key == pygame.K_UP:
            key_down = 0 
        if event.key == pygame.K_DOWN:
            key_down = 0 

    game_display_array[:,:,0] = ((DEM) - min_ele) / (max_ele - min_ele) * 255
    game_display_array[:,:,1] = ((DEM) - min_ele) / (max_ele - min_ele) * 255
    game_display_array[:,:,2] = ((DEM) - min_ele) / (max_ele - min_ele) * 255

    #DEM surface
    array_to_surface = pygame.surfarray.make_surface(np.rot90(np.rot90(game_display_array)))
    scaled_surface = pygame.transform.scale(array_to_surface,(int(res_width * scale),int(res_height * scale)))
    plot_to_surface = pygame.surfarray.make_surface(plot_2)
    gameDisplay.blit(scaled_surface,(0,0))
    gameDisplay.blit(plot_to_surface,(res_width * scale,0))
    
    pygame.draw.line(gameDisplay,(255,0,0),[x_start*scale,y_start*scale],[x_end*scale,y_end*scale]) 
    gameDisplay.blit(text1,(x_start*scale,y_start*scale))
    gameDisplay.blit(text2,(x_end*scale,y_end*scale))
    gameDisplay.blit(text_year,(0.01 * res_width*scale,0.01 * res_height*scale))
    
    #update screen
    pygame.display.update()
    clock.tick(f_rate)

#unintialize and quit pygame
pygame.quit()
quit()



###Extras

#bilinear interpolation
            
##            SOC_arr[ijk,:] = SOC[int(x_lower + 0.5),int(y_lower + 0.5),:]
##            
##            y_upper = y_lower + 1
##            mean_1 = (DEM[x_lower,y_upper] - DEM[x_lower,y_lower]) * (y_arr[ijk] - float(y_lower)) + DEM[x_lower,y_lower]
##            mean_2 = (DEM[x_upper,y_upper] - DEM[x_upper,y_lower]) * (y_arr[ijk] - float(y_lower)) + DEM[x_upper,y_lower]
##            mean_3 = (mean_2 - mean_1) * (x_arr[ijk] - float(x_lower)) + mean_1
##            eta_arr[ijk] = mean_3
##            
##            mean_ini_1 = (DEM_INI[x_lower,y_upper] - DEM_INI[x_lower,y_lower]) * (y_arr[ijk] - float(y_lower)) + DEM_INI[x_lower,y_lower]
##            mean_ini_2 = (DEM_INI[x_upper,y_upper] - DEM_INI[x_upper,y_lower]) * (y_arr[ijk] - float(y_lower)) + DEM_INI[x_upper,y_lower]
##            mean_ini_3 = (mean_ini_2 - mean_ini_1) * (x_arr[ijk] - float(x_lower)) + mean_ini_1
##            eta_ini_arr[ijk] = mean_ini_3
