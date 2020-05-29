import sys
import os
import importlib
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LightSource
import matplotlib.colors as colors
from scipy.ndimage import gaussian_filter
ls = LightSource(azdeg=315, altdeg=15)
hillshade_boolean = 0
plot_real_data = 1
real_data_file = 'D:\\SOC_LEM\\SOC_LEM\\SOC_LEM\\input\\willis_2_soci.txt'
sigma_value = 1

cmap = matplotlib.cm.gray
##parent = 'C:\\Users\\jeffs\\Desktop\\SOC_LEM_results'#laptop
##parent = 'C:\\Users\\jkwang\\Desktop\\SOC_LEM_results'#umass desktop
parent = 'D:\\SOC_LEM\\SOC_LEM_results'#desktop

if os.path.isdir(parent+'\\results\\plots')==False:
    os.makedirs(parent+'\\results\\plots')

with open(parent+'\\results\\'+'parameters_and_IC.asc') as file:
    lines = [next(file) for x in range(10)]

def read_parameters(field,lines):
    value = None
    for line in lines:
        if line.startswith(field):
            end_of_string = line[len(field):]
            value = end_of_string.strip()
    return value

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

new_cmap = truncate_colormap(matplotlib.cm.BrBG_r, 0.7, 1.0)

runname = read_parameters('runname',lines)
cellsx = int(read_parameters('cellsx',lines))
cellsy = int(read_parameters('cellsy',lines))
dx = float(read_parameters('dx',lines))
dy = float(read_parameters('dy',lines))
T = float(read_parameters('T',lines))
dt= float(read_parameters('dt',lines))
K = float(read_parameters('K',lines))
D = float(read_parameters('D',lines))
U = float(read_parameters('U',lines))
x_plot = np.linspace(0,dx*(float(cellsx-1.)),cellsx)
y_plot = np.linspace(0,dy*(float(cellsy-1.)),cellsy)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${}\times 10^{{{}}}$'.format(a, b)

#plot_function
def plot(x_plot,y_plot,filename,slabel,log_boolean,temp_cmap):
    s = np.load(parent+'\\results\\'+filename)
    data_holder = np.zeros((cellsx,cellsy))
    if log_boolean == 1:
        s[s==0] = np.min(s[s!=0])
    for t in range (0,s.shape[0]):
        if log_boolean == 1:
            data_holder = np.log10(s[int(t),:].reshape(cellsx,cellsy))
        elif log_boolean == 0:
            data_holder =s[int(t),:].reshape(cellsx,cellsy)
            
        fig = plt.figure(t,figsize = (10.,8.),facecolor='white')
        ax = fig.add_axes([0.0, 0.075, .85, .85])
        cax = fig.add_axes([0.8, 0.075, 0.02, 0.85])
        im = ax.imshow(np.fliplr(np.rot90(np.rot90(data_holder))),extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=temp_cmap,interpolation='none')
        fig.colorbar(im, cax=cax, label = slabel, format=ticker.FuncFormatter(fmt),fraction=0.046, pad=0.04)
        ax.set_xlim(0,y_plot[-1])
        ax.set_ylim(0,x_plot[-1])
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        plt.savefig(parent+'\\results\\plots\\'+filename[:-4] + '_'+ '%06d' % t + '.png',dpi=300)
        plt.clf()
        plt.close('all')
         
def plot_hillshade(x_plot,y_plot,filename_eta,slabel,temp_cmap):
    s_eta = np.load(parent+'\\results\\'+filename_eta)
    data_holder_eta = np.zeros((cellsx,cellsy))
        
    for t in range (0,s_eta.shape[0]):
        data_holder_eta = s_eta[int(t),:].reshape(cellsx,cellsy)
                
        fig = plt.figure(t,figsize = (10.,8.),facecolor='white')
        ax = fig.add_axes([0.0, 0.075, .85, .85])
        cax = fig.add_axes([0.8, 0.075, 0.02, 0.85])
        
        im1 = ax.imshow(np.fliplr(np.rot90(np.rot90(data_holder_eta))),extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=temp_cmap,interpolation='none')
        fig.colorbar(im1, cax=cax, label = slabel, format=ticker.FuncFormatter(fmt),fraction=0.046, pad=0.04)
        hill = ls.hillshade(np.fliplr(np.rot90(np.rot90(data_holder_eta))),vert_exag=1,dx=dx,dy=dx)
        im2 = ax.imshow(hill,extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=cmap,alpha=0.5)
        ax.set_xlim(0,y_plot[-1])
        ax.set_ylim(0,x_plot[-1])
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        plt.savefig(parent+'\\results\\plots\\hillshade_'+ '%06d' % t + '.png',dpi=300)
        plt.clf()
        plt.close('all')
        
def plot_SOC(x_plot,y_plot,filename_eta,filename_SOC,filename_SOC_real,slabel,temp_cmap,max_soci, min_soci):
    s_eta = np.load(parent+'\\results\\'+filename_eta)
    s_SOC = np.load(parent+'\\results\\'+filename_SOC)
    data_holder_eta = np.zeros((cellsx,cellsy))
    data_holder_SOC = np.zeros((cellsx,cellsy))
    data_holder_SOC_real = np.loadtxt(filename_SOC_real,skiprows=6)
    data_holder_SOC_real = gaussian_filter(data_holder_SOC_real, sigma=sigma_value,mode='nearest')
        
    for t in range (0,s_eta.shape[0]):
        data_holder_eta = s_eta[int(t),:].reshape(cellsx,cellsy)
        data_holder_SOC = s_SOC[int(t),:].reshape(cellsx,cellsy)
                
        fig = plt.figure(t,figsize = (10.,8.),facecolor='white')
        ax = fig.add_axes([0.0, 0.075, .85, .85])
        cax = fig.add_axes([0.8, 0.075, 0.02, 0.85])
        
        im1 = ax.imshow(np.fliplr(np.rot90(np.rot90(data_holder_SOC))),\
                        extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],\
                        cmap=temp_cmap,interpolation='none',vmin=min_soci,vmax=max_soci)
        fig.colorbar(im1, cax=cax, label = slabel, format=ticker.FuncFormatter(fmt),fraction=0.046, pad=0.04)
##        hill = ls.hillshade(np.fliplr(np.rot90(np.rot90(data_holder_eta))),vert_exag=1,dx=dx,dy=dx)
##        im2 = ax.imshow(hill,extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=cmap,alpha=0.25)
        ax.set_xlim(0,y_plot[-1])
        ax.set_ylim(0,x_plot[-1])
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        plt.savefig(parent+'\\results\\plots\\'+filename_SOC[:-4] + '_'+ '%06d' % t + '.png',dpi=300)
        plt.clf()
        plt.close('all')
                        
        fig2 = plt.figure(t,figsize = (10.,8.),facecolor='white')
        ax2 = fig2.add_axes([0.0, 0.075, .85, .85])
        cax2 = fig2.add_axes([0.8, 0.075, 0.02, 0.85])
        max_soci_diff = np.max(np.abs(np.fliplr(np.rot90(np.rot90(data_holder_SOC))) - np.flipud(np.fliplr(np.rot90(np.rot90(data_holder_SOC_real))))))
        im2 = ax2.imshow(np.fliplr(np.rot90(np.rot90(data_holder_SOC))) - np.flipud(np.fliplr(np.rot90(np.rot90(data_holder_SOC_real)))),\
                        extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],\
                        cmap=matplotlib.cm.bwr,interpolation='none',vmin=-max_soci_diff,vmax=max_soci_diff)
        fig2.colorbar(im2, cax=cax2, label = slabel, format=ticker.FuncFormatter(fmt),fraction=0.046, pad=0.04)
##        hill = ls.hillshade(np.fliplr(np.rot90(np.rot90(data_holder_eta))),vert_exag=1,dx=dx,dy=dx)
##        im2 = ax.imshow(hill,extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=cmap,alpha=0.25)
        ax2.set_xlim(0,y_plot[-1])
        ax2.set_ylim(0,x_plot[-1])
        ax2.set_xlabel('x [m]')
        ax2.set_ylabel('y [m]')
        plt.savefig(parent+'\\results\\plots\\'+filename_SOC[:-4] + '_difference_'+ '%06d' % t + '.png',dpi=300)
        plt.clf()
        plt.close('all')
        
def plot_SOC_real(x_plot,y_plot,filename_eta,filename_SOC,slabel,temp_cmap,max_soci, min_soci,gauss):
    s_eta = np.load(parent+'\\results\\'+filename_eta)
    data_holder_eta = np.zeros((cellsx,cellsy))
    data_holder_SOC = np.loadtxt(filename_SOC,skiprows=6)
    
    t = s_eta.shape[0] - 1
    data_holder_eta = s_eta[int(t),:].reshape(cellsx,cellsy)
    if gauss == 1:
        data_holder_SOC = gaussian_filter(data_holder_SOC, sigma=sigma_value,mode='nearest')
    fig = plt.figure(t,figsize = (10.,8.),facecolor='white')
    ax = fig.add_axes([0.0, 0.075, .85, .85])
    cax = fig.add_axes([0.8, 0.075, 0.02, 0.85])
    
    im1 = ax.imshow(np.flipud(np.fliplr(np.rot90(np.rot90(data_holder_SOC)))),\
                    extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],\
                        cmap=temp_cmap,interpolation='none')#,vmin=min_soci,vmax=max_soci)
    fig.colorbar(im1, cax=cax, label = slabel, format=ticker.FuncFormatter(fmt),fraction=0.046, pad=0.04)
##    hill = ls.hillshade(np.fliplr(np.rot90(np.rot90(data_holder_eta))),vert_exag=1,dx=dx,dy=dx)
##    im2 = ax.imshow(hill,extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=cmap,alpha=0.25)
    ax.set_xlim(0,y_plot[-1])
    ax.set_ylim(0,x_plot[-1])
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    if gauss == 0:
        plt.savefig(parent+'\\results\\plots\\real_data.png',dpi=300)
    else:
        plt.savefig(parent+'\\results\\plots\\real_data_gauss.png',dpi=300)
    plt.clf()
    plt.close('all')
        
def plot_SOC_cali(filename_SOC,filename_SOC_real):
    s_SOC = np.load(parent+'\\results\\'+filename_SOC)
    data_holder_SOC_real = np.flipud(np.fliplr(np.rot90(np.rot90(np.loadtxt(filename_SOC_real,skiprows=6)))))
    data_holder_SOC_real = gaussian_filter(data_holder_SOC_real, sigma=sigma_value,mode='nearest')
    
    for t in range (0,s_SOC.shape[0]):
        data_holder_SOC = np.zeros((cellsx,cellsy))
        data_holder_SOC = np.fliplr(np.rot90(np.rot90(s_SOC[int(t),:].reshape(cellsx,cellsy))))
        text_holder = np.zeros((cellsx*cellsy,2))
        text_int = 0
        for x in range(0,cellsx):
            for y in range(0,cellsy):
                text_holder[text_int,0] = data_holder_SOC[x,y]
                text_holder[text_int,1] = data_holder_SOC_real[x,y]
                text_int += 1
        np.savetxt(parent+'\\results\\plots\\soci_data_'+str(t)+'.txt',text_holder)
            
def plot_difference(x_plot,y_plot,filename,slabel,log_boolean,temp_cmap):
    s = np.load(parent+'\\results\\'+filename)
    data_holder = np.zeros((cellsx,cellsy))
    data_holder_ini = np.zeros((cellsx,cellsy))
    if log_boolean == 1:
        s[s==0] = np.min(s[s!=0])
    
    for t in range (0,s.shape[0]):
        data_holder = s[int(t),:].reshape(cellsx,cellsy) - s[0,:].reshape(cellsx,cellsy)
        fig = plt.figure(t,figsize = (10.,8.),facecolor='white')
        ax = fig.add_axes([0.0, 0.075, .85, .85])
        cax = fig.add_axes([0.8, 0.075, 0.02, 0.85])
        im = ax.imshow(np.fliplr(np.rot90(np.rot90(data_holder))),extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=temp_cmap,interpolation='none',vmin=-1.0,vmax=1.0)
        fig.colorbar(im, cax=cax, label = slabel, format=ticker.FuncFormatter(fmt),fraction=0.046, pad=0.04)
        ax.set_xlim(0,y_plot[-1])
        ax.set_ylim(0,x_plot[-1])
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        plt.savefig(parent+'\\results\\plots\\'+filename[:-4] + '_difference_'+ '%06d' % t + '.png',dpi=300)
        plt.clf()
        plt.close('all')

max_soci, min_soci = 8.0 , 2.0
if plot_real_data == 1:
##    plot_SOC_real(x_plot,y_plot,'elevation.npy',real_data_file, r'SOC [% by wt.]', matplotlib.cm.viridis,max_soci, min_soci)#,new_cmap)
    plot_SOC_real(x_plot,y_plot,'elevation.npy',real_data_file, r'SOCI [-]', matplotlib.cm.viridis,max_soci, min_soci,0)#,new_cmap)
    plot_SOC_real(x_plot,y_plot,'elevation.npy',real_data_file, r'SOCI [-]', matplotlib.cm.viridis,max_soci, min_soci,1)#,new_cmap)
    
    plot_SOC_cali('soc.npy',real_data_file)

print ('Plotting Elevation...')
plot(x_plot,y_plot,'elevation.npy', r'$\eta$ [$m$]',0,cmap)
print ('Plotting Hillsahde...')
plot_hillshade(x_plot,y_plot,'elevation.npy', r'$\eta$ [$m$]',cmap)
print ('Plotting Drainage Area...')
plot(x_plot,y_plot,'area.npy', r'$log(A)$ [$-$]',1,cmap)
print ('Plotting SOC...')
##plot_SOC(x_plot,y_plot,'elevation.npy','soc.npy', r'SOC [% by wt.]', matplotlib.cm.viridis,max_soci, min_soci)#,new_cmap)
plot_SOC(x_plot,y_plot,'elevation.npy','soc.npy',real_data_file, r'SOCI [-]', matplotlib.cm.viridis,max_soci, min_soci)#,new_cmap)
print ('Plotting Difference...')
plot_difference(x_plot,y_plot,'elevation.npy', r'$\Delta\eta$ [$m$]',0,cmap)
print ('renaming...')
os.rename(parent+'\\results',parent+'\\'+runname)

