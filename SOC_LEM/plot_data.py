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
ls = LightSource(azdeg=315, altdeg=15)
hillshade_boolean = 1

cmap = matplotlib.cm.gray
parent = 'C:\\Users\\jeffs\\Desktop\\SOC_LEM_results'#os.getcwd()

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
        
def plot_SOC(x_plot,y_plot,filename_eta,filename_SOC,slabel,temp_cmap):
    s_eta = np.load(parent+'\\results\\'+filename_eta)
    s_SOC = np.load(parent+'\\results\\'+filename_SOC)
    data_holder_eta = np.zeros((cellsx,cellsy))
    data_holder_SOC = np.zeros((cellsx,cellsy))
        
    for t in range (0,s_eta.shape[0]):
        data_holder_eta = s_eta[int(t),:].reshape(cellsx,cellsy)
        data_holder_SOC = s_SOC[int(t),:].reshape(cellsx,cellsy)
                
        fig = plt.figure(t,figsize = (10.,8.),facecolor='white')
        ax = fig.add_axes([0.0, 0.075, .85, .85])
        cax = fig.add_axes([0.8, 0.075, 0.02, 0.85])
        
        im1 = ax.imshow(np.fliplr(np.rot90(np.rot90(data_holder_SOC))),extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=temp_cmap,interpolation='none',vmin=0.0,vmax=0.035)
        fig.colorbar(im1, cax=cax, label = slabel, format=ticker.FuncFormatter(fmt),fraction=0.046, pad=0.04)
        hill = ls.hillshade(np.fliplr(np.rot90(np.rot90(data_holder_eta))),vert_exag=1,dx=dx,dy=dx)
        im2 = ax.imshow(hill,extent=[y_plot[0],y_plot[-1],x_plot[0],x_plot[-1]],cmap=cmap,alpha=0.25)
        ax.set_xlim(0,y_plot[-1])
        ax.set_ylim(0,x_plot[-1])
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        plt.savefig(parent+'\\results\\plots\\'+filename_SOC[:-4] + '_'+ '%06d' % t + '.png',dpi=300)
        plt.clf()
        plt.close('all')
        
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

print ('Plotting Elevation...')
plot(x_plot,y_plot,'elevation.npy', r'$\eta$ [$m$]',0,cmap)
print ('Plotting Hillsahde...')
plot_hillshade(x_plot,y_plot,'elevation.npy', r'$\eta$ [$m$]',cmap)
print ('Plotting Drainage Area...')
plot(x_plot,y_plot,'area.npy', r'$log(A)$ [$-$]',1,cmap)
print ('Plotting SOC...')
plot_SOC(x_plot,y_plot,'elevation.npy','soc.npy', r'SOC [% by wt.]',new_cmap)
print ('Plotting Difference...')
plot_difference(x_plot,y_plot,'elevation.npy', r'$\Delta\eta$ [$m$]',0,cmap)
print ('renaming...')
os.rename(parent+'\\results',parent+'\\'+runname)

