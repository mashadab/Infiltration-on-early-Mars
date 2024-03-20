#Supporting tools
#Mohammad Afzal Shadab and Marc Hesse
#Date modified: 08/05/21

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import scipy.optimize as opt
from matplotlib import cm         #color map
plt.rcParams.update({'font.family': "Serif"})
plt.rcParams.update({'font.size': 22})
#For plotting in a different window:  %matplotlib qt 

#############################################################
#Color palette
#############################################################

brown  = [181/255 , 101/255, 29/255]
red    = [255/255 ,255/255 ,255/255 ]
blue   = [ 30/255 ,144/255 , 255/255 ]
green  = [  0/255 , 166/255 ,  81/255]
orange = [247/255 , 148/255 ,  30/255]
purple = [102/255 ,  45/255 , 145/255]
brown  = [155/255 ,  118/255 ,  83/255]
tan    = [199/255 , 178/255 , 153/255]
gray   = [100/255 , 100/255 , 100/255]

#############################################################
#Classes used
#############################################################

class Grid:
    def __init__(self):
        self.xmin = []
        self.xmax = []
        self.Nx = []
        self.N = []
        
#############################################################
#Flux used in present work: f_Cm*f_Cn
#############################################################

#Dimensionless saturated flux, f/f_fc
def f_Cm(phi,m,phi_L):
    fC = np.zeros_like(phi)        
    fC = phi**m / phi_L**m           #Power law porosity
    return fC

#Relative permeability of water, k_rw
def f_Cn(theta,phi,s_gr,s_wr,n):
    fC = np.zeros_like(phi)
    fC = ((theta/phi-s_wr)/(1-s_gr-s_wr))**n    #Power law rel perm    
    return fC


#############################################################
#Modified van-Genuchten model: Used in Hydrus
#############################################################

#For more information about the parameters visit
#http://www.pc-progress.com/en/OnlineHelp/HYDRUS3/Hydrus.html?WaterFlowParameters.html

Se = lambda  theta  ,theta_s,theta_r: (theta - theta_r)/(theta_s - theta_r)  #effective water content [-]
Se_k=lambda  theta_k,theta_s,theta_r: (theta_k - theta_r)/(theta_s - theta_r)#modified van-Genuchten model parameter [c-]
def Kr(Kk,Ks,theta,theta_r,theta_s,theta_k,theta_m,theta_a,m):
            return Kk/Ks*(Se(theta,theta_s,theta_r)/Se_k(theta_k,theta_s,theta_r))*\
            ((F(theta_r,theta_m,theta_a,m)-F(theta,theta_m,theta_a,m))/(F(theta_r,theta_m,theta_a,m)-F(theta_k,theta_m,theta_a,m)))**2 #relative hydraulic conductivity [-]
F  = lambda  theta,theta_m,theta_a,m:       (1-((theta-theta_a)/(theta_m-theta_a))**(1/m))**m #modified van-Genuchten model parameter [-]


#############################################################
#Plotting the time sequenced images
#############################################################

def time_sequenced_figure(zbottom,zsurface,xlim_min,xlim_max,t_interest):

    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    ax1 = fig.add_subplot(1, 6, 1)
    ax2 = fig.add_subplot(1, 6, 2)
    ax3 = fig.add_subplot(1, 6, 3)
    ax4 = fig.add_subplot(1, 6, 4)
    ax5 = fig.add_subplot(1, 6, 5)
    ax6 = fig.add_subplot(1, 6, 6)
    
    ax1.set_ylabel(r'Dimensionless depth $z/z_0$', fontsize='medium')
    ax1.set_ylim([zbottom,zsurface])
    ax1.set_xlim([xlim_min,xlim_max])
    
    ax2.set_xlim([xlim_min,xlim_max])
    ax2.set_ylim([zbottom,zsurface])
    ax2.axes.yaxis.set_visible(False)
    
    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", bottom=False, left=False)
    
    ax3.set_xlim([xlim_min,xlim_max])
    ax3.axes.yaxis.set_visible(False)
    ax3.set_ylim([zbottom,zsurface])
    
    ax4.set_xlim([xlim_min,xlim_max])
    ax4.axes.yaxis.set_visible(False)
    ax4.set_ylim([zbottom,zsurface])
    
    ax5.set_xlim([xlim_min,xlim_max])
    ax5.axes.yaxis.set_visible(False)
    ax5.set_ylim([zbottom,zsurface])
    
    ax6.set_xlim([xlim_min,xlim_max])
    ax6.axes.yaxis.set_visible(False)
    ax6.set_ylim([zbottom,zsurface])
    
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    
    plt.yticks(fontsize='medium')
    plt.xticks(fontsize='medium')
    ax1.set_title(r'''$t'=$%.3f'''%t_interest[0], fontsize='medium')
    ax2.set_title(r'%.3f'%t_interest[1], fontsize='medium')
    ax4.set_title(r'%.3f'%t_interest[3], fontsize='medium')
    ax6.set_title(r'%.3f'%t_interest[5], fontsize='medium')
    plt.subplots_adjust(wspace=0.25, hspace=0)
    
    return fig,([ax1,ax2,ax3,ax4,ax5,ax6]) 





