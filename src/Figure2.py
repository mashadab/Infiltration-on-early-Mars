#Infiltration in an initially-dry soil with power-law porosity decay
#Mohammad Afzal Shadab and Marc Hesse
#Date modified: 07/17/23

from supporting_tools import *

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
plt.rcParams.update({'font.family': "Serif"})

#Mars parameters
qp = 10# 4e-6 # [m/year] precipitation  1e-6, 2e-6, 4e-6, 5e-6
yr2s = 60**2*24*365.25 # seconds per year
g_mars = 3.711 # [m/s^2] grav. acceleration on Mars
rho = 1e3 # [kg/m^3] desity of water
mu = 1e-3 # [Pa s] water viscosity
m_exp = 2.5422 #power law exponent: variation in porosity with depth
n_exp = 3*m_exp#power law exponent: variation in permeability with depth

qp    = qp / yr2s

################################################
#Dimensional values
################################################
fc_dim= 1e-6 #infiltration capacity [m/s]
z0_dim = 10e3 #bedrock depth [m]
day2s = 24 * 60 * 60
yr2s = 24 * 60 * 60 * 365.25
################################################

## Finding the hydraulic conductivity K0
dmax  = 10e3 # [km] max. aquifer depth
k0_MI = 10**(-12.65) # [m^2] ref. permeability at 1km depth Manning and Ingebritsen (1999)
alpha = 3.2
k_MI  = lambda d: k0_MI*d**(-alpha)
k0    = lambda n: k0_MI/(dmax-1e3)**n
K0_func = lambda n: k0(n)*rho*g_mars/mu  #Hydraulic conductivity

K0 = K0_func(n_exp) #Basal hydraulic conductivity [m/s]
Ks = K0_func(n_exp)*dmax**n_exp  #Surface hydraulic conductivity [m/s]


#parameters
simulation_name = f'Mars_power_law_analytical-inst-ponding'
m = 3 #Cozeny-Karman coefficient for numerator K = K0 (1-phi_i)^m
n = 2 #Corey-Brooks coefficient krw = krw0 * sw^n
s_wr = 0.0 #Residual water saturation
s_gr = 0.0 #Residual gas saturation
phi_0 = 0.5#Porosity at the surface
z0 = 0.06   #point specified for power law porosity variation
p  = 2.5422#index specified for power law porosity variation 
R  = 1#qp/Ks #0.16  #dimensionless rainfall rates

#spatial discretization
zbottom = 0.06
Nz   = 1000
zc   = np.linspace(0,zbottom,Nz) #depth array
phi  = np.transpose([phi_0*(1-zc/z0)**p])#porosity vector: upper layer

#temporal discretization
tmax = 0.05   #time scaling with respect to z0/fc
Nt   = 20000 #time steps
t = np.linspace(0,tmax,Nt+1) #time array

############################################################# 
#Functions for two-layered soils
#############################################################

#Flux in the saturated region qs
def qs(zu,zl):
    return ((m*p - 1)*(zu/z0-zl/z0))/((1-zu/z0)**(1-m*p)-(1-zl/z0)**(1-m*p))
    
#Inverting for moisture content at the surface from rainfall
def func_R2phil_U(theta_0):
    return (f_Cm(np.array([phi_0]),m,phi_0)[0]*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))[0] - R 
theta_0 = opt.fsolve(func_R2phil_U,0.01)[0]  
fcbyR = 1/f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n)[0]

    
#Stage 3: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
def rhs_stage3(t, y): 
    return [ 0, qs(y[0],y[1])/(phi_0*(1-y[1]/z0)**p*(1-s_wr-s_gr))] 

#time stamps of interest: [beginning, stage1, stage2: just after saturation (ts), stage2, stage2: ponding time (tp), stage3]
t_interest = np.array([0,0.01, 0.02, 0.03, 0.04,0.05])
    
zs = 0 #dimensionless depth of complete saturation
ts = 0 #dimensionless time of saturation
tnew = 0
tnew = 0
res = solve_ivp(rhs_stage3, (0, tmax), [0,1e-10],t_eval=t_interest)
tp = 0

print("############################################################# \n")
print("10 (top) Volume fraction with depth")
print("############################################################# \n")

# First set up the figure, the axis
fig,([ax1,ax2,ax3,ax4,ax5,ax6]) = time_sequenced_figure(zbottom,0,0,1,t_interest)

#time_sequenced_figure_dim(zbottom,zsurface,xlim_min,xlim_max,t_interest,fc_dim,z0_dim)

ytop = res.y[0,0]
ybot = res.y[0,0]

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,0]] = s_wr

S_w_analy_int_combined = S_w_analy_int

ax1.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax1.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax1.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')
    
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,1]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,1]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax2.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax2.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax2.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')


S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,2]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,2]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax3.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax3.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax3.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,3]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,3]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax4.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax4.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax4.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,4]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,4]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax5.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax5.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax5.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,5]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,5]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax6.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax6.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax6.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

ax3.set_title(r'''%.3f'''%t_interest[2], fontsize='medium')
ax5.set_title(r'''%.3f'''%t_interest[4], fontsize='medium')
ax1.legend(loc='lower left', shadow=False, fontsize='medium')
plt.xlabel("Volume fractions $\phi$", fontsize='medium')
plt.savefig(f"../Figures/{simulation_name}_CL{theta_0}.pdf")


print("############################################################# \n")
print("10 (top) Volume fraction with depth dimensional")
print("############################################################# \n")

# First set up the figure, the axis
fig,([ax1,ax2,ax3,ax4,ax5,ax6]) = time_sequenced_figure_dim(zbottom,0,0,1,t_interest,fc_dim,z0_dim)

ytop = res.y[0,0]
ybot = res.y[0,0]

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,0]] = s_wr

S_w_analy_int_combined = S_w_analy_int

ax1new = ax1.twinx()
ax1new.set_ylim([0.06,0.0])

ax1new.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{g}$')
ax1new.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{w}$')
ax1new.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{s}$')

ax1new.set_yticks([])
    
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,1]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,1]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax2.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax2.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax2.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')


S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,2]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,2]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax3.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax3.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax3.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,3]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,3]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax4.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax4.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax4.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,4]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,4]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax5.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax5.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax5.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,5]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,5]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax6.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax6.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax6.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

ax6new = ax6.twinx()
ax6new.set_ylim([0.06,0.0])
ax6new.set_ylabel(r"Dimensionless depth, $z/z_0$")

ax1.set_yticklabels(z0_dim*(ax6.get_yticks()))
ax1.set_ylabel(r"Depth, $z$ [m]")

t_dim = (z0_dim /fc_dim * t_interest)/yr2s

ax3.set_title(r"%.2f" "\n" r"%.1f" %(t_interest[2],t_dim[2]), fontsize='medium')
ax5.set_title(r"%.2f" "\n" r"%.1f" %(t_interest[4],t_dim[4]), fontsize='medium')

ax1.legend(loc='lower left', shadow=False, fontsize='medium')
plt.xlabel("Volume fractions, $\phi$", fontsize='medium')
plt.savefig(f"../Figures/{simulation_name}_CL{theta_0}_dimensional.pdf")

print("############################################################# \n")
print("10 (bottom) Saturation with depth")
print("############################################################# \n")

fig,([ax1,ax2,ax3,ax4,ax5,ax6]) = time_sequenced_figure(zbottom,0,-0.05,1.05,t_interest)
ax1.plot(S_w_analy_int_combined[:,0],zc , c = 'k',linestyle='--')
ax2.plot(S_w_analy_int_combined[:,1],zc , c = 'k',linestyle='--')
ax3.plot(S_w_analy_int_combined[:,2],zc , c = 'k',linestyle='--')
ax4.plot(S_w_analy_int_combined[:,3],zc , c = 'k',linestyle='--')
ax5.plot(S_w_analy_int_combined[:,4],zc , c = 'k',linestyle='--')
ax6.plot(S_w_analy_int_combined[:,5],zc , c = 'k',linestyle='--')

ax3.set_title(r'''%.3f'''%(t_interest[2]+0.005), fontsize='medium')
ax5.set_title(r'''%.3f'''%t_interest[4], fontsize='medium')
plt.xlabel("Water saturation $s_w$", fontsize='medium')
plt.savefig(f"../Figures/swvsZpanelshock_{simulation_name}_phi0{theta_0}.pdf")

print("############################################################# \n")
print("11 Dimensionless infiltration rate with dimensionless time")
print("############################################################# \n")

print('R/fc:',R,''', tp':''',tp,''', zs':''',z0,''', ts':''',ts,'\n')

      

fig = plt.figure(figsize=(8,8) , dpi=100)
t0 = np.linspace(0,tp,1000)
qs0 = (f_Cm(np.array([phi_0]),m,phi_0)*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))* np.ones_like(t0)
t1 = np.linspace(tp,tmax*10,10000)
res= solve_ivp(rhs_stage3, (t1[0],t1[-1]), [0, 1e-10],t_eval=t1)
y  = res.y
qs1 = qs(y[0],y[1])

T  = np.concatenate([t0,t1])
QS = np.concatenate([qs0,qs1])
plot = plt.plot(T,QS/f_Cm(np.array([phi_0]),m,phi_0),c='r')
plt.ylabel(r'''$I(t')/f_c$''', fontsize='medium')
plt.xlabel(r'''$t'$''', fontsize='medium')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
plt.xlim([np.min(T), tmax])
plt.ylim([0,round(np.max(QS),1)+0.03])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f"../Figures/ICvsT_phi0array_{simulation_name}_phi0{phi_0}.pdf")




fig = plt.figure(figsize=(8,8) , dpi=100)
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx().twiny()
ax2.plot(T,QS/f_Cm(np.array([phi_0]),m,phi_0),c='r')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
ax2.set_xlim([np.min(T), tmax])
ax2.set_ylim([0,1.1])

def tick_function(X):
    
    V = (X * z0_dim / fc_dim)/yr2s
    return ["%.1f" % z for z in V]

def ytick_function(Y):
    
    V = (Y * 1e6*fc_dim)
    return ["%.1f" % z for z in V]

ax1.set_xticklabels(tick_function(ax2.get_xticks()))
ax1.set_xlabel(r"time, $t$ [years]")

ax1.set_yticklabels(ytick_function(ax2.get_yticks()))
ax1.set_ylabel(r"Infiltration rate, $I(t) \times 10^{6}$ [m/s]")

ax1.set_xlim(np.array([np.min(T), 5]))
ax1.set_ylim(np.array([0,11]))

plt.tight_layout()#pad=0.4, w_pad=0.5, h_pad=1.0)
ax2.set_ylabel(r"$I(t')/f_c$")
ax2.set_xlabel(r"$t'$")

plt.savefig(f"../Figures/dimensional_ICvsT_phi0array_{simulation_name}_phi0{phi_0}.pdf")