#Infiltration in a two-layered initially-dry soil (Combined)
#Mohammad Afzal Shadab and Marc Hesse
#Date modified: 08/05/21

from supporting_tools import *

#parameters
simulation_name = f'two_layered'
m = 3      #Cozeny-Karman coefficient for numerator K = K0 (1-phi_i)^m
n = 2      #Corey-Brooks coefficient krw = krw0 * se(sw)^n
s_wr = 0.0 #Residual water saturation
s_gr = 0.0 #Residual gas saturation
phi_U = 0.4; phi_L = 0.2 #Porosity of upper and lower layer
theta_L = s_gr*phi_L #setting the saturation in domain
zsurface = 0    #dimensionless surface height
z0 = 1.0         #dimensionless location of jump
R = 0.999 #Array of dimensionless rainfall rates


################################################
#Dimensional values
################################################
fc_dim= 1e-5 #infiltration capacity [m/s]
z0_dim = 10 #dune depth [m]
day2s = 24 * 60 * 60
yr2s = 24 * 60 * 60 * 365.25
################################################


#spatial discretization
zbottom = 2
Nz = 1000
zc   = np.linspace(zsurface,zbottom,Nz) #depth array
phi  = phi_U*np.ones((len(zc),1)) #porosity vector: upper layer
phi[zc>z0] = phi_L                #porosity vector: upper layer

#temporal discretization
tmax = 20   #time scaling with respect to z0/fc
Nt   = 20000#time steps
t = np.linspace(0,tmax,Nt+1) #time array
#time stamps of interest: [beginning, stage1, stage2: just after saturation (ts), stage2, stage2: ponding time (tp), stage3]
t_interest = np.array([0,0.2,0.6249999999999999 +0.005,0.6,0.8,1.0])

############################################################# 
#Functions for two-layered soils
#############################################################

#Flux in the saturated region qs
def qs(x):
    return (x-1)/(x/f_Cm(np.array([phi_U]),m,phi_U) - 1/f_Cm(np.array([phi_L]),m,phi_U))

#To convert dimensionless rainfall to moisture content at the surface
def funtheta_L2phil_U(theta_U):
    return (f_Cm(np.array([phi_U]),m,phi_U)[0]*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0] - R 
theta_U = opt.fsolve(funtheta_L2phil_U,0.01)[0]   

#Stage 2: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
def rhs_stage2(t, y): 
    return [((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0]/(phi_U*(1-s_gr)-theta_U), \
            ((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-s_gr)-theta_L)]

#Stage 3: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
def rhs_stage3(t, y): 
    return [0, \
            ((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-s_gr)-theta_L)]


print("############################################################# \n")
print("2 (top) Volume fraction with depth")
print("############################################################# \n")

# First set up the figure, the axis
fig,([ax1,ax2,ax3,ax4,ax5,ax6]) = time_sequenced_figure(zbottom,zsurface,0,1,t_interest)

#Relation (32) - (35) to find the shock speed ratio k
aa = phi_U/phi_L*(1-theta_U/phi_U-s_gr)/(1-theta_L/phi_L-s_gr)*(1- f_Cm(np.array([phi_L]),m,phi_U) * f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n) / f_Cm(np.array([phi_U]),m,phi_U))
bb =-phi_U/phi_L*(1-theta_U/phi_U-s_gr)/(1-theta_L/phi_L-s_gr)*(1- f_Cm(np.array([phi_L]),m,phi_U) * f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n) / f_Cm(np.array([phi_L]),m,phi_U)) - (1 - f_Cm(np.array([phi_U]),m,phi_U) * f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n) / f_Cm(np.array([phi_U]),m,phi_U))
cc = 1 - f_Cm(np.array([phi_U]),m,phi_U) * f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n) / f_Cm(np.array([phi_L]),m,phi_U)
k = (-bb-np.sqrt(bb**2-4*aa*cc))/(2*aa)

#Calculating the shock speeds
s_U_analy =((qs(k) - f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0]/(phi_U*(1-s_gr-theta_U/phi_U)))
s_L_analy =((qs(k) - f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-theta_L/phi_L-s_gr)))

Sf = ((f_Cm(np.array([phi_U]),m,phi_U)[0]*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0] / (phi_U*(theta_U/phi_U - s_wr))) #initial front speed
ts = (z0 - zsurface)/Sf             #dimensionless time of saturation
tp = ts + (zsurface-z0)/s_U_analy   #dimensionless time of ponding

xf = zsurface + Sf*t_interest[0]
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=xf] = theta_U/phi_U
S_w_analy_int[zc>=xf] = theta_L/phi_L

S_w_analy_int_combined = []
S_w_analy_int_combined = S_w_analy_int.copy()

ax1.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{g}$')
ax1.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{w}$')
ax1.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{s}$')

xf = zsurface + Sf*t_interest[1]
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=xf] = theta_U/phi_U
S_w_analy_int[zc>=xf] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax2.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_g$')
ax2.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_w$')
ax2.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_s$')

res = solve_ivp(rhs_stage2, (0, 0.005), [-1e-14,1e-15],t_eval=[0.005])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]+z0] = theta_U/phi_U
S_w_analy_int[zc>=res.y[1,0]+z0] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax3.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{g}$')
ax3.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{w}$')
ax3.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{s}$')

ytop = -1
ybot = 0

res = solve_ivp(rhs_stage3, (0, t_interest[3]-tp), [ytop,ybot],t_eval=[t_interest[3]-tp])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]+z0] = theta_U/phi_U
S_w_analy_int[zc>=res.y[1,0]+z0] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax4.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_g$')
ax4.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_w$')
ax4.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_s$')

ytop = -1
ybot = 0

res = solve_ivp(rhs_stage3, (0, t_interest[4]-tp), [ytop,ybot],t_eval=[t_interest[4]-tp])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]+z0] = theta_U/phi_U
S_w_analy_int[zc>=res.y[1,0]+z0] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax5.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax5.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax5.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

ytop = -1
ybot = 0

res = solve_ivp(rhs_stage3, (0, t_interest[5]-tp), [ytop,ybot],t_eval=[t_interest[5]-tp])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]+z0] = theta_U/phi_U
S_w_analy_int[zc>=res.y[1,0]+z0] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax6.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax6.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax6.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

ax3.set_title(r'''%.2f'''%(ts+0.005), fontsize='medium')
ax5.set_title(r'''%.2f'''%t_interest[4], fontsize='medium')
ax1.legend(loc='lower left', shadow=False, fontsize='medium')

plt.xlabel("Volume fractions $\phi$", fontsize='medium')
plt.savefig(f"../Figures/{simulation_name}_CL{theta_U}CR{theta_L}_phiL{phi_U}phiR{phi_L}.pdf")



print("############################################################# \n")
print("2 (top) Volume fraction with depth with dimensions")
print("############################################################# \n")

# First set up the figure, the axis 
fig,([ax1,ax2,ax3,ax4,ax5,ax6]) = time_sequenced_figure_dim_days(zbottom,zsurface,0,1,t_interest,fc_dim,z0_dim)

#Relation (32) - (35) to find the shock speed ratio k
aa = phi_U/phi_L*(1-theta_U/phi_U-s_gr)/(1-theta_L/phi_L-s_gr)*(1- f_Cm(np.array([phi_L]),m,phi_U) * f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n) / f_Cm(np.array([phi_U]),m,phi_U))
bb =-phi_U/phi_L*(1-theta_U/phi_U-s_gr)/(1-theta_L/phi_L-s_gr)*(1- f_Cm(np.array([phi_L]),m,phi_U) * f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n) / f_Cm(np.array([phi_L]),m,phi_U)) - (1 - f_Cm(np.array([phi_U]),m,phi_U) * f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n) / f_Cm(np.array([phi_U]),m,phi_U))
cc = 1 - f_Cm(np.array([phi_U]),m,phi_U) * f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n) / f_Cm(np.array([phi_L]),m,phi_U)
k = (-bb-np.sqrt(bb**2-4*aa*cc))/(2*aa)

#Calculating the shock speeds
s_U_analy =((qs(k) - f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0]/(phi_U*(1-s_gr-theta_U/phi_U)))
s_L_analy =((qs(k) - f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-theta_L/phi_L-s_gr)))

Sf = ((f_Cm(np.array([phi_U]),m,phi_U)[0]*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0] / (phi_U*(theta_U/phi_U - s_wr))) #initial front speed
ts = (z0 - zsurface)/Sf             #dimensionless time of saturation
tp = ts + (zsurface-z0)/s_U_analy   #dimensionless time of ponding

xf = zsurface + Sf*t_interest[0]
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=xf] = theta_U/phi_U
S_w_analy_int[zc>=xf] = theta_L/phi_L

S_w_analy_int_combined = []
S_w_analy_int_combined = S_w_analy_int.copy()

ax1new = ax1.twinx()
ax1new.set_ylim([2,0])
ax1new.set_yticks([])

ax1.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{g}$')
ax1.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{w}$')
ax1.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{s}$')

xf = zsurface + Sf*t_interest[1]
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=xf] = theta_U/phi_U
S_w_analy_int[zc>=xf] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax2.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_g$')
ax2.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_w$')
ax2.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_s$')

res = solve_ivp(rhs_stage2, (0, 0.005), [-1e-14,1e-15],t_eval=[0.005])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]+z0] = theta_U/phi_U
S_w_analy_int[zc>=res.y[1,0]+z0] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax3.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{g}$')
ax3.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{w}$')
ax3.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{s}$')

ytop = -1
ybot = 0

res = solve_ivp(rhs_stage3, (0, t_interest[3]-tp), [ytop,ybot],t_eval=[t_interest[3]-tp])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]+z0] = theta_U/phi_U
S_w_analy_int[zc>=res.y[1,0]+z0] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax4.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_g$')
ax4.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_w$')
ax4.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_s$')

ytop = -1
ybot = 0

res = solve_ivp(rhs_stage3, (0, t_interest[4]-tp), [ytop,ybot],t_eval=[t_interest[4]-tp])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]+z0] = theta_U/phi_U
S_w_analy_int[zc>=res.y[1,0]+z0] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax5.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax5.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax5.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

ytop = -1
ybot = 0

res = solve_ivp(rhs_stage3, (0, t_interest[5]-tp), [ytop,ybot],t_eval=[t_interest[5]-tp])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]+z0] = theta_U/phi_U
S_w_analy_int[zc>=res.y[1,0]+z0] = theta_L/phi_L

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax6.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax6.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax6.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

ax6new = ax6.twinx()
ax6new.set_ylim([2,0])
ax6new.set_ylabel(r"Dimensionless depth, $z/z_0$")

ax1.set_yticklabels(z0_dim*(ax6.get_yticks()))
ax1.set_ylabel(r"Depth, $z$ [m]")

t_dim = (z0_dim /fc_dim * t_interest)/day2s
ts_dim = (z0_dim /fc_dim * ts)/day2s 

ax3.set_title(r"%.2f" "\n" r"%.1f" %(ts,ts_dim), fontsize='medium')
ax5.set_title(r"%.2f" "\n" r"%.1f" %(t_interest[4],t_dim[4]), fontsize='medium')

ax1.legend(loc='lower left', shadow=False, fontsize='medium')

plt.xlabel("Volume fractions $\phi$", fontsize='medium')
plt.savefig(f"../Figures/{simulation_name}_CL{theta_U}CR{theta_L}_phiL{phi_U}phiR{phi_L}.pdf")


print("############################################################# \n")

print("############################################################# \n")
print("2 (bottom) Saturation with depth")
print("############################################################# \n")

fig,([ax1,ax2,ax3,ax4,ax5,ax6]) = time_sequenced_figure(zbottom,zsurface,-0.05,1.05,t_interest)
ax1.plot(S_w_analy_int_combined[:,0],zc , c = 'k',linestyle='--')
ax2.plot(S_w_analy_int_combined[:,1],zc , c = 'k',linestyle='--')
ax3.plot(S_w_analy_int_combined[:,2],zc , c = 'k',linestyle='--')
ax4.plot(S_w_analy_int_combined[:,3],zc , c = 'k',linestyle='--')
ax5.plot(S_w_analy_int_combined[:,4],zc , c = 'k',linestyle='--')
ax6.plot(S_w_analy_int_combined[:,5],zc , c = 'k',linestyle='--')

ax3.set_title(r'''%.2f ($t_s'$)'''%(ts+0.005), fontsize='medium')
ax5.set_title(r'''%.2f '''%t_interest[4], fontsize='medium')
plt.xlabel("Water saturation $s_w$", fontsize='medium')
plt.savefig(f"../Figures/swvsZpanelshock_{simulation_name}_CL{theta_U}CR{theta_L}_phiL{phi_U}phiR{phi_L}.pdf")


print("############################################################# \n")
print("3 Dimensionless infiltration rate with dimensionless time")
print("############################################################# \n")

print('R/fc:',R,''', tp':''',tp,''', zs':''',z0,''', ts':''',ts,'\n')
      
res = solve_ivp(rhs_stage2, (0, tp-ts), [-1e-14,1e-15],t_eval=[tp-ts])
t0 = np.linspace(0,ts,1000)
qs0 = (f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))* np.ones_like(t0)
t1 = np.linspace(ts,tp,1000)
qs1 = (f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))*np.ones_like(t1)
t2 = np.linspace(tp,10,10000)
res= solve_ivp(rhs_stage3, (t2[0],t2[-1]), [ytop,ybot],t_eval=t2)
y  = res.y
qs2 = (y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U) - y[1]/f_Cm(np.array([phi_L]),m,phi_U))
T  = np.concatenate([t0,t1,t2])
QS = np.concatenate([qs0,qs1,qs2])

fig = plt.figure(figsize=(8,8) , dpi=100)
plot = plt.plot(T,QS/f_Cm(np.array([phi_U]),m,phi_U),c='r')
plt.hlines(f_Cm(np.array([phi_L]),m,phi_U), np.min(T), np.max(T), colors=gray, linestyles='--')
plt.vlines(ts, 0, qs0[0]/f_Cm(np.array([phi_U]),m,phi_U), colors=gray, linestyles='--')
plt.vlines(tp, 0, qs0[0]/f_Cm(np.array([phi_U]),m,phi_U), colors=gray, linestyles='--')
plt.hlines(f_Cm(np.array([phi_L]),m,phi_U), np.min(T), np.max(T), colors=gray, linestyles='--')
plt.ylabel(r'''$I(t')/f_c$''', fontsize='medium')
plt.xlabel(r'''$t'$''', fontsize='medium')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
plt.xlim([np.min(T), np.max(T)])
plt.ylim([0,round(np.max(QS),1)+0.1])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f"../Figures/ICvsT_shock_{simulation_name}_CL{theta_U}CR{theta_L}_phiL{phi_U}phiR{phi_L}.pdf")







fig = plt.figure(figsize=(8,8) , dpi=100)
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx().twiny()
ax2.plot(T,QS/f_Cm(np.array([phi_U]),m,phi_U),c='r')
ax2.hlines(f_Cm(np.array([phi_L]),m,phi_U), np.min(T), np.max(T), colors=gray, linestyles='--')
ax2.vlines(ts, 0, qs0[0]/f_Cm(np.array([phi_U]),m,phi_U), colors=gray, linestyles='--')
ax2.vlines(tp, 0, qs0[0]/f_Cm(np.array([phi_U]),m,phi_U), colors=gray, linestyles='--')
ax2.hlines(f_Cm(np.array([phi_L]),m,phi_U), np.min(T), np.max(T), colors=gray, linestyles='--')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
ax2.set_xlim([np.min(T), 2])
ax2.set_ylim([0,round(np.max(QS),1)+0.1])

def tick_function(X):
    
    V = (X * z0_dim / fc_dim)/day2s
    return ["%.1f" % z for z in V]

def ytick_function(Y):
    
    V = (Y * 1e6*fc_dim)
    return ["%.1f" % z for z in V]

ax1.set_xticklabels(tick_function(ax2.get_xticks()))
ax1.set_xlabel(r"time, $t$ [days]")

ax1.set_yticklabels(ytick_function(ax2.get_yticks()))
ax1.set_ylabel(r"Infiltration rate, $I(t) \times 10^{6}$ [m/s]")

ax1.set_xlim(np.array([np.min(T), 2]))
ax1.set_ylim(np.array([0,round(np.max(QS),1)+0.1]))

plt.tight_layout()#pad=0.4, w_pad=0.5, h_pad=1.0)
ax2.set_ylabel(r"$I(t')/f_c$")
ax2.set_xlabel(r"$t'$")

plt.savefig(f"../Figures/newICvsT_shock_{simulation_name}_CL{theta_U}CR{theta_L}_phiL{phi_U}phiR{phi_L}.pdf")



fc_dim= 1e-5 #infiltration capacity [m/s]
z0_dim = 10 #dune depth [m]





















