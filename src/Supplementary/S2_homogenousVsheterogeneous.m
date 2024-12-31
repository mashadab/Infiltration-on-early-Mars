%%
%Calculating the mean
phi_surf = 0.5
yr2s = 60^2*24*365.25; % seconds per year
R    = 3389508; % [m] Mars' mean radius
g_mars = 3.711; % [m/s^2] grav. acceleration on Mars
rho = 1e3; % [kg/m^3] desity of water
mu = 1e-3; % [Pa s] water viscosity
theta_p_main = 0%deg2rad(60); % [rad] angle for beginning of the precipitation region 
theta_b_main = pi-acos(1/3); % [rad] angle for dichotomy boundary 
%%
m_exp = 2.5422;
n_exp = 3*m_exp;
%% Finding the hydraulic conductivity K0
dmax  = 10e3; % [km] max. aquifer depth
k0_MI = 10^(-12.65); % [m^2] ref. permeability at 1km depth Manning and Ingebritsen (1999)
%k0_MI_avg = @(n) k0_MI * dmax^n / ((dmax-1e3).^n*(n+1)) %To calculate the
%average k0_MI for homogeneous aquifer
%homogeneous equivalent of k0_MI
alpha = 3.2;
k_MI  = @(d) k0_MI*d.^(-alpha);
k0    = @(n) k0_MI/(dmax-1e3)^n;
K0_func = @(n) k0(n)*rho*g_mars/mu;

K0 = K0_func(n_exp);

phi0  = phi_surf / dmax.^m_exp; %base porosity dimensional
phi_mean = phi0*dmax.^(m_exp + 1) / ((m_exp+1)*dmax);

k_mean= k0(n_exp)*dmax.^(n_exp + 1) / ((n_exp+1)*dmax);

K0_mean = k_mean*rho*g_mars/mu;

%Tiled plot
set(0,'defaultAxesFontSize',35)
fig1d = figure('position',[100 100 250 600])
ylim([-9 1])
zz = linspace(0,dmax,10000)
y = (zz/dmax).^m_exp;
t = tiledlayout(1,1);
ax1 = axes(t);
plot(y,zz/1e3-9,'-r');
xline(phi_mean/phi_surf,'r--');
%ax1.XColor = 'r';
%xlabel(ax1,'$\frac{\phi}{\phi_s}$ [-]','Interpreter','latex','Color','r')
xlabel(ax1,'$\frac{K}{K_s}$ [-]','Interpreter','latex','Color','b')

ylabel('Elevation, $z-9$ [km]','Interpreter','latex')
ylim([-9 1])
ax2 = axes(t);
r = (zz/dmax).^n_exp;
plot(r,zz/1e3-9,'-k');
hold on
xline(k_mean/(k0(n_exp)*dmax.^n_exp),'k--');
hold on
plot(r,zz/1e3-9,'-b');
xline(k_mean/(k0(n_exp)*dmax.^n_exp),'b--');
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
%ax2.XColor = 'b';
ax1.Box = 'off';
ax2.Box = 'off';
%set(gca,'ytick',[])
legend('Actual','Mean','Location','southeast','color','none','FontSize',20)
% set(fig1d,'Units','Inches');
% pos = get(fig1d,'Position');
% set(fig1d,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
ylim([-9 1])
print(fig1d,'kandphi.pdf','-dpdf','-r0')

%averaging and flux
z0 = 10e3;
phi0 = 0.5;
z  = linspace(0,z0,1000);
%%vertical averaging
I_v = 1./z .* (z0 - (z0 - z).*(1-z/z0).^(n_exp))/(n_exp + 1);

%%harmonic averaging
I_h = z ./ (-z0./(n_exp - 1) + (z0 - z).*(1-z/z0).^(-n_exp)./(n_exp - 1));

set(0,'defaultAxesFontSize',35)
set(0,'defaultAxesFontName','Times')
fig2d = figure('position',[100 100 350 600])
plot(I_v,z/1e3,'r-',linewidth=4);
hold on
plot(I_h,z/1e3,'b-',linewidth=4);
set(gca, 'YDir','reverse')
xline(I_v(end),'k--',linewidth=4)
xlabel('${I}/{K_0}$ [-]','Interpreter','latex')
ylabel('Depth, $z$ [km]','Interpreter','latex')
legend('Arithmetic mean', 'Harmonic mean','Average','Location','southeast','color','none','FontSize',20,'Interpreter','latex')
print(fig2d,'arithmetic_vs_harmonic_conductivity.pdf','-dpdf','-r0')

%%Porosity decay
phi_array = (1-z/z0).^m_exp;
phi_avg = 1./(m_exp+1);

set(0,'defaultAxesFontSize',35)
set(0,'defaultAxesFontName','Times')
fig2d = figure('position',[100 100 350 600])
plot(I_v(end)/phi_avg*ones(1000,1),z/1e3,'r-',linewidth=4); %homo
hold on
plot(I_h./phi_array,z/1e3,'b-',linewidth=4); %hetero
set(gca, 'YDir','reverse')
xlabel('$dz/dt$ [-]','Interpreter','latex')
ylabel('Depth, $z$ [km]','Interpreter','latex')
legend('Arithmetic mean', 'Harmonic mean','Average','Location','southeast','color','none','FontSize',20,'Interpreter','latex')
print(fig2d,'arithmetic_vs_harmonic_shock_speed.pdf','-dpdf','-r0')



%%Porosity decay
phi_array = (1-z/z0).^m_exp;
phi_avg = 1./(m_exp+1);

set(0,'defaultAxesFontSize',35)
set(0,'defaultAxesFontName','Times')
fig2d = figure('position',[100 100 350 600])
plot(phi_avg*ones(1000,1),z/1e3,'r-',linewidth=4); %homo
hold on
plot(phi_array,z/1e3,'b-',linewidth=4); %hetero
set(gca, 'YDir','reverse')
xlabel('$\phi/\phi_s$ [-]','Interpreter','latex')
ylabel('Depth, $z$ [km]','Interpreter','latex')
legend('Average', 'Harmonic mean','Location','southeast','color','none','FontSize',20,'Interpreter','latex')
print(fig2d,'arithmetic_vs_harmonic_porosity.pdf','-dpdf','-r0')

