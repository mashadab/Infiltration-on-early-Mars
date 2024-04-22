%Set default values
set(groot,'defaultAxesFontName','Times')
set(groot,'defaultAxesFontSize',20)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot, 'DefaultFigureVisible', 'on');

%Load file
%Loading topography file
load('./Data/Mars_topography.mat')
%Loading groundwater table height
load('./Data/thirty_to_negthirty_Arabia_shore_10_microns_per_year_recharge.mat')

%Fixing the contour plot
mars_topo1 = mars_topo;
mars_topo2 = mars_topo;
mars_topo(1:720,:) = mars_topo1(721:end,:);
mars_topo(721:end,:) = mars_topo1(1:720,:);
mars_topo = flipud(mars_topo);

% Need to use Xc Yc in data
zm1 = zm;
zm2 = zm;
zm2(1:150,:) = zm1(151:end,:);
zm2(151:end,:) = zm1(1:150,:);
zm2 = flipud(zm2);

%subtract Mars topo from zm 
mars_topo = interp2(Theta,Phi,mars_topo,Xc,Yc) ; %Mars topography at a certain point
Arabia = find(mars_topo>=2080.0 & mars_topo<2100.0); %Arabia shoreline (m)
zm_subtract = zm2-mars_topo;   %Difference in the topo and GW elevation
mars_topoplot = mars_topo;
mars_topoplot(mars_topoplot>9e3) = nan;

%Power law decay
g_mars = 3.711; % [m/s^2] grav. acceleration on Mars
rho = 1e3; % [kg/m^3] desity of water
mu = 1e-3; % [Pa s] water viscosity
n = 2; %Corey-Brooks coefficient krw = krw0 * sw^n
s_wr = 0.0; %Residual water saturation
s_gr = 0.0; %Residual gas saturation
z0 = 0.06  ; %point specified for power law porosity variation
p  = 2.5422; %index specified for power law porosity variation 
m_exp = p; %This is for Mars work: depth variation
m = 3; %This is for Mars work: Cozeny-Karman coefficient for numerator K = K0 (1-phi_i)^m
n_exp = m; %This is for Mars work: depth variation
yr2s = 365.25*60*60*24;

%Calculate porosity at the surface phi0
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
phi_surf = 0.5;
phi_top_surf = phi_surf;
K0 = K0_func(n_exp);
phi0  = phi_surf / dmax.^m_exp; %base porosity dimensional for elevation dependence


%What is the zbase for topo case?
zbase = - 9e3; %elevation of the base is at -9 km
phi_surf = phi0.*((mars_topo - zbase)).^m_exp ; %Porosity at surface everywhere 

phi0 = phi_surf; %surface porosity for depth dependence

%Calculate ponding time everywhere
%Equivalent Homogeneous 
zm_subtract = zm2-mars_topoplot;
zm_subtract(zm_subtract>0) =  nan;

t_base = (0-zm_subtract)/K*(phi_top_surf/(p+1)*(1-s_gr-s_wr)) ;  %time_to_reach_the_aquifer


%Vertical variation
n= 2; m=3; p=2.5422;
tbase = 0;
z0_calc = (mars_topo - zbase+1e-14);
Ratio = (mars_topoplot - zm2)./z0_calc;
fc_basalt = 1e-6 %Infiltration capacity for basalt [m/s]
%Removing case where table is above the ground
Ratio(Ratio<0) = nan;
t_inf_vertical_var = (1-Ratio.^(((-m*p + n*p + n))/n)).*(n/(-m*p + n*p + n)).*phi_top_surf.*z0_calc.*(1-s_wr-s_gr)/fc_basalt  ; 
t_inf_vertical_var(isinf(t_inf_vertical_var)) = nan;


%New calculation plot

hh= figure('Position', [1, 1, 500, 800]);
h = tiledlayout(3,1, 'TileSpacing', 'none', 'Padding', 'none');
nexttile
%contourf(Yc,Xc,zm_subtract); %Plotting subtracted topo; Yc is longitude, Xc is latitute
%Calculate ponding time everywhere
%zm_subtract

contourf(-180+rad2deg(Yc),-90+rad2deg(Xc),zm2 - mars_topoplot,50,'edgecolor','none'); %GW - topo in m
hold on
shoreline = table2array(readtable('/Users/afzal-admin/Documents/Research/mars-project/Afzal/Infiltration/Infiltration/Mars-map/Shorelines.csv'));
plot(shoreline(:,1),shoreline(:,2),'k.','MarkerSize',6);
%plot(-90+rad2deg(Phi(Arabia)),-180+rad2deg(Theta(Arabia)),'k.','MarkerSize',10)
yline(-45,'k--','Linewidth',3);
yline(45,'k--','Linewidth',3);
pbaspect([1.8 1 1])
%set(gca,'ColorScale','log')
ylabel('Latitude [$$^\textrm{o}$$]');
c1 = colorbar()
c1.Label.String = 'G.W. - Topo. elev. [m]';
colormap(flipud(turbo))
nexttile
contourf(-180+rad2deg(Yc),-90+rad2deg(Xc),log10(t_base/yr2s),50,'edgecolor','none'); %time in log of years
hold on
plot(shoreline(:,1),shoreline(:,2),'k.','MarkerSize',6);
pbaspect([1.8 1 1])        
colormap(flipud(turbo))
caxis([1.5 3]);
%set(gca,'ColorScale','log')
c1 = colorbar()
c1.Label.String = 'Infiltration time [log_{10}(years)]';
ylabel('Latitude [$$^\textrm{o}$$]');

nexttile
contourf(-180+rad2deg(Yc),-90+rad2deg(Xc),log10(t_inf_vertical_var/yr2s),50,'edgecolor','none'); %time in log of years
pbaspect([1.8 1 1])        
colormap(flipud(turbo))
caxis([1.5 3]);
hold on
plot(shoreline(:,1),shoreline(:,2),'k.','MarkerSize',6);
%set(gca,'ColorScale','log')
c2 = colorbar()
c2.Label.String = 'Infiltration time [log_{10}(years)]';
ylabel('Latitude [$$^\textrm{o}$$]');
xlabel('Longitude [$$^\textrm{o}$$]');
set(gcf,'PaperType','A4')
saveas(h,'./res_fig_final.pdf');
print(hh,'-painters','-opengl', '-r1000','../Figures/Figure3.pdf','-dpdf','-fillpage');


%GEL calculation
GEL_angle = 45; %Band angle degrees
diff_new = zm2 - mars_topoplot;
diff_new(diff_new>0)= 0;
diff_new(isnan(diff_new)) = 0;
diff_new(abs(-90+rad2deg(Xc))>GEL_angle) = 0;
diff_no_GW = zbase - mars_topoplot;
diff_no_GW(abs(-90+rad2deg(Xc))>GEL_angle) = 0;
diff_no_GW(diff_new>0)= 0;
diff_no_GW(isnan(diff_no_GW)) = 0;

figure
contourf(-180+rad2deg(Yc),-90+rad2deg(Xc),diff_new,50,'edgecolor','none'); %GW - topo in m
colorbar

figure
contourf(-180+rad2deg(Yc),-90+rad2deg(Xc),diff_no_GW,50,'edgecolor','none'); %no GW - topo in m
colorbar

r = 3389.5; SA_Mars = 4.* pi.*r.^2; %radius  [km] and SA [km^2] of Mars
XXc_plus = Xc + (Xc(1,2)-Xc(1,1))/2; XXc_minus = Xc - (Xc(1,2)-Xc(1,1))/2; %longitude [radians]
YYc_plus = Yc + (Yc(2,1)-Yc(1,1))/2; YYc_minus = Yc - (Yc(2,1)-Yc(1,1))/2; %latitude [radians]
Grid_area = pi * r.^2 .* (sin(XXc_plus) - sin(XXc_minus)).* (YYc_plus - YYc_minus); %area of the grid [m^2]
Vol_GW = sum(Grid_area.*(-diff_new).*phi_top_surf/(p+1)*(1-s_gr-s_wr),'all'); GEL_GW = Vol_GW./SA_Mars %GEL without Groundwater table [m]
Vol_no_GW = sum(Grid_area.*(-diff_no_GW).*phi_top_surf/(p+1)*(1-s_gr-s_wr),'all'); GEL_no_GW = Vol_no_GW./SA_Mars %GEL without Groundwater table [m]
