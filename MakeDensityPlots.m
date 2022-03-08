geofile = '/Volumes/Samsung_T5/FromThilo/NISKINEthilo_4km_negLons.nc';
historyfile = '/Volumes/Samsung_T5/FromThilo/niskine_his_00040.nc';
history2file = '/Volumes/Samsung_T5/FromThilo/niskine_his2_00040.nc';

% Read from geometry file
x_rho = ncread(geofile,'x_rho');
y_rho = ncread(geofile,'y_rho');
h = ncread(geofile,'h');

% Read constant values from history file
Vtransform = ncread(historyfile,'Vtransform');
hc = ncread(historyfile,'hc');
sc_r = ncread(historyfile,'s_rho');
Cs_r = ncread(historyfile,'Cs_r');
sc_w = ncread(historyfile,'s_w');
Cs_w = ncread(historyfile,'Cs_w');

tindex = 12;
igrid = 1; % density points
idims = 0; % Fortran based ordering.

% Read time-varying values from history 2 file
zeta = nc_read(history2file,'zeta',tindex,NaN);

Zr = DepthsFromVars(Vtransform, hc, sc_r, Cs_r, sc_w, Cs_w, h, igrid, idims, zeta);

% Read time-varying values from history file
T = nc_read(historyfile,'temp',tindex,NaN);
S = nc_read(historyfile,'salt',tindex,NaN);

[pden,den] = DensityFromTSZ(T,S,Zr);
rho = nc_read(historyfile,'rho',tindex,NaN);

x = repmat(squeeze(x_rho(:,100,:)),1,size(Zr,3));

figure
subplot(2,1,1)
pcolor(x/1e3,squeeze(Zr(:,100,:)),squeeze(pden(:,100,:))), shading interp
colormap(flip(colormap('parula')));
ylabel('depth')
xlabel('x (km)')
colorbar('eastoutside')


y = repmat(shiftdim(y_rho(175,:,:),1),1,size(Zr,3));
subplot(2,1,2)
pcolor(y/1e3,squeeze(Zr(175,:,:)),squeeze(pden(175,:,:))), shading interp
ylabel('depth')
xlabel('y (km)')
colorbar('eastoutside')