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

% if you pass an empty zeta ([]), then it will ignore the contribution from
% the free-surface. I think this is reasonable when computing a mean
% density?
Zr = DepthsFromVars(Vtransform, hc, sc_r, Cs_r, sc_w, Cs_w, h, igrid, idims, zeta);

% Read time-varying values from history file
T = nc_read(historyfile,'temp',tindex,NaN);
S = nc_read(historyfile,'salt',tindex,NaN);

% The density save in the output is in-situ density, but we need potential
% density.
% rho = nc_read(historyfile,'rho',tindex,NaN);
pden = DensityFromTSZ(T,S,Zr);

iX = 200;
iY = 125;

x = repmat(squeeze(x_rho(:,iY)),1,size(Zr,3));
y = repmat(shiftdim(y_rho(iX,:),1),1,size(Zr,3));

figure
subplot(2,3,[1 4])
pcolor(x_rho/1e3,y_rho/1e3,pden(:,:,end)), shading interp
hold on
plot(x_rho(iX,:)/1e3,y_rho(iX,:)/1e3,'LineWidth',2,'Color','black')
plot(x_rho(:,iY)/1e3,y_rho(:,iY)/1e3,'LineWidth',2,'Color','black')
colormap(flip(colormap('parula')));
xlabel('x (km)')
ylabel('y (km)')
caxis([1025.5 1028])

subplot(2,3,[2 3])
pcolor(x/1e3,squeeze(Zr(:,iY,:)),squeeze(pden(:,iY,:))), shading interp
ylabel('depth')
xlabel('x (km)')
colorbar('eastoutside')
caxis([1025.5 1028])


subplot(2,3,[5 6])
pcolor(y/1e3,squeeze(Zr(iX,:,:)),squeeze(pden(iX,:,:))), shading interp
ylabel('depth')
xlabel('y (km)')
colorbar('eastoutside')
caxis([1025.5 1028])