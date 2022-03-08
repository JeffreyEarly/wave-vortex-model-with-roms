% Just playing with the data to try to figure out what is going on
geofile = '/Volumes/Samsung_T5/FromThilo/NISKINEthilo_4km_negLons.nc';
historyfile = '/Volumes/Samsung_T5/FromThilo/niskine_his_00040.nc';
history2file = '/Volumes/Samsung_T5/FromThilo/niskine_his2_00040.nc';

% because zeta (the free surface height) is in his2, but s_rho is in his,
% this function cannot be used to automatically grab all the grid
% information. Annoyingly.
Gout = get_roms_grid(geofile, historyfile);

Gout = get_roms_grid(geofile, {historyfile,history2file});

% specifying t=0 will ignore zeta
z = depths(historyfile,geofile,1,0,0);

figure, pcolor(Gout.x_rho,Gout.y_rho,z), shading interp
figure, pcolor(Gout.x_rho,Gout.y_rho,Gout.z_r), shading interp

figure, pcolor(Gout.x_rho,Gout.y_rho,z(:,:,50)), shading interp
% figure, plot(Gout.x_rho(:,1),squeeze(z(:,1,:)))


tindex = 1;
Zr = depths(historyfile,geofile,1,0,0);
T = nc_read(historyfile,'temp',tindex);
S = nc_read(historyfile,'salt',tindex);
[pden,den] = densityFromTSZ(T,S,Zr);
rho = nc_read(historyfile,'rho',tindex);

figure, pcolor(Gout.x_rho,Gout.y_rho,rho(:,:,50)), shading interp

% https://www.myroms.org/forum/viewtopic.php?p=23209&hilit=potential+density#p23209
% potential density can be computed from eos.m

% F=plot_field(Gout, historyfile, 'rho', 1, 1);


% 
% grd = roms_get_grid(geofile);
% roms_zview(historyfile,'temp',1,-10,roms_get_grid)