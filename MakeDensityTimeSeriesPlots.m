geofile = '/Volumes/Samsung_T5/FromThilo/NISKINEthilo_4km_negLons.nc';
UniqueHistoryFiles=dir('/Volumes/Samsung_T5/FromThilo/niskine_his_*.nc');
UniqueHistory2Files=dir('/Volumes/Samsung_T5/FromThilo/niskine_his2_*.nc');

% Read from geometry file
x_rho = ncread(geofile,'x_rho');
y_rho = ncread(geofile,'y_rho');
h = ncread(geofile,'h');

% Read constant values from first history file
historyfile = fullfile(UniqueHistoryFiles(1).folder,UniqueHistoryFiles(1).name);
Vtransform = ncread(historyfile,'Vtransform');
hc = ncread(historyfile,'hc');
sc_r = ncread(historyfile,'s_rho');
Cs_r = ncread(historyfile,'Cs_r');
sc_w = ncread(historyfile,'s_w');
Cs_w = ncread(historyfile,'Cs_w');

% read in all the time points we have available across files
t = [];
for iFile=1:length(UniqueHistoryFiles)
    t = cat(1,t,ncread(fullfile(UniqueHistoryFiles(iFile).folder,UniqueHistoryFiles(iFile).name),'ocean_time'));
end

% if you pass an empty zeta ([]), then it will ignore the contribution from
% the free-surface. I think this is reasonable when computing a mean
% density?
igrid = 1; % density points
idims = 0; % Fortran based ordering.
zMean = DepthsFromVars(Vtransform, hc, sc_r, Cs_r, sc_w, Cs_w, h, igrid, idims, []);
iX = 200;
iY = 125;

rhoTimeSeries = zeros(size(zMean,3),length(t));
zTimeSeries = zeros(size(zMean,3),length(t));

iTime = 1;
for iFile=1:length(UniqueHistoryFiles)
    historyfile = fullfile(UniqueHistoryFiles(iFile).folder,UniqueHistoryFiles(iFile).name);
    history2file = fullfile(UniqueHistory2Files(iFile).folder,UniqueHistory2Files(iFile).name);
    tfile = ncread(historyfile,'ocean_time');

    for tindex=1:length(tfile)
        % Read time-varying values from history file
        T = nc_read(historyfile,'temp',tindex,NaN);
        S = nc_read(historyfile,'salt',tindex,NaN);

        zeta = nc_read(history2file,'zeta',tindex,NaN);
        Zr = DepthsFromVars(Vtransform, hc, sc_r, Cs_r, sc_w, Cs_w, h, igrid, idims, zeta);

        % The density save in the output is in-situ density, but we need potential
        % density.
        % rho = nc_read(historyfile,'rho',tindex,NaN);
        pden = DensityFromTSZ(T,S,Zr);

        zTimeSeries(:,iTime) = Zr(iX,iY,:);
        rhoTimeSeries(:,iTime) = pden(iX,iY,:);

        iTime = iTime + 1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make a few figures
% 

rhobar = mean(rhoTimeSeries,2);
rhoprime = rhoTimeSeries - rhobar;

z = squeeze(zMean(iX,iY,:));
figure, pcolor(t,z,rhoprime), shading flat

im = InternalModesSpectral(rhobar,z,z,45);
im.upperBoundary = UpperBoundary.freeSurface;
[F,G,h] = im.ModesAtFrequency(0);

a = G(:,1:35)\rhoprime(:,1);
figure, plot(rhoprime(:,1),z), hold on, plot(G(:,1:35)*a,z)