WaveVortexModel with ROMS
==========================

March 9th, 2022
---------------

Thilo,

The overall goal here is to apply the wave-vortex model to your ROMS output. The projection requires Fourier transforms in the horizontal, and a projection onto the internal modes in the vertical. To compute the vertical modes we need a mean density at each location. In practice, the modal projection is not so sensitive that we need a mean density at *each* location, so we can easily bin nearby locations if that is computationally helpful.

To compute the mean density from ROMS output I used the [scripts available from svn at myroms.org](https://www.myroms.org/wiki/Matlab_Processing_Scripts). Those scripts assume that all the output data will be in a single netcdf file, so it is hard to use them directly for the NISKINE run, which separates the interior output from the surface output. So I copied the relevant chunks and made my own functions, ```DensityFromTSZ``` and ```DepthsFromVars```.

If you run my script ```MakeDensitySpatialSlicesPlots.m``` it should show you two planar slices of density from a given time point and file, using the above two functions. Please note that ```ncread``` is the built-in matlab command and ```nc_read``` is the ROMS version which does a few extra things.

The script ```MakeDensityTimeSeriesPlots.m``` is similar, but now grabs a whole time-series of density data (as a function of time and depth) for a given location.

Finally, the script ```MakeMeanDensity.m``` just reads the whole folder of ROMS output and tries to make a mean density profile and write it out.

So, if you could try to use that ```MakeMeanDensity.m``` script to create some kind of reasonable mean, that would be awesome. I don't think you need to grab every time point from every file, as I have done in the example. I think that if you simply ensemble together enough profiles over a long enough time period, it should work fine. The time period should be long enough that there is some geostrophic variation, not just internal wave.

After this is done, we then need to actually do the projection! There will a few details to sort out there, as well.
