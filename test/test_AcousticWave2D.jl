using PyPlot, Seismic
nz = 201; nx = 201; dz = 5.; dx = 5.; dt = 5e-4; f0 = 30.0
v  = 2500. * ones(nz, nx); v[101:end,:] = 4000.0
# receivers location
ix = collect(1:2:nx); iz = ones(Int64, length(ix)); pos = hcat(iz, ix);
# set up SparseMatrixCSC for computing derivative in spatial
fidMtx = AcousticSetup(v, dz, dx, dt, f0, iflag=2)
# source location
isx = 101; isz = 1;
shot = SeisAcousticWave(fidMtx, pos, isz, isx, f0, dt, tmax=2.0)
SeisPlot(shot.d)
