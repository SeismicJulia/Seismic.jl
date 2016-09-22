using PyPlot, Seismic

download("http://seismic.physics.ualberta.ca/data/marmvel.bin","marmvel.bin");
n1 = 751; n2 = 2301; dz =4.0; dx = 4.0;
v = reshape(read(open("marmvel.bin", "r"), Float32, n1*n2), n1, n2)
v = convert(Array{Float64,2}, 1000.*v);
vp = v[400:end, 1:500]
(nz, nx) = size(vp)
# receivers location
ix = collect(1:5:nx); iz = ones(Int64, length(ix)); pos = hcat(iz, ix);
# set up SparseMatrixCSC for computing derivative in spatial
f0 = 30.0; dt = 3e-4;
fidMtx = AcousticSetup(vp, dz, dx, dt, f0, iflag=2);
# source location
isx = 251; isz = 1;
shot = SeisAcousticWave(fidMtx, pos, isz, isx, f0, dt, tmax=1.0, print_flag=true, interval=500)
SeisPlot(shot.d, name="ShotGather.pdf")
SeisPlot(vp, name = "vp.pdf")
