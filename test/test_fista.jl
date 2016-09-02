# generate synthetic data from Marmousi model
using PyPlot, Seismic
download("http://seismic.physics.ualberta.ca/data/marmvel.bin","marmvel.bin");
n1 = 751; n2 = 2301; dz =4.0; dx = 4.0;
path = join([pwd() "/marmvel.bin"])
v = reshape(read(open(path, "r"), Float32, n1*n2), n1, n2)
rm(path)
v = convert(Array{Float64,2}, 1000.*v);
vp = v[300:end-100, 1:200]
dt = 0.002; f0=20.0;
d  = SeisConv(vp, dx, dt, f0)

# deconvolution
operators = [conv1]
parameters = [Dict(:w=>Ricker(f0, dt))]
(nt, nx) = size(d)
tmp = randn(nt, nx)

# estimate maximum eigenvalue
lambda = power_method(tmp, operators, parameters)
lambda = ceil(lambda)
mu = 0.1
(m, J) = FISTA(d, operators, parameters, mu, lambda)

# Plotting
SeisPlot(d, title="Input" , dx=dx, dy=dt, xlabel="Offset (m)", ylabel="Time (s)", cmap="seismic")
SeisPlot(m, title="Output", dx=dx, dy=dt, xlabel="Offset (m)", ylabel="Time (s)", cmap="seismic")
