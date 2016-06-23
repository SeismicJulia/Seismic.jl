# Multiple supresion using the parabolic Radon transform. Synthetic example.

using PyPlot, Seismic

close("all")

# 1-  Make an ideal Radon Gather 
source = Ricker()
ns = length(source)
np = 100
nt = 600
nw = nt
nh = 80
dt = 0.004
h = collect(linspace(0.0, 1000.0, nh))
p = collect(linspace(-0.04, 2.2, np))
m = zeros(Float64, nt, np)
for is = 1:ns
    m[160+is,  3] =  source[is]
    m[ 60+is, 44] = -source[is]
    m[280+is, 44] = -source[is]
    m[240+is,  3] = -source[is]
end
href = maximum(abs(h))

# 2- Model the data using the forward Radon transform 
param = Dict(:order=>"parab", :dt=>dt, :p=>p, :h=>h, :href=>href, :flow=>0.0,
             :fhigh=>125.0)
d = SeisRadonFreqFor(m, nt; param...)

# 3- Recover Radon gather via LS inversion
m = SeisRadonFreqInv(d; param..., mu=0.001)

# 4- Filter primaries and keep multiples in Radon gather  
mf = copy(m);
pcut = 0.3
icut = find(p .< pcut)
mf[:, icut] = 0.0

# 5- Model multiples by forward modelling
d_mult = SeisRadonFreqFor(mf, nt; param...)
 
# 6- Substract multiples from data to obtain primaries
d_prim = d - d_mult 
dp = p[2] - p[1]
dh = h[2] - h[1]

# 7- Plotting  
figure(1, figsize=(12,6));
subplot(141)
SeisPlot(d, title="Data", xlabel="Offset [m]", ylabel="Time [s]",
         vmin=-0.8, vmax = 0.8, ox=h[1], dx=dh, dy=dt, fignum=1)
subplot(142)
SeisPlot(m, title="Radon gather", xlabel="Residual moveout [s]",
         ox=p[1], dx=dp, dy=dt, fignum=1)
subplot(143)
SeisPlot(d_mult, title="Multiples", xlabel="Offset [m]",
         vmin=-0.8, vmax = 0.8, ox=h[1], dx=dh, dy=dt, fignum=1)
subplot(144)
SeisPlot(d_prim, title="Primaries", xlabel="Offset [m]",
         vmin=-0.8, vmax = 0.8, ox=h[1], dx=dh, dy=dt, fignum=1)
tight_layout()
