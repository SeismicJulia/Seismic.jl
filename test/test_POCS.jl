using Seismic
using Base.Test

d,e = SeisLinearEvents()
ddec = SeisDecimate(d,perc=50)

dpocs = SeisPOCS(ddec,dt=0.004,fmax=100,Niter=100,p=1.5)

# test that quality factor is greater than 10 Decibels
quality_factor = 10*log10(norm(d[:],2)/norm(dpocs[:]-d[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 5

